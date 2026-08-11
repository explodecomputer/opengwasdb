"""Integrated Hybrid build — one read per study, routing on-panel → dense fill,
off-panel → ragged overflow (ADR 0026, PRD "Generation").

This is a **thin integration layer** (PRD "Component reuse"): it composes the
dense two-pass builder's liftover, key-lookup, disk-spill, band-write and
fork-pool machinery (``opengwasdb.layouts.dense.build_vcf``) and the ragged
``RaggedCSRWriter`` — the only new logic is per-variant on-panel/off-panel
routing in the single Pass 2 read that both components share.

Like the dense builder, association streaming and the union-variant pass go
through a ``SourceReader`` resolved from each row's ``source_reader_capability``
(issue #20) rather than importing ``opengwasdb.build.vcf_source`` directly.
"""

from __future__ import annotations

import json
import logging
import shutil
import tempfile
import time
from concurrent.futures import as_completed
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from opengwasdb.build.liftover import LiftoverFailureError, build_liftover_lookup
from opengwasdb.layouts.dense.build import AnalysisMetadata
from opengwasdb.layouts.dense.build_vcf import (
    _RESOLVE_BATCH,
    _alid_sort_key,
    _apply_se_divisor,
    _create_dense_zarr,
    _encode_variant_keys,
    _fork_pool,
    _log_progress,
    _read_manifest,
    _write_dense_bands,
    _write_index,
)
from opengwasdb.layouts.dense.constants import (
    DEFAULT_CHUNK_SHAPE,
    DEFAULT_COMPRESSOR,
    DEFAULT_DTYPE,
)
from opengwasdb.layouts.dense.top_hits import write_top_hit_indexes
from opengwasdb.layouts.hybrid.layout import (
    DENSE_SUBDIR,
    dense_component_path,
    dense_to_shared_path,
)
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRWriter
from opengwasdb.model.enums import (
    AssociationCoverage,
    CompletionState,
    PrimaryStorageLayout,
    StoredEffectScale,
)
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.readers.gwas_vcf import GWAS_VCF_CAPABILITY
from opengwasdb.readers.registry import resolve_reader
from opengwasdb.variants import CanonicalVariant, write_variant_axis

log = logging.getLogger(__name__)

__all__ = [
    "build_hybrid_from_vcf_manifest",
    "HybridBuildResult",
    "read_reference_panel_alids",
    "LiftoverFailureError",
]


@dataclass(frozen=True)
class HybridBuildResult:
    output_path: Path
    n_variants: int          # shared union table
    n_analyses: int
    n_panel: int             # Dense Component rows
    n_off_panel: int         # off-panel observed variants (overflow variant axis)
    n_overflow: int          # overflow associations


# ── Reference panel input ────────────────────────────────────────────────────


def read_reference_panel_alids(panel_path: str | Path) -> set[str]:
    """Read the reference-panel variant set as canonical hg38 ALIDs.

    Accepts either a plain-text file of one ``chr:pos:a1:a2`` ALID per line, or a
    Store Variant Table (``*.tsv.gz`` with an ``alid`` column). The panel defines
    the Dense Component axis (the imputable set).
    """
    path = Path(panel_path)
    name = str(path).lower()
    alids: set[str] = set()
    if name.endswith((".tsv.gz", ".tsv.bgz", ".bgz")):
        import pysam

        with pysam.BGZFile(str(path), "r") as handle:  # type: ignore[call-arg]
            header: list[str] | None = None
            for raw in handle:
                line = raw.decode("utf-8").rstrip("\n")
                if line.startswith("#"):
                    header = line.lstrip("#").split("\t")
                    continue
                fields = line.split("\t")
                if header is not None and "alid" in header:
                    alids.add(fields[header.index("alid")])
                elif len(fields) >= 6:
                    alids.add(fields[5])
        return alids
    opener = open
    if name.endswith(".gz"):
        import gzip

        opener = gzip.open  # type: ignore[assignment]
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            token = line.strip().split()[0] if line.strip() else ""
            if token and not token.startswith("#") and token.count(":") == 3:
                alids.add(token)
    return alids


# ── Pass 2 fork-inherited routing lookup (see dense.build_vcf for the rationale) ─
_pass2_keys_sorted: np.ndarray | None = None
_pass2_targets_sorted: np.ndarray | None = None   # int64: dense row (panel) or shared idx (off)
_pass2_ispanel_sorted: np.ndarray | None = None   # bool: True = on-panel target
_pass2_spill_dir: Path | None = None


def _build_routing_index(
    hg19_lookup: dict[tuple[str, int, str, str], str],
    dense_row: dict[str, int],
    shared_index: dict[str, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compose the liftover + partition into a single fork-safe sorted lookup.

    Returns ``(keys_sorted, targets_sorted, ispanel_sorted)``: for every hg19
    source variant that lifts to a stored variant, its byte key, the target index
    (dense row when on-panel, shared variant_index when off-panel), and whether
    it is on-panel. Workers binary-search this once per association.
    """
    chroms: list[str] = []
    poss: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    targets: list[int] = []
    ispanel: list[bool] = []
    for (chrom, pos, ref, alt), alid in hg19_lookup.items():
        row = dense_row.get(alid)
        if row is not None:
            targets.append(row)
            ispanel.append(True)
        else:
            sidx = shared_index.get(alid)
            if sidx is None:
                continue
            targets.append(sidx)
            ispanel.append(False)
        chroms.append(chrom)
        poss.append(pos)
        refs.append(ref)
        alts.append(alt)
    keys = _encode_variant_keys(chroms, poss, refs, alts)
    targets_arr = np.array(targets, dtype=np.int64)
    ispanel_arr = np.array(ispanel, dtype=bool)
    order = np.argsort(keys, kind="stable")
    return keys[order], targets_arr[order], ispanel_arr[order]


def _dedup_last_wins(
    target: np.ndarray, z: np.ndarray, se: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Keep the last stream occurrence per target index (matches dense semantics)."""
    if len(target) == 0:
        return target, z, se
    _, first_in_rev = np.unique(target[::-1], return_index=True)
    keep = np.sort(len(target) - 1 - first_in_rev)
    return target[keep], z[keep], se[keep]


def _resolve_column_hybrid(
    file_path: str,
    keys_sorted: np.ndarray,
    targets_sorted: np.ndarray,
    ispanel_sorted: np.ndarray,
    se_divisor: float = 1.0,
    *,
    capability: str = GWAS_VCF_CAPABILITY,
    stored_effect_scale: str = StoredEffectScale.SD.value,
) -> tuple[
    tuple[np.ndarray, np.ndarray, np.ndarray],
    tuple[np.ndarray, np.ndarray, np.ndarray],
]:
    """Stream one study once, routing each association to the dense fill (on-panel)
    or the ragged overflow (off-panel). Returns ``(dense, overflow)`` where each is
    ``(index int64, z f32, se f32)`` deduped last-wins by index.

    ``capability`` resolves a ``SourceReader`` (issue #20) rather than this
    module streaming a VCF itself; ``stored_effect_scale`` is required to
    construct one (see ``dense.build_vcf._resolve_column``).

    ``se_divisor`` divides every returned ``se`` value, dense and overflow alike
    (continuous-trait phenotype-SD standardisation, issue #18): a study's SD
    rescaling applies uniformly regardless of which component an association
    routes to. Defaults to 1.0 (no-op).
    """
    reader = resolve_reader(capability, file_path, StoredEffectScale(stored_effect_scale))
    d_idx: list[np.ndarray] = []
    d_z: list[np.ndarray] = []
    d_se: list[np.ndarray] = []
    o_idx: list[np.ndarray] = []
    o_z: list[np.ndarray] = []
    o_se: list[np.ndarray] = []

    chroms: list[str] = []
    poss: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    zs: list[float] = []
    ses: list[float] = []

    def _flush() -> None:
        if not zs:
            return
        query = _encode_variant_keys(chroms, poss, refs, alts)
        idx = np.searchsorted(keys_sorted, query)
        idx_clip = np.minimum(idx, len(keys_sorted) - 1)
        matched = keys_sorted[idx_clip] == query
        tgt = targets_sorted[idx_clip[matched]]
        panel = ispanel_sorted[idx_clip[matched]]
        z_arr = np.array(zs, dtype=np.float32)[matched]
        se_arr = np.array(ses, dtype=np.float32)[matched]
        if panel.any():
            d_idx.append(tgt[panel])
            d_z.append(z_arr[panel])
            d_se.append(se_arr[panel])
        off = ~panel
        if off.any():
            o_idx.append(tgt[off])
            o_z.append(z_arr[off])
            o_se.append(se_arr[off])
        for lst in (chroms, poss, refs, alts, zs, ses):
            lst.clear()

    for assoc in reader.stream_associations():
        chroms.append(assoc.chromosome)
        poss.append(assoc.position)
        refs.append(assoc.ref)
        alts.append(assoc.alt)
        zs.append(assoc.z)
        ses.append(assoc.se)
        if len(zs) >= _RESOLVE_BATCH:
            _flush()
    _flush()

    def _assemble(
        parts_i: list[np.ndarray], parts_z: list[np.ndarray], parts_se: list[np.ndarray]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not parts_i:
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float32),
                np.empty(0, dtype=np.float32),
            )
        idx, z, se = _dedup_last_wins(
            np.concatenate(parts_i), np.concatenate(parts_z), np.concatenate(parts_se)
        )
        return idx, z, _apply_se_divisor(se, se_divisor)

    return _assemble(d_idx, d_z, d_se), _assemble(o_idx, o_z, o_se)


def _spill_hybrid_column(
    spill_dir: Path,
    col_idx: int,
    dense: tuple[np.ndarray, np.ndarray, np.ndarray],
    overflow: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    """Spill one resolved study column: dense rows to ``{col}.npz`` (the layout the
    dense band-writer consumes) and overflow to ``{col}.ovf.npz``."""
    d_rows, d_z, d_se = dense
    o_idx, o_z, o_se = overflow
    for suffix, arrs in (
        ("", {"rows": d_rows, "z": d_z, "se": d_se}),
        (".ovf", {"variant_index": o_idx, "z": o_z, "se": o_se}),
    ):
        final = spill_dir / f"{col_idx}{suffix}.npz"
        tmp = spill_dir / f"{col_idx}{suffix}.tmp.npz"
        np.savez(tmp, **arrs)
        tmp.replace(final)


def _pass2_worker(task: tuple[int, str, float, str, str]) -> int:
    assert _pass2_keys_sorted is not None
    assert _pass2_targets_sorted is not None
    assert _pass2_ispanel_sorted is not None
    assert _pass2_spill_dir is not None
    col_idx, file_path, se_divisor, capability, stored_effect_scale = task
    dense, overflow = _resolve_column_hybrid(
        file_path,
        _pass2_keys_sorted,
        _pass2_targets_sorted,
        _pass2_ispanel_sorted,
        se_divisor,
        capability=capability,
        stored_effect_scale=stored_effect_scale,
    )
    _spill_hybrid_column(_pass2_spill_dir, col_idx, dense, overflow)
    return col_idx


def _assemble_overflow_csr(
    spill_dir: Path, n_analyses: int
) -> RaggedCSRWriter:
    """Assemble the overflow CSR from per-column ``.ovf.npz`` spills, in analysis
    order (so CSR offsets align with analysis_index)."""
    csr = RaggedCSRWriter()
    for col in range(n_analyses):
        path = spill_dir / f"{col}.ovf.npz"
        if not path.exists():
            csr.add_analysis(
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float16),
                np.empty(0, dtype=np.float16),
            )
            continue
        with np.load(path) as data:
            vi = data["variant_index"].astype(np.int32)
            z = data["z"].astype(np.float16)
            se = data["se"].astype(np.float16)
        # Sort by variant_index for consistent within-analysis ordering (matches
        # the ragged BESD builder and lets top-hit CSR cross-validation searchsort).
        order = np.argsort(vi, kind="stable")
        csr.add_analysis(vi[order], z[order], se[order])
        path.unlink()
    return csr


# ── Public build entry point ─────────────────────────────────────────────────


def build_hybrid_from_vcf_manifest(
    manifest_path: str | Path,
    output_path: str | Path,
    *,
    reference_panel: str | Path,
    chain_file: str | Path | None = None,
    store_id: str,
    release_id: str,
    liftover_failure_threshold: float = 0.01,
    chunk_shape: tuple[int, int] = DEFAULT_CHUNK_SHAPE,
    dtype: str = DEFAULT_DTYPE,
    overwrite: bool = False,
    n_workers: int = 1,
) -> HybridBuildResult:
    """Build a Hybrid store from a manifest of GWAS-VCF files and a reference panel.

    The Dense Component axis is exactly ``reference_panel`` (hg38 ALIDs). Each
    study is read **once**: on-panel associations fill the nested Dense Component,
    off-panel associations go to the Ragged Overflow. VCFs are assumed hg19 and
    lifted to hg38 inline (like the dense builder).
    """
    manifest_rows = _read_manifest(manifest_path)
    if not manifest_rows:
        raise ValueError(f"manifest {manifest_path} contains no rows")

    out = Path(output_path)
    if out.exists():
        if not overwrite:
            raise FileExistsError(f"output path already exists: {out}")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    dense_dir = dense_component_path(out)
    dense_dir.mkdir()

    panel_alids = read_reference_panel_alids(reference_panel)
    if not panel_alids:
        raise ValueError(f"reference panel {reference_panel} contained no ALIDs")
    log.info("Reference panel: %d variants", len(panel_alids))

    # ── Pass 1: union of source variants + liftover ──────────────────────────
    log.info("Pass 1: collecting source variants from %d VCFs (serial)", len(manifest_rows))
    hg19_tuples: set[tuple[str, int, str, str]] = set()
    t0 = time.monotonic()
    for i, row in enumerate(manifest_rows):
        reader = resolve_reader(
            row.source_reader_capability, row.file_path, StoredEffectScale(row.stored_effect_scale)
        )
        hg19_tuples.update(reader.stream_variants())
        _log_progress("Pass 1", i + 1, len(manifest_rows), t0,
                      f"{len(hg19_tuples)} unique variants", every=250)
    log.info("Running liftover hg19 → hg38 (%d variants)", len(hg19_tuples))
    hg19_lookup = build_liftover_lookup(
        hg19_tuples, from_build="hg19", to_build="hg38",
        failure_threshold=liftover_failure_threshold, chain_file=chain_file,
    )
    log.info("Liftover complete: %d variants mapped", len(hg19_lookup))

    # Partition observed hg38 ALIDs into on-panel and off-panel.
    observed_alids = set(hg19_lookup.values())
    off_panel_alids = sorted(observed_alids - panel_alids, key=_alid_sort_key)
    panel_sorted = sorted(panel_alids, key=_alid_sort_key)
    shared_sorted = sorted(panel_alids | set(off_panel_alids), key=_alid_sort_key)

    dense_row = {alid: i for i, alid in enumerate(panel_sorted)}
    shared_index = {alid: i for i, alid in enumerate(shared_sorted)}
    n_panel = len(panel_sorted)
    n_off_panel = len(off_panel_alids)
    n_shared = len(shared_sorted)
    n_analyses = len(manifest_rows)
    analysis_index = {row.trait_id: i for i, row in enumerate(manifest_rows)}
    log.info(
        "Partition: %d panel (dense), %d off-panel (overflow), %d shared variants, %d analyses",
        n_panel, n_off_panel, n_shared, n_analyses,
    )

    # Provenance: source (hg19) ALID each hg38 row was lifted from (None on collision).
    hg38_to_source: dict[str, str | None] = {}
    for (chrom, pos, ref, alt), hg38_alid in hg19_lookup.items():
        a1, a2 = sorted((ref, alt))
        origin = f"{chrom}:{pos}:{a1}:{a2}"
        if hg38_alid in hg38_to_source:
            if hg38_to_source[hg38_alid] != origin:
                hg38_to_source[hg38_alid] = None
        else:
            hg38_to_source[hg38_alid] = origin

    # stored_effect_scale comes from the manifest, not the VCF header (issue
    # #17 -- the ieu-a-7 fix: the source header is not authoritative for
    # effect scale).
    analyses = [
        AnalysisMetadata(
            analysis_id=row.trait_id,
            phenotype_id=row.trait_id,
            phenotype_label=row.trait_name,
            analysis_label=row.trait_id,
            stored_effect_scale=row.stored_effect_scale,
        )
        for row in manifest_rows
    ]

    # ── Write the Dense Component skeleton (a valid dense store) ──────────────
    _write_index(dense_dir, panel_sorted, analyses, chunk_shape, dtype)
    _write_variant_table(dense_dir, panel_sorted, hg38_to_source)
    effective_chunks = _create_dense_zarr(dense_dir, n_panel, n_analyses, chunk_shape, dtype)

    # dense row -> shared variant_index (ascending — panel keeps genomic order).
    dense_to_shared = np.array(
        [shared_index[alid] for alid in panel_sorted], dtype=np.int32
    )
    np.save(dense_to_shared_path(out), dense_to_shared)

    # ── Pass 2: route each study once (dense spill + overflow spill) ──────────
    log.info("Pass 2: routing %d analyses (n_workers=%d)", n_analyses, n_workers)
    keys_sorted, targets_sorted, ispanel_sorted = _build_routing_index(
        hg19_lookup, dense_row, shared_index
    )
    del hg19_lookup
    t2 = time.monotonic()
    spill_dir = Path(tempfile.mkdtemp(prefix=f".{out.name}.hybridspill.", dir=out.parent))
    try:
        if n_workers <= 1:
            for i, row in enumerate(manifest_rows):
                col = analysis_index[row.trait_id]
                dense, overflow = _resolve_column_hybrid(
                    row.file_path,
                    keys_sorted,
                    targets_sorted,
                    ispanel_sorted,
                    row.se_divisor,
                    capability=row.source_reader_capability,
                    stored_effect_scale=row.stored_effect_scale,
                )
                _spill_hybrid_column(spill_dir, col, dense, overflow)
                _log_progress("Pass 2", i + 1, n_analyses, t2, f"last: {row.trait_id}", every=25)
        else:
            global _pass2_keys_sorted, _pass2_targets_sorted, _pass2_ispanel_sorted
            global _pass2_spill_dir
            _pass2_keys_sorted = keys_sorted
            _pass2_targets_sorted = targets_sorted
            _pass2_ispanel_sorted = ispanel_sorted
            _pass2_spill_dir = spill_dir
            id_by_col = {analysis_index[row.trait_id]: row.trait_id for row in manifest_rows}
            try:
                with _fork_pool(n_workers) as pool:
                    tasks = [
                        (
                            analysis_index[row.trait_id],
                            row.file_path,
                            row.se_divisor,
                            row.source_reader_capability,
                            row.stored_effect_scale,
                        )
                        for row in manifest_rows
                    ]
                    futures = [pool.submit(_pass2_worker, t) for t in tasks]
                    for i, fut in enumerate(as_completed(futures)):
                        col = fut.result()
                        _log_progress("Pass 2", i + 1, n_analyses, t2,
                                      f"last: {id_by_col[col]}", every=25)
            finally:
                _pass2_keys_sorted = None
                _pass2_targets_sorted = None
                _pass2_ispanel_sorted = None
                _pass2_spill_dir = None

        # Dense band-write (reuses the dense builder's band-streamer + top-hit harvest).
        all_rows, all_cols, all_z, all_se = _write_dense_bands(
            dense_dir, spill_dir, n_panel, n_analyses, effective_chunks, dtype, t2
        )
        log.info("Writing Dense Component top-hit index (%d candidates)", len(all_rows))
        write_top_hit_indexes(dense_dir, all_rows, all_cols, all_z, all_se)

        # Overflow CSR assembly (reuses RaggedCSRWriter).
        log.info("Assembling Ragged Overflow CSR from %d columns", n_analyses)
        csr = _assemble_overflow_csr(spill_dir, n_analyses)
    finally:
        shutil.rmtree(spill_dir, ignore_errors=True)

    csr.flush(out)
    n_overflow = csr.n_associations
    log.info("Building Ragged Overflow top-hit index")
    build_ragged_top_hit_indexes(out)

    # ── Shared union table + shared index + manifests ────────────────────────
    _write_index(out, shared_sorted, analyses, chunk_shape, dtype)
    _write_variant_table(out, shared_sorted, hg38_to_source)

    _write_dense_manifest(dense_dir, store_id, release_id, n_panel, n_analyses, chain_file, dtype)
    _write_hybrid_manifest(
        out, store_id, release_id, n_shared, n_analyses, n_panel, n_off_panel,
        n_overflow, chain_file, chunk_shape, dtype,
    )

    log.info(
        "Hybrid build complete: %d shared variants (%d panel + %d off-panel), "
        "%d analyses, %d overflow associations",
        n_shared, n_panel, n_off_panel, n_analyses, n_overflow,
    )
    return HybridBuildResult(
        output_path=out, n_variants=n_shared, n_analyses=n_analyses,
        n_panel=n_panel, n_off_panel=n_off_panel, n_overflow=n_overflow,
    )


def _write_variant_table(
    store_path: Path, alids: list[str], hg38_to_source: dict[str, str | None]
) -> None:
    canonical = [
        CanonicalVariant(
            chromosome=chrom, position=int(pos), effect_allele=a1, other_allele=a2
        )
        for alid in alids
        for chrom, pos, a1, a2 in [alid.split(":")]
    ]
    source_alids = [hg38_to_source.get(alid) for alid in alids]
    write_variant_axis(store_path, canonical, {}, source_alids)


def _write_dense_manifest(
    dense_dir: Path,
    store_id: str,
    release_id: str,
    n_variants: int,
    n_analyses: int,
    chain_file: str | Path | None,
    dtype: str,
) -> None:
    manifest = StoreManifest(
        store_id=f"{store_id}-dense",
        release_id=release_id,
        format_version="0.1",
        primary_layout=PrimaryStorageLayout.DENSE,
        association_coverage=AssociationCoverage.FULL,
        completion_state=CompletionState.OBSERVED_ONLY,
        reference_assembly="GRCh38",
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            "builder": "opengwasdb.v0.1_hybrid_dense_component",
            "chain_file": str(chain_file) if chain_file else "pyliftover_builtin_hg19_hg38",
            "n_variants": n_variants,
            "n_analyses": n_analyses,
            "dense": {"statistic_arrays": ["z", "se"], "dtype": dtype},
        },
    )
    (dense_dir / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_hybrid_manifest(
    out: Path,
    store_id: str,
    release_id: str,
    n_variants: int,
    n_analyses: int,
    n_panel: int,
    n_off_panel: int,
    n_overflow: int,
    chain_file: str | Path | None,
    chunk_shape: tuple[int, int],
    dtype: str,
) -> None:
    manifest = StoreManifest(
        store_id=store_id,
        release_id=release_id,
        format_version="0.1",
        primary_layout=PrimaryStorageLayout.HYBRID,
        association_coverage=AssociationCoverage.FULL,
        completion_state=CompletionState.OBSERVED_ONLY,
        reference_assembly="GRCh38",
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            "builder": "opengwasdb.v0.1_hybrid_vcf",
            "chain_file": str(chain_file) if chain_file else "pyliftover_builtin_hg19_hg38",
            "n_variants": n_variants,
            "n_analyses": n_analyses,
            "hybrid": {
                "dense_component": DENSE_SUBDIR,
                "n_panel": n_panel,
                "n_off_panel": n_off_panel,
                "n_overflow_associations": n_overflow,
                "dtype": dtype,
                "chunk_shape": list(chunk_shape),
                "compressor": DEFAULT_COMPRESSOR,
            },
        },
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
