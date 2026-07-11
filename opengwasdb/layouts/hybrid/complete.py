"""Hybrid Reference Completion — impute only the Dense Component (ADR 0026, issue 058).

The Dense Component's axis is the LD reference panel, so completing it is exactly
the dense completion problem: this reuses ``complete_dense_store`` **unchanged** on
``<store>/dense``. The Ragged Overflow Component is off-panel (no LD structure) and
is left observed-only — its associations are copied through untouched.

Because dense completion may extend the panel axis (if the LD panel is a superset
of the build panel), the shared union table and the ``dense_to_shared`` map are
rebuilt from the *completed* Dense Component's axis, and the overflow
``variant_index`` values are remapped onto the rebuilt shared index space. When the
build panel already equals the LD panel (the intended case) nothing is added, so
the overflow indices are unchanged.
"""

from __future__ import annotations

import json
import logging
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from opengwasdb.index import connect
from opengwasdb.layouts.dense.build import AnalysisMetadata
from opengwasdb.layouts.dense.build_vcf import _alid_sort_key, _write_index
from opengwasdb.layouts.dense.complete import complete_dense_store
from opengwasdb.layouts.dense.constants import DEFAULT_COMPRESSOR, DEFAULT_DTYPE
from opengwasdb.layouts.hybrid.build import _write_variant_table
from opengwasdb.layouts.hybrid.layout import (
    DENSE_SUBDIR,
    dense_component_path,
    dense_to_shared_path,
)
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader, RaggedCSRWriter
from opengwasdb.model.enums import (
    AssociationCoverage,
    CompletionState,
    PrimaryStorageLayout,
)
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.variants import VariantAxis

log = logging.getLogger(__name__)

__all__ = ["complete_hybrid_store", "HybridCompletionResult"]


@dataclass(frozen=True)
class HybridCompletionResult:
    output_path: Path
    n_variants: int
    n_analyses: int
    n_panel: int
    n_off_panel: int
    n_overflow: int
    n_imputed: int


def complete_hybrid_store(
    source_path: str | Path,
    dest_path: str | Path,
    ld_dir: str | Path,
    *,
    ancestry: str = "EUR",
    min_cor: float = 0.7,
    thresh: float = 0.9,
    release_id: str | None = None,
    n_workers: int = 1,
    overwrite: bool = False,
) -> HybridCompletionResult:
    """Produce a Reference-Completed Hybrid store from an Observed-Only Hybrid source."""
    src = Path(source_path)
    dst = Path(dest_path)

    src_manifest = StoreManifest.load(src)
    if src_manifest.primary_layout is not PrimaryStorageLayout.HYBRID:
        raise ValueError(
            f"source store is not Hybrid (primary_layout={src_manifest.primary_layout})"
        )
    if src_manifest.completion_state is not CompletionState.OBSERVED_ONLY:
        raise ValueError(
            "source hybrid store is not Observed-Only "
            f"(completion_state={src_manifest.completion_state})"
        )

    if dst.exists():
        if not overwrite:
            raise FileExistsError(f"Destination already exists: {dst}. Use overwrite=True.")
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    # Per-Analysis ancestry-match filter (ADR 0028): if the store carries an
    # ancestry sidecar, only impute Analyses whose Assigned Ancestry matches the
    # applied panel; the rest are carried through observed-only. Absent sidecar =
    # impute everything (no behaviour change).
    from opengwasdb.ancestry.store import (
        read_ancestry_provenance,
        read_ancestry_sidecar,
        write_ancestry_provenance,
        write_ancestry_sidecar,
    )

    src_ancestry = read_ancestry_sidecar(src)
    impute_ids: set[str] | None = None
    if src_ancestry:
        matched = {r["trait_id"] for r in src_ancestry if r["assigned_ancestry"] == ancestry}
        impute_ids = matched
        log.info(
            "Ancestry-matched completion: %d/%d analyses match panel ancestry %s",
            len(matched), len(src_ancestry), ancestry,
        )

    try:
        # ── 1. Complete the Dense Component (dense pipeline, unchanged) ────────
        log.info("Completing Dense Component via the dense reference-completion pipeline")
        dense_result = complete_dense_store(
            dense_component_path(src),
            dense_component_path(dst),
            ld_dir,
            ancestry=ancestry,
            min_cor=min_cor,
            thresh=thresh,
            release_id=release_id,
            n_workers=n_workers,
            overwrite=True,
            impute_analysis_ids=impute_ids,
        )

        # ── 2. Rebuild the shared union table from the completed dense axis ────
        dense_axis = VariantAxis(dense_component_path(dst))
        dense_records = dense_axis.all()
        dense_axis.close()
        dense_alids = [r.alid for r in dense_records]  # dense row order

        src_csr = RaggedCSRReader(src)
        offsets = src_csr._offsets[:]
        n_analyses = len(offsets) - 1
        src_vi = src_csr._variant_index[:]
        src_z = src_csr._z[:]
        src_se = src_csr._se[:]

        src_shared_axis = VariantAxis(src)
        vi_to_record = src_shared_axis.by_indices(np.unique(src_vi).tolist())
        src_shared_all = src_shared_axis.all()
        src_shared_axis.close()
        source_alid_by_alid = {r.alid: r.source_alid for r in src_shared_all}

        overflow_alids = np.array(
            [vi_to_record[int(v)].alid for v in src_vi], dtype=object
        ) if len(src_vi) else np.empty(0, dtype=object)

        union = sorted(set(dense_alids) | set(overflow_alids.tolist()), key=_alid_sort_key)
        new_shared_index = {alid: i for i, alid in enumerate(union)}
        n_shared = len(union)
        n_panel = len(dense_alids)
        n_off_panel = n_shared - n_panel

        # ── 3. Write shared table + index at the destination ──────────────────
        analyses = _read_analyses(dense_component_path(dst))
        _write_index(dst, union, analyses, _chunk_shape(src_manifest), DEFAULT_DTYPE)
        source_by_alid = {a: source_alid_by_alid.get(a) for a in union}
        _write_variant_table(dst, union, source_by_alid)

        # ── 4. dense_to_shared map for the completed axis ─────────────────────
        dense_to_shared = np.array([new_shared_index[a] for a in dense_alids], dtype=np.int32)
        np.save(dense_to_shared_path(dst), dense_to_shared)

        # ── 5. Rebuild the overflow CSR with remapped shared indices ──────────
        csr = RaggedCSRWriter()
        for ai in range(n_analyses):
            s, e = int(offsets[ai]), int(offsets[ai + 1])
            if s == e:
                csr.add_analysis(
                    np.empty(0, dtype=np.int32),
                    np.empty(0, dtype=np.float16),
                    np.empty(0, dtype=np.float16),
                )
                continue
            new_vi = np.array(
                [new_shared_index[vi_to_record[int(v)].alid] for v in src_vi[s:e]],
                dtype=np.int32,
            )
            z = src_z[s:e].astype(np.float16)
            se = src_se[s:e].astype(np.float16)
            order = np.argsort(new_vi, kind="stable")
            csr.add_analysis(new_vi[order], z[order], se[order])
        csr.flush(dst)
        build_ragged_top_hit_indexes(dst)

        # ── 6. Hybrid manifest (reference-completed) ──────────────────────────
        new_release = release_id or f"{src_manifest.release_id}-completed"
        _write_completed_manifest(
            dst, src_manifest, new_release, n_shared, n_analyses, n_panel, n_off_panel,
            csr.n_associations, dense_result.n_imputed,
        )

        # ── 7. Carry the ancestry sidecar forward, recording completed_against ─
        if src_ancestry:
            analyses_ancestry = [(r["trait_id"], r["assigned_ancestry"]) for r in src_ancestry]
            completed_against = {
                r["trait_id"]: (ancestry if r["assigned_ancestry"] == ancestry else "")
                for r in src_ancestry
            }
            write_ancestry_sidecar(dst, analyses_ancestry, completed_against=completed_against)
            prov = read_ancestry_provenance(src)
            if prov:
                prov_n = prov.get("n_analyses")
                write_ancestry_provenance(
                    dst,
                    catalogue_version=str(prov.get("catalogue_version", "")),
                    subset_filter=str(prov.get("subset_filter", "")),
                    ancestry_reference_version=str(prov.get("ancestry_reference_version", "")),
                    n_analyses=prov_n if isinstance(prov_n, int) else len(src_ancestry),
                )

        log.info(
            "Hybrid completion complete: %d shared variants (%d panel + %d off-panel), "
            "%d imputed dense cells, %d overflow associations (observed-only)",
            n_shared, n_panel, n_off_panel, dense_result.n_imputed, csr.n_associations,
        )
        return HybridCompletionResult(
            output_path=dst, n_variants=n_shared, n_analyses=n_analyses,
            n_panel=n_panel, n_off_panel=n_off_panel, n_overflow=csr.n_associations,
            n_imputed=dense_result.n_imputed,
        )
    except Exception:
        shutil.rmtree(dst, ignore_errors=True)
        raise


def _chunk_shape(manifest: StoreManifest) -> tuple[int, int]:
    cs = manifest.provenance.get("hybrid", {}).get("chunk_shape")
    if cs and len(cs) == 2:
        return (int(cs[0]), int(cs[1]))
    from opengwasdb.layouts.dense.constants import DEFAULT_CHUNK_SHAPE

    return DEFAULT_CHUNK_SHAPE


def _read_analyses(dense_dir: Path) -> list[AnalysisMetadata]:
    with connect(dense_dir / "index.sqlite") as conn:
        rows = conn.execute(
            "SELECT analysis_index, analysis_id, phenotype_id, phenotype_label, "
            "analysis_label, stored_effect_scale FROM analyses ORDER BY analysis_index"
        ).fetchall()
    return [
        AnalysisMetadata(
            analysis_id=r["analysis_id"],
            phenotype_id=r["phenotype_id"],
            phenotype_label=r["phenotype_label"],
            analysis_label=r["analysis_label"],
            stored_effect_scale=r["stored_effect_scale"],
        )
        for r in rows
    ]


def _write_completed_manifest(
    dst: Path,
    src_manifest: StoreManifest,
    release_id: str,
    n_variants: int,
    n_analyses: int,
    n_panel: int,
    n_off_panel: int,
    n_overflow: int,
    n_imputed: int,
) -> None:
    manifest = StoreManifest(
        store_id=src_manifest.store_id,
        release_id=release_id,
        format_version=src_manifest.format_version,
        primary_layout=PrimaryStorageLayout.HYBRID,
        association_coverage=AssociationCoverage.FULL,
        completion_state=CompletionState.REFERENCE_COMPLETED,
        reference_assembly=src_manifest.reference_assembly,
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            **src_manifest.provenance,
            "source_release_id": src_manifest.release_id,
            "hybrid": {
                **src_manifest.provenance.get("hybrid", {}),
                "dense_component": DENSE_SUBDIR,
                "n_panel": n_panel,
                "n_off_panel": n_off_panel,
                "n_overflow_associations": n_overflow,
                "compressor": DEFAULT_COMPRESSOR,
            },
            "n_variants": n_variants,
            "n_analyses": n_analyses,
            "completion": {
                "component_completed": "dense",
                "overflow_completion_state": "observed_only",
                "n_imputed_dense": n_imputed,
            },
        },
    )
    (dst / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
