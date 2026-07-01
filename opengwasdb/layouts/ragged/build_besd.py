"""Build a Ragged Observed-Only Store from BESD files."""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from opengwasdb.layouts.ragged.besd_reader import BESDReader, read_epi, read_esi
from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRWriter
from opengwasdb.model.enums import (
    AssociationCoverage,
    CompletionState,
    PrimaryStorageLayout,
)
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.traits.axis import TraitRecord, write_traits_axis
from opengwasdb.variants.axis import (
    VARIANT_AXIS_FORMAT,
    VARIANT_TABLE_FILENAME,
    VARIANT_TABIX_FILENAME,
    write_variant_axis,
)
from opengwasdb.build.liftover import build_liftover_lookup
from opengwasdb.variants.normalise import (
    CanonicalVariant,
    VariantNormalisationError,
    chromosome_sort_key,
    normalise_chromosome,
    orient_to_canonical,
)

_FORMAT_VERSION = "0.1"
_REFERENCE_ASSEMBLY = "GRCh38"


@dataclass(frozen=True)
class RaggedBuildResult:
    output_path: Path
    n_variants: int
    n_analyses: int
    n_associations: int


def build_ragged_from_besd(
    besd_prefix: str | Path,
    output_path: str | Path,
    *,
    store_id: str,
    release_id: str,
    tissue: str | None = None,
    source_build: str = "hg38",
    overwrite: bool = False,
) -> RaggedBuildResult:
    """Build a Ragged Observed-Only Store from BESD files.

    besd_prefix: path without extension (.esi, .epi, .besd are appended).
    source_build: genome assembly of the input BESD ("hg38" or "hg19").
    When source_build is "hg19", SNP coordinates are lifted over to hg38 inline.
    """
    prefix = Path(besd_prefix)
    out = Path(output_path)

    if out.exists() and not overwrite:
        raise FileExistsError(f"Store already exists: {out}. Use overwrite=True.")
    if out.exists() and overwrite:
        import shutil
        shutil.rmtree(out)
    out.mkdir(parents=True)

    # ── 1. Read ESI / EPI ────────────────────────────────────────────────────
    print(f"Reading ESI: {prefix}.esi")
    snps = read_esi(f"{prefix}.esi")
    print(f"Reading EPI: {prefix}.epi")
    probes = read_epi(f"{prefix}.epi")
    print(f"Loaded {len(snps)} SNPs and {len(probes)} probes")

    # ── 2. Optionally liftover ESI coordinates hg19 → hg38 ──────────────────
    # Maps (chr, hg19_bp, a1, a2) → lifted hg38 (chr, bp) or None if failed.
    _lifted: dict[int, tuple[str, int]] | None = None
    _normalised_source = source_build.lower().strip()
    if _normalised_source not in ("hg38", "grch38", "b38", "38"):
        print(f"Lifting over {source_build} → hg38 ...")
        lo_input = [
            (s.chromosome, s.bp, s.a1 or "A", s.a2 or "T")
            for s in snps if s.a1 and s.a2
        ]
        lo_lookup = build_liftover_lookup(
            lo_input, from_build=source_build, to_build="hg38",
        )
        # Rebuild a map: esi_row_idx → (hg38_chr, hg38_bp)
        _lifted = {}
        for s in snps:
            if not s.a1 or not s.a2:
                continue
            alid = lo_lookup.get((s.chromosome, s.bp, s.a1, s.a2))
            if alid is None:
                continue
            parts = alid.split(":")
            _lifted[s.row_idx] = (parts[0], int(parts[1]))
        n_lifted = len(_lifted)
        n_failed = len([s for s in snps if s.a1 and s.a2]) - n_lifted
        print(f"Liftover {source_build}→hg38: {n_failed}/{len(snps)} variants failed")

    # ── 3. Build canonical variant list from ESI ─────────────────────────────
    # esi_row_idx → (variant_index, flipped)
    esi_to_variant: dict[int, tuple[int, bool]] = {}
    variants: list[CanonicalVariant] = []
    rsid_by_alid: dict[str, str] = {}

    # Collect valid variants, deduplicate by ALID, sort by chr/pos
    alid_to_idx: dict[str, int] = {}
    candidate: list[tuple[tuple, int, bool, str | None]] = []  # (sort_key, esi_idx, flipped, rsid)

    for snp in snps:
        if snp.a1 is None or snp.a2 is None:
            continue
        # Use lifted coordinates when liftover was performed
        if _lifted is not None:
            pos_info = _lifted.get(snp.row_idx)
            if pos_info is None:
                continue
            chrom, bp = pos_info
        else:
            chrom, bp = snp.chromosome, snp.bp
        try:
            ori = orient_to_canonical(chrom, bp, snp.a1, snp.a2)
        except VariantNormalisationError:
            continue
        sort_key = (chromosome_sort_key(ori.variant.chromosome), ori.variant.position)
        rsid = snp.snp_id if snp.snp_id.startswith("rs") else None
        candidate.append((sort_key, snp.row_idx, ori.variant, ori.flipped, rsid))

    # Sort by genomic position so variant_index is position-ordered
    candidate.sort(key=lambda x: x[0])

    for sort_key, esi_row_idx, variant, flipped, rsid in candidate:
        alid = variant.alid
        if alid in alid_to_idx:
            # Duplicate ALID — map to the same variant_index
            existing_idx = alid_to_idx[alid]
            esi_to_variant[esi_row_idx] = (existing_idx, flipped)
        else:
            variant_index = len(variants)
            alid_to_idx[alid] = variant_index
            variants.append(variant)
            esi_to_variant[esi_row_idx] = (variant_index, flipped)
            if rsid:
                rsid_by_alid[alid] = rsid

    print(f"Canonical variants: {len(variants)} (from {len(snps)} ESI entries)")

    # ── 3. Write variant axis ────────────────────────────────────────────────
    print("Writing variants.tsv.gz ...")
    write_variant_axis(out, variants, rsid_by_alid)

    # ── 4. Build and write traits.tsv.gz ─────────────────────────────────────
    trait_records: list[TraitRecord] = []
    for probe in probes:
        try:
            probe_chr = normalise_chromosome(probe.chromosome)
        except VariantNormalisationError:
            probe_chr = None

        analysis_id = probe.probe_id
        if tissue:
            analysis_id = f"{probe.probe_id}::{tissue}"

        trait_records.append(TraitRecord(
            analysis_index=probe.row_idx,
            analysis_id=analysis_id,
            probe_id=probe.probe_id,
            n=None,
            probe_chr=probe_chr,
            probe_bp=probe.probe_bp if probe.probe_bp > 0 else None,
            gene_id=probe.probe_id if probe.probe_id.startswith("ENSG") else None,
            gene_name=probe.gene,
            tissue=tissue,
            context=None,
        ))

    print("Writing traits.tsv.gz ...")
    write_traits_axis(out, trait_records)

    # ── 5. Write index.sqlite ────────────────────────────────────────────────
    print("Writing index.sqlite ...")
    _write_index_sqlite(out, trait_records)

    # ── 6. Stream BESD associations into zarr CSR ────────────────────────────
    print(f"Reading BESD: {prefix}.besd")
    besd = BESDReader(f"{prefix}.besd", len(probes))
    print(f"BESD format: SPARSE_FILE_TYPE_{besd.format_type}")

    csr = RaggedCSRWriter()
    skipped_probes = 0

    for probe in probes:
        raw_snp_idx, betas, ses = besd.get_probe_associations(probe.row_idx)

        if len(raw_snp_idx) == 0:
            csr.add_analysis(
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float16),
                np.empty(0, dtype=np.float16),
            )
            continue

        # Map ESI row indices → variant indices, apply orientation flips
        vi_list: list[int] = []
        z_list: list[float] = []
        se_list: list[float] = []

        for esi_idx, beta, se in zip(
            raw_snp_idx.tolist(), betas.tolist(), ses.tolist()
        ):
            mapping = esi_to_variant.get(int(esi_idx))
            if mapping is None:
                continue
            variant_index, flipped = mapping
            if se <= 0 or not np.isfinite(beta) or not np.isfinite(se):
                continue
            z = beta / se
            if flipped:
                z = -z
            vi_list.append(variant_index)
            z_list.append(z)
            se_list.append(se)

        if vi_list:
            # Sort by variant_index for consistent ordering within each analysis
            order = np.argsort(vi_list)
            csr.add_analysis(
                np.array(vi_list, dtype=np.int32)[order],
                np.array(z_list, dtype=np.float16)[order],
                np.array(se_list, dtype=np.float16)[order],
            )
        else:
            skipped_probes += 1
            csr.add_analysis(
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float16),
                np.empty(0, dtype=np.float16),
            )

        if (probe.row_idx + 1) % 1000 == 0:
            print(f"  Processed {probe.row_idx + 1} / {len(probes)} probes")

    print(f"Flushing zarr CSR ({csr.n_associations:,} associations) ...")
    csr.flush(out)

    if skipped_probes:
        print(f"  {skipped_probes} probes had no valid associations after filtering")

    # ── 7. Build top-hit indexes ─────────────────────────────────────────────
    print("Building top-hit indexes ...")
    build_ragged_top_hit_indexes(out)

    # ── 8. Write manifest.json ───────────────────────────────────────────────
    _write_manifest(
        out, store_id, release_id,
        n_variants=len(variants),
        n_analyses=len(probes),
        n_associations=csr.n_associations,
        besd_prefix=str(prefix),
        source_build=source_build,
    )

    result = RaggedBuildResult(
        output_path=out,
        n_variants=len(variants),
        n_analyses=len(probes),
        n_associations=csr.n_associations,
    )
    print(
        f"Build complete: {result.n_variants:,} variants, "
        f"{result.n_analyses:,} analyses, "
        f"{result.n_associations:,} associations"
    )
    return result


def _write_manifest(
    out: Path,
    store_id: str,
    release_id: str,
    *,
    n_variants: int,
    n_analyses: int,
    n_associations: int,
    besd_prefix: str,
    source_build: str = "hg38",
) -> None:
    manifest = StoreManifest(
        store_id=store_id,
        release_id=release_id,
        format_version=_FORMAT_VERSION,
        primary_layout=PrimaryStorageLayout.RAGGED,
        association_coverage=AssociationCoverage.CIS_AND_SIGNALS,
        completion_state=CompletionState.OBSERVED_ONLY,
        reference_assembly=_REFERENCE_ASSEMBLY,
        created_at=datetime.now(UTC).isoformat(),
        provenance={
            "builder": "opengwasdb.v0.1_ragged_observed_besd",
            "source_besd_prefix": besd_prefix,
            "source_build": source_build,
            "n_variants": n_variants,
            "n_analyses": n_analyses,
            "n_associations": n_associations,
            "ragged": {
                "statistic_arrays": ["z", "se"],
                "dtype": "float16",
                "variant_axis": {
                    "format": VARIANT_AXIS_FORMAT,
                    "table": VARIANT_TABLE_FILENAME,
                    "tabix_index": VARIANT_TABIX_FILENAME,
                },
                "traits_axis": {
                    "format": "tabix_tsv_v1",
                },
            },
        },
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_index_sqlite(out: Path, trait_records: list[TraitRecord]) -> None:
    db_path = out / "index.sqlite"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE analyses (
            analysis_index INTEGER PRIMARY KEY,
            analysis_id    TEXT NOT NULL,
            probe_id       TEXT NOT NULL,
            gene_id        TEXT,
            gene_name      TEXT,
            tissue         TEXT,
            context        TEXT,
            probe_chr      TEXT,
            probe_bp       INTEGER,
            n              INTEGER
        )
    """)
    cursor.execute("CREATE INDEX idx_analyses_probe_id ON analyses(probe_id)")
    cursor.execute("CREATE INDEX idx_analyses_gene_id  ON analyses(gene_id)")
    cursor.execute("CREATE INDEX idx_analyses_probe_loc ON analyses(probe_chr, probe_bp)")

    for rec in trait_records:
        cursor.execute("""
            INSERT INTO analyses
              (analysis_index, analysis_id, probe_id, gene_id, gene_name,
               tissue, context, probe_chr, probe_bp, n)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            rec.analysis_index, rec.analysis_id, rec.probe_id,
            rec.gene_id, rec.gene_name, rec.tissue, rec.context,
            rec.probe_chr, rec.probe_bp, rec.n,
        ))

    conn.commit()
    conn.close()
