"""Build a Ragged Observed-Only Store from filtered GWAS-SSF files.

Adapter mirroring :func:`build_ragged_from_besd`, but the source is one filtered
GWAS-SSF ``.tsv.gz`` per analysis (as produced by the opengwasdb-stores
download+filter step) plus an analyses manifest carrying per-analysis metadata
(trait id, gene, trait position, N, tissue/context, MHC flag).

The ragged storage, axes, top-hit indexes and manifest are all reused unchanged;
the only new code is the read + group + variant-index-mapping glue.
"""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from opengwasdb.layouts.ragged.top_hits import build_ragged_top_hit_indexes
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRWriter
from opengwasdb.model.enums import (
    AssociationCoverage,
    CompletionState,
    PrimaryStorageLayout,
    StoredEffectScale,
)
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.traits.axis import TraitRecord, write_traits_axis
from opengwasdb.variants.axis import (
    VARIANT_AXIS_FORMAT,
    VARIANT_TABIX_FILENAME,
    VARIANT_TABLE_FILENAME,
    write_variant_axis,
)
from opengwasdb.variants.normalise import (
    CanonicalVariant,
    VariantNormalisationError,
    chromosome_sort_key,
    orient_to_canonical,
)

_FORMAT_VERSION = "0.1"
_REFERENCE_ASSEMBLY = "GRCh38"
_MISSING = {"", ".", "NA", "NaN", "nan", "None"}


@dataclass(frozen=True)
class RaggedBuildResult:
    output_path: Path
    n_variants: int
    n_analyses: int
    n_associations: int


@dataclass
class AnalyteInput:
    analysis_index: int
    analysis_id: str
    trait_id: str
    gene_id: str | None
    gene_name: str | None
    trait_chr: str | None
    trait_bp: int | None
    n: int | None
    tissue: str | None
    context: str | None
    mhc: bool
    filtered_path: Path


def _opt(value: str | None) -> str | None:
    if value is None or value.strip() in _MISSING:
        return None
    return value.strip()


def _read_manifest(manifest_path: str | Path, filtered_dir: str | Path) -> list[AnalyteInput]:
    rows: list[AnalyteInput] = []
    with open(manifest_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bp = _opt(r.get("trait_bp"))
            n = _opt(r.get("n"))
            rows.append(AnalyteInput(
                analysis_index=int(r["analysis_index"]),
                analysis_id=r["analysis_id"],
                trait_id=r["trait_id"],
                gene_id=_opt(r.get("gene_id")),
                gene_name=_opt(r.get("gene_name")),
                trait_chr=_opt(r.get("trait_chr")),
                trait_bp=int(bp) if bp else None,
                n=int(n) if n else None,
                tissue=_opt(r.get("tissue")),
                context=_opt(r.get("context")),
                mhc=str(r.get("mhc", "")).strip().upper() in {"TRUE", "1", "YES"},
                filtered_path=Path(filtered_dir) / r["filtered_file"],
            ))
    rows.sort(key=lambda a: a.analysis_index)
    # analysis_index must be a dense 0..n-1 sequence for CSR offset alignment.
    for expected, a in enumerate(rows):
        if a.analysis_index != expected:
            raise ValueError(
                f"analysis_index must be 0..n-1 in order; got {a.analysis_index} "
                f"at position {expected}"
            )
    return rows


def _read_filtered(path: Path) -> Iterator[tuple[CanonicalVariant, float, float, str | None]]:
    """Yield (CanonicalVariant, z, se, rsid) for each usable row of a filtered file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                ori = orient_to_canonical(
                    row["chromosome"], row["base_pair_location"],
                    row["effect_allele"], row["other_allele"],
                )
            except (VariantNormalisationError, KeyError):
                continue
            try:
                se = float(row.get("standard_error", ""))
                beta = float(row.get("beta", ""))
            except (TypeError, ValueError):
                continue
            if not (se > 0) or not np.isfinite(se) or not np.isfinite(beta):
                continue
            z = beta / se
            if ori.flipped:
                z = -z
            # GWAS-SSF harmonised files carry a dedicated `rsid` column; `variant_id`
            # is the harmonised (non-rs) identifier and only usable as a fallback
            # when `rsid` is absent (opengwasdb-stores' download-filter step selects
            # both columns when present).
            rsid = _opt(row.get("rsid")) or _opt(row.get("variant_id"))
            yield ori.variant, z, se, rsid


def build_ragged_from_ssf(
    manifest_path: str | Path,
    filtered_dir: str | Path,
    output_path: str | Path,
    *,
    store_id: str,
    release_id: str,
    stored_effect_scale: str = StoredEffectScale.SD.value,
    overwrite: bool = False,
) -> RaggedBuildResult:
    """Build a Ragged Observed-Only Store from filtered GWAS-SSF inputs."""
    try:
        StoredEffectScale(stored_effect_scale)
    except ValueError as exc:
        allowed = [member.value for member in StoredEffectScale]
        raise ValueError(
            f"invalid stored_effect_scale {stored_effect_scale!r}; expected one of {allowed}"
        ) from exc

    out = Path(output_path)
    if out.exists() and not overwrite:
        raise FileExistsError(f"Store already exists: {out}. Use overwrite=True.")
    if out.exists() and overwrite:
        import shutil
        shutil.rmtree(out)
    out.mkdir(parents=True)

    analytes = _read_manifest(manifest_path, filtered_dir)
    print(f"Manifest: {len(analytes)} analyses")

    # ── Pass 1: read every filtered file; collect per-analysis (alid,z,se) and
    #            the global set of canonical variants ─────────────────────────
    per_analysis: list[list[tuple[str, float, float]]] = []
    alid_variant: dict[str, CanonicalVariant] = {}
    rsid_by_alid: dict[str, str] = {}
    for a in analytes:
        recs: list[tuple[str, float, float]] = []
        for variant, z, se, rsid in _read_filtered(a.filtered_path):
            alid = variant.alid
            recs.append((alid, z, se))
            if alid not in alid_variant:
                alid_variant[alid] = variant
                if rsid and rsid.startswith("rs"):
                    rsid_by_alid[alid] = rsid
        per_analysis.append(recs)
        print(f"  {a.gene_name or a.analysis_id}: {len(recs):,} associations")

    # ── Variant axis: sort unique variants by (chr,pos), assign variant_index ─
    variants = sorted(
        alid_variant.values(),
        key=lambda v: (chromosome_sort_key(v.chromosome), v.position),
    )
    alid_to_idx = {v.alid: i for i, v in enumerate(variants)}
    print(f"Canonical variants: {len(variants):,}")
    write_variant_axis(out, variants, rsid_by_alid)

    # ── Traits axis + SQLite index (builder schema) ──────────────────────────
    trait_records = [
        TraitRecord(
            analysis_index=a.analysis_index, analysis_id=a.analysis_id, trait_id=a.trait_id,
            n=a.n, trait_chr=a.trait_chr, trait_bp=a.trait_bp,
            gene_id=a.gene_id, gene_name=a.gene_name, tissue=a.tissue, context=a.context,
        )
        for a in analytes
    ]
    write_traits_axis(out, trait_records)
    _write_index_sqlite(out, trait_records)

    # ── Stream per-analysis associations into the CSR store ──────────────────
    csr = RaggedCSRWriter()
    for recs in per_analysis:
        if not recs:
            csr.add_analysis(
                np.empty(0, np.int32), np.empty(0, np.float16), np.empty(0, np.float16)
            )
            continue
        vi = np.fromiter(
            (alid_to_idx[alid] for alid, _, _ in recs), dtype=np.int32, count=len(recs)
        )
        z_arr = np.fromiter((zz for _, zz, _ in recs), dtype=np.float16, count=len(recs))
        se_arr = np.fromiter((ss for _, _, ss in recs), dtype=np.float16, count=len(recs))
        order = np.argsort(vi, kind="stable")
        csr.add_analysis(vi[order], z_arr[order], se_arr[order])
    csr.flush(out)

    build_ragged_top_hit_indexes(out)

    _write_manifest(
        out, store_id, release_id,
        n_variants=len(variants), n_analyses=len(analytes),
        n_associations=csr.n_associations, manifest_path=str(manifest_path),
        stored_effect_scale=stored_effect_scale,
        mhc_analyses=[a.analysis_id for a in analytes if a.mhc],
    )

    result = RaggedBuildResult(out, len(variants), len(analytes), csr.n_associations)
    print(
        f"Build complete: {result.n_variants:,} variants, "
        f"{result.n_analyses:,} analyses, {result.n_associations:,} associations"
    )
    return result


def _write_manifest(
    out: Path, store_id: str, release_id: str, *,
    n_variants: int, n_analyses: int, n_associations: int,
    manifest_path: str, stored_effect_scale: str, mhc_analyses: list[str],
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
            "builder": "opengwasdb.v0.1_ragged_observed_ssf",
            "source_manifest": manifest_path,
            "stored_effect_scale": stored_effect_scale,
            "mhc_analyses": mhc_analyses,
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
                "traits_axis": {"format": "tabix_tsv_v1"},
            },
        },
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_index_sqlite(out: Path, trait_records: list[TraitRecord]) -> None:
    # Same schema as build_besd._write_index_sqlite (kept compatible with the
    # completion pipeline and query paths).
    conn = sqlite3.connect(str(out / "index.sqlite"))
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE analyses (
            analysis_index INTEGER PRIMARY KEY,
            analysis_id    TEXT NOT NULL,
            trait_id       TEXT NOT NULL,
            gene_id        TEXT,
            gene_name      TEXT,
            tissue         TEXT,
            context        TEXT,
            trait_chr      TEXT,
            trait_bp       INTEGER,
            n              INTEGER
        )
    """)
    cur.execute("CREATE INDEX idx_analyses_trait_id  ON analyses(trait_id)")
    cur.execute("CREATE INDEX idx_analyses_gene_id   ON analyses(gene_id)")
    cur.execute("CREATE INDEX idx_analyses_trait_loc ON analyses(trait_chr, trait_bp)")
    for r in trait_records:
        cur.execute(
            "INSERT INTO analyses (analysis_index, analysis_id, trait_id, gene_id, "
            "gene_name, tissue, context, trait_chr, trait_bp, n) VALUES (?,?,?,?,?,?,?,?,?,?)",
            (r.analysis_index, r.analysis_id, r.trait_id, r.gene_id, r.gene_name,
             r.tissue, r.context, r.trait_chr, r.trait_bp, r.n),
        )
    conn.commit()
    conn.close()
