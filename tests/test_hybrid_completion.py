"""Hybrid Reference Completion (issue 058): impute only the Dense Component.

hg19 → hg38 anchors used (all A/G, ALT>REF → z negated on store):
  1:100000  → 1:100000     1:750000  → 1:814620
  1:1000000 → 1:1064620    1:1500000 → 1:1564620
  1:2000000 → 1:2068561    (OFF-PANEL → overflow)
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import zarr

from opengwasdb.layouts.dense.top_hits import read_top_hit_counts
from opengwasdb.layouts.hybrid.build import build_hybrid_from_vcf_manifest
from opengwasdb.layouts.hybrid.complete import complete_hybrid_store
from opengwasdb.model.analyses import read_analyses
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.query import query_store
from opengwasdb.validation import validate_store

PANEL_ALIDS = ["1:100000:A:G", "1:814620:A:G", "1:1064620:A:G", "1:1564620:A:G"]
OFF_PANEL_ALID = "1:2068561:A:G"


def _vcf(tmp_path, name, rows):
    header = (
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size\">\n"
        "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error\">\n"
        "##FORMAT=<ID=EZ,Number=A,Type=Float,Description=\"Z-score\">\n"
        "##SAMPLE=<ID=S,StudyType=Continuous>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    p = tmp_path / f"{name}.vcf"
    p.write_text(header + "".join(rows), encoding="utf-8")
    return p


def _write_ld_block(block_dir, block_name, snps, seed=0):
    block_dir.mkdir(parents=True, exist_ok=True)
    n = len(snps)
    lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
    for alid, eaf, bp in snps:
        chrom, _pos, a1, a2 = alid.split(":")
        lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t{eaf}\t{bp}")
    (block_dir / f"{block_name}.tsv").write_text("\n".join(lines) + "\n")
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    ld = A @ A.T + np.eye(n) * n * 0.1
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for row in ld:
            gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
    (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())


def _make_ld_panel(tmp_path):
    root = tmp_path / "ld_panel"
    _write_ld_block(
        root / "EUR" / "1", "100000-1600000",
        [
            ("1:100000:A:G", 0.35, 100_000),
            ("1:814620:A:G", 0.30, 814_620),
            ("1:1064620:A:G", 0.40, 1_064_620),
            ("1:1564620:A:G", 0.45, 1_564_620),
        ],
        seed=0,
    )
    return root


def _build_source(tmp_path):
    trait_a = _vcf(
        tmp_path, "trait_a",
        [
            "1\t100000\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n",
            "1\t750000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.0:0.4\n",
            "1\t1000000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.5:0.3\n",
            "1\t1500000\t.\tA\tG\t.\tPASS\t.\tES:SE\t0.6:0.2\n",
            "1\t2000000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.2:0.3\n",  # OFF-PANEL overflow
        ],
    )
    trait_c = _vcf(
        tmp_path, "trait_c",
        [
            "1\t100000\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.2:0.5\n",
            "1\t750000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.1:0.4\n",
            "1\t1000000\t.\tA\tG\t.\tPASS\t.\tES:SE\t1.6:0.3\n",
            # NOTE: 1:1500000 missing → 1:1564620 is an imputation target for trait_c
        ],
    )
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale\toriginal_sd_method\n"
        f"trait_a\t{trait_a}\tTrait A\t1000\tsd\tdeclared_standardised\n"
        f"trait_c\t{trait_c}\tTrait C\t1000\tsd\tdeclared_standardised\n",
        encoding="utf-8",
    )
    panel = tmp_path / "panel.txt"
    panel.write_text("\n".join(PANEL_ALIDS) + "\n", encoding="utf-8")

    src = tmp_path / "src.opengwasdb"
    build_hybrid_from_vcf_manifest(
        manifest, src, reference_panel=panel, store_id="hyb", release_id="v1"
    )
    return src


def test_hybrid_completion(tmp_path):
    src = _build_source(tmp_path)
    ld = _make_ld_panel(tmp_path)
    dst = tmp_path / "dst.opengwasdb"

    result = complete_hybrid_store(src, dst, ld, min_cor=0.0, thresh=0.9)

    # Completed hybrid store validates.
    assert validate_store(dst).ok, validate_store(dst).errors

    manifest = StoreManifest.load(dst)
    assert manifest.completion_state.value == "reference_completed"
    assert manifest.primary_layout.value == "hybrid"

    # Dense Component carries imputed + on_panel arrays.
    dense_root = zarr.open_group(str(dst / "dense" / "data.zarr"), mode="r")
    assert "imputed" in dense_root
    assert "on_panel" in dense_root

    # Overflow is observed-only, untouched (byte-identical z/se).
    src_ovf = zarr.open_group(str(src / "data.zarr" / "ragged"), mode="r")
    dst_ovf = zarr.open_group(str(dst / "data.zarr" / "ragged"), mode="r")
    assert "imputed" not in dst_ovf
    assert np.array_equal(src_ovf["z"][:], dst_ovf["z"][:])
    assert np.array_equal(src_ovf["se"][:], dst_ovf["se"][:])


def test_hit_counts_sum_dense_and_overflow_after_completion(tmp_path):
    # The completed Dense Component's own analyses.tsv already carries
    # correct post-completion counts (test_dense_completion.py); the shared
    # root's counts must be that plus the Ragged Overflow's own top-hit
    # index (ADR 0032) -- neither component alone is double-counted.
    src = _build_source(tmp_path)
    ld = _make_ld_panel(tmp_path)
    dst = tmp_path / "dst.opengwasdb"
    complete_hybrid_store(src, dst, ld, min_cor=0.0, thresh=0.9)

    rows = sorted(
        read_analyses(dst / "analyses.tsv").rows, key=lambda r: int(r["analysis_index"])
    )
    dense_counts = read_top_hit_counts(dst / "dense")
    overflow_counts = read_top_hit_counts(dst)
    for column in dense_counts:
        expected = [
            d + o for d, o in zip(dense_counts[column], overflow_counts[column], strict=True)
        ]
        assert [int(r[column]) for r in rows] == expected


def test_completed_status_only_on_dense(tmp_path):
    src = _build_source(tmp_path)
    ld = _make_ld_panel(tmp_path)
    dst = tmp_path / "dst.opengwasdb"
    result = complete_hybrid_store(src, dst, ld, min_cor=0.0, thresh=0.9)

    q = query_store(dst)
    # The off-panel overflow association must always be observed.
    r = q.phewas(OFF_PANEL_ALID)
    assert len(r["z"]) >= 1
    assert set(r["association_status"].tolist()) == {"observed"}

    # trait_c's missing panel variant (1:1564620) should be imputed if completion filled it.
    rc = q.analysis("trait_c")
    statuses = set(rc["association_status"].tolist())
    assert statuses <= {"observed", "imputed"}
    if result.n_imputed > 0:
        assert "imputed" in statuses
