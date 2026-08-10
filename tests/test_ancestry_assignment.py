"""Tracer tests for ancestry assignment (issue 062).

Exercises the whole chain on tiny synthetic data: reference loader → AF extraction
(oriented, palindromic-filtered) → NNLS mixture → multi-gate rule → Catalogue TSV.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from opengwasdb.ancestry import (
    AncestryReference,
    Gates,
    assign_ancestry,
    assign_from_vcf,
    catalogue_fieldnames,
    extract_af_at_sites,
    is_palindromic,
    load_reference,
    write_catalogue,
)
from opengwasdb.ancestry.catalogue import BUILD_COLUMNS, CatalogueRow
from opengwasdb.ancestry.extract import write_regions_file
from opengwasdb.layouts.dense.build_vcf import _read_manifest

# Fine groups (two EUR subgroups + one AFR + one EAS) and their super-populations.
GROUPS = ["United Kingdom", "Finland", "Africa (West)", "Asia (East)"]
GROUP_TO_SUPERPOP = {
    "United Kingdom": "EUR",
    "Finland": "EUR",
    "Africa (West)": "AFR",
    "Asia (East)": "EAS",
}
N_VARIANTS = 40


def _reference_freqs() -> np.ndarray:
    """Distinct, well-separated frequency profiles per group (identifiable fit)."""
    rng = np.random.default_rng(20260710)
    return rng.uniform(0.05, 0.95, size=(N_VARIANTS, len(GROUPS)))


def _write_reference(tmp_path: Path, freqs: np.ndarray) -> tuple[Path, Path, list[str]]:
    """Write a tiny reference panel + group map; return (freqs, groups, alids)."""
    alids = [f"1:{1000 + i}:A:C" for i in range(N_VARIANTS)]
    freqs_path = tmp_path / "ref_freqs.tsv"
    header = ["alid", "chromosome", "position", "effect_allele", "other_allele", "rsid", *GROUPS]
    lines = ["\t".join(header)]
    for i, alid in enumerate(alids):
        chrom, pos, a1, a2 = alid.split(":")
        cells = [alid, chrom, pos, a1, a2, f"rs{i}", *[f"{f:.6g}" for f in freqs[i]]]
        lines.append("\t".join(cells))
    freqs_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    groups_path = tmp_path / "ancestry_groups.tsv"
    glines = ["group\tsuper_pop"] + [f"{g}\t{GROUP_TO_SUPERPOP[g]}" for g in GROUPS]
    groups_path.write_text("\n".join(glines) + "\n", encoding="utf-8")
    return freqs_path, groups_path, alids


@pytest.fixture
def reference(tmp_path) -> AncestryReference:
    freqs = _reference_freqs()
    freqs_path, groups_path, _alids = _write_reference(tmp_path, freqs)
    # maf_floor=0: keep every synthetic variant regardless of MAF.
    return load_reference(freqs_path, groups_path, maf_floor=0.0)


# --- reference loader ------------------------------------------------------


def test_reference_loader_shape_and_superpops(reference):
    assert reference.n_variants == N_VARIANTS
    assert reference.groups == GROUPS
    assert reference.superpops == ["AFR", "EAS", "EUR"]
    # aggregate() sums fine proportions into super-pops in superpops order.
    fine = np.array([0.3, 0.2, 0.4, 0.1])  # UK, Finland, AfrWest, AsiaEast
    agg = reference.aggregate(fine)
    assert agg[reference.superpops.index("EUR")] == pytest.approx(0.5)
    assert agg[reference.superpops.index("AFR")] == pytest.approx(0.4)
    assert agg[reference.superpops.index("EAS")] == pytest.approx(0.1)


def test_reference_maf_floor_drops_monomorphic(tmp_path):
    freqs = _reference_freqs()
    freqs[0, :] = 0.001  # near-monomorphic in every group
    freqs_path, groups_path, _ = _write_reference(tmp_path, freqs)
    ref = load_reference(freqs_path, groups_path, maf_floor=0.01)
    assert ref.n_variants == N_VARIANTS - 1


# --- AF extraction: orientation & palindromic exclusion --------------------


def _af_vcf(tmp_path: Path, rows: list[str]) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
        "##SAMPLE=<ID=S1,StudyType=Continuous>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    path = tmp_path / "study.vcf"
    path.write_text(header + "".join(rows), encoding="utf-8")
    return path


def test_extract_orients_and_excludes_palindromic(tmp_path):
    # Variant 1: REF=C ALT=A → A1=A=ALT (not flipped) → af stays 0.30.
    # Variant 2: REF=A ALT=G → A1=A=REF (ALT flipped) → af = 1 - 0.30 = 0.70.
    # Variant 3: REF=A ALT=T → palindromic → excluded.
    vcf = _af_vcf(
        tmp_path,
        [
            "1\t1000\t.\tC\tA\t.\tPASS\t.\tAF\t0.30\n",
            "1\t1001\t.\tA\tG\t.\tPASS\t.\tAF\t0.30\n",
            "1\t1002\t.\tA\tT\t.\tPASS\t.\tAF\t0.30\n",
        ],
    )
    wanted = {"1:1000:A:C", "1:1001:A:G", "1:1002:A:T"}
    afs = extract_af_at_sites(vcf, wanted)
    assert afs["1:1000:A:C"] == pytest.approx(0.30)
    assert afs["1:1001:A:G"] == pytest.approx(0.70)
    assert "1:1002:A:T" not in afs  # palindromic dropped


def test_is_palindromic():
    assert is_palindromic("A", "T") and is_palindromic("C", "G")
    assert not is_palindromic("A", "C") and not is_palindromic("A", "G")


class _FakeLiftover:
    """Stand-in for pyliftover.LiftOver: shifts chr1 positions by +64620."""

    def convert_coordinate(self, chrom, pos0):  # noqa: D401
        if chrom != "chr1":
            return []
        return [("chr1", pos0 + 64620, "+", 0)]


def test_extract_applies_liftover(tmp_path):
    # Study is on an older build; the reference ALIDs are on the lifted build.
    vcf = _af_vcf(
        tmp_path,
        ["1\t1000\t.\tC\tA\t.\tPASS\t.\tAF\t0.25\n"],  # pos 1000 → 65620 after lift
    )
    wanted = {"1:65620:A:C"}  # canonical A1=A (ALT), lifted position
    afs = extract_af_at_sites(vcf, wanted, liftover=_FakeLiftover())
    assert set(afs) == {"1:65620:A:C"}
    assert afs["1:65620:A:C"] == pytest.approx(0.25)  # unflipped (A is canonical A1)


# --- mixture + gates -------------------------------------------------------


def _study_af(reference, weights: dict[str, float]) -> dict[str, float]:
    """Synthesize per-variant AF as a group mixture of reference frequencies."""
    b = np.zeros(reference.n_variants)
    for group, w in weights.items():
        b += w * reference.freqs[:, reference.groups.index(group)]
    return dict(zip(reference.alids.tolist(), b.tolist(), strict=True))


def _gates() -> Gates:
    return Gates(tau=0.90, delta=0.20, n_min=10, residual_max=0.05)


def test_clean_single_ancestry_is_assigned(reference):
    study = _study_af(reference, {"United Kingdom": 1.0})
    result = assign_ancestry(study, reference, _gates())
    assert result.assigned_ancestry == "EUR"
    assert result.gate_reason == "ok"
    assert result.dominant_proportion > 0.9
    assert result.residual < 0.01


def test_admixed_is_unassigned(reference):
    study = _study_af(reference, {"United Kingdom": 0.5, "Africa (West)": 0.5})
    result = assign_ancestry(study, reference, _gates())
    assert result.assigned_ancestry is None
    assert result.gate_reason in {"proportion", "margin"}


def test_low_overlap_is_unassigned(reference):
    study = _study_af(reference, {"United Kingdom": 1.0})
    few = dict(list(study.items())[:5])  # 5 < n_min=10
    result = assign_ancestry(few, reference, _gates())
    assert result.assigned_ancestry is None
    assert result.gate_reason == "overlap"
    assert result.af_overlap == 5


def test_corrupt_af_is_unassigned_by_residual(reference):
    # Mis-orient half the sites: no single mixture fits → large residual.
    study = _study_af(reference, {"United Kingdom": 1.0})
    for i, alid in enumerate(list(study)):
        if i % 2 == 0:
            study[alid] = 1.0 - study[alid]
    result = assign_ancestry(study, reference, _gates())
    assert result.assigned_ancestry is None
    assert result.gate_reason == "residual"


# --- end-to-end via VCF ----------------------------------------------------


def _bgzip_index(vcf: Path) -> Path:
    """Compress + tabix-index a plain VCF via bcftools (what -R index-jumps need)."""
    import subprocess

    out = vcf.with_suffix(".vcf.gz")
    subprocess.run(
        ["bcftools", "view", str(vcf), "-Oz", "-o", str(out), "--write-index=tbi"],
        check=True,
        capture_output=True,
    )
    return out


def test_assign_from_vcf_end_to_end(tmp_path, reference):
    study = _study_af(reference, {"United Kingdom": 1.0})
    rows = []
    for alid, af in study.items():
        _chrom, pos, _a1, _a2 = alid.split(":")
        # REF=C ALT=A → canonical A1=A=ALT, AF unflipped.
        rows.append(f"1\t{pos}\t.\tC\tA\t.\tPASS\t.\tAF\t{af:.6g}\n")
    vcf = _bgzip_index(_af_vcf(tmp_path, rows))
    regions = write_regions_file(reference.index.keys(), tmp_path / "regions.txt")
    result = assign_from_vcf(vcf, reference, _gates(), regions_file=regions)
    assert result.assigned_ancestry == "EUR"
    assert result.af_overlap == N_VARIANTS


# --- catalogue writer: superset of the build manifest ----------------------


def test_catalogue_is_manifest_superset(tmp_path, reference):
    study = _study_af(reference, {"United Kingdom": 1.0})
    assigned = assign_ancestry(study, reference, _gates())
    admixed = assign_ancestry(
        _study_af(reference, {"United Kingdom": 0.5, "Asia (East)": 0.5}), reference, _gates()
    )
    rows = [
        CatalogueRow("t1", "/data/t1.vcf.gz", "Trait One", 1000, "European", assigned),
        CatalogueRow("t2", "/data/t2.vcf.gz", "Trait Two", 2000, "Mixed", admixed),
    ]
    path = write_catalogue(
        tmp_path / "catalogue.tsv",
        rows,
        reference.superpops,
        catalogue_version="cat-v1",
        ancestry_reference_version="prive2022-hg38",
    )

    # Header carries build columns first, then annotations + version stamps.
    header = path.read_text().splitlines()[0].split("\t")
    assert header[:4] == BUILD_COLUMNS
    assert "assigned_ancestry" in header and "catalogue_version" in header
    assert catalogue_fieldnames(reference.superpops) == header

    # The Catalogue carries BUILD_COLUMNS (trait_id/file_path/trait_name/n)
    # plus its ancestry annotations, but is not on its own a complete build
    # manifest as of issue #17: stored_effect_scale is a genuinely separate
    # build input ancestry assignment never needs (it may run before a
    # study's effect scale is even resolved), so _read_manifest correctly
    # rejects a Catalogue file until something adds that column
    # (opengwasdb.ancestry.subset does this when bridging into an actual
    # build).
    import csv

    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["trait_id"] for r in rows] == ["t1", "t2"]
    assert rows[0]["file_path"] == "/data/t1.vcf.gz"
    assert rows[0]["n"] == "1000"
    with pytest.raises(ValueError, match="stored_effect_scale"):
        _read_manifest(path)

    # Parked (non-EUR/Unassigned) analyses are present and labelled, not dropped.
    import csv

    with open(path, newline="", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh, delimiter="\t"))
    assert records[0]["assigned_ancestry"] == "EUR"
    assert records[1]["assigned_ancestry"] == "Unassigned"
    assert records[0]["catalogue_version"] == "cat-v1"
