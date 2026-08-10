"""Build a EUR Hybrid store from an Analysis Catalogue subset (issue 065).

The Catalogue is a superset of the build manifest, so selecting EUR is a pure row
filter. The unchanged build-hybrid consumes it; the store then carries per-Analysis
Assigned Ancestry + Catalogue provenance in a sidecar. Non-EUR/Unassigned Analyses
are absent from the store (still parked in the Catalogue).
"""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from opengwasdb.ancestry.store import read_ancestry_map, read_ancestry_provenance
from opengwasdb.ancestry.subset import build_hybrid_from_catalogue, subset_catalogue
from opengwasdb.cli.main import app
from opengwasdb.validation import validate_store


def _store_analysis_ids(store: Path) -> set[str]:
    """Read the store's analysis ids from its SQLite index."""
    import sqlite3

    conn = sqlite3.connect(store / "index.sqlite")
    try:
        return {row[0] for row in conn.execute("SELECT analysis_id FROM analyses")}
    finally:
        conn.close()

# hg19 coordinates that liftover to the hybrid fixture's hg38 ALIDs.
HG19_POS_1 = 100_000
HG19_POS_2 = 1_000_000
HG19_POS_3 = 1_500_000
HG38_ALID_1 = "1:100000:A:G"  # on-panel → Dense Component
HG38_ALID_3 = "1:1564620:A:G"  # on-panel → Dense Component


def _vcf(tmp_path: Path, name: str) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size">\n'
        '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">\n'
        '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score">\n'
        "##SAMPLE=<ID=STUDY1,StudyType=Continuous>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTUDY1\n"
    )
    rows = (
        f"1\t{HG19_POS_1}\t.\tA\tG\t.\tPASS\t.\tES:SE\t2.0:0.5\n"
        f"1\t{HG19_POS_2}\t.\tC\tT\t.\tPASS\t.\tES:SE\t1.5:0.3\n"
        f"1\t{HG19_POS_3}\t.\tG\tA\t.\tPASS\t.\tES:SE\t0.6:0.2\n"
    )
    path = tmp_path / f"{name}.vcf"
    path.write_text(header + rows, encoding="utf-8")
    return path


def _panel(tmp_path: Path) -> Path:
    panel = tmp_path / "panel.txt"
    panel.write_text(f"{HG38_ALID_1}\n{HG38_ALID_3}\n", encoding="utf-8")
    return panel


def _catalogue(tmp_path: Path) -> Path:
    """A Catalogue with two EUR, one AFR, and one Unassigned Analysis."""
    vcfs = {name: _vcf(tmp_path, name) for name in ("eur_a", "eur_b", "afr_c", "unassigned_d")}
    header = [
        "trait_id",
        "file_path",
        "trait_name",
        "n",
        "assigned_ancestry",
        "reported_population",
        "catalogue_version",
        "ancestry_reference_version",
    ]
    ver = ["cat-v1", "prive2022-hg38"]
    rows = [
        ["eur_a", str(vcfs["eur_a"]), "EUR A", "1000", "EUR", "European", *ver],
        ["eur_b", str(vcfs["eur_b"]), "EUR B", "1200", "EUR", "European", *ver],
        ["afr_c", str(vcfs["afr_c"]), "AFR C", "900", "AFR", "African", *ver],
        ["und_d", str(vcfs["unassigned_d"]), "Und", "800", "Unassigned", "Mixed", *ver],
    ]
    path = tmp_path / "catalogue.tsv"
    path.write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n", encoding="utf-8"
    )
    return path


def test_subset_catalogue_is_pure_row_filter(tmp_path):
    catalogue = _catalogue(tmp_path)
    manifest = tmp_path / "eur_manifest.tsv"
    result = subset_catalogue(catalogue, manifest, ancestry="EUR", stored_effect_scale="sd")
    assert result.n_total == 4 and result.n_kept == 2
    assert result.subset_filter == "assigned_ancestry == EUR"
    assert result.catalogue_version == "cat-v1"
    # The manifest is still a Catalogue superset: same header, EUR rows only.
    header = manifest.read_text().splitlines()[0].split("\t")
    assert header[0] == "trait_id" and "assigned_ancestry" in header
    from opengwasdb.layouts.dense.build_vcf import _read_manifest

    loaded = _read_manifest(manifest)
    assert [m.trait_id for m in loaded] == ["eur_a", "eur_b"]


def test_build_eur_hybrid_records_provenance_and_validates(tmp_path):
    catalogue = _catalogue(tmp_path)
    store = tmp_path / "eur.opengwasdb"
    subset = build_hybrid_from_catalogue(
        catalogue,
        store,
        reference_panel=_panel(tmp_path),
        store_id="eur-store",
        release_id="v1",
        stored_effect_scale="sd",
        ancestry="EUR",
    )
    assert subset.n_kept == 2

    # The store validates.
    result = validate_store(store)
    assert result.ok, result.errors

    # Only EUR-assigned Analyses are present in the store.
    assert _store_analysis_ids(store) == {"eur_a", "eur_b"}

    # Per-Analysis Assigned Ancestry sidecar.
    ancestry_map = read_ancestry_map(store)
    assert ancestry_map == {"eur_a": "EUR", "eur_b": "EUR"}

    # Store-level provenance links back to the Catalogue.
    prov = read_ancestry_provenance(store)
    assert prov["catalogue_version"] == "cat-v1"
    assert prov["subset_filter"] == "assigned_ancestry == EUR"
    assert prov["ancestry_reference_version"] == "prive2022-hg38"
    assert prov["n_analyses"] == 2


def test_build_hybrid_from_catalogue_cli(tmp_path):
    catalogue = _catalogue(tmp_path)
    store = tmp_path / "eur_cli.opengwasdb"
    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "build-hybrid-from-catalogue",
            str(catalogue),
            str(store),
            "--reference-panel",
            str(_panel(tmp_path)),
            "--store-id",
            "eur-store",
            "--release-id",
            "v1",
            "--stored-effect-scale",
            "sd",
        ],
    )
    assert result.exit_code == 0, result.output
    import json

    summary = json.loads(result.output.strip().splitlines()[-1])
    assert summary["ancestry"] == "EUR"
    assert summary["n_kept"] == 2 and summary["n_total"] == 4
    assert read_ancestry_map(store) == {"eur_a": "EUR", "eur_b": "EUR"}
