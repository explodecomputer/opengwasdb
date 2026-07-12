"""Tests for store routing: reported-fallback + coverage gate (rule 2)."""

from __future__ import annotations

import csv

from opengwasdb.ancestry.routing import (
    ROUTING_COLUMNS,
    finalize_catalogue,
    reported_to_superpop,
    resolve_routing,
)

GW = dict(total_variants=6_000_000, n_autosomes=22, frac_largest_chrom=0.09)  # genome-wide


def test_reported_to_superpop():
    assert reported_to_superpop("European") == "EUR"
    assert reported_to_superpop("European (Sardinian)") == "EUR"
    assert reported_to_superpop("East Asian") == "EAS"
    assert reported_to_superpop("African American") == "AFR"
    assert reported_to_superpop("Hispanic or Latin American") == "AMR"
    assert reported_to_superpop("Iranian") == "MID"
    assert reported_to_superpop("Indian") == "SAS"
    assert reported_to_superpop("Mixed") is None
    assert reported_to_superpop("") is None


def test_af_assigned_takes_precedence_over_reported():
    r = resolve_routing("EUR", "Mixed", **GW)
    assert r.routing_ancestry == "EUR"
    assert r.routing_source == "af-assigned"
    assert r.store_eligible


def test_reported_fallback_when_unassigned():
    r = resolve_routing("Unassigned", "East Asian", **GW)
    assert r.routing_ancestry == "EAS"
    assert r.routing_source == "reported-fallback"
    assert r.store_eligible


def test_unassigned_and_mixed_is_unroutable():
    r = resolve_routing("Unassigned", "Mixed", **GW)
    assert r.routing_ancestry == ""
    assert r.routing_source == ""
    assert not r.store_eligible


def test_low_coverage_keeps_ancestry_but_ineligible():
    # Valid AF ancestry, but an array-based study (<500k) → flagged, not dropped.
    r = resolve_routing(
        "EUR", "European", total_variants=191_000, n_autosomes=22, frac_largest_chrom=0.08
    )
    assert r.routing_ancestry == "EUR"  # ancestry retained
    assert r.low_coverage
    assert not r.store_eligible


def test_missing_autosome_is_low_coverage():
    r = resolve_routing("EUR", "European", total_variants=6_000_000, n_autosomes=21,
                        frac_largest_chrom=0.09)
    assert not r.genome_distributed and r.low_coverage and not r.store_eligible


def test_single_chromosome_concentration_is_low_coverage():
    r = resolve_routing("EUR", "European", total_variants=6_000_000, n_autosomes=22,
                        frac_largest_chrom=0.5)  # half the variants on one chrom
    assert r.low_coverage and not r.store_eligible


def test_finalize_catalogue_augments_and_tallies(tmp_path):
    cat = tmp_path / "catalogue.tsv"
    cat.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tassigned_ancestry\treported_population\n"
        "s1\t/d/s1.vcf.gz\tS1\t100\tEUR\tMixed\n"            # af-assigned, GW → keep
        "s2\t/d/s2.vcf.gz\tS2\t100\tUnassigned\tEuropean\n"  # fallback EUR, GW → keep
        "s3\t/d/s3.vcf.gz\tS3\t100\tUnassigned\tMixed\n"     # unroutable → drop (ancestry)
        "s4\t/d/s4.vcf.gz\tS4\t100\tEUR\tEuropean\n",        # af EUR but low cov → drop (coverage)
        encoding="utf-8",
    )
    cov = tmp_path / "coverage.tsv"
    cov.write_text(
        "trait_id\ttotal_variants\tn_autosomes\tfrac_largest_chrom\n"
        "s1\t6000000\t22\t0.09\n"
        "s2\t6000000\t22\t0.09\n"
        "s3\t6000000\t22\t0.09\n"
        "s4\t191000\t22\t0.08\n",
        encoding="utf-8",
    )
    out = tmp_path / "routed.tsv"
    tally = finalize_catalogue(cat, cov, out)
    assert tally == {"kept": 2, "dropped_ancestry": 1, "dropped_coverage": 1,
                     "reported_fallback": 1}

    rows = {r["trait_id"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    for col in ROUTING_COLUMNS:
        assert col in rows["s1"]
    assert rows["s1"]["routing_ancestry"] == "EUR" and rows["s1"]["store_eligible"] == "True"
    assert rows["s2"]["routing_source"] == "reported-fallback"
    assert rows["s3"]["routing_ancestry"] == "" and rows["s3"]["store_eligible"] == "False"
    assert rows["s4"]["routing_ancestry"] == "EUR" and rows["s4"]["low_coverage"] == "True"
    assert rows["s4"]["store_eligible"] == "False"

    # Still a manifest superset: the unchanged reader consumes it.
    from opengwasdb.layouts.dense.build_vcf import _read_manifest
    assert {m.trait_id for m in _read_manifest(out)} == {"s1", "s2", "s3", "s4"}
