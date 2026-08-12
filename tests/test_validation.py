from __future__ import annotations

import sqlite3

import numpy as np
import zarr

from opengwasdb.layouts.dense.rho import build_dense_rho
from opengwasdb.model.analyses import read_analyses, write_analyses
from opengwasdb.validation import validate_store


def test_validator_rejects_missing_required_array(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    del root["se"]

    result = validate_store(dense_store_path)

    assert not result.ok
    assert "missing data.zarr/se" in result.errors


def test_validator_rejects_negative_se(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["se"][0, 0] = -0.1

    result = validate_store(dense_store_path)

    assert not result.ok
    assert "se contains negative finite values" in result.errors


def test_validator_rejects_inconsistent_missingness(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["se"][0, 0] = float("nan")

    result = validate_store(dense_store_path)

    assert not result.ok
    assert "z and se missingness is inconsistent" in result.errors


def test_validator_rejects_invalid_stored_effect_scale(dense_store_path):
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        if row["analysis_id"] == "a1":
            row["stored_effect_scale"] = "kg"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("invalid stored_effect_scale" in error for error in result.errors)


def test_validator_rejects_inconsistent_top_hit_index(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["top_hits"]["p_5e_08"]["z"][0] = 0.0

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("top-hit index p_5e_08" in error for error in result.errors)


def test_validator_rejects_invalid_top_hit_offsets(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="r+")
    offsets = root["top_hits/p_5e_08/analysis_offsets"][:]
    offsets[-1] -= 1
    root["top_hits/p_5e_08/analysis_offsets"][:] = offsets

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("invalid analysis offsets" in error for error in result.errors)


def test_validator_rejects_missing_variant_axis_file(dense_store_path):
    (dense_store_path / "variants.tsv.gz").unlink()

    result = validate_store(dense_store_path)

    assert not result.ok
    assert "missing variants.tsv.gz" in result.errors


def test_validator_rejects_bad_variant_offset(dense_store_path):
    offsets_path = dense_store_path / "variant_offsets.npy"
    offsets = np.load(offsets_path)
    offsets[1] = offsets[0]
    np.save(offsets_path, offsets)

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("variant offset for row 1" in error for error in result.errors)


def test_validator_rejects_leftover_sqlite_analyses_table(dense_store_path):
    # ADR 0030/issue #22: analyses.tsv is the sole source of truth for
    # Analytical Metadata; a store carrying a leftover SQLite `analyses`
    # table (store-format spec §20) is invalid even though nothing else
    # about it is malformed.
    connection = sqlite3.connect(dense_store_path / "index.sqlite")
    try:
        connection.execute("CREATE TABLE analyses (analysis_index INTEGER PRIMARY KEY)")
        connection.commit()
    finally:
        connection.close()

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("still contains an analyses table" in error for error in result.errors)


def test_validator_accepts_store_without_sqlite_analyses_table(dense_store_path):
    result = validate_store(dense_store_path)

    assert result.ok, result.errors


def test_validator_rejects_invalid_original_sd_method(dense_store_path):
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        if row["analysis_id"] == "a1":
            row["original_sd_method"] = "bogus_method"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("invalid original_sd_method" in error for error in result.errors)


def test_validator_accepts_valid_original_sd_method(dense_store_path):
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        row["original_sd_method"] = "declared_standardised"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert result.ok, result.errors


def test_validator_rejects_invalid_ancestry_assignment_method(dense_store_path):
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        if row["analysis_id"] == "a2":
            row["ancestry_assignment_method"] = "bogus_method"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("invalid ancestry_assignment_method" in error for error in result.errors)


def test_validator_accepts_valid_ancestry_assignment_method(dense_store_path):
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        row["ancestry_assignment_method"] = "af_assigned"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert result.ok, result.errors


def test_validator_rejects_analyses_tsv_missing_referenced_analysis_index(dense_store_path):
    # store-format spec §20: analyses.tsv must cover every analysis_index
    # referenced by the store's positional indexing -- a gap (here, a2's
    # index jumps from 1 to 2, leaving 1 uncovered) is the one relationship
    # SQLite cannot enforce as a foreign key, since analyses.tsv is a
    # separate file.
    analyses_path = dense_store_path / "analyses.tsv"
    table = read_analyses(analyses_path)
    for row in table.rows:
        if row["analysis_id"] == "a2":
            row["analysis_index"] = "2"
    write_analyses(analyses_path, table)

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any(
        "analysis_index values are not exactly the range" in error for error in result.errors
    )


def test_validator_accepts_analyses_tsv_covering_every_analysis_index(dense_store_path):
    result = validate_store(dense_store_path)

    assert result.ok, result.errors


# ── Rho Matrix (issue 049) ────────────────────────────────────────────────────
# window_bp=1 keeps every variant (100/200/300bp spacing); z_thresh=10.0 makes
# every finite z "null" -- the fixture then has exactly one shared-null variant
# (rs1, the only row both a1 and a2 have data for) for the single a1/a2 pair.


def _build_fixture_rho(dense_store_path, min_nulls: int = 1) -> None:
    build_dense_rho(dense_store_path, window_bp=1, z_thresh=10.0, min_nulls=min_nulls, n_workers=1)


def test_validator_accepts_store_without_rho_group(dense_store_path):
    result = validate_store(dense_store_path)

    assert result.ok, result.errors


def test_validator_accepts_well_formed_rho_matrix(dense_store_path):
    _build_fixture_rho(dense_store_path)

    result = validate_store(dense_store_path)

    assert result.ok, result.errors


def test_validator_rejects_rho_wrong_packed_length(dense_store_path):
    _build_fixture_rho(dense_store_path)
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    del root["rho"]["rho"]
    root["rho"].create_dataset("rho", data=np.array([0.1, 0.2], dtype="float16"))

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("rho array has" in error for error in result.errors)


def test_validator_rejects_rho_finite_with_insufficient_support(dense_store_path):
    _build_fixture_rho(dense_store_path)
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["rho"]["n_null"][:] = 0  # below min_nulls=1, but rho stays finite

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("rho finiteness is inconsistent" in error for error in result.errors)


def test_validator_rejects_rho_out_of_range(dense_store_path):
    _build_fixture_rho(dense_store_path)
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["rho"]["rho"][:] = np.array([1.5], dtype="float16")

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("rho contains values outside" in error for error in result.errors)


def test_validator_rejects_rho_missing_provenance_attr(dense_store_path):
    _build_fixture_rho(dense_store_path)
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    del root["rho"].attrs["window_bp"]

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("missing provenance attrs" in error for error in result.errors)


def test_validator_rejects_rho_variant_index_length_mismatch(dense_store_path):
    _build_fixture_rho(dense_store_path)
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["rho"].attrs["n_variants_used"] = 999

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("rho variant_index has" in error for error in result.errors)
