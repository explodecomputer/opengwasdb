from __future__ import annotations

import sqlite3

import numpy as np
import zarr

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
    with sqlite3.connect(dense_store_path / "index.sqlite") as connection:
        connection.execute(
            "UPDATE analyses SET stored_effect_scale = 'kg' WHERE analysis_id = 'a1'"
        )
        connection.commit()

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("invalid stored_effect_scale" in error for error in result.errors)


def test_validator_rejects_inconsistent_top_hit_index(dense_store_path):
    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="a")
    root["top_hits"]["p_5e_08"]["z"][0] = 0.0

    result = validate_store(dense_store_path)

    assert not result.ok
    assert any("top-hit index p_5e_08" in error for error in result.errors)


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
