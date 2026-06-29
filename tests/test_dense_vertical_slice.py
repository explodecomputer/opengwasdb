from __future__ import annotations

import math

import numpy as np
import pytest
import zarr

from opengwasdb.index import connect, get_metadata
from opengwasdb.query import query_store
from opengwasdb.validation import validate_store


def test_dense_build_writes_standard_envelope_and_metadata(dense_store_path):
    assert (dense_store_path / "manifest.json").exists()
    assert (dense_store_path / "index.sqlite").exists()
    assert (dense_store_path / "data.zarr").exists()

    result = validate_store(dense_store_path)
    assert result.ok, result.errors

    root = zarr.open_group(str(dense_store_path / "data.zarr"), mode="r")
    assert root["z"].shape == (3, 2)
    assert root["se"].shape == (3, 2)
    assert root["z"].chunks == (1000, 1000)

    z = root["z"][:].astype("float32")
    se = root["se"][:].astype("float32")
    assert z[0, 0] == 2.0
    assert z[0, 1] == 6.0
    assert z[1, 0] == -3.0
    assert math.isnan(float(z[1, 1]))
    assert math.isnan(float(se[2, 0]))
    assert z[2, 1] == 6.0

    with connect(dense_store_path / "index.sqlite") as connection:
        assert get_metadata(connection, "n_variants") == 3
        assert get_metadata(connection, "n_analyses") == 2
        dense_meta = get_metadata(connection, "dense")
    assert dense_meta["chunk_shape"] == [1000, 1000]
    assert dense_meta["compressor"]["cname"] == "zstd"
    assert dense_meta["compressor"]["shuffle"] == "bitshuffle"


def test_query_facade_supports_variant_range_analysis_phewas_and_top_hits(dense_store_path):
    query = query_store(dense_store_path)

    by_alias = query.variant("rs1")
    assert [row.analysis_id for row in by_alias] == ["a1", "a2"]
    assert by_alias[1].beta == pytest.approx(1.2, rel=5e-4)
    assert by_alias[1].stored_effect_scale == "log_or"
    assert by_alias[1].p_value < 5e-8

    range_results = query.range("1", 1, 250)
    assert [(row.alid, row.analysis_id) for row in range_results] == [
        ("1:100:A:G", "a1"),
        ("1:100:A:G", "a2"),
        ("1:200:C:T", "a1"),
    ]

    a1 = query.analysis("a1")
    assert [(row.alid, row.z) for row in a1] == [("1:100:A:G", 2.0), ("1:200:C:T", -3.0)]

    phewas = query.phewas("1:100:A:G")
    assert [row.analysis_id for row in phewas] == ["a1", "a2"]

    lookup = query.lookup(["1:100:A:G", "1:300:A:G"], ["a1", "a2"])
    assert [(row.alid, row.analysis_id) for row in lookup] == [
        ("1:100:A:G", "a1"),
        ("1:100:A:G", "a2"),
        ("1:300:A:G", "a2"),
    ]

    top_hits = query.top_hits(threshold=5e-8)
    assert [(row.alid, row.analysis_id, row.z) for row in top_hits] == [
        ("1:100:A:G", "a2", 6.0),
        ("1:300:A:G", "a2", 6.0),
    ]


def test_query_excludes_missing_dense_cells(dense_store_path):
    query = query_store(dense_store_path)

    results = query.range("1", 150, 350)

    assert [(row.alid, row.analysis_id) for row in results] == [
        ("1:200:C:T", "a1"),
        ("1:300:A:G", "a2"),
    ]
    assert all(np.isfinite(row.z) and np.isfinite(row.se) for row in results)
