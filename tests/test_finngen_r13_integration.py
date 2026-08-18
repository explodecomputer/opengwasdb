"""Opt-in integration proof against one complete public FinnGen R13 endpoint.

Run with ``OPENGWASDB_RUN_FINNGEN_R13_INTEGRATION=1 pixi run -e dev pytest
-q tests/test_finngen_r13_integration.py``. The source download is roughly
760 MB, so it is deliberately excluded from the default unit-test run.
"""

from __future__ import annotations

import os
import urllib.request

import pytest

from opengwasdb.layouts.dense.build_vcf import build_dense_from_vcf_manifest
from opengwasdb.query import query_store
from opengwasdb.readers import FINNGEN_R13_CAPABILITY
from opengwasdb.validation import validate_store

ENDPOINT_URL = (
    "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
    "finngen_R13_AB1_ACTINOMYCOSIS.gz"
)


@pytest.mark.skipif(
    os.environ.get("OPENGWASDB_RUN_FINNGEN_R13_INTEGRATION") != "1",
    reason="set OPENGWASDB_RUN_FINNGEN_R13_INTEGRATION=1 for the 760 MB endpoint build",
)
def test_real_finngen_r13_endpoint_builds_and_queries(tmp_path):
    source = tmp_path / "finngen_R13_AB1_ACTINOMYCOSIS.gz"
    urllib.request.urlretrieve(ENDPOINT_URL, source)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "trait_id\tfile_path\ttrait_name\tn\tstored_effect_scale"
        "\toriginal_sd_method\toriginal_sd\tsource_reader_capability\tsource_assembly\n"
        f"finngen-r13-actinomycosis\t{source}\tActinomycosis\t433889\tlog_or"
        f"\tbinary_trait\t\t{FINNGEN_R13_CAPABILITY}\tGRCh38\n",
        encoding="utf-8",
    )
    store_path = tmp_path / "finngen-r13-actinomycosis.opengwasdb"

    build_dense_from_vcf_manifest(
        manifest, store_path, store_id="finngen-r13-integration", release_id="r13"
    )

    validation = validate_store(store_path)
    assert validation.ok, validation.errors
    query = query_store(store_path)
    result = query.lookup(["1:13668:A:G"], ["finngen-r13-actinomycosis"])
    query.close()
    assert result["z"][0] == pytest.approx(1.99175 / 1.46977, rel=5e-3)
