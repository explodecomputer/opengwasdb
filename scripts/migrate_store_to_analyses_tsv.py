#!/usr/bin/env python3
"""Migrate a Store Release from the pre-issue-#22 layout to `analyses.tsv`
(ADR 0030, issue #24).

Reads a store's old SQLite `analyses` table (dense/hybrid stores only -- Ragged
stores keep their own divergent `analyses` schema, out of scope per ADR 0030)
and its `ancestry.tsv` sidecar, when present, derives what it can from them,
writes `analyses.tsv` + `overview.html`, drops the SQLite `analyses` table, and
folds `ancestry_provenance.json` (when present) into `manifest.json`.

Top-Hit Counts (`n_hits_5e8`/`n_hits_5e6`/`n_hits_5e4`, ADR 0032) are computed
for real from the store's existing top-hit index -- that index predates
`analyses.tsv` entirely, so unlike the fields below it is genuinely recoverable,
not something this script's "don't guess" rule applies to.

Fields the old layout never recorded are written honestly, not guessed:

  * `original_sd_method` becomes `unavailable` (`OriginalSdMethod.UNAVAILABLE`)
    -- phenotype-SD provenance (issue #18) postdates every old-layout store,
    so there is no source to derive it from. `unavailable` is a real vocabulary
    member for this column, unlike the others below.
  * `sample_size_kind`/`sample_size_scope`/`sample_size`, `n_cases`/`n_controls`,
    `original_effect_scale`, and `original_sd_dispersion` stay blank. None of
    them have an `unavailable` vocabulary member (only `original_sd_method`
    and `ancestry_assignment_method` are controlled vocabularies with such a
    value), and writing an arbitrary non-vocabulary placeholder into a
    validated column would be worse than an honest blank.
  * `ancestry_assignment_method` is `af_assigned` when the sidecar recorded an
    Assigned Ancestry (every pre-#22 sidecar was written by the AF-based NNLS
    mixture fit), else blank -- not `unassigned`, which means "assignment was
    attempted and the gates rejected it," a different state from "the old
    layout recorded nothing for this Analysis."

Migration is one-way: this script has no reverse mode, and the library carries
no runtime shim for reading the old layout afterward.

Usage:
  uv run python scripts/migrate_store_to_analyses_tsv.py /path/to/store.opengwasdb
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path

from opengwasdb.layouts.dense.build import AnalysisMetadata, add_hit_counts, write_analyses_tsv
from opengwasdb.model.enums import OriginalSdMethod
from opengwasdb.store.open import open_store

_ANCESTRY_SIDECAR = "ancestry.tsv"
_ANCESTRY_PROVENANCE = "ancestry_provenance.json"

# Issue #16 renamed StoredEffectScale.SD_UNITS ("sd_units") to
# StoredEffectScale.SD ("sd") without a migration path -- every store built
# before that rename (including every real UKB-B production store) still
# carries the old spelling in its SQLite `analyses` table. This is a known
# historical rename, not a guess: translate it so a migrated store's
# stored_effect_scale passes the current controlled vocabulary.
_LEGACY_STORED_EFFECT_SCALE = {"sd_units": "sd"}


def read_old_layout_analysis_ids(store_path: str | Path) -> list[str]:
    """`analysis_id`s from the old SQLite `analyses` table, in `analysis_index`
    order -- the pre-migration snapshot a caller diffs against
    `migrate_store`'s return value to verify the migration (issue #24 AC4)."""
    rows = _read_old_analyses(Path(store_path) / "index.sqlite")
    return [str(row["analysis_id"]) for row in rows]


def _read_old_analyses(index_path: Path) -> list[dict[str, object]]:
    connection = sqlite3.connect(str(index_path))
    connection.row_factory = sqlite3.Row
    try:
        columns = {row[1] for row in connection.execute("PRAGMA table_info(analyses)")}
        if not columns:
            raise ValueError(
                f"{index_path} has no analyses table -- nothing to migrate "
                "(already migrated, or not a dense/hybrid store)"
            )
        rows = connection.execute("SELECT * FROM analyses ORDER BY analysis_index").fetchall()
        return [dict(row) for row in rows]
    finally:
        connection.close()


def _read_ancestry_sidecar(store_path: Path) -> dict[str, dict[str, str]]:
    path = store_path / _ANCESTRY_SIDECAR
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["trait_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _drop_analyses_table(index_path: Path) -> None:
    connection = sqlite3.connect(str(index_path))
    try:
        connection.execute("DROP TABLE analyses")
        connection.commit()
    finally:
        connection.close()


def _fold_ancestry_provenance(store_path: Path) -> None:
    provenance_path = store_path / _ANCESTRY_PROVENANCE
    if not provenance_path.exists():
        return
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    open_store(store_path).amend_provenance({"ancestry": provenance})
    provenance_path.unlink()


def migrate_store(store_path: str | Path) -> list[str]:
    """Migrate one store in place. Returns the migrated `analysis_id`s, in
    `analysis_index` order -- compare against `read_old_layout_analysis_ids`
    (captured before calling this) to verify the migration preserved the same
    Analyses (issue #24 AC1/AC4)."""
    store_path = Path(store_path)
    index_path = store_path / "index.sqlite"
    old_rows = _read_old_analyses(index_path)
    ancestry_by_id = _read_ancestry_sidecar(store_path)

    analyses: list[AnalysisMetadata] = []
    for row in old_rows:
        analysis_id = str(row["analysis_id"])
        ancestry_row = ancestry_by_id.get(analysis_id, {})
        assigned_ancestry = ancestry_row.get("assigned_ancestry", "")
        old_scale = str(row["stored_effect_scale"])
        analyses.append(
            AnalysisMetadata(
                analysis_id=analysis_id,
                phenotype_id=row.get("phenotype_id") or None,  # type: ignore[arg-type]
                phenotype_label=row.get("phenotype_label") or None,  # type: ignore[arg-type]
                analysis_label=row.get("analysis_label") or None,  # type: ignore[arg-type]
                stored_effect_scale=_LEGACY_STORED_EFFECT_SCALE.get(old_scale, old_scale),
                assigned_ancestry=assigned_ancestry,
                ancestry_assignment_method="af_assigned" if assigned_ancestry else "",
                completed_against=ancestry_row.get("completed_against", ""),
                completion_n_missing_total=(
                    str(row["n_missing_off_panel"]) if "n_missing_off_panel" in row else ""
                ),
                original_sd_method=OriginalSdMethod.UNAVAILABLE.value,
            )
        )

    write_analyses_tsv(store_path, add_hit_counts(store_path, analyses))
    _drop_analyses_table(index_path)
    _fold_ancestry_provenance(store_path)
    (store_path / _ANCESTRY_SIDECAR).unlink(missing_ok=True)

    return [a.analysis_id for a in analyses]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store_path", type=Path, help="Store Release directory to migrate in place")
    args = parser.parse_args()

    before = read_old_layout_analysis_ids(args.store_path)
    after = migrate_store(args.store_path)
    assert before == after, "migration changed the Analysis set -- this is a bug"
    print(f"Migrated {len(after)} analyses at {args.store_path}")


if __name__ == "__main__":
    main()
