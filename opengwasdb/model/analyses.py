"""Analytical Metadata contract for `analyses.tsv` (store-format spec §7a).

`analyses.tsv` spans three column classes (opengwasdb-stores
docs/release-metadata-schema.md, "Column classes"): a *shared core* owned by
`opengwasdb` that carries interpretation-bearing Analysis metadata and appears
in both a registry release manifest and a built store; *registry-only*
build-input columns (source file location, checksum, licence, inclusion
reason, ...) that a manifest carries but a built store does not; and
*store-only* columns (Reference Completion quality rollups, ...) that cannot
exist before a build has run. A manifest and a built store's `analyses.tsv`
are therefore superset-related, not identical: this module's reader accepts
any column beyond the ones it knows about and carries it through unchanged in
`AnalysesTable.rows`, so a caller reading a manifest never has to strip
registry-only columns first, and a caller reading a built store never has to
strip store-only ones.

This is the schema `opengwasdb-stores#51` needs to validate its emitted
manifests against, replacing that repository's hand-rolled, per-generator
column-name checks (issue #16). Controlled-vocabulary literals here match the
registry's already-shipped vocabulary (`opengwasdb-stores`
docs/release-metadata-schema.md) rather than this package's earlier, unshipped
draft spelling -- see `opengwasdb.model.enums`.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from opengwasdb.model.enums import (
    AncestryAssignmentMethod,
    OriginalSdMethod,
    SampleSizeKind,
    SampleSizeScope,
    StoredEffectScale,
)

ANCESTRY_PROP_PREFIX = "ancestry_prop_"

# Columns interpretation-bearing enough that opengwasdb owns their meaning and
# vocabulary, present in both a registry manifest and a built store's
# analyses.tsv (store-format spec §7a). `ancestry_prop_<population>` columns
# are dynamic (one per reference population) and matched via
# ANCESTRY_PROP_PREFIX rather than listed here.
SHARED_CORE_COLUMNS: tuple[str, ...] = (
    "analysis_index",
    "analysis_id",
    "phenotype_id",
    "phenotype_label",
    "analysis_label",
    "stored_effect_scale",
    "assigned_ancestry",
    "ancestry_assignment_method",
    "sample_size_kind",
    "sample_size_scope",
    "sample_size",
    "n_cases",
    "n_controls",
    "original_effect_scale",
    "original_sd",
    "original_sd_method",
    "original_sd_dispersion",
)

# Top-Hit Count columns (ADR 0032, store-format spec §7a): one persisted
# per-Analysis count per threshold tier in
# opengwasdb.layouts.dense.constants.TOP_HIT_THRESHOLDS (5e-8/5e-6/5e-4, in
# that order). Not imported from there to keep this module layout-independent
# -- it is also the schema opengwasdb-stores validates manifests against.
TOP_HIT_COUNT_COLUMNS: tuple[str, ...] = ("n_hits_5e8", "n_hits_5e6", "n_hits_5e4")

# Produced during or after the build; must not be required of a release
# manifest, which is generated before a build runs.
STORE_ONLY_COLUMNS: tuple[str, ...] = (
    "completed_against",
    "completion_median_pearson_r",
    "completion_n_imputed_total",
    "completion_n_missing_total",
    *TOP_HIT_COUNT_COLUMNS,
)

# A reference listing of known registry build-input columns (issue #9's
# worked examples), for callers that want to distinguish "a column the
# registry is known to use" from "a column nobody has documented yet."
# classify_column() itself does NOT consult this tuple: every column outside
# SHARED_CORE_COLUMNS and STORE_ONLY_COLUMNS classifies as registry-only
# regardless of whether it appears here, so an undocumented registry
# build-input column still passes through correctly (the superset property).
REGISTRY_ONLY_COLUMNS: tuple[str, ...] = (
    "source_analysis_id",
    "source_label",
    "trait_ontology_name",
    "trait_ontology_id",
    "source_file",
    "source_bundle_id",
    "checksum",
    "checksum_algorithm",
    "size_bytes",
    "source_genome_build",
    "license",
    "publication_doi",
    "publication_pmid",
    "consortium",
    "source_ancestry_label",
    "analysis_group_id",
    "inclusion_reason",
    "exclude_from_build",
)

# The subset of SHARED_CORE_COLUMNS that must be present -- as columns, not
# necessarily non-blank per-row values -- in any analyses.tsv, manifest or
# store. Excludes columns opengwasdb only knows how to populate once a build
# has assigned ancestry/phenotype identity (phenotype_id, phenotype_label,
# analysis_label, assigned_ancestry, ancestry_prop_*, analysis_index) and
# columns that are legitimately blank pending resolution (original_sd,
# original_sd_dispersion, n_cases, n_controls -- the latter two required only
# conditionally, see _validate_case_control_counts below).
REQUIRED_COLUMNS: tuple[str, ...] = (
    "analysis_id",
    "stored_effect_scale",
    "sample_size_kind",
    "sample_size_scope",
    "sample_size",
    "original_effect_scale",
    "original_sd_method",
    "ancestry_assignment_method",
)

_VOCABULARIES: dict[str, type[StrEnum]] = {
    "stored_effect_scale": StoredEffectScale,
    "sample_size_kind": SampleSizeKind,
    "sample_size_scope": SampleSizeScope,
    "original_sd_method": OriginalSdMethod,
    "ancestry_assignment_method": AncestryAssignmentMethod,
}

_CASE_CONTROL_SCALE_VALUES = {StoredEffectScale.LOG_OR.value, StoredEffectScale.LOG_HAZARD.value}


class ColumnClass(StrEnum):
    """Which of the three column classes a column belongs to."""

    SHARED_CORE = "shared_core"
    REGISTRY_ONLY = "registry_only"
    STORE_ONLY = "store_only"


def classify_column(name: str) -> ColumnClass:
    """Classify `name` as shared-core, registry-only, or store-only.

    A column absent from every list above still classifies as registry-only:
    the superset property means a manifest may legitimately carry a
    build-input column this module has never been told about.
    """
    if name in SHARED_CORE_COLUMNS or name.startswith(ANCESTRY_PROP_PREFIX):
        return ColumnClass.SHARED_CORE
    if name in STORE_ONLY_COLUMNS:
        return ColumnClass.STORE_ONLY
    return ColumnClass.REGISTRY_ONLY


@dataclass(frozen=True)
class AnalysesTable:
    """A parsed `analyses.tsv`: column order plus raw string rows.

    Rows are kept as raw strings, like `opengwasdb.ancestry.catalogue`'s
    Catalogue reader -- `analyses.tsv` is inherently a text format, and typed
    access belongs to callers that know which columns they need, not to this
    module.
    """

    fieldnames: tuple[str, ...]
    rows: tuple[dict[str, str], ...]


def read_analyses(path: str | Path) -> AnalysesTable:
    """Read an `analyses.tsv` file, preserving column order and every field.

    Columns beyond the known contract (registry-only build inputs, or
    anything else) are read without error and carried through in each row's
    dict unchanged -- the superset property a manifest and a built store's
    analyses.tsv are related by.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        rows = tuple(dict(row) for row in reader)
    return AnalysesTable(fieldnames=fieldnames, rows=rows)


def write_analyses(path: str | Path, table: AnalysesTable) -> Path:
    """Write `table` back to a TSV file, in its own column order."""
    out_path = Path(path)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table.fieldnames), delimiter="\t")
        writer.writeheader()
        for row in table.rows:
            writer.writerow(row)
    return out_path


def validate_analyses(table: AnalysesTable) -> list[str]:
    """Validate `table` against the analyses.tsv contract.

    Returns a list of error strings (empty means valid), following this
    package's existing convention (`opengwasdb.validation.validate`) of
    accumulating actionable errors rather than raising on the first problem.
    """
    errors: list[str] = []
    fieldnames = set(table.fieldnames)

    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        errors.append(f"analyses.tsv is missing required column(s): {', '.join(missing)}")

    for row in table.rows:
        analysis_id = row.get("analysis_id") or "<unknown analysis_id>"
        # A required column present in the header but blank on this row is
        # equally a missing-required-column failure from the row's point of
        # view -- it just cannot be caught by the header-level check above.
        for column in REQUIRED_COLUMNS:
            if column in fieldnames and not row.get(column, ""):
                errors.append(
                    f"analysis {analysis_id!r} has no value for required column {column!r}"
                )
        for column, vocabulary in _VOCABULARIES.items():
            if column not in fieldnames:
                continue
            value = row.get(column, "")
            if not value:
                continue
            try:
                vocabulary(value)
            except ValueError:
                allowed = [member.value for member in vocabulary]
                errors.append(
                    f"analysis {analysis_id!r} has invalid {column} {value!r}; "
                    f"expected one of {allowed}"
                )
        _validate_case_control_counts(row, analysis_id, fieldnames, errors)

    return errors


def _validate_case_control_counts(
    row: dict[str, str],
    analysis_id: str,
    fieldnames: set[str],
    errors: list[str],
) -> None:
    if row.get("stored_effect_scale", "") not in _CASE_CONTROL_SCALE_VALUES:
        return
    scale = row["stored_effect_scale"]
    for column in ("n_cases", "n_controls"):
        if column in fieldnames and row.get(column, ""):
            continue
        errors.append(f"analysis {analysis_id!r} has stored_effect_scale={scale!r} but no {column}")


def to_json_schema() -> dict[str, Any]:
    """A JSON Schema description of the analyses.tsv contract.

    The machine-readable artifact issue #16 asks for, so a consumer outside
    this package (`opengwasdb-stores`) can validate emitted manifests without
    vendoring or hand-copying the vocabulary. `additionalProperties: true`
    encodes the superset property directly: registry-only, store-only, and
    genuinely unrecognised columns are all schema-valid.
    """
    properties: dict[str, Any] = {}
    for column in (*SHARED_CORE_COLUMNS, *STORE_ONLY_COLUMNS):
        prop: dict[str, Any] = {"type": "string"}
        vocabulary = _VOCABULARIES.get(column)
        if vocabulary is not None:
            prop["enum"] = [member.value for member in vocabulary]
        properties[column] = prop
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "OpenGWASDB analyses.tsv",
        "description": "Analytical Metadata contract, store-format spec §7a.",
        "type": "object",
        "properties": properties,
        "required": list(REQUIRED_COLUMNS),
        "additionalProperties": True,
    }
