"""Subset the Analysis Catalogue to one ancestry and build its Hybrid store.

A build manifest is a pure **row filter** over the Catalogue (ADR 0027): because
the Catalogue is a superset of the build manifest's ``trait_id``/``file_path``/
``trait_name``/``n`` columns, filtering to ``assigned_ancestry == EUR`` yields
those columns already in the shape ``build-hybrid`` reads. Two columns the
Catalogue does not carry are ``stored_effect_scale`` (issue #17) and
``original_sd_method``/``original_sd`` (issue #18): ancestry assignment runs on
allele frequencies alone and may run before a study's effect scale or
phenotype SD is even resolved, so they stay genuinely separate build inputs --
supplied by the caller here and stamped onto the filtered rows before writing
the manifest, not something ancestry assignment should ever need to know
about. After the build, the store's per-Analysis Assigned Ancestry and
Catalogue provenance are recorded in a sidecar (``opengwasdb.ancestry.store``).
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from opengwasdb.ancestry.store import (
    write_ancestry_provenance,
    write_ancestry_sidecar,
)

log = logging.getLogger(__name__)

DEFAULT_ANCESTRY = "EUR"


@dataclass(frozen=True)
class SubsetResult:
    manifest_path: Path
    ancestry: str
    subset_filter: str
    n_total: int
    n_kept: int
    catalogue_version: str
    ancestry_reference_version: str
    analyses: list[tuple[str, str]]  # (trait_id, assigned_ancestry)


def subset_catalogue(
    catalogue_path: str | Path,
    manifest_path: str | Path,
    *,
    stored_effect_scale: str,
    original_sd_method: str,
    ancestry: str = DEFAULT_ANCESTRY,
    original_sd: str | float | None = None,
) -> SubsetResult:
    """Row-filter the Catalogue to ``assigned_ancestry == ancestry``.

    Writes the filtered rows verbatim (all Catalogue columns), plus
    ``stored_effect_scale`` and ``original_sd_method``/``original_sd`` stamped
    onto every row -- the columns `opengwasdb.layouts.dense.build_vcf._read_manifest`
    now requires that the Catalogue itself never carries (issues #17, #18).
    ``original_sd`` applies uniformly to every kept Analysis, like
    ``stored_effect_scale``; leave it ``None`` for ``original_sd_method``
    values that carry no SD magnitude (``declared_standardised``,
    ``binary_trait``). Returns the subset provenance.
    """
    catalogue_path = Path(catalogue_path)
    manifest_path = Path(manifest_path)
    with open(catalogue_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "assigned_ancestry" not in fieldnames:
            raise ValueError(f"catalogue lacks assigned_ancestry column: {catalogue_path}")
        rows = list(reader)

    kept = [r for r in rows if r.get("assigned_ancestry") == ancestry]
    for row in kept:
        row["stored_effect_scale"] = stored_effect_scale
        row["original_sd_method"] = original_sd_method
        row["original_sd"] = "" if original_sd is None else str(original_sd)
    out_fieldnames = [*fieldnames, "stored_effect_scale", "original_sd_method", "original_sd"]
    with open(manifest_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)

    catalogue_version = _single_value(kept or rows, "catalogue_version")
    reference_version = _single_value(kept or rows, "ancestry_reference_version")
    subset_filter = f"assigned_ancestry == {ancestry}"
    log.info(
        "Catalogue subset %s: %d/%d analyses → %s",
        subset_filter, len(kept), len(rows), manifest_path,
    )
    return SubsetResult(
        manifest_path=manifest_path,
        ancestry=ancestry,
        subset_filter=subset_filter,
        n_total=len(rows),
        n_kept=len(kept),
        catalogue_version=catalogue_version,
        ancestry_reference_version=reference_version,
        analyses=[(r["trait_id"], r["assigned_ancestry"]) for r in kept],
    )


def record_store_ancestry(store_path: str | Path, subset: SubsetResult) -> None:
    """Write the ancestry sidecar + provenance into a freshly built store."""
    write_ancestry_sidecar(store_path, subset.analyses)
    write_ancestry_provenance(
        store_path,
        catalogue_version=subset.catalogue_version,
        subset_filter=subset.subset_filter,
        ancestry_reference_version=subset.ancestry_reference_version,
        n_analyses=subset.n_kept,
    )


def build_hybrid_from_catalogue(
    catalogue_path: str | Path,
    output_path: str | Path,
    *,
    reference_panel: str | Path,
    store_id: str,
    release_id: str,
    stored_effect_scale: str,
    original_sd_method: str,
    ancestry: str = DEFAULT_ANCESTRY,
    original_sd: str | float | None = None,
    overwrite: bool = False,
    n_workers: int = 1,
    chunk_shape: tuple[int, int] | None = None,
    manifest_path: str | Path | None = None,
) -> SubsetResult:
    """Subset the Catalogue to ``ancestry`` and build a Hybrid store from it.

    ``stored_effect_scale`` and ``original_sd_method``/``original_sd`` apply
    uniformly to every kept Analysis (issues #17, #18) -- this entry point is
    for a single Source Collection release where that is expected to hold; a
    release mixing scales or SD methods per Analysis needs a manifest built
    some other way.

    Uses the unchanged ``build_hybrid_from_vcf_manifest`` on the row-filtered
    manifest, then records per-Analysis Assigned Ancestry + Catalogue provenance
    in the store sidecar.
    """
    # Imported here to keep the ancestry package importable without the heavy
    # build stack (and to avoid a build → ancestry → build import cycle).
    from opengwasdb.layouts.dense.constants import DEFAULT_CHUNK_SHAPE
    from opengwasdb.layouts.hybrid.build import build_hybrid_from_vcf_manifest

    output_path = Path(output_path)
    if manifest_path is None:
        manifest_path = output_path.with_suffix(".manifest.tsv")
    manifest_path = Path(manifest_path)

    subset = subset_catalogue(
        catalogue_path,
        manifest_path,
        ancestry=ancestry,
        stored_effect_scale=stored_effect_scale,
        original_sd_method=original_sd_method,
        original_sd=original_sd,
    )
    build_hybrid_from_vcf_manifest(
        subset.manifest_path,
        output_path,
        reference_panel=reference_panel,
        store_id=store_id,
        release_id=release_id,
        overwrite=overwrite,
        n_workers=n_workers,
        chunk_shape=chunk_shape or DEFAULT_CHUNK_SHAPE,
    )
    record_store_ancestry(output_path, subset)
    return subset


def _single_value(rows: list[dict[str, str]], column: str) -> str:
    """The common value of ``column`` across rows (empty when absent/mixed)."""
    values = {r.get(column, "") for r in rows}
    return next(iter(values)) if len(values) == 1 else ""
