"""Structural validation shared by every layout's Top-Hit Index."""

from __future__ import annotations

from typing import Any

import numpy as np


def validate_group_structure(
    key: str,
    group: Any,
    n_analyses: int,
    *,
    imputed_required: bool,
    check_cell_uniqueness: bool = False,
) -> tuple[list[str], dict[str, np.ndarray] | None]:
    """Checks common to every layout's threshold group: required fields
    present, array lengths consistent, ``analysis_offsets`` valid, each
    Analysis's slice in genomic (variant_index-ascending) order, and the
    ``imputed`` column (when present) holding only 0/1 -- required when the
    underlying store is Reference-Completed.

    Returns ``(errors, arrays)``. ``arrays`` is ``None`` if any structural
    check failed (nothing further should be validated against this group);
    otherwise a dict of the loaded columns (``rows``, ``cols``, ``z_values``,
    ``abs_z``, ``offsets``, ``imputed_values``) for the caller's own
    layout-specific cross-validation (streaming the Dense matrix; sampling
    the Ragged CSR).
    """
    errors: list[str] = []
    required = {
        "analysis_offsets", "variant_index", "analysis_index", "abs_z", "z", "se", "p_value",
    }
    missing = sorted(required.difference(group.keys()))
    if missing:
        errors.append(f"top-hit index {key} is missing {', '.join(missing)}")
        return errors, None

    rows = group["variant_index"][:].astype(np.int64)
    cols = group["analysis_index"][:].astype(np.int64)
    z_values = group["z"][:].astype("float32")
    abs_z = group["abs_z"][:].astype("float32")
    offsets = group["analysis_offsets"][:].astype(np.int64)
    imputed_values = group["imputed"][:].astype(np.uint8) if "imputed" in group else None

    if imputed_required and imputed_values is None:
        errors.append(f"top-hit index {key} is missing imputed completion status")
        return errors, None

    if (
        len(rows) != len(cols)
        or len(rows) != len(z_values)
        or len(rows) != len(abs_z)
        or len(rows) != len(group["se"])
        or len(rows) != len(group["p_value"])
        or (imputed_values is not None and len(rows) != len(imputed_values))
    ):
        errors.append(f"top-hit index {key} has inconsistent array lengths")
        return errors, None

    if (
        len(offsets) != n_analyses + 1
        or len(offsets) == 0
        or offsets[0] != 0
        or np.any(offsets[:-1] > offsets[1:])
        or offsets[-1] != len(rows)
    ):
        errors.append(f"top-hit index {key} has invalid analysis offsets")
        return errors, None

    ordered = True
    for analysis_index in range(n_analyses):
        start, stop = int(offsets[analysis_index]), int(offsets[analysis_index + 1])
        if np.any(cols[start:stop] != analysis_index):
            ordered = False
            break
        if stop - start > 1 and np.any(rows[start : stop - 1] >= rows[start + 1 : stop]):
            ordered = False
            break
    if not ordered:
        errors.append(f"top-hit index {key} has incorrect or non-genomic analysis slices")
        return errors, None

    if imputed_values is not None and not np.all((imputed_values == 0) | (imputed_values == 1)):
        errors.append(f"top-hit index {key} has invalid imputed completion status")
        return errors, None

    if check_cell_uniqueness and len(rows):
        flat = rows * n_analyses + cols
        if len(np.unique(flat)) != len(flat):
            errors.append(f"top-hit index {key} does not match stored z values")
            return errors, None

    return errors, {
        "rows": rows,
        "cols": cols,
        "z_values": z_values,
        "abs_z": abs_z,
        "offsets": offsets,
        "imputed_values": imputed_values,
    }
