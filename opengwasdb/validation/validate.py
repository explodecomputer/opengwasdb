"""Store validation."""

from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import zarr

from opengwasdb.index import connect, count_rows
from opengwasdb.layouts.dense.top_hits import threshold_key
from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.stats import p_value_from_z


@dataclass(frozen=True)
class ValidationResult:
    """Validation outcome with actionable error strings."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def validate_store(path: str | Path) -> ValidationResult:
    """Validate a v0.1 Store Release directory."""

    store_path = Path(path)
    errors: list[str] = []
    manifest = _load_manifest(store_path, errors)
    if manifest is None:
        return ValidationResult(errors=errors)
    index_path = store_path / "index.sqlite"
    data_path = store_path / "data.zarr"
    if not index_path.exists():
        errors.append("missing index.sqlite")
    if not data_path.exists():
        errors.append("missing data.zarr")
    if errors:
        return ValidationResult(errors=errors)

    try:
        with connect(index_path) as connection:
            _validate_sqlite(connection, errors)
            n_variants = count_rows(connection, "variants")
            n_analyses = count_rows(connection, "analyses")
            root = zarr.open_group(str(data_path), mode="r")
            _validate_dense_arrays(root, n_variants, n_analyses, errors)
            if not errors:
                _validate_top_hits(root, errors)
    except Exception as exc:  # noqa: BLE001 - validators should report actionable failures
        errors.append(f"validation failed: {exc}")
    return ValidationResult(errors=errors)


def _load_manifest(store_path: Path, errors: list[str]) -> StoreManifest | None:
    try:
        return StoreManifest.load(store_path)
    except KeyError as exc:
        errors.append(f"manifest missing required field: {exc.args[0]}")
    except ValueError as exc:
        errors.append(f"manifest has invalid enum value: {exc}")
    except FileNotFoundError:
        errors.append("missing manifest.json")
    except Exception as exc:  # noqa: BLE001
        errors.append(f"manifest is malformed: {exc}")
    return None


def _validate_sqlite(connection: sqlite3.Connection, errors: list[str]) -> None:
    duplicates = connection.execute(
        """
        SELECT alid, COUNT(*) AS n
        FROM variants
        GROUP BY alid
        HAVING n > 1
        """
    ).fetchall()
    if duplicates:
        errors.append("duplicate canonical variants in variants table")
    rows = connection.execute("SELECT analysis_id, stored_effect_scale FROM analyses").fetchall()
    for row in rows:
        try:
            StoredEffectScale(row["stored_effect_scale"])
        except ValueError:
            errors.append(
                f"analysis {row['analysis_id']} has invalid stored_effect_scale "
                f"{row['stored_effect_scale']!r}"
            )


def _validate_dense_arrays(
    root: Any,
    n_variants: int,
    n_analyses: int,
    errors: list[str],
) -> None:
    for name in ("z", "se"):
        if name not in root:
            errors.append(f"missing data.zarr/{name}")
    if errors:
        return
    z = root["z"][:].astype("float32")
    se = root["se"][:].astype("float32")
    expected_shape = (n_variants, n_analyses)
    if tuple(z.shape) != expected_shape:
        errors.append(f"z shape {tuple(z.shape)} does not match {expected_shape}")
    if tuple(se.shape) != expected_shape:
        errors.append(f"se shape {tuple(se.shape)} does not match {expected_shape}")
    if np.any(np.isfinite(se) & (se < 0)):
        errors.append("se contains negative finite values")
    if np.any(np.isnan(z) != np.isnan(se)):
        errors.append("z and se missingness is inconsistent")


def _validate_top_hits(root: Any, errors: list[str]) -> None:
    if "top_hits" not in root:
        return
    top = root["top_hits"]
    z_matrix = root["z"][:].astype("float32")
    for key in top:
        group = top[key]
        threshold = float(group.attrs.get("threshold", key.replace("p_", "").replace("_", "-")))
        rows = group["variant_index"][:]
        cols = group["analysis_index"][:]
        z_values = group["z"][:].astype("float32")
        abs_z = group["abs_z"][:].astype("float32")
        observed = set()
        if len(rows) != len(cols) or len(rows) != len(z_values) or len(rows) != len(abs_z):
            errors.append(f"top-hit index {key} has inconsistent array lengths")
            continue
        previous_abs = math.inf
        for row, col, z, indexed_abs in zip(rows, cols, z_values, abs_z, strict=True):
            stored_z = float(z_matrix[int(row), int(col)])
            z_matches = np.isclose(stored_z, float(z), rtol=1e-3, atol=1e-3)
            if not math.isfinite(stored_z) or not z_matches:
                errors.append(f"top-hit index {key} contains z value inconsistent with z array")
                break
            if p_value_from_z(stored_z) > threshold:
                errors.append(f"top-hit index {key} contains association above threshold")
                break
            if float(indexed_abs) > previous_abs:
                errors.append(f"top-hit index {key} is not ranked by descending significance")
                break
            previous_abs = float(indexed_abs)
            observed.add((int(row), int(col)))
        expected = {
            (int(row), int(col))
            for row, col in zip(*np.where(np.isfinite(z_matrix)), strict=True)
            if p_value_from_z(float(z_matrix[int(row), int(col)])) <= threshold
        }
        if observed != expected:
            errors.append(f"top-hit index {key} does not match stored z values")


def default_top_hit_key(threshold: float = 5e-8) -> str:
    return threshold_key(threshold)
