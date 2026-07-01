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
from opengwasdb.model.enums import PrimaryStorageLayout, StoredEffectScale
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.stats import p_value_from_z
from opengwasdb.traits.axis import traits_table_path, traits_tabix_path
from opengwasdb.variants import (
    VariantAxis,
    variant_alid_bytes_path,
    variant_alid_rows_path,
    variant_offsets_path,
    variant_tabix_path,
    variant_table_path,
)


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

    if manifest.primary_layout is PrimaryStorageLayout.RAGGED:
        return _validate_ragged_store(store_path, errors)
    return _validate_dense_store(store_path, errors)


def _validate_dense_store(store_path: Path, errors: list[str]) -> ValidationResult:
    index_path = store_path / "index.sqlite"
    data_path = store_path / "data.zarr"
    variants_path = variant_table_path(store_path)
    tabix_path = variant_tabix_path(store_path)
    offsets_path = variant_offsets_path(store_path)
    if not index_path.exists():
        errors.append("missing index.sqlite")
    if not data_path.exists():
        errors.append("missing data.zarr")
    if not variants_path.exists():
        errors.append("missing variants.tsv.gz")
    if not tabix_path.exists():
        errors.append("missing variants.tsv.gz.tbi")
    if not offsets_path.exists():
        errors.append("missing variant_offsets.npy")
    alid_bytes_path = variant_alid_bytes_path(store_path)
    alid_rows_path = variant_alid_rows_path(store_path)
    if not alid_bytes_path.exists():
        errors.append(
            "missing variant_alid_bytes.npy — rebuild the store to generate the ALID search index"
        )
    if not alid_rows_path.exists():
        errors.append(
            "missing variant_alid_rows.npy — rebuild the store to generate the ALID search index"
        )
    if errors:
        return ValidationResult(errors=errors)

    try:
        with connect(index_path) as connection:
            variant_axis = VariantAxis(store_path, connection)
            try:
                n_variants = _validate_variant_axis(variant_axis, errors)
                _validate_sqlite(connection, n_variants, errors)
            finally:
                variant_axis.close()
            n_analyses = count_rows(connection, "analyses")
            root = zarr.open_group(str(data_path), mode="r")
            _validate_dense_arrays(root, n_variants, n_analyses, errors)
            if not errors:
                _validate_top_hits(root, errors)
    except Exception as exc:  # noqa: BLE001 - validators should report actionable failures
        errors.append(f"validation failed: {exc}")
    return ValidationResult(errors=errors)


def _validate_ragged_store(store_path: Path, errors: list[str]) -> ValidationResult:
    index_path = store_path / "index.sqlite"
    data_path = store_path / "data.zarr"
    ragged_path = data_path / "ragged"

    for label, p in [
        ("index.sqlite", index_path),
        ("data.zarr", data_path),
        ("data.zarr/ragged", ragged_path),
        ("variants.tsv.gz", variant_table_path(store_path)),
        ("variants.tsv.gz.tbi", variant_tabix_path(store_path)),
        ("variant_alid_bytes.npy", variant_alid_bytes_path(store_path)),
        ("variant_alid_rows.npy", variant_alid_rows_path(store_path)),
        ("traits.tsv.gz", traits_table_path(store_path)),
    ]:
        if not p.exists():
            errors.append(f"missing {label}")
    if errors:
        return ValidationResult(errors=errors)

    try:
        root = zarr.open_group(str(ragged_path), mode="r")
        for name in ("offsets", "variant_index", "z", "se"):
            if name not in root:
                errors.append(f"missing data.zarr/ragged/{name}")
        if errors:
            return ValidationResult(errors=errors)

        offsets = root["offsets"][:]
        n_assoc = int(offsets[-1])
        for name in ("variant_index", "z", "se"):
            if len(root[name]) != n_assoc:
                errors.append(
                    f"data.zarr/ragged/{name} has {len(root[name])} entries "
                    f"but offsets imply {n_assoc}"
                )
        se_vals = root["se"][:].astype("float32")
        if np.any(np.isfinite(se_vals) & (se_vals < 0)):
            errors.append("se contains negative finite values")

        with sqlite3.connect(str(index_path)) as conn:
            tables = {r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
            if "analyses" not in tables:
                errors.append("index.sqlite is missing the analyses table")
            else:
                n_analyses_db = conn.execute("SELECT COUNT(*) FROM analyses").fetchone()[0]
                n_analyses_csr = len(offsets) - 1
                if n_analyses_db != n_analyses_csr:
                    errors.append(
                        f"analyses table has {n_analyses_db} rows but "
                        f"zarr CSR offsets imply {n_analyses_csr} analyses"
                    )
    except Exception as exc:  # noqa: BLE001
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


def _validate_sqlite(
    connection: sqlite3.Connection,
    n_variants: int,
    errors: list[str],
) -> None:
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
    alias_rows = connection.execute("SELECT alias, variant_index FROM variant_aliases").fetchall()
    for row in alias_rows:
        variant_index = int(row["variant_index"])
        if variant_index < 0 or variant_index >= n_variants:
            errors.append(f"alias {row['alias']!r} points to missing variant {variant_index}")
    rows = connection.execute("SELECT analysis_id, stored_effect_scale FROM analyses").fetchall()
    for row in rows:
        try:
            StoredEffectScale(row["stored_effect_scale"])
        except ValueError:
            errors.append(
                f"analysis {row['analysis_id']} has invalid stored_effect_scale "
                f"{row['stored_effect_scale']!r}"
            )


def _validate_variant_axis(variant_axis: VariantAxis, errors: list[str]) -> int:
    records = variant_axis.all()
    if variant_axis.n_variants != len(records):
        errors.append(
            f"variant_offsets.npy has {variant_axis.n_variants} rows but "
            f"variants.tsv.gz has {len(records)} rows"
        )
    if variant_axis._alid_bytes is not None and len(variant_axis._alid_bytes) != len(records):
        errors.append(
            f"variant_alid_bytes.npy has {len(variant_axis._alid_bytes)} entries but "
            f"variants.tsv.gz has {len(records)} rows"
        )
    seen_alids: set[str] = set()
    for expected_index, record in enumerate(records):
        if record.variant_index != expected_index:
            errors.append(
                f"variant table row {expected_index} has variant_index {record.variant_index}"
            )
            break
        if record.alid in seen_alids:
            errors.append("duplicate canonical variants in variant table")
            break
        seen_alids.add(record.alid)
    for expected_index in _representative_variant_indices(len(records)):
        record = records[expected_index]
        try:
            offset_record = variant_axis.by_index(expected_index)
        except ValueError as exc:
            errors.append(str(exc))
            break
        if offset_record != record:
            errors.append(f"variant offset for row {expected_index} points to a different row")
            break
        fetched = variant_axis.range(record.chromosome, record.position, record.position)
        if record not in fetched:
            errors.append(f"tabix index cannot fetch variant {record.alid}")
            break
    return len(records)


def _representative_variant_indices(n_variants: int) -> list[int]:
    if n_variants <= 0:
        return []
    if n_variants <= 1000:
        return list(range(n_variants))
    anchors = {0, n_variants // 2, n_variants - 1}
    step = max(1, n_variants // 997)
    anchors.update(range(0, n_variants, step))
    return sorted(index for index in anchors if 0 <= index < n_variants)


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
