"""Store validation."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import zarr

from opengwasdb.index import connect, count_rows
from opengwasdb.layouts.dense.top_hits import threshold_key, z_critical
from opengwasdb.layouts.ragged.zarr_csr import RaggedCSRReader
from opengwasdb.model.enums import CompletionState, PrimaryStorageLayout, StoredEffectScale
from opengwasdb.model.manifest import StoreManifest
from opengwasdb.stats import p_value_from_z
from opengwasdb.traits.axis import traits_table_path
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
        return _validate_ragged_store(store_path, manifest, errors)
    return _validate_dense_store(store_path, manifest, errors)


def _validate_dense_store(
    store_path: Path, manifest: StoreManifest, errors: list[str]
) -> ValidationResult:
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
            imputed_arr = None
            on_panel_arr = None
            if manifest.completion_state is CompletionState.REFERENCE_COMPLETED:
                imputed_arr, on_panel_arr = _validate_completion_metadata(
                    root, connection, n_variants, n_analyses, errors
                )
            if not errors:
                _validate_dense_arrays(
                    root, n_variants, n_analyses, errors, imputed_arr, on_panel_arr
                )
            if not errors:
                _validate_top_hits(root, errors)
    except Exception as exc:  # noqa: BLE001 - validators should report actionable failures
        errors.append(f"validation failed: {exc}")
    return ValidationResult(errors=errors)


def _validate_completion_metadata(
    root: Any,
    connection: sqlite3.Connection,
    n_variants: int,
    n_analyses: int,
    errors: list[str],
) -> tuple[Any, np.ndarray | None]:
    """Validate the non-matrix completion metadata for a Reference-Completed
    Dense store and hand back the arrays the streamed band pass needs.

    The matrix-touching checks (imputed values are 0/1, imputed cells have
    finite z/se, off-panel rows are never imputed) are *not* done here — they
    are folded into ``_validate_dense_arrays``' row-band loop so the full
    z/se/imputed matrices are never resident (issue 045). This function only
    runs the cheap structural/metadata checks and returns the ``imputed`` zarr
    array (lazy) and the fully-loaded 1-D ``on_panel`` array (~n_variants bytes)
    for that loop. Returns ``(None, None)`` if the arrays are missing/malformed,
    which the caller treats as a hard error and skips the band pass.
    """
    for name in ("imputed", "on_panel"):
        if name not in root:
            errors.append(f"reference-completed store is missing data.zarr/{name}")
    if errors:
        return None, None

    imputed_arr = root["imputed"]
    expected_shape = (n_variants, n_analyses)
    if tuple(imputed_arr.shape) != expected_shape:
        errors.append(
            f"data.zarr/imputed shape {tuple(imputed_arr.shape)} does not match {expected_shape}"
        )

    on_panel = root["on_panel"][:]
    if len(on_panel) != n_variants:
        errors.append(f"data.zarr/on_panel has {len(on_panel)} entries but expected {n_variants}")
    if not np.all((on_panel == 0) | (on_panel == 1)):
        errors.append("data.zarr/on_panel contains values other than 0 and 1")

    tables = {r[0] for r in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()}
    if "completion_quality" not in tables:
        errors.append(
            "index.sqlite is missing the completion_quality table for a reference-completed store"
        )
    else:
        cols = {
            r[1] for r in connection.execute("PRAGMA table_info(completion_quality)").fetchall()
        }
        required = {"analysis_index", "block_id", "pearson_r", "n_imputed", "n_missing"}
        missing_cols = required - cols
        if missing_cols:
            errors.append(
                f"completion_quality table is missing columns: {', '.join(sorted(missing_cols))}"
            )

    analyses_cols = {r[1] for r in connection.execute("PRAGMA table_info(analyses)").fetchall()}
    if "n_missing_off_panel" not in analyses_cols:
        errors.append(
            "analyses table is missing n_missing_off_panel for a reference-completed store"
        )

    # Only hand the arrays to the band pass if they are well-formed; a shape
    # mismatch above already recorded a hard error, so the caller skips the pass.
    if errors:
        return None, None
    return imputed_arr, on_panel


def _validate_ragged_store(
    store_path: Path,
    manifest: StoreManifest,
    errors: list[str],
) -> ValidationResult:
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

        data_root = zarr.open_group(str(data_path), mode="r")

        # Reference-completed stores: validate imputed array and quality table.
        if manifest.completion_state is CompletionState.REFERENCE_COMPLETED:
            _validate_ragged_completion(ragged_path, index_path, n_assoc, errors)

        if not errors and "top_hits" in data_root:
            _validate_ragged_top_hits(store_path, data_root, errors)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"validation failed: {exc}")
    return ValidationResult(errors=errors)


def _validate_ragged_completion(
    ragged_path: Path,
    index_path: Path,
    n_assoc: int,
    errors: list[str],
) -> None:
    """Validate the imputed mask and completion_quality table in a Reference-Completed store."""
    root = zarr.open_group(str(ragged_path), mode="r")
    if "imputed" not in root:
        errors.append("reference-completed store is missing data.zarr/ragged/imputed")
        return

    imp = root["imputed"][:]
    if len(imp) != n_assoc:
        errors.append(
            f"data.zarr/ragged/imputed has {len(imp)} entries but offsets imply {n_assoc}"
        )
        return

    # imputed values must be 0 or 1
    if not np.all((imp == 0) | (imp == 1)):
        errors.append("data.zarr/ragged/imputed contains values other than 0 and 1")

    # Where imputed=1: z and se must both be finite
    z_vals = root["z"][:].astype("float32")
    se_vals = root["se"][:].astype("float32")
    imp_mask = imp == 1
    if imp_mask.any():
        if not np.all(np.isfinite(z_vals[imp_mask])):
            errors.append("imputed=1 rows have NaN z-scores")
        if not np.all(np.isfinite(se_vals[imp_mask])):
            errors.append("imputed=1 rows have NaN se values")

    # Where z is NaN: se must also be NaN and imputed must be 0
    nan_z = ~np.isfinite(z_vals)
    if nan_z.any():
        if not np.all(~np.isfinite(se_vals[nan_z])):
            errors.append("NaN z rows have finite se values (inconsistent missingness)")
        if np.any(imp[nan_z] == 1):
            errors.append("NaN z rows have imputed=1 (inconsistent)")

    with sqlite3.connect(str(index_path)) as conn:
        tables = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        if "completion_quality" not in tables:
            errors.append(
                "index.sqlite is missing the completion_quality table "
                "for a reference-completed store"
            )
        else:
            cols = {r[1] for r in conn.execute("PRAGMA table_info(completion_quality)").fetchall()}
            required = {"analysis_index", "block_id", "pearson_r", "n_imputed", "n_missing"}
            missing_cols = required - cols
            if missing_cols:
                errors.append(
                    "completion_quality table is missing columns: "
                    f"{', '.join(sorted(missing_cols))}"
                )


_TOP_HIT_SAMPLE_SIZE = 1000  # cross-validate this many sampled entries against the CSR


def _validate_ragged_top_hits(
    store_path: Path,
    data_root: Any,
    errors: list[str],
) -> None:
    """Validate the top-hit index groups in a ragged zarr store.

    Structural checks (lengths, sort order) are exhaustive.
    CSR cross-validation is sampled to keep validation O(1) for large stores.
    """
    top = data_root["top_hits"]
    csr = RaggedCSRReader(store_path)
    offsets = csr._offsets[:]
    vi_all = csr._variant_index[:].astype(np.int32)
    z_all = csr._z[:].astype(np.float32)

    for key in top:
        group = top[key]
        threshold = float(group.attrs.get("threshold", 0))
        vis = group["variant_index"][:].astype(np.int32)
        ais = group["analysis_index"][:].astype(np.int32)
        zs = group["z"][:].astype(np.float32)
        abs_zs = group["abs_z"][:].astype(np.float32)

        if not (len(vis) == len(ais) == len(zs) == len(abs_zs)):
            errors.append(f"top-hit index {key} has inconsistent array lengths")
            continue

        # Exhaustive: check sort order is non-increasing |z|
        if len(abs_zs) > 1 and np.any(np.diff(abs_zs) > 1e-4):
            errors.append(f"top-hit index {key} is not ranked by descending |z|")

        # Exhaustive: all abs_z must match |z|
        if not np.allclose(abs_zs, np.abs(zs), rtol=1e-3, atol=1e-3):
            errors.append(f"top-hit index {key} abs_z inconsistent with z")

        # Sampled cross-validation against the CSR
        n = len(vis)
        if n == 0:
            continue
        rng = np.random.default_rng(0)
        sample = rng.choice(n, size=min(_TOP_HIT_SAMPLE_SIZE, n), replace=False)
        for idx in sample.tolist():
            ai = int(ais[idx])
            vi = int(vis[idx])
            start = int(offsets[ai])
            end = int(offsets[ai + 1])
            pos = start + int(np.searchsorted(vi_all[start:end], vi))
            if pos >= end or int(vi_all[pos]) != vi:
                errors.append(f"top-hit index {key} references missing association")
                break
            stored_z = float(z_all[pos])
            if not np.isclose(stored_z, float(zs[idx]), rtol=1e-3, atol=1e-3):
                errors.append(f"top-hit index {key} z value inconsistent with CSR")
                break
            if p_value_from_z(stored_z) > threshold:
                errors.append(f"top-hit index {key} contains association above threshold")
                break


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


# Row-band height for streaming the dense matrix during validation, so no check
# ever materialises the full (n_variants × n_analyses) array (issue 045).
_VALIDATE_BAND_ROWS = 250_000


def _validate_dense_arrays(
    root: Any,
    n_variants: int,
    n_analyses: int,
    errors: list[str],
    imputed_arr: Any = None,
    on_panel: np.ndarray | None = None,
) -> None:
    """Stream the dense z/se (and, for a completed store, imputed) matrices in
    row-bands and run every matrix-touching content check in one pass.

    ``imputed_arr``/``on_panel`` are supplied only for Reference-Completed
    stores (from ``_validate_completion_metadata``); when present, the imputed
    checks that used to load the whole matrix are folded into this same band
    loop so peak memory stays at one band rather than the full matrices
    (issue 045).
    """
    for name in ("z", "se"):
        if name not in root:
            errors.append(f"missing data.zarr/{name}")
    if errors:
        return
    z_arr = root["z"]
    se_arr = root["se"]
    expected_shape = (n_variants, n_analyses)
    if tuple(z_arr.shape) != expected_shape:
        errors.append(f"z shape {tuple(z_arr.shape)} does not match {expected_shape}")
    if tuple(se_arr.shape) != expected_shape:
        errors.append(f"se shape {tuple(se_arr.shape)} does not match {expected_shape}")
    if errors:
        return

    # Stream in row-bands; the finite/sign checks work on float16 directly, so
    # there is no full-matrix load and no float32 upcast.
    neg_se = False
    missingness = False
    imp_not_binary = False
    imp_nan_z = False
    imp_nan_se = False
    off_panel_imputed = False
    for r0 in range(0, n_variants, _VALIDATE_BAND_ROWS):
        r1 = min(r0 + _VALIDATE_BAND_ROWS, n_variants)
        z = z_arr[r0:r1]
        se = se_arr[r0:r1]
        if not neg_se and np.any(np.isfinite(se) & (se < 0)):
            neg_se = True
        if not missingness and np.any(np.isnan(z) != np.isnan(se)):
            missingness = True
        if imputed_arr is not None:
            imp = imputed_arr[r0:r1]
            if not imp_not_binary and not np.all((imp == 0) | (imp == 1)):
                imp_not_binary = True
            imp_mask = imp == 1
            if imp_mask.any():
                if not imp_nan_z and not np.all(np.isfinite(z[imp_mask])):
                    imp_nan_z = True
                if not imp_nan_se and not np.all(np.isfinite(se[imp_mask])):
                    imp_nan_se = True
                # Off-panel rows can never be imputed (no LD structure).
                if not off_panel_imputed:
                    off_band = on_panel[r0:r1] == 0
                    if off_band.any() and np.any(imp[off_band] == 1):
                        off_panel_imputed = True
    if neg_se:
        errors.append("se contains negative finite values")
    if missingness:
        errors.append("z and se missingness is inconsistent")
    if imp_not_binary:
        errors.append("data.zarr/imputed contains values other than 0 and 1")
    if imp_nan_z:
        errors.append("imputed=1 cells have NaN z-scores")
    if imp_nan_se:
        errors.append("imputed=1 cells have NaN se values")
    if off_panel_imputed:
        errors.append(
            "off-panel (on_panel=0) rows have imputed=1 cells — off-panel is never imputable"
        )


def _validate_top_hits(root: Any, errors: list[str]) -> None:
    if "top_hits" not in root:
        return
    top = root["top_hits"]
    z_arr = root["z"]
    n_variants, n_analyses = int(z_arr.shape[0]), int(z_arr.shape[1])

    # Load each threshold group's (comparatively small) index arrays and run the
    # checks that don't need the matrix: array-length consistency, descending-|z|
    # ranking, and cell uniqueness. Keep the survivors for the streamed matrix pass.
    groups: dict[str, dict[str, Any]] = {}
    for key in top:
        group = top[key]
        threshold = float(group.attrs.get("threshold", key.replace("p_", "").replace("_", "-")))
        rows = group["variant_index"][:].astype(np.int64)
        cols = group["analysis_index"][:].astype(np.int64)
        z_values = group["z"][:].astype("float32")
        abs_z = group["abs_z"][:].astype("float32")
        if len(rows) != len(cols) or len(rows) != len(z_values) or len(rows) != len(abs_z):
            errors.append(f"top-hit index {key} has inconsistent array lengths")
            continue
        if len(abs_z) > 1 and np.any(abs_z[:-1] < abs_z[1:]):
            errors.append(f"top-hit index {key} is not ranked by descending significance")
            continue
        if len(rows):
            flat = rows * n_analyses + cols
            if len(np.unique(flat)) != len(flat):
                errors.append(f"top-hit index {key} does not match stored z values")
                continue
        groups[key] = {
            "z_crit": z_critical(threshold),
            "rows": rows,
            "cols": cols,
            "z_values": z_values,
            "n": len(rows),
            "pass_count": 0,
            "consistent": True,
        }
    if not groups:
        return

    # Single streamed pass over the matrix in row-bands: per band, count cells that
    # pass each threshold (completeness) and verify every index cell in that band
    # (matrix z finite, matches the stored index z, and itself clears the cutoff).
    for r0 in range(0, n_variants, _VALIDATE_BAND_ROWS):
        r1 = min(r0 + _VALIDATE_BAND_ROWS, n_variants)
        z_band = z_arr[r0:r1].astype("float32")
        abs_band = np.abs(z_band)
        for g in groups.values():
            g["pass_count"] += int(np.count_nonzero(abs_band >= g["z_crit"]))
            in_band = (g["rows"] >= r0) & (g["rows"] < r1)
            if not np.any(in_band):
                continue
            gathered = z_band[g["rows"][in_band] - r0, g["cols"][in_band]]
            matches = np.isclose(gathered, g["z_values"][in_band], rtol=1e-3, atol=1e-3)
            if not (
                np.all(np.isfinite(gathered))
                and np.all(matches)
                and np.all(np.abs(gathered) >= g["z_crit"])
            ):
                g["consistent"] = False

    for key, g in groups.items():
        if not g["consistent"]:
            errors.append(f"top-hit index {key} contains z value inconsistent with z array")
        elif g["n"] != g["pass_count"]:
            errors.append(f"top-hit index {key} does not match stored z values")


def default_top_hit_key(threshold: float = 5e-8) -> str:
    return threshold_key(threshold)
