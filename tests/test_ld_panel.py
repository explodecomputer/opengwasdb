"""Unit tests for the LD panel loader (opengwasdb#10).

A block must load when it can produce eigenvectors by either route (a stored
eigendecomposition or a raw LD matrix) and be skipped only when neither
artifact is present. A stored eigendecomposition that falls short of a
requested truncation threshold must be reported, not silently absorbed. See
ADR 0031 and store-format spec §13.1.
"""

from __future__ import annotations

import gzip
import io
import logging
from pathlib import Path

import numpy as np
import pytest

from opengwasdb.completion.ld_panel import (
    canonical_panel_alid,
    load_block,
    load_ld_eigenvectors,
    load_ld_matrix,
)

_LOGGER_NAME = "opengwasdb.completion.ld_panel"

SNPS = [
    ("1:100:A:G", 0.3, 100),
    ("1:200:A:C", 0.4, 200),
    ("1:300:C:T", 0.5, 300),
    ("1:400:A:T", 0.2, 400),
]


def _write_tsv(block_dir: Path, block_name: str, snps: list[tuple[str, float, int]]) -> None:
    block_dir.mkdir(parents=True, exist_ok=True)
    lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
    for alid, eaf, bp in snps:
        chrom, _pos, a1, a2 = alid.split(":")
        lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t{eaf}\t{bp}")
    (block_dir / f"{block_name}.tsv").write_text("\n".join(lines) + "\n")


def _write_matrix_text(block_dir: Path, block_name: str, text: str) -> None:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(text.encode())
    (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())


def _write_matrix(block_dir: Path, block_name: str, ld: np.ndarray) -> None:
    text = "\n".join("\t".join(f"{v:.6f}" for v in row) for row in ld) + "\n"
    _write_matrix_text(block_dir, block_name, text)


def _write_npz(
    block_dir: Path, block_name: str, vals: np.ndarray, vecs: np.ndarray, k: int
) -> None:
    np.savez_compressed(block_dir / f"{block_name}.ldeig.npz", values=vals, vectors=vecs[:, :k])


def _random_ld(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n))
    return a @ a.T + np.eye(n) * n * 0.1


def _eigh_desc(ld: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(ld)
    return vals[::-1], vecs[:, ::-1]


class TestLoadBlockArtifactRequirement:
    def test_npz_only_loads_and_yields_eigenvectors(self, tmp_path):
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", SNPS)
        vals, vecs = _eigh_desc(_random_ld(len(SNPS), seed=0))
        _write_npz(block_dir, "blk", vals, vecs, k=len(SNPS))

        block = load_block(block_dir / "blk.tsv")
        assert block is not None
        assert block.ld_path is None
        assert block.ldeig_npz_path is not None

        eigenvalues, eigenvectors = load_ld_eigenvectors(block, thresh=0.9)
        assert eigenvalues.shape[0] > 0
        assert eigenvectors.shape[0] == len(SNPS)

    def test_matrix_present_but_npz_preferred(self, tmp_path):
        """A block with all three artifacts loads via the npz — proven by
        corrupting the matrix file so using it would raise."""
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", SNPS)
        vals, vecs = _eigh_desc(_random_ld(len(SNPS), seed=0))
        _write_npz(block_dir, "blk", vals, vecs, k=len(SNPS))
        _write_matrix_text(block_dir, "blk", "not,a,valid,ld,matrix\n")

        block = load_block(block_dir / "blk.tsv")
        assert block is not None
        assert block.ld_path is not None
        assert block.ldeig_npz_path is not None

        eigenvalues, eigenvectors = load_ld_eigenvectors(block, thresh=0.9)
        assert eigenvectors.shape[0] == len(SNPS)

    def test_no_artifacts_returns_none(self, tmp_path):
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", SNPS)
        assert load_block(block_dir / "blk.tsv") is None

    def test_load_ld_matrix_raises_when_ld_path_none(self, tmp_path):
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", SNPS)
        vals, vecs = _eigh_desc(_random_ld(len(SNPS), seed=0))
        _write_npz(block_dir, "blk", vals, vecs, k=len(SNPS))

        block = load_block(block_dir / "blk.tsv")
        assert block is not None and block.ld_path is None

        with pytest.raises(RuntimeError, match="no LD matrix"):
            load_ld_matrix(block, [0, 1])


class TestTruncationShortfall:
    def test_underresolved_npz_logs_warning_and_returns_available(self, tmp_path, caplog):
        n = 10
        snps = [(f"1:{100 * i}:A:G", 0.3, 100 * i) for i in range(1, n + 1)]
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", snps)
        vals, vecs = _eigh_desc(_random_ld(n, seed=42))
        # Deliberately store only 1 component so a high threshold cannot be met.
        _write_npz(block_dir, "blk", vals, vecs, k=1)

        block = load_block(block_dir / "blk.tsv")
        assert block is not None

        with caplog.at_level(logging.WARNING, logger=_LOGGER_NAME):
            eigenvalues, eigenvectors = load_ld_eigenvectors(block, thresh=0.99)

        # Truncated to what's stored, not raised.
        assert eigenvectors.shape[1] == 1
        assert eigenvalues.shape[0] == 1

        assert block.block_id in caplog.text
        assert "thresh=0.990" in caplog.text

    def test_fully_resolved_npz_logs_nothing(self, tmp_path, caplog):
        vals, vecs = _eigh_desc(_random_ld(len(SNPS), seed=0))
        block_dir = tmp_path / "EUR" / "1"
        _write_tsv(block_dir, "blk", SNPS)
        _write_npz(block_dir, "blk", vals, vecs, k=len(SNPS))

        block = load_block(block_dir / "blk.tsv")
        assert block is not None

        with caplog.at_level(logging.WARNING, logger=_LOGGER_NAME):
            load_ld_eigenvectors(block, thresh=0.5)

        assert caplog.text == ""


class TestCanonicalPanelAlid:
    def test_both_conventions_resolve_to_same_alid(self):
        legacy = canonical_panel_alid("1:100_A_G")
        canonical = canonical_panel_alid("1:100:A:G")
        assert legacy == canonical == "1:100:A:G"

    def test_chr_prefix_is_stripped(self):
        assert canonical_panel_alid("chr1:100:A:G") == "1:100:A:G"

    def test_allele_orientation_is_canonicalised(self):
        forward = canonical_panel_alid("1:100:A:G")
        reversed_order = canonical_panel_alid("1:100:G:A")
        assert forward == reversed_order == "1:100:A:G"

    def test_unparseable_id_returns_none_and_logs(self, caplog):
        with caplog.at_level(logging.WARNING, logger=_LOGGER_NAME):
            result = canonical_panel_alid("not-an-id")
        assert result is None
        assert "not-an-id" in caplog.text

    def test_invalid_position_returns_none_and_logs(self, caplog):
        with caplog.at_level(logging.WARNING, logger=_LOGGER_NAME):
            result = canonical_panel_alid("1:notanumber:A:G")
        assert result is None
        assert "1:notanumber:A:G" in caplog.text
