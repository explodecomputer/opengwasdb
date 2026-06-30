"""Batch hg19 → hg38 liftover utility for the dense build pipeline.

Converts a collection of (bare_chrom, pos, ref, alt) tuples in hg19
coordinates to hg38 ALIDs in a single pass using one LiftOver object.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path

log = logging.getLogger(__name__)

_BUILD_ALIASES: dict[str, str] = {
    "hg19": "hg19", "grch37": "hg19", "b37": "hg19", "37": "hg19",
    "hg38": "hg38", "grch38": "hg38", "b38": "hg38", "38": "hg38",
}


class LiftoverFailureError(RuntimeError):
    """Raised when the liftover failure rate exceeds the configured threshold."""


def _normalise_build(build: str) -> str:
    canonical = _BUILD_ALIASES.get(build.lower().strip())
    if canonical is None:
        raise ValueError(
            f"Unknown genome build {build!r}. "
            "Use hg19/hg38 or aliases GRCh37/GRCh38, b37/b38, 37/38."
        )
    return canonical


def build_liftover_lookup(
    variants: Iterable[tuple[str, int, str, str]],
    from_build: str = "hg19",
    to_build: str = "hg38",
    failure_threshold: float = 0.01,
    chain_file: str | Path | None = None,
) -> dict[tuple[str, int, str, str], str]:
    """Convert hg19 (bare_chrom, pos, ref, alt) tuples to hg38 ALIDs.

    Returns a dict mapping each input tuple to its hg38 ALID
    (``chr:pos:A1:A2`` where A1 = alphabetically first allele).  Variants
    that fail liftover are omitted; a warning is logged.

    When the failure rate exceeds ``failure_threshold``, ``LiftoverFailureError``
    is raised before any output is returned.

    A single ``LiftOver`` object is created for the entire batch (not per-variant).

    Parameters
    ----------
    variants:
        Iterable of ``(bare_chrom, pos, ref, alt)`` tuples.  Both bare (``1``)
        and chr-prefixed (``chr1``) chrom forms are accepted.
    from_build / to_build:
        Genome build identifiers.  Ignored when ``chain_file`` is provided.
    failure_threshold:
        Maximum allowed fraction of failures before raising.  Default 0.01 (1%).
    chain_file:
        Path to a chain file.  When provided, ``LiftOver`` is initialised from
        the file rather than by downloading by build name.
    """
    try:
        from pyliftover import LiftOver  # type: ignore[import-untyped]
    except ImportError as exc:
        raise ImportError(
            "pyliftover is required for coordinate liftover: pip install pyliftover"
        ) from exc

    if chain_file is not None:
        lo = LiftOver(str(chain_file))
    else:
        lo = LiftOver(_normalise_build(from_build), _normalise_build(to_build))

    variant_list = list(variants)
    total = len(variant_list)
    if total == 0:
        return {}

    result: dict[tuple[str, int, str, str], str] = {}
    n_fail = 0

    for chrom, pos, ref, alt in variant_list:
        bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
        mapped = lo.convert_coordinate(f"chr{bare}", pos - 1)  # 1-based → 0-based
        if not mapped:
            n_fail += 1
            continue
        new_chrom_full: str = mapped[0][0]
        new_pos = int(mapped[0][1]) + 1  # 0-based → 1-based
        new_bare = new_chrom_full[3:] if new_chrom_full.lower().startswith("chr") else new_chrom_full  # noqa: E501
        a1 = min(ref, alt)
        a2 = max(ref, alt)
        result[(chrom, pos, ref, alt)] = f"{new_bare}:{new_pos}:{a1}:{a2}"

    if n_fail:
        rate = n_fail / total
        log.warning(
            "Liftover %s→%s: %d/%d variants failed (%.1f%%)",
            from_build, to_build, n_fail, total, rate * 100,
        )
        if rate > failure_threshold:
            raise LiftoverFailureError(
                f"Liftover failure rate {rate:.1%} ({n_fail}/{total}) exceeds "
                f"threshold {failure_threshold:.1%}"
            )

    return result
