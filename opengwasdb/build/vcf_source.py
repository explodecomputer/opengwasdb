"""Streaming GWAS-VCF reader using bcftools subprocesses.

Follows the GWAS-VCF / GWAS-SSF spec (Lyon et al. 2021).
All orientation is normalised to canonical ALID convention: A1 = alphabetically
first allele.  Z-scores are negated when the VCF effect allele (ALT) is not A1.
"""

from __future__ import annotations

import logging
import math
import re
import shutil
import subprocess
from collections.abc import Iterator
from pathlib import Path

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.variants.normalise import normalise_chromosome

log = logging.getLogger(__name__)

_STUDY_TYPE_RE = re.compile(r"StudyType=([^,>\s]+)")


def _require_bcftools() -> str:
    path = shutil.which("bcftools")
    if path is None:
        raise RuntimeError(
            "bcftools not found in PATH — install via conda: conda install -c bioconda bcftools"
        )
    return path


def _infer_study_type(path: str, header: str) -> StoredEffectScale:
    match = _STUDY_TYPE_RE.search(header)
    if match is None:
        raise ValueError(f"StudyType not found in ##SAMPLE header of {path}")
    study_type = match.group(1).strip()
    if study_type == "CaseControl":
        return StoredEffectScale.LOG_OR
    if study_type == "Continuous":
        return StoredEffectScale.SD_UNITS
    raise ValueError(f"Unrecognised StudyType {study_type!r} in {path}")


def read_vcf_study_type(path: str | Path) -> StoredEffectScale:
    """Return the StoredEffectScale inferred from a GWAS-VCF header.

    Reads only the header via bcftools view -h (fast).
    Raises ValueError if StudyType is absent or unrecognised.
    """
    bcftools = _require_bcftools()
    result = subprocess.run(
        [bcftools, "view", "-h", str(path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return _infer_study_type(str(path), result.stdout)


def stream_vcf_variants(path: str | Path) -> Iterator[tuple[str, int, str, str]]:
    """Yield (bare_chrom, pos, ref, alt) for every biallelic record.

    Multi-allelic records (comma in ALT) are skipped silently.
    CHROM is normalised to bare form (no chr prefix).
    """
    bcftools = _require_bcftools()
    proc = subprocess.Popen(
        [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        for line in proc.stdout:  # type: ignore[union-attr]
            line = line.rstrip("\n")
            if not line:
                continue
            chrom_raw, pos_str, ref, alt = line.split("\t")
            if "," in alt:
                continue
            yield normalise_chromosome(chrom_raw), int(pos_str), ref, alt
    finally:
        proc.stdout.close()  # type: ignore[union-attr]
        proc.wait()


def stream_vcf_associations(
    path: str | Path,
) -> Iterator[tuple[str, int, str, str, float, float, StoredEffectScale]]:
    """Yield (bare_chrom, pos, ref, alt, z, se, stored_effect_scale) for each biallelic record.

    z is oriented to canonical ALID convention: A1 = min(ref, alt).  When the
    VCF effect allele (ALT) is not A1, z is negated.  SE is always positive.

    Records with SE ≤ 0, non-finite z, or all EZ/ES/SE missing are skipped.
    """
    bcftools = _require_bcftools()
    stored_effect_scale = read_vcf_study_type(path)
    proc = subprocess.Popen(
        [
            bcftools,
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t[%EZ]\t[%ES]\t[%SE]\n",
            str(path),
        ],
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        for line in proc.stdout:  # type: ignore[union-attr]
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 7:
                continue
            chrom_raw, pos_str, ref, alt, ez_str, es_str, se_str = parts
            if "," in alt:
                continue

            se = _parse_float(se_str)
            if se is None or se <= 0:
                continue

            z = _derive_z(ez_str, es_str, se)
            if z is None or not math.isfinite(z):
                continue

            if alt > ref:
                z = -z

            yield normalise_chromosome(chrom_raw), int(pos_str), ref, alt, z, se, stored_effect_scale
    finally:
        proc.stdout.close()  # type: ignore[union-attr]
        proc.wait()


def _parse_float(s: str) -> float | None:
    if s in {".", ""}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _derive_z(ez_str: str, es_str: str, se: float) -> float | None:
    """Prefer EZ; fall back to ES/SE."""
    ez = _parse_float(ez_str)
    if ez is not None and math.isfinite(ez):
        return ez
    es = _parse_float(es_str)
    if es is not None and se > 0:
        return es / se
    return None
