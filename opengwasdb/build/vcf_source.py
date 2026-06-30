"""Streaming GWAS-VCF reader for the two-pass dense build pipeline.

Follows the GWAS-VCF / GWAS-SSF spec (Lyon et al. 2021).
All orientation is normalised to canonical ALID convention: A1 = alphabetically
first allele.  Z-scores are negated when the VCF effect allele (ALT) is not A1.
"""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Iterator
from pathlib import Path

from opengwasdb.model.enums import StoredEffectScale
from opengwasdb.variants.normalise import normalise_chromosome

log = logging.getLogger(__name__)

_STUDY_TYPE_RE = re.compile(r"StudyType=([^,>\s]+)")


def _get_cyvcf2():
    try:
        import cyvcf2  # type: ignore[import-untyped]

        return cyvcf2
    except ImportError as exc:
        raise ImportError(
            "cyvcf2 is required for GWAS-VCF reading: "
            "conda install -c bioconda cyvcf2"
        ) from exc


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

    Opens only the header (fast).  Raises ValueError if StudyType is absent or
    unrecognised.
    """
    cyvcf2 = _get_cyvcf2()
    vcf = cyvcf2.VCF(str(path))
    header = vcf.raw_header
    vcf.close()
    return _infer_study_type(str(path), header)


def stream_vcf_variants(path: str | Path) -> Iterator[tuple[str, int, str, str]]:
    """Yield (bare_chrom, pos, ref, alt) for every biallelic record.

    Multi-allelic records (more than one ALT) are skipped silently.
    CHROM is normalised to bare form (no chr prefix).
    """
    cyvcf2 = _get_cyvcf2()
    vcf = cyvcf2.VCF(str(path))
    try:
        for rec in vcf:
            if len(rec.ALT) != 1:
                continue
            chrom = normalise_chromosome(rec.CHROM)
            yield chrom, rec.POS, rec.REF, rec.ALT[0]
    finally:
        vcf.close()


def stream_vcf_associations(
    path: str | Path,
) -> Iterator[tuple[str, int, str, str, float, float, StoredEffectScale]]:
    """Yield (bare_chrom, pos, ref, alt, z, se, stored_effect_scale) for each biallelic record.

    z is oriented to canonical ALID convention: A1 = min(ref, alt).  When the
    VCF effect allele (ALT) is not A1, z is negated.  SE is always taken
    positive from FORMAT/SE and is never sign-flipped.

    Records with SE ≤ 0, non-finite z, or EZ/ES/SE fields absent are skipped.
    """
    cyvcf2 = _get_cyvcf2()
    vcf_path = str(path)
    vcf = cyvcf2.VCF(vcf_path)
    try:
        stored_effect_scale = _infer_study_type(vcf_path, vcf.raw_header)
        for rec in vcf:
            if len(rec.ALT) != 1:
                continue
            ref = rec.REF
            alt = rec.ALT[0]
            chrom = normalise_chromosome(rec.CHROM)

            se = _extract_se(rec)
            if se is None or se <= 0:
                continue

            z = _extract_z(rec, se)
            if z is None or not math.isfinite(z):
                continue

            # Orient z to canonical ALID convention (A1 = alphabetically first)
            if alt > ref:
                z = -z

            yield chrom, rec.POS, ref, alt, z, se, stored_effect_scale
    finally:
        vcf.close()


def _extract_z(rec, se: float) -> float | None:
    """Prefer FORMAT/EZ; fall back to FORMAT/ES / FORMAT/SE."""
    try:
        ez = rec.format("EZ")
        if ez is not None:
            v = float(ez[0][0])
            if math.isfinite(v):
                return v
    except (TypeError, IndexError, ValueError):
        pass
    try:
        es = rec.format("ES")
        if es is not None and se > 0:
            return float(es[0][0]) / se
    except (TypeError, IndexError, ValueError):
        pass
    return None


def _extract_se(rec) -> float | None:
    try:
        se_field = rec.format("SE")
        if se_field is not None:
            v = float(se_field[0][0])
            if v > 0:
                return v
    except (TypeError, IndexError, ValueError):
        pass
    return None
