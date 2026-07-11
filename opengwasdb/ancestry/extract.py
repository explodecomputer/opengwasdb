"""Allele-frequency extraction from GWAS-VCF at reference sites (ADR 0028).

Pulls ``FORMAT/AF`` from a GWAS-VCF, orients each frequency to the canonical A1
(the same A1 = min(ref, alt) convention and EAF flip used by the stores), drops
strand-ambiguous palindromic variants, and keeps only sites present in the
Ancestry Reference Panel. Extraction can be restricted to the reference sites via
a bcftools ``-R`` regions file so a full genome scan is never performed per study.
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterable, KeysView
from pathlib import Path

from opengwasdb.build.vcf_source import _require_bcftools
from opengwasdb.variants.normalise import VariantNormalisationError, orient_to_canonical

# Strand-ambiguous allele pairs: unalignable without strand info, so excluded.
_PALINDROMES = frozenset({frozenset({"A", "T"}), frozenset({"C", "G"})})


def is_palindromic(ref: str, alt: str) -> bool:
    """True for strand-ambiguous A/T or C/G variants."""
    return frozenset({ref.upper(), alt.upper()}) in _PALINDROMES


def extract_af_at_sites(
    vcf_path: str | Path,
    wanted_alids: Iterable[str],
    *,
    regions_file: str | Path | None = None,
    region: str | None = None,
    liftover: object | None = None,
    exclude_palindromic: bool = True,
) -> dict[str, float]:
    """Return ``{canonical_alid: A1-oriented AF}`` for reference sites in the VCF.

    ``wanted_alids`` restricts the result to the reference panel. ``regions_file``,
    when given, is passed to ``bcftools query -R`` so only those regions are read
    (targeted, not a full scan). ``region`` is a bcftools ``-r`` region string
    (e.g. ``"1"``) for cheap chromosome subsetting. ``liftover`` is a pyliftover
    ``LiftOver`` object (see :func:`load_liftover`): when the study is on a
    different assembly than the reference (GWAS-VCF is typically GRCh37, the
    reference GRCh38), each variant is lifted before its canonical ALID is formed,
    so study and reference join in one coordinate system. Palindromic variants are
    excluded by default; records with missing/invalid AF, unliftable coordinates,
    or non-canonical alleles are skipped.
    """
    # Accept containers that already support O(1) membership (a reference index dict
    # or its keys) without rebuilding a multi-million-element set per call.
    wanted = (
        wanted_alids
        if isinstance(wanted_alids, (set, frozenset, dict, KeysView))
        else set(wanted_alids)
    )
    bcftools = _require_bcftools()
    cmd = [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n"]
    if regions_file is not None:
        cmd += ["-R", str(regions_file)]
    if region is not None:
        cmd += ["-r", region]
    cmd.append(str(vcf_path))

    out: dict[str, float] = {}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    try:
        for line in proc.stdout:  # type: ignore[union-attr]
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 5:
                continue
            chrom_raw, pos_str, ref, alt, af_str = parts
            if "," in alt:
                continue
            if exclude_palindromic and is_palindromic(ref, alt):
                continue
            af = _parse_af(af_str)
            if af is None:
                continue
            chrom, pos = chrom_raw, pos_str
            if liftover is not None:
                lifted = _lift(liftover, chrom_raw, pos_str)
                if lifted is None:
                    continue
                chrom, pos = lifted
            try:
                orientation = orient_to_canonical(chrom, pos, alt, ref)
            except VariantNormalisationError:
                continue
            alid = orientation.variant.alid
            if alid not in wanted:
                continue
            # AF is the frequency of ALT (the VCF effect allele). ``flipped`` means
            # ALT is not the canonical A1, so orient the frequency to A1.
            out[alid] = (1.0 - af) if orientation.flipped else af
    finally:
        proc.stdout.close()  # type: ignore[union-attr]
        proc.wait()
    return out


def load_liftover(
    from_build: str = "hg19", to_build: str = "hg38", chain_file: str | None = None
) -> object:
    """Load a pyliftover ``LiftOver`` for orienting study AF onto the reference build."""
    from pyliftover import LiftOver  # type: ignore[import-untyped]

    return LiftOver(chain_file) if chain_file else LiftOver(from_build, to_build)


def _lift(liftover: object, chrom: str, pos_str: str) -> tuple[str, str] | None:
    """Lift ``(chrom, 1-based pos)`` to the target build; None if unmapped/non-standard."""
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    try:
        mapped = liftover.convert_coordinate(f"chr{bare}", int(pos_str) - 1)  # type: ignore[attr-defined]
    except (ValueError, KeyError):
        return None
    if not mapped:
        return None
    new_chrom = mapped[0][0]
    new_chrom = new_chrom[3:] if new_chrom.lower().startswith("chr") else new_chrom
    return new_chrom, str(int(mapped[0][1]) + 1)


def write_regions_file(alids: Iterable[str], path: str | Path) -> Path:
    """Write a bcftools ``-R`` regions file (``CHROM\tPOS``) from canonical ALIDs.

    Positions are de-duplicated and sorted so ``bcftools query -R`` reads each
    site once in coordinate order.
    """
    seen: set[tuple[str, int]] = set()
    for alid in alids:
        chrom, pos_str, _a1, _a2 = alid.split(":", 3)
        seen.add((chrom, int(pos_str)))

    def _key(cp: tuple[str, int]) -> tuple[int, str, int]:
        chrom, pos = cp
        return (int(chrom) if chrom.isdigit() else 1000, chrom, pos)

    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        for chrom, pos in sorted(seen, key=_key):
            fh.write(f"{chrom}\t{pos}\n")
    return path


def _parse_af(value: str) -> float | None:
    if value in {".", ""}:
        return None
    try:
        af = float(value)
    except ValueError:
        return None
    if not (0.0 <= af <= 1.0):
        return None
    return af
