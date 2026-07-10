"""Allele-frequency extraction from GWAS-VCF at reference sites (ADR 0028).

Pulls ``FORMAT/AF`` from a GWAS-VCF, orients each frequency to the canonical A1
(the same A1 = min(ref, alt) convention and EAF flip used by the stores), drops
strand-ambiguous palindromic variants, and keeps only sites present in the
Ancestry Reference Panel. Extraction can be restricted to the reference sites via
a bcftools ``-R`` regions file so a full genome scan is never performed per study.
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterable
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
    exclude_palindromic: bool = True,
) -> dict[str, float]:
    """Return ``{canonical_alid: A1-oriented AF}`` for reference sites in the VCF.

    ``wanted_alids`` restricts the result to the reference panel. ``regions_file``,
    when given, is passed to ``bcftools query -R`` so only those regions are read
    (targeted, not a full scan). Palindromic variants are excluded by default;
    records with missing/invalid AF or non-canonical alleles are skipped.
    """
    wanted = wanted_alids if isinstance(wanted_alids, (set, frozenset)) else set(wanted_alids)
    bcftools = _require_bcftools()
    cmd = [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n"]
    if regions_file is not None:
        cmd += ["-R", str(regions_file)]
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
            try:
                orientation = orient_to_canonical(chrom_raw, pos_str, alt, ref)
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
