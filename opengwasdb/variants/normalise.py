"""Canonical variant identity and association orientation."""

from __future__ import annotations

from dataclasses import dataclass

VALID_BASES = frozenset("ACGT")


class VariantNormalisationError(ValueError):
    """Raised when a source row cannot be represented as a canonical variant."""


@dataclass(frozen=True)
class CanonicalVariant:
    """Store-local canonical variant identity before variant-index assignment."""

    chromosome: str
    position: int
    effect_allele: str
    other_allele: str

    @property
    def alid(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.effect_allele}:{self.other_allele}"


@dataclass(frozen=True)
class Orientation:
    """How a source association maps to canonical ALID orientation."""

    variant: CanonicalVariant
    flipped: bool


def normalise_chromosome(chromosome: str) -> str:
    """Normalise chromosome labels without changing assembly coordinates."""

    chrom = str(chromosome).strip()
    if not chrom:
        raise VariantNormalisationError("chromosome is missing")
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper() if chrom.upper() in {"X", "Y", "MT", "M"} else chrom


def normalise_allele(allele: str) -> str:
    """Upper-case and validate a simple v0.1 allele string."""

    value = str(allele).strip().upper()
    if not value:
        raise VariantNormalisationError("allele is missing")
    if any(base not in VALID_BASES for base in value):
        raise VariantNormalisationError(f"unsupported allele {value!r}")
    return value


def orient_to_canonical(
    chromosome: str,
    position: int | str,
    source_effect_allele: str,
    source_other_allele: str,
) -> Orientation:
    """Create an ALID with alphabetically first allele as canonical effect allele.

    If the source effect allele is not the canonical A1, signed statistics must
    be negated by the caller.
    """

    chrom = normalise_chromosome(chromosome)
    try:
        pos = int(position)
    except (TypeError, ValueError) as exc:
        raise VariantNormalisationError(f"invalid position {position!r}") from exc
    if pos <= 0:
        raise VariantNormalisationError(f"invalid position {position!r}")

    effect = normalise_allele(source_effect_allele)
    other = normalise_allele(source_other_allele)
    if effect == other:
        raise VariantNormalisationError("effect and other alleles are identical")

    a1, a2 = sorted((effect, other))
    variant = CanonicalVariant(chromosome=chrom, position=pos, effect_allele=a1, other_allele=a2)
    return Orientation(variant=variant, flipped=(effect != a1))


def chromosome_sort_key(chromosome: str) -> tuple[int, str]:
    """Sort chromosomes in natural human order where possible."""

    chrom = normalise_chromosome(chromosome)
    if chrom.isdigit():
        return (int(chrom), chrom)
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (special.get(chrom.upper(), 1000), chrom)
