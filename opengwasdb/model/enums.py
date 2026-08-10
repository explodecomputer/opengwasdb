"""Controlled vocabularies for the store contract."""

from enum import StrEnum


class PrimaryStorageLayout(StrEnum):
    DENSE = "dense"
    RAGGED = "ragged"
    HYBRID = "hybrid"


class AssociationCoverage(StrEnum):
    FULL = "full"
    CIS_AND_SIGNALS = "cis_and_signals"


class CompletionState(StrEnum):
    OBSERVED_ONLY = "observed_only"
    REFERENCE_COMPLETED = "reference_completed"


class StoredEffectScale(StrEnum):
    SD = "sd"
    LOG_OR = "log_or"
    LOG_HAZARD = "log_hazard"


class SampleSizeKind(StrEnum):
    TOTAL = "total"
    CASE_CONTROL = "case_control"
    EFFECTIVE = "effective"
    VARIANT_LEVEL = "variant_level"


class SampleSizeScope(StrEnum):
    ANALYSIS_LEVEL = "analysis_level"
    VARIANT_LEVEL = "variant_level"


# Analytical Metadata vocabularies (store-format spec §7a). Literal values match
# the registry's already-shipped analyses.tsv vocabulary (opengwasdb-stores
# docs/release-metadata-schema.md), adopted as canonical here rather than
# reworked registry-side, since opengwasdb had not yet implemented analyses.tsv
# I/O when the two vocabularies diverged.
class AncestryAssignmentMethod(StrEnum):
    AF_ASSIGNED = "af_assigned"
    SOURCE_FALLBACK = "source_fallback"
    SOURCE_TRUSTED_NO_AF = "source_trusted_no_af"
    UNASSIGNED = "unassigned"


# Priority order documented in ADR-0029.
class OriginalSdMethod(StrEnum):
    DECLARED_STANDARDISED = "declared_standardised"
    SOURCE_PROVIDED = "source_provided"
    ESTIMATED_FROM_SOURCE_MAF = "estimated_from_source_maf"
    ESTIMATED_FROM_REFERENCE_MAF = "estimated_from_reference_maf"
    ESTIMATED_FROM_BETA_DISTRIBUTION = "estimated_from_beta_distribution"
    BINARY_TRAIT = "binary_trait"
    UNAVAILABLE = "unavailable"


class EafScope(StrEnum):
    ABSENT = "absent"
    VARIANT = "variant"
    ASSOCIATION = "association"


class InfoScope(StrEnum):
    ABSENT = "absent"
    VARIANT = "variant"
    ASSOCIATION = "association"


class AssociationStatus(StrEnum):
    MISSING = "missing"
    OBSERVED = "observed"
    IMPUTED = "imputed"

