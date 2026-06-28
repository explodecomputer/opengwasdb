"""Controlled vocabularies for the store contract."""

from enum import StrEnum


class PrimaryStorageLayout(StrEnum):
    DENSE = "dense"
    RAGGED = "ragged"


class AssociationCoverage(StrEnum):
    FULL = "full"
    CIS_AND_SIGNALS = "cis_and_signals"


class CompletionState(StrEnum):
    OBSERVED_ONLY = "observed_only"
    REFERENCE_COMPLETED = "reference_completed"


class StoredEffectScale(StrEnum):
    SD_UNITS = "sd_units"
    LOG_OR = "log_or"
    LOG_HAZARD = "log_hazard"


class SampleSizeKind(StrEnum):
    PARTICIPANTS = "participants"
    CASE_CONTROL = "case_control"
    EFFECTIVE = "effective"
    UNKNOWN = "unknown"


class SampleSizeScope(StrEnum):
    ANALYSIS = "analysis"
    VARIANT = "variant"
    NONE = "none"


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

