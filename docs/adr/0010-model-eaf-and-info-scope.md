# Model EAF and INFO scope explicitly

Effect allele frequency and imputation INFO use explicit scope rather than assuming they belong to the canonical variant table. Dense Stores may store one value per variant when genuinely shared, or one value per association when values differ across Analyses; Ragged Stores use association scope because retained associations are Analysis-specific sequences.

EAF and INFO are optional metadata and are not used to reconstruct beta, SE, Z, or p-value. Variant-scoped EAF or INFO is allowed only when the builder can establish that one value is genuinely shared; differing values are represented at association scope or omitted, not averaged.

For imputed associations, stored EAF comes from the LD Reference Panel rather than being inferred from the source study. In v1, EAF provenance is inferred from Association Status: observed associations use source EAF and imputed associations use reference-panel EAF.

