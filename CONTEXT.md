# OpenGWASDB Domain Glossary

OpenGWASDB stores, validates, and queries large collections of GWAS and QTL summary statistics as standalone releases. The domain language separates source observations, physical storage layout, reference completion, and query behaviour.

## Language

**Analysis**:
A particular statistical analysis of one **Trait**, producing associations between that trait and variants. One Trait may have many Analyses that differ by cohort, ancestry, model, sample subset, or meta-analysis.

**Trait**:
A measured or derived outcome, such as disease status, LDL cholesterol, gene expression, methylation, or protein abundance. _Avoid_: Phenotype, except as a synonym in user-facing explanatory text.

**OpenGWASDB Store**:
A self-contained logical distribution unit containing one or more Analyses and everything required to interpret and query them. A Store has a stable identity and one **Primary Storage Layout**.

**Store Release**:
An immutable, self-identifying published version of a Store. A release records Store identity, release identity, format version, creation time, completion state, and provenance so downloaded or mirrored copies remain interpretable without a catalogue service.

**Format Version**:
The version of the OpenGWASDB storage contract required to interpret a Store Release. It describes representation compatibility, not the version of the biological data.

**Primary Storage Layout**:
The main physical organisation of associations within a Store. Dense and Ragged are alternative primary layouts behind the same metadata, identity, validation, and query concepts.

**Dense Layout**:
A Primary Storage Layout in which Analyses share a variant axis and associations occupy cells in a variant-by-Analysis matrix. For Reference-Completed Dense stores, the dense variant axis is the **Reference Variant Set** and observed off-panel variants are retained in **Ragged Overflow**.

**Ragged Layout**:
A Primary Storage Layout in which each Analysis has its own sequence of retained associations referencing a Store-wide variant table. In Reference-Completed Ragged stores, observed, imputed, and missing reference variants for a completed region belong to the same Analysis association sequence.

**Ragged Overflow**:
An optional ragged component of a Dense Reference-Completed Store Release that preserves observed associations outside the **Reference Variant Set**. The dense component remains the primary query grid.

**Query Component**:
The part of a multi-component Store Release from which a returned association was read, such as Dense Grid or Ragged Overflow. Off-panel observed variants inside a Ragged primary layout's completed region are ordinary observed associations, not a separate Query Component.

**Association Coverage**:
The guarantee a Store Release makes about which source associations it retains, independently of Primary Storage Layout. Full Coverage retains every usable source association after normalisation and quality control; Cis-and-Signals Coverage retains complete cis regions plus selected significant and suggestive trans associations.

**Completion State**:
Whether a Store Release contains only source-observed associations or has also been completed against a reference variant set. Observed-Only and Reference-Completed releases share query concepts but differ in whether imputed associations may be present.

**LD Reference Panel**:
The declared ancestry-specific LD resource used for reference completion. It defines the **Reference Variant Set** and provides LD information used to infer imputed associations.

**Reference Variant Set**:
The canonical variant set defined by an LD Reference Panel for a Reference-Completed release. Reference completion attempts to provide associations on this set, subject to missingness where imputation fails or is out of scope.

**Reference Panel Membership**:
Store variant metadata indicating whether a variant belongs to the Reference Variant Set for a Reference-Completed release. It distinguishes reference-panel variants from observed off-panel variants.

**Reference Completion Region**:
A genomic interval within which a Ragged Cis-and-Signals Store attempts reference completion. A completed region contains every Reference Variant Set variant inside its boundary plus observed off-panel variants in the same boundary; singleton suggestive associations are not expanded.

**Reference Completion Method**:
The release-level algorithm and parameters used to infer imputed Z and SE values from observed associations and the LD Reference Panel. It is recorded as provenance for a Reference-Completed release rather than repeated per association.

**Reference Completion Quality**:
A quality summary for imputed associations at LD-block-by-Analysis granularity. It describes confidence in a block of imputed values rather than attaching quality metadata to every association.

**Observed Association**:
An association whose statistics come from the source dataset after OpenGWASDB normalisation. Observed associations are the authoritative basis for a Store Release.

**Imputed Association**:
An association whose Z and SE were inferred during reference completion rather than reported by the source dataset. Imputed associations must remain distinguishable from observed associations.

**Association Status**:
The origin state of an association in a Reference-Completed Store Release: Missing, Observed, or Imputed. Missing means the association was not observed and reference completion did not impute it.

**Observed-Only Query**:
A query mode that excludes Imputed Associations and returns only source-observed results. Reference-Completed releases include imputed associations by default, but callers may request observed-only results.

**Top-Hit Query**:
A query that returns associations ranked by statistical significance. Dense, Ragged, Ragged Overflow, observed, and imputed results have equal priority; ranking is determined by significance, not by storage component or association status.

**Top-Hit Index**:
A layout-specific acceleration structure that supports Top-Hit Queries using the Store's shared significance thresholds and result contract. Dense and Ragged components may encode the index differently but expose the same query semantics.

**Variant Identity**:
The canonical within-Store variant key is ALID, `chr:pos:A1:A2`, where A1 is alphabetically first and A2 is the other allele after trimming and left alignment. Cross-Store identity is the pair (**Reference Assembly**, ALID).

**Store Variant Table**:
The Store-wide union of canonical variants referenced by its Analyses. Each variant occurs once in a Store Release, independently of how many Analyses report it.

**Variant Index**:
A compact Store-local reference to a row in the Store Variant Table. It has no identity or stability guarantee outside its Store Release.

**Reference Assembly**:
The genome assembly to which every variant coordinate in a Store Release refers, such as GRCh37 or GRCh38. Each Store Release declares exactly one Reference Assembly.

**Stored Effect Scale**:
The controlled scale in which an Analysis stores beta and SE: SD Units for continuous traits, Log Odds for binary traits, or Log Hazard for survival traits.

**Original Effect Scale**:
The source-reported or source-measurement scale of beta before OpenGWASDB normalisation, such as kg/m², mmol/L, log-odds, or source-specific free text. OpenGWASDB records this as provenance rather than forcing it into a strict ontology.

**Phenotype Standard Deviation**:
The Analysis-level scale factor used to convert linear continuous effects from Original Effect Scale to SD Units. Its value and provenance are part of Analysis metadata when conversion is meaningful.

**Z-Score**:
`z = beta / se`. Z is signed, carries effect direction, and is invariant to simple rescaling of beta and SE.

**Standard Error**:
The non-negative uncertainty of beta on the same Stored Effect Scale. OpenGWASDB stores SE directly as part of the canonical statistic pair.

**Sample Size**:
The source-reported number or effective number of participants contributing to an Analysis or to one of its variant associations. OpenGWASDB does not present a value inferred from effect statistics as an observed participant count.

**Sample Size Kind**:
The interpretation of a Sample Size: Participants, Case-Control, Effective, or Unknown.

**Sample Size Scope**:
Where a Sample Size applies: Analysis when one value applies throughout, Variant when values may differ between associations, or None when sample size is unknown.

**Effect Allele Frequency**:
The frequency of the effect allele associated with an association. For observed associations EAF comes from the source dataset when available; for imputed associations EAF comes from the LD Reference Panel when stored.

**EAF Scope**:
Where an Effect Allele Frequency value applies: Variant when one value is shared for a Store variant, Association when values may differ by Analysis-variant association, or Absent when EAF is not stored.

**Imputation INFO**:
The source-reported imputation quality or information score for a variant association. INFO is optional association metadata and is not required to reconstruct beta, SE, Z, or p-value.

**INFO Scope**:
Where an Imputation INFO value applies: Variant when one value is shared for a Store variant, Association when values may differ by Analysis-variant association, or Absent when INFO is not stored.

## Example dialogue

Developer: "This UK Biobank biomarker batch has 1,000 Analyses sharing a variant axis. Should it be Dense?"

Domain expert: "Yes. Build an Observed-Only Dense Store first. Later, Reference Completion can create a new release using the LD Reference Panel as the dense axis."

Developer: "A queried variant was not source-observed but was imputed during completion. Should it appear by default?"

Domain expert: "Yes. Reference-Completed queries include imputed associations by default, but the row must expose Association Status so users can request observed-only results when needed."

Developer: "For a molecular QTL store, do suggestive trans hits get expanded to reference-panel regions?"

Domain expert: "No. Complete cis regions and significant trans regions within their existing boundaries. Suggestive singleton associations remain singletons."

