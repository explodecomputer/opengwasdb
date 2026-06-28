# Separate layout, coverage, and completion

OpenGWASDB has two Primary Storage Layouts: Dense matrices for Analyses sharing a variant axis and Ragged sequences for Analyses with differing variant sets. Full versus Cis-and-Signals data is modelled separately as Association Coverage, and Observed-Only versus Reference-Completed is modelled separately as Completion State.

Keeping these axes separate allows full piecemeal GWAS, dense multi-phenotype batches, and filtered molecular QTL data to share concepts without conflating physical layout with retained-association guarantees. A Dense Reference-Completed release may also include Ragged Overflow for observed variants outside the reference panel without changing the dense grid into a different primary layout.

