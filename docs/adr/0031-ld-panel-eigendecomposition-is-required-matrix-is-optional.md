# LD Reference Panel loading requires an eigendecomposition; the LD matrix becomes optional

`opengwasdb/completion/ld_panel.py` obtains eigenvectors for a block through
one entry point, `load_ld_eigenvectors()`, which reads the per-block
`.ldeig.npz` eigendecomposition cache and falls back to deriving it from the
full `.unphased.vcor1.gz` LD matrix only when that cache is absent. The
imputation kernel (`opengwasdb/completion/impute.py`) consumes eigenvectors
alone; eigenvalues are used only to count components. The LD matrix is
therefore dead weight on the consumption path — yet `load_block()` has always
treated it as mandatory, returning `None` for any block whose matrix file is
missing.

That mismatch became load-bearing once `opengwasdb-stores` began generating
AFR, EAS and SAS LD Reference Panels (ADR 0020 there) from a cohort roughly
two orders of magnitude smaller than the UK Biobank cohort behind the
production EUR panel. Those panels are eigendecomposition-only by design — no
matrix is ever written. Under the matrix-mandatory rule, every block of every
one of those panels would silently vanish from completion with no error.

A related silent failure sits in the same loader: the number of components
needed to reach a requested cumulative-variance threshold is computed from the
full stored eigenvalue spectrum, then capped at however many eigenvectors are
actually stored, with no signal when the cap binds. The historical fixed
component count (250) was calibrated on the ~50,000-individual EUR cohort;
backfilling the 12 EUR blocks that had never been decomposed showed that even
at that sample size, reaching 99% cumulative variance took 250–516
components. Panels built from cohorts roughly 100x smaller have flatter
eigenvalue spectra and will need more components still, so an unreported cap
binds hardest on exactly the panels whose imputation quality is already most
fragile.

## Decision

An LD Reference Panel block's required storage contract is a **variant table
plus an eigendecomposition** — eigenvalues and retained eigenvectors. The LD
matrix becomes an **optional, legacy artifact**: it may still be used to
derive a missing decomposition when no `.ldeig.npz` cache exists (preserving
the production EUR panel's ability to load without a full backfill), but its
absence no longer disqualifies a block from loading. A block is admitted when
it can produce eigenvectors by either route, and skipped only when it can do
so by neither.

When the stored eigendecomposition cannot satisfy a requested truncation
threshold, that shortfall is reported (logged, naming the block, the
requested threshold, and the achieved cumulative variance) rather than
absorbed — completion proceeds with the available components instead of
aborting, since one marginal block failing a genome-wide run is a worse
failure mode than a recorded shortfall.

Identifier normalisation that already existed but lived duplicated and
private inside the Dense completion layout — handling both the legacy
`chr:pos_ref_alt` convention and canonical `chr:pos:a1:a2`, canonicalising
allele orientation — moves into the panel module as `canonical_panel_alid()`,
the single implementation both completion layouts (and any future consumer)
share. The panel module's separate `match_variants()` helper, which compared
raw identifiers without normalisation and matched nothing against any real
panel, is deleted rather than repaired; it was already dead code (imported by
the Ragged completion layout but never called).

Governs `docs/spec/store-format.md` §13.1. Implements `opengwasdb#10`.

## Considered options

- **Require both a matrix and an eigendecomposition unconditionally.** Rejected:
  this is the status quo, and it makes the eigendecomposition-only panels
  `opengwasdb-stores` already generates unusable without a pointless
  from-scratch matrix rebuild that the imputation path never reads anyway.
- **Drop matrix support entirely, require only the eigendecomposition.**
  Rejected: the production EUR panel predates the eigendecomposition-only
  storage contract and, absent a full backfill, some of its blocks may only
  ever carry a matrix. Keeping the matrix as a fallback costs nothing at the
  loader and avoids a forced migration.
- **Raise on truncation shortfall instead of reporting it.** Rejected: a
  single under-resolved block would abort an entire genome-wide completion
  run over a shortfall that is, at worst, a quality degradation for that one
  block — worse than continuing with a visible warning.

## Consequences

- `LDBlock.ld_path` changes from `Path` to `Path | None`; any code reading it
  directly (rather than through `load_ld_eigenvectors()`/`load_ld_matrix()`)
  must handle the `None` case.
- `load_ld_matrix()` and the internal matrix-fallback path raise a clear
  `RuntimeError` when called on a block with no matrix, rather than failing
  opaquely inside a file read.
- Existing EUR-panel behaviour is unchanged: blocks with a matrix and no
  cached decomposition still derive one via `ld_pca()`; blocks with both
  still prefer the cached decomposition, as before.
- A panel generator may now ship eigendecomposition-only blocks with no
  operational cost at the consumer; reclaiming disk space by deleting the
  EUR panel's now-unnecessary matrices remains a separate, unforced decision.
