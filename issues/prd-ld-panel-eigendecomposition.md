# Consume LD Reference Panels as eigendecompositions

## Problem Statement

Reference completion can only use an LD Reference Panel that ships full LD
matrices, even though it never reads them.

Tracing the path: both completion layouts obtain eigenvectors through a single
panel-module entry point, which returns from the per-block eigendecomposition
cache and only falls back to reading the LD matrix when that cache is missing.
Downstream, the imputation kernel uses eigenvalues solely to count components —
the imputation itself consumes eigenvectors alone. The LD matrix is dead weight in
the consumption path.

Despite that, the panel loader treats the matrix as mandatory: a block whose LD
matrix file is absent is silently skipped, returning nothing. So a panel built as
eigendecompositions — which is what the panels being generated for AFR, EAS and
SAS will be, and what ADR 0031 establishes as the storage contract — cannot be
loaded at all. Every one of its blocks would vanish from completion with no error.

Two further defects sit in the same module:

A block's stored eigenvector count silently bounds the requested truncation. The
loader computes how many components are needed to reach the requested cumulative
variance, then caps that at however many are stored, and returns. If the cap binds,
imputation proceeds against less variance than it asked for and nothing reports it.
The historical fixed count of 250 was calibrated on a cohort of ~50,000
individuals, whose eigenvalue spectrum is sharply concentrated. Backfilling the 12
EUR blocks that lacked decompositions showed that reaching 99% variance took
250–516 components even at that sample size. Panels built from ~10² fewer
individuals have flatter spectra and will need more, so the cap will bind hardest
on the panels whose imputation quality is already most fragile.

Finally, the module exports a variant-matching helper that does not work against
any real panel. It compares panel SNP identifiers verbatim against store ALIDs
after stripping only a `chr` prefix, so it returns zero matches against the
production EUR panel, whose identifiers use a different separator convention.
Dense completion works only because it carries its own identifier normalisation
that handles both conventions. The broken helper is imported by the ragged
completion path but never called — so it is currently inert, and is exactly the
kind of apparent entry point someone wiring up the new panels would reach for.

## Solution

Make the eigendecomposition the artifact the panel loader requires, and the LD
matrix an optional extra.

A block loads when it has a variant table and an eigendecomposition. A block that
additionally carries an LD matrix still loads, and the matrix may still be used to
derive a missing decomposition — but its absence is no longer disqualifying. This
makes eigendecomposition-only panels loadable, which is the precondition for the
AFR/EAS/SAS panels being generated in the registry.

When a block's stored components cannot satisfy the requested truncation, that is
surfaced rather than absorbed, so a panel that under-resolves a block is
discoverable instead of quietly degrading imputation.

The broken matcher is removed, and identifier normalisation — which already exists
and works, but lives in one layout — becomes the panel module's responsibility, so
every consumer shares one correct implementation rather than one correct and one
broken.

Governed by ADR 0031 and store-format spec §13.1.

## User Stories

1. As a build engineer completing an AFR store, I want the panel loader to accept eigendecomposition-only panels, so that panels built under the current storage contract can be used at all.
2. As a build engineer, I want a block missing its LD matrix to load normally when it has an eigendecomposition, so that panel storage decisions do not silently remove blocks from completion.
3. As a build engineer, I want a panel that mixes blocks with and without LD matrices to load uniformly, so that partially migrated panels behave predictably.
4. As a build engineer, I want the existing EUR panel to keep working unchanged, so that this change carries no regression for the ancestry already in production.
5. As a build engineer, I want to be told when a block's stored components cannot reach the requested variance, so that I learn about under-resolved blocks at completion time rather than inferring them from poor results.
6. As a build engineer, I want that signal to identify the block, the variance requested, and the variance actually achieved, so that I can decide whether to regenerate the panel or accept the shortfall.
7. As a build engineer, I want under-resolution to be reported without aborting the run, so that one marginal block does not fail a genome-wide completion.
8. As an analyst, I want imputation quality to reflect the data rather than an undisclosed truncation, so that completion quality figures mean what they appear to mean.
9. As a developer, I want one implementation of panel-identifier normalisation, so that a second consumer cannot be written against a broken one.
10. As a developer, I want the panel module to expose no matcher that fails against real panels, so that the obvious entry point is the correct entry point.
11. As a developer, I want identifier normalisation to accept both the legacy separator convention and canonical ALIDs, so that the EUR panel and the new panels are both consumable.
12. As a developer, I want allele orientation canonicalised during matching, so that a panel variant recorded in the opposite orientation cannot silently invert an imputed effect direction.
13. As a developer wiring up a new panel, I want a variant identifier that fails to parse to be visible rather than dropped, so that a format mistake does not present as reduced coverage.
14. As a developer, I want the ragged completion path to stop importing a helper it does not use, so that the dependency graph reflects reality.
15. As a developer, I want the panel module's documented layout to state which files are required and which are optional, so that panel generators know what they must produce.
16. As a maintainer, I want the panel loader's behaviour to match the store-format spec's LD representation section, so that the specification remains trustworthy.
17. As a maintainer of the registry's panel generation pipeline, I want a clear statement of the minimum artifact set a panel must ship, so that the generator and the consumer agree.
18. As a reviewer, I want the number of components stored and the variance achieved to be readable per block, so that I can judge whether a panel is adequately resolved without re-running completion.
19. As a developer, I want panel loading tested against a fixture with no LD matrix, so that the core behaviour of this change is locked in.
20. As a developer, I want panel loading tested against a fixture whose stored components are deliberately insufficient, so that the reporting path is exercised rather than assumed.

## Implementation Decisions

**The panel module is the deep module here.** It already encapsulates panel
layout, block discovery, variant-table parsing, eigendecomposition loading and
truncation behind a narrow interface that the two completion layouts consume. This
work deepens it — identifier normalisation moves in, the required-artifact rule
changes — without widening its interface. Completion layouts should continue to
see the same entry points.

**Required artifacts become the variant table plus the eigendecomposition.** The
LD matrix moves from mandatory to optional. Block discovery is driven by the
variant table, as now; a block is admitted when it can produce eigenvectors by
either route, and skipped only when it can produce them by neither. The existing
fallback that derives a decomposition from a matrix is retained for panels that
predate the storage contract.

**Truncation shortfall is reported, not raised.** When the components needed to
reach the requested cumulative variance exceed those stored, the loader returns
what it has and emits a warning naming the block, the requested threshold, and the
achieved proportion. Completion continues. Aborting would make one marginal block
fail a genome-wide run, which is a worse failure mode than a recorded shortfall.

**Identifier normalisation moves into the panel module.** The working
implementation currently living in the dense completion layout — which accepts
both the legacy `chr:pos_ref_alt` convention and canonical `chr:pos:a1:a2`, and
canonicalises allele orientation — becomes the panel module's single
implementation, with the dense layout delegating to it. The broken matcher is
deleted rather than repaired, and the ragged layout's unused import of it removed.

**No new panel format is introduced.** The eigendecomposition file keeps its
existing structure (all eigenvalues, plus the retained eigenvectors), so panels
already on disk remain valid. This PRD changes what is *required*, not what is
*written*.

**No changes to the imputation kernel.** It already consumes eigenvectors and uses
eigenvalues only for a component count; nothing about this work alters the
statistics.

## Testing Decisions

Tests target observable behaviour of the panel module — which blocks load, what
eigenvectors come back, what is reported — rather than internal call structure.
The existing completion tests are the model: they construct a small synthetic
panel on disk from a fixture helper, then assert on completion results. That
helper currently writes a variant table and an LD matrix per block, and will need
to write eigendecompositions too, which makes it the natural place to express the
new required-artifact rule.

**Panel loading is the priority.** A block with a variant table and an
eigendecomposition but no LD matrix must load and yield usable eigenvectors — this
is the single behaviour the AFR/EAS/SAS panels depend on, and it currently fails.
A block with all three must load identically, proving no regression for the EUR
panel. A block with neither an eigendecomposition nor a matrix must be skipped, as
now.

**Truncation reporting.** A fixture whose stored eigenvector count is deliberately
too small to reach the requested variance must produce a warning identifying the
block and the achieved variance, and must still return the components it has.
Asserting on the warning matters: the entire point is converting a silent loss
into a visible one, so a test that only checks the returned array would pass
against the current broken behaviour.

**Identifier normalisation.** Both conventions resolve to the same canonical ALID;
allele orientation is canonicalised so a reversed panel entry matches its store
variant rather than being dropped; an unparseable identifier is handled
predictably. These are pure-function tests and cheap.

**Completion end to end.** The existing dense and ragged completion tests must
pass unchanged — that is the regression signal for identifier normalisation
moving modules. One additional case builds a completed store from an
eigendecomposition-only panel and asserts imputed cells appear, proving the
change works through the full path and not just at the loader.

## Out of Scope

- **Generating the AFR, EAS and SAS panels.** That is the registry's work; see its
  companion PRD. This PRD only makes such panels consumable.
- **Changing the imputation method or its quality gate.** The statistics are
  untouched.
- **Removing LD matrices from the existing EUR panel.** The matrices become
  unnecessary, not forbidden; reclaiming that disk space is a separate operational
  decision.
- **Recording panel provenance inside the store.** Which panel a release was
  completed against is already manifest provenance; surfacing panel sample size or
  variant-inclusion thresholds to query consumers is not addressed here.
- **Shrinkage or other regularised LD estimation.** Named in ADR 0020 as the escape
  hatch if low-sample-size panels prove inadequate; it would change the completion
  method and is not part of this work.

## Further Notes

The 12 EUR blocks that previously lacked eigendecompositions have been backfilled,
so the production EUR panel is complete at 1,357 blocks and safe to load once the
LD-matrix requirement is dropped. Without that backfill, this change would have
silently removed those blocks from completion — the inverse of the bug it fixes.

The backfill also produced the empirical case for variance-driven component
counts: at ~50,000 individuals those blocks needed 250–516 components to reach 99%
of eigenvalue mass. A panel built from several hundred individuals will have a
flatter spectrum, so any fixed count inherited from the EUR pipeline would bind.
That is why the shortfall must be reported rather than absorbed.
