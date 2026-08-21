# Statistic array encodings: fixed-point z, residual-coded EAF, reference EAF for imputed cells

Supersedes ADR 0036's Decision 3 (`float32` for `eaf`) and refines its Decision 5
(Reference Completion). Leaves ADR 0036's semantics intact: EAF is still per
(variant, Analysis), still oriented to the stored effect allele, still declared
per Analysis by `eaf_scope`. Only the physical encoding changes — and, unlike
ADR 0036, this **is** a breaking format change.

## Context

ADR 0036 stored EAF as a `float32` plane parallel to `z`/`se`. That roughly
doubled a store's statistic bytes, which is a large price for a column that is
annotation rather than the finding. Reviewing that cost surfaced three
questions, each answered by measurement on this project's own pilot data
rather than from theory. All figures below are **compressed** bytes per stored
cell under the store's own codec (Blosc zstd level 3 + bitshuffle), because
raw size is not what a Store Release occupies.

### `float16` is the wrong encoding for `z`, and always was

`z` is bounded (|z| ≤ 48.4 across every pilot) and needs *uniform* precision,
because p-value error scales with `z · Δz`. `float16` gives the opposite:
precision degrades exactly where p-values are steepest.

| z | float16 step | stored p uncertain by |
|---|---|---|
| 5.0 | 0.0039 | 1.02× |
| 30.0 | 0.0156 | 1.60× |
| 47.8 | 0.0312 | **4.45×** |

The largest |z| observed across the FinnGen, metabolome and pQTL pilots is
48.4, and the codebase already documents a FADS1/FADS2 hit at 47.8 — so this
is not a theoretical tail, it is the top of the existing data.

| encoding | B/cell | max \|Δz\| | worst p error |
|---|---|---|---|
| float16 | 1.85 | 0.0125 | 1.8× |
| **int16, scale 1/512** | **1.47** | 0.0010 | **1.02×** |
| int8 | 0.47 | 0.25 | 14,500× |

Fixed-point int16 is smaller *and* 16× more accurate. Smaller because
bitshuffle groups the near-constant high bytes of a bounded quantity, where
float16's exponent bits churn for the small values that dominate.

### `float16` is the *right* encoding for `se`

The mirror image, and worth stating so it is not "fixed" later by analogy:
`se` spans 3.2 decades and needs relative precision, which is what a
floating-point exponent already provides.

| encoding | B/cell | max relative error |
|---|---|---|
| **float16** | **1.75** | 4.9e-4 |
| int16 log-se | 1.84 | 1.1e-4 |
| int16 uniform | 1.36 | 4.1e-1 |

### EAF is highly redundant, but only against the right baseline

Cross-analysis EAF agreement, as the log-residual against a candidate
baseline:

| store | baseline | residual sd | p99 | max |
|---|---|---|---|---|
| FinnGen, 8 endpoints (one cohort) | within-store median | **0.0121** | 0.059 | 0.63 |
| FinnGen | EUR reference panel | 0.436 | 1.82 | **8.07** |
| GWAS Catalog, 7 EUR studies | within-store median | **0.0876** | 0.154 | 4.22 |

A within-store baseline leaves residuals small enough to code in 8 bits. The
reference panel does not: a max log-residual of 8.07 is a **3000× discrepancy
in MAF**, the Finnish bottleneck showing up exactly as it should. This is the
central negative result — reference-panel EAF is not a proxy for a cohort's
own allele frequency, and must never be substituted for a missing one.

Cross-*study* residuals (GWAS Catalog, one ancestry) are 7× wider than
cross-*endpoint* residuals within one cohort, which is what separates the
Hybrid case from Dense/Ragged — not the layout itself, but the cohort
heterogeneity the layout happens to correlate with.

### The measurement found a live data defect

While computing the above, GCST003566 — one of the ten studies in the
`gwas-catalog-eur-hybrid` pilot — turned out to report
`effect_allele_frequency` against the *other* allele: r = −0.999 against every
other study and against the reference panel. Nothing in the pipeline could
have caught it, because until ADR 0036 no store retained EAF at all. See #115.

## Decision

### 1. `z` is `int16` fixed-point at scale 1/512

Range |z| ≤ 64, uniform step 0.00195, p accurate to 1.02% throughout. A build
encountering |z| > 64 fails loudly rather than clipping: silently flattening
the most significant association in a store is precisely the failure mode this
project has spent a release stage eliminating.

`se` stays `float16`. Both choices follow from the shape of the quantity, not
from a general preference for integers.

### 2. EAF is a per-variant baseline plus a per-cell `int8` log-residual

- `eaf_baseline` — one `float32` per variant, the within-store representative
  observed EAF. Amortises to nothing: 0.42 B/cell across 8 Analyses, 0.014
  across 1000.
- `eaf` — one `int8` per cell, the quantised log-ratio to that baseline.

Two codes are reserved: one for "this Analysis has no EAF here", one for "out
of range — exact value in a sparse side table". 253 levels remain for the
residual.

The builder **measures** the residual spread and selects the range, rather
than inferring it from the layout. Nothing stops a Dense manifest spanning
several cohorts, and a store that did would silently clip against a range
chosen on the assumption it did not.

| case | range | B/cell | clipped | p99 error |
|---|---|---|---|---|
| one cohort (Dense/Ragged) | ±0.5 | 0.14 | 1 in 2M | 0.18% |
| many studies (Hybrid) | ±1.0 | 0.52 | 0.11% | 0.39% |

### 3. `se` may be coded as an `int8` residual, but only where every Analysis has EAF

`log(se) ≈ a + b·log(2f(1−f))` fits real data at R² = 0.992 with b = −0.52
against a theoretical −0.5 — `opengwasdb.completion.block` already relies on
this relationship to derive imputed SE. Coding the residual costs 0.86 B/cell
instead of 1.75, at 0.12% median / 0.24% p99 round-trip error, far below the
sampling error of the SE estimate itself.

It requires EAF **in the same cell**, and a zarr array has one dtype, so a
single EAF-less Analysis forces the whole array back to `float16`. This is why
the Hybrid case differs from Dense/Ragged.

### 4. Reference-panel EAF is stored once per variant, for imputed cells only

An imputed cell's EAF *is* the panel's, identical for every Analysis imputed
at that variant, so it is a per-variant constant: `eaf_reference`, one
`float32` per variant, ~0 B/cell. Which cells it describes is already
recorded, by `association_status` (ADR 0013), exactly as spec §9 specifies.

Observed cells whose source reported no EAF stay NaN. They do **not** fall
back to the panel value — see the 3000× result above.

This supersedes ADR 0036's reason for deferring imputed EAF to #113. That
deferral assumed the value had to travel through the completion checkpoint
shards; as a per-variant constant read straight from the panel, it does not.

Constraint: one panel per completed store. True of
`complete_{dense,ragged,hybrid}_store` today, but it becomes load-bearing
here, so it is recorded in `manifest.json` rather than left implicit.

### 5. Scenarios

| | z | se | eaf | total | vs today |
|---|---|---|---|---|---|
| **A** Dense/Ragged, no EAF | int16 1.47 | f16 1.75 | — | **3.22** | −11% |
| **B** Dense/Ragged, EAF | int16 1.47 | int8 0.86 | 0.14 | **2.47** | −31% |
| **C** Hybrid, mixed | int16 1.47 | f16 1.75 | 0.52 | **3.74** | +4% |
| **D** + Reference-Completed | | | +ref ~0.00 | | |

Today's baseline is 3.60 B/cell carrying no EAF at all. Every scenario stores
strictly more information than that; two of the three are smaller.

### 6. Allele-flipped EAF is rejected at build time

Each Analysis's A1-oriented EAF is correlated against the reference panel;
`r < 0` fails the build. Measured separation is unambiguous, and population
bottlenecking does not confuse it — Finnish frequencies differ from EUR by up
to 3000× in magnitude but not in direction.

| source | r vs EUR panel |
|---|---|
| GCST003566 | −0.9992 |
| GCST005076 | +0.9996 |
| FinnGen (bottlenecked isolate) | +0.9954 |
| Metabolome | +0.9996 |

## Consequences

- **This is a breaking format change**, unlike ADR 0036. It needs a
  `format_version` bump and therefore depends on #112 settling the
  compatibility and migration policy first.
- **Existing stores must be rebuilt** to gain any of it. All four pilots need
  rebuilding regardless, since #83/#106/#109 were all build-time defects.
- **Arrays stop being independently interpretable.** An `int8` `eaf` plane is
  meaningless without `eaf_baseline`; a residual-coded `se` is meaningless
  without `eaf` and the per-Analysis fit coefficients. That is a real cost
  against a format whose stated goal is self-containment, and it is accepted
  deliberately, not chosen on byte count. It also argues for landing the three
  changes in order of increasing coupling — `z` first, `se` last.
- **Decode is no longer free.** Reads gain a vectorised transform and, for
  `se`, a second array plus per-Analysis coefficients. Negligible per query,
  but it removes the option of handing a caller a raw mmap'd slice.
- **`EafScope.VARIANT` stays reserved and unimplemented.** The measurements
  rule out collapsing per-Analysis EAF to one value per variant: no variant in
  either single-cohort pilot has identical EAF across analyses.
