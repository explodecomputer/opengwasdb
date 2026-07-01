# ukb-b Dense Store

A Dense Observed-Only Store built from all 2514 UKB-b traits in the OpenGWAS IGD
collection.

## Source data

| Property | Value |
|---|---|
| Source directory | `/local-scratch/data/opengwas/igd/ukb-b-*/` |
| VCF files | `ukb-b-XXXXX.vcf.gz` (bgzipped, tabix-indexed) |
| Input assembly | GRCh37 / hg19 |
| Output assembly | GRCh38 / hg38 (liftover applied inline at build time) |
| Traits | 2514 |
| Variants (expected) | ~9.8M (genome-wide UKB imputed array) |

## Output

| Property | Value |
|---|---|
| Store path | `/local-scratch/data/opengwas/opengwasdb/ukb-b.opengwasdb` |
| Store ID | `ukb-b` |
| Release ID | `dense-observed-hg38-v1` |
| Format | Dense Observed-Only (zarr + tabix variant axis) |

## Build

The manifest (`manifest.tsv`) was generated from the JSON metadata files in each
trait directory.  To regenerate it:

```bash
python3 - <<'EOF'
import json
from pathlib import Path

igd = Path('/local-scratch/data/opengwas/igd')
dirs = sorted(d for d in igd.iterdir() if d.name.startswith('ukb-b'))
rows = ['trait_id\tfile_path\ttrait_name\tn']
for d in dirs:
    trait_id = d.name
    vcf = d / f'{trait_id}.vcf.gz'
    if not vcf.exists():
        continue
    meta = json.loads((d / f'{trait_id}.json').read_text())
    rows.append(f'{trait_id}\t{vcf}\t{meta.get("trait", trait_id)}\t{meta.get("sample_size", 0) or 0}')
Path('data/ukb-b/manifest.tsv').write_text('\n'.join(rows) + '\n')
print(f'Wrote {len(rows)-1} rows')
EOF
```

To run the full build:

```bash
bash data/ukb-b/build.sh
```

Expected build time: ~17 hours single-threaded (same two-pass pipeline as the
100-trait chr1 benchmark, scaled to ~9.8M variants × 2514 analyses).

Expected store size: ~45 GB (extrapolated from 1774 MB at 1000 analyses × chr1,
scaled by 2514/1000 analyses × 12.9 chromosomes, with expected compression
improvement at larger chunk fill).

## Query performance (projected from chr1 × 1000 benchmark)

| Pattern | chr1 × 1000 | Full genome × 2514 (estimate) |
|---|---|---|
| regional (1 Mb window) | 130ms | ~130ms — unchanged, tabix window |
| phewas (1 variant) | 7.8ms | ~10ms — O(log n) ALID lookup |
| random_lookup (100 × 10) | 650ms | ~800ms — scattered zarr rows |
| tophits | 8.9ms | ~50–200ms — scales with hit count |
| bulk (1 analysis) | 5547ms | ~70s — reads all 9.8M rows |

Latency-sensitive queries (regional, phewas, random_lookup) are essentially
unaffected by genome size.  Bulk queries are O(n_variants) and will be slow
for full-genome; consider whether that pattern is required before optimising.

## Validation

After the build completes:

```bash
conda run -n snakemake opengwasdb validate \
  /local-scratch/data/opengwas/opengwasdb/ukb-b.opengwasdb
```
