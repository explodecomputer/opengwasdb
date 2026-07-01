#!/usr/bin/env bash
# Build the ukb-b Dense Observed-Only Store from all 2514 UKB-b GWAS-VCFs.
#
# Run from the repo root:
#   bash data/ukb-b/build.sh
#
# Requires the snakemake conda environment with opengwasdb installed.
# Expected runtime: ~17 hours.  Output: ~45 GB.

set -euo pipefail

MANIFEST="$(dirname "$0")/manifest.tsv"
OUTPUT="/local-scratch/data/opengwas/opengwasdb/ukb-b.opengwasdb"
STORE_ID="ukb-b"
RELEASE_ID="dense-observed-hg38-v1"
LOG="$(dirname "$0")/build.log"

echo "Build started: $(date)" | tee -a "$LOG"
echo "  Manifest: $MANIFEST" | tee -a "$LOG"
echo "  Output:   $OUTPUT" | tee -a "$LOG"

conda run -n snakemake opengwasdb build-dense-vcf \
  "$MANIFEST" \
  "$OUTPUT" \
  --store-id "$STORE_ID" \
  --release-id "$RELEASE_ID" \
  2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Validating store..." | tee -a "$LOG"
conda run -n snakemake opengwasdb validate "$OUTPUT" 2>&1 | tee -a "$LOG"

echo "Build complete: $(date)" | tee -a "$LOG"
