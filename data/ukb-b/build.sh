#!/usr/bin/env bash
# Build the ukb-b Dense Observed-Only Store from all 2514 UKB-b GWAS-VCFs.
#
# Run from anywhere, ideally inside a persistent tmux/screen session:
#   bash data/ukb-b/build.sh
#
# Runs through this repo's own uv-managed environment (opengwasdb installed
# editable via `uv sync`) — no separate conda env required.
#
# Runtime: unknown precisely, but Pass 1 and Pass 2 are now parallelised
# across analyses with a fork-based process pool (N_WORKERS below) — each of
# the 2514 files is independent work, so this should scale close to linearly
# with worker count until I/O or the box's core count becomes the limit.
# Pass 1 (union variant set + liftover) took ~4h15m single-threaded on the
# first attempt; Pass 2 (filling the association matrix) is the expensive
# part — scaling from the chr1x1000 benchmark suggested 30-40+ hours
# single-threaded. Both passes log progress (elapsed + ETA) to build.log —
# watch it to see real throughput rather than guessing.
# Output: ~45-86 GB (see README.md for the estimate breakdown).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MANIFEST="$(dirname "$0")/manifest.tsv"
OUTPUT="/local-scratch/data/opengwas/opengwasdb/ukb-b.opengwasdb"
STORE_ID="ukb-b"
RELEASE_ID="dense-observed-hg38-v1"
LOG="$(dirname "$0")/build.log"
N_WORKERS="${N_WORKERS:-32}"

if [ -f "$LOG" ]; then
  mv "$LOG" "$LOG.$(date +%Y%m%dT%H%M%S).bak"
fi

echo "Build started: $(date)" | tee -a "$LOG"
echo "  Manifest:  $MANIFEST" | tee -a "$LOG"
echo "  Output:    $OUTPUT" | tee -a "$LOG"
echo "  N_WORKERS: $N_WORKERS" | tee -a "$LOG"

uv run --project "$REPO_ROOT" opengwasdb build-dense-vcf \
  "$MANIFEST" \
  "$OUTPUT" \
  --store-id "$STORE_ID" \
  --release-id "$RELEASE_ID" \
  --overwrite \
  --n-workers "$N_WORKERS" \
  2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Validating store..." | tee -a "$LOG"
uv run --project "$REPO_ROOT" opengwasdb validate "$OUTPUT" 2>&1 | tee -a "$LOG"

echo "Build complete: $(date)" | tee -a "$LOG"
