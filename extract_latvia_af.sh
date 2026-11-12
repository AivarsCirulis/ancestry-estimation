#!/bin/bash
###############################################################################
# Compute the ALT allele frequency at every site across ALL samples in
# Latvia's national reference VCF (a single population -- unlike AADR/AC82,
# there's no per-sample population map here; every sample in this VCF is
# "Latvian" by construction).
#
# Uses bcftools' own AF computation (same as every other AF-extraction step
# in this pipeline), so ALT-allele orientation is guaranteed correct by
# construction -- AF is computed relative to the VCF's own ALT column,
# matching the VAR_ID=chrom_pos_ref_alt convention used everywhere else.
#
# Note on missingness: this VCF has already had merge-induced missingness
# (sites called in one batch but absent from others during multi-batch
# merging, which really means "reference" there, not "unknown") explicitly
# set to 0/0 upstream ("mis2ref") -- so no special missing-data handling is
# needed here beyond what bcftools does by default.
#
# Run format_latvia_af.py afterward to turn this into the
# CHROM/START/STOP/VAR_ID/Latvian wide format the final merge step expects.
#
# Usage:
#   ./extract_latvia_af.sh <latvia.vcf.gz> <out.tsv>
###############################################################################
set -euo pipefail

VCF="$1"
OUT_TSV="$2"

BCFTOOLS="${BCFTOOLS:-bcftools}"

log() { printf "[%(%F %T)T] %s\n" -1 "$*"; }

if [[ ! -f "$VCF" ]]; then
  echo "FATAL: VCF not found: $VCF" >&2
  exit 1
fi

N_SAMPLES=$($BCFTOOLS query -l "$VCF" | wc -l)
log "VCF has $N_SAMPLES samples -- computing AF across all of them as the Latvian reference population."
if [[ "$N_SAMPLES" -eq 0 ]]; then
  echo "FATAL: zero samples found in $VCF (bcftools query -l returned nothing)." >&2
  exit 1
fi

if [[ -s "$OUT_TSV" ]]; then
  log "reuse $OUT_TSV (already exists and is non-empty)"
  exit 0
fi

log "Running bcftools +fill-tags -t AF | bcftools query ... (this may take a while for a WGS/large-panel VCF)"
$BCFTOOLS +fill-tags "$VCF" -- -t AF 2>/dev/null \
  | $BCFTOOLS query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
  > "$OUT_TSV"

N_LINES=$(wc -l < "$OUT_TSV")
log "[DONE] wrote $N_LINES lines to $OUT_TSV"
log "[DONE] Next: python format_latvia_af.py $OUT_TSV latvia_allele_freqs_hg38"
