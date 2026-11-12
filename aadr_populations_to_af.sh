#!/bin/bash
###############################################################################
# For each kept AADR population (per aadr_population_sample_map.DRAFT.tsv),
# subset the AADR hg38 VCF to that population's samples, compute the ALT
# allele frequency at every site with bcftools +fill-tags, and extract a
# simple CHROM/POS/REF/ALT/AF table -- one pop_<name>.tsv per population.
#
# This deliberately uses bcftools' own AF computation (same tool already
# trusted throughout this project) rather than hand-rolled genotype-dosage
# counting, so the ALT-allele orientation is guaranteed correct by
# construction (bcftools computes AF relative to the VCF's own ALT column,
# exactly matching the VAR_ID=chrom_pos_ref_alt convention used everywhere
# else in this pipeline).
#
# Run pivot_aadr_af.py afterward to combine all pop_*.tsv files into one
# wide aeon -a format table (aadr_allele_freqs_hg38.tsv.gz).
#
# Usage:
#   ./aadr_populations_to_af.sh <aadr.hg38.snps.sorted.vcf.gz> <sample_map.tsv> <outdir>
###############################################################################
set -euo pipefail

VCF="$1"
SAMPLE_MAP="$2"   # columns: sample  final_population  (with header)
OUTDIR="$3"

BCFTOOLS="${BCFTOOLS:-bcftools}"
mkdir -p "$OUTDIR"

log() { printf "[%(%F %T)T] %s\n" -1 "$*"; }

# --- Validate sample IDs against what's actually in the VCF header, rather
# than silently dropping/ignoring mismatches inside bcftools. A naming
# mismatch here (e.g. sample_map uses "I0001.AG" but the VCF header uses a
# different ID string) would otherwise silently produce empty/wrong
# per-population subsets.
VCF_SAMPLES="$OUTDIR/.vcf_samples.txt"
$BCFTOOLS query -l "$VCF" | sort > "$VCF_SAMPLES"

MAP_SAMPLES="$OUTDIR/.map_samples.txt"
tail -n +2 "$SAMPLE_MAP" | cut -f1 | sort > "$MAP_SAMPLES"

N_MAP=$(wc -l < "$MAP_SAMPLES")
N_MATCHED=$(comm -12 "$VCF_SAMPLES" "$MAP_SAMPLES" | wc -l)
log "Sample map has $N_MAP samples; $N_MATCHED of those are present in the VCF header."
if [[ "$N_MATCHED" -eq 0 ]]; then
  echo "FATAL: zero overlap between sample_map sample IDs and the VCF's own sample IDs." >&2
  echo "First 5 sample_map IDs: $(head -5 "$MAP_SAMPLES" | tr '\n' ' ')" >&2
  echo "First 5 VCF IDs:        $(head -5 "$VCF_SAMPLES" | tr '\n' ' ')" >&2
  echo "Fix the ID mismatch before proceeding -- do not rerun with --force-samples until you understand why they differ." >&2
  exit 1
fi
if [[ "$N_MATCHED" -lt "$N_MAP" ]]; then
  log "WARNING: $((N_MAP - N_MATCHED)) sample_map samples were NOT found in the VCF and will be silently excluded by bcftools. Investigate if this number is large."
fi

# --- One sample list + bcftools pass per final population.
POPULATIONS=$(tail -n +2 "$SAMPLE_MAP" | cut -f2 | sort -u)
N_POPS=$(echo "$POPULATIONS" | wc -l)
log "Processing $N_POPS final populations..."

i=0
for pop in $POPULATIONS; do
  i=$((i + 1))
  OUT_TSV="$OUTDIR/pop_${pop}.tsv"
  if [[ -s "$OUT_TSV" ]]; then
    log "[$i/$N_POPS] reuse $OUT_TSV"
    continue
  fi
  SAMPLE_LIST="$OUTDIR/.samples_${pop}.txt"
  awk -F'\t' -v p="$pop" 'NR>1 && $2==p {print $1}' "$SAMPLE_MAP" > "$SAMPLE_LIST"
  N_SAMPLES=$(wc -l < "$SAMPLE_LIST")
  log "[$i/$N_POPS] $pop: $N_SAMPLES samples -> $OUT_TSV"
  $BCFTOOLS view -S "$SAMPLE_LIST" "$VCF" 2>/dev/null \
    | $BCFTOOLS +fill-tags -- -t AF 2>/dev/null \
    | $BCFTOOLS query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    > "$OUT_TSV"
  rm -f "$SAMPLE_LIST"
done

log "[DONE] Per-population AF tables written to $OUTDIR/pop_*.tsv"
log "[DONE] Next: python pivot_aadr_af.py $OUTDIR <out_prefix>"