#!/bin/bash
###############################################################################
# Split the full Latvia cohort VCF (already subset down to the reference's
# ~36,188 loci) into N smaller per-batch VCFs by sample, so aeon can be run
# on each batch separately instead of all 2191 samples at once.
#
# Why: aeon's genotype-extraction code (genotype_from_vcf.py) builds a
# Python object (list-of-tuples) structure sized n_samples x n_loci and
# does 3 non-vectorized pandas .apply() passes PER SAMPLE -- this scales
# very badly with sample count and is almost certainly why the full
# 2191-sample run gets OOM-killed even against the already loci-subsetted
# VCF. Since aeon's ancestry estimate for one sample doesn't depend on
# which OTHER samples are in the same query VCF, splitting into batches
# and merging results afterward is mathematically identical to one big run.
#
# Usage:
#   ./split_samples_into_batches.sh <subset.vcf.gz> <batch_size> <outdir>
#
# Produces <outdir>/batch_001.vcf.gz, batch_002.vcf.gz, ... plus an index
# file <outdir>/batch_sample_lists/batch_NNN.samples.txt for each.
###############################################################################
set -euo pipefail

VCF="$1"
BATCH_SIZE="$2"
OUTDIR="$3"

BCFTOOLS="${BCFTOOLS:-bcftools}"
mkdir -p "$OUTDIR" "$OUTDIR/batch_sample_lists"

log() { printf "[%(%F %T)T] %s\n" -1 "$*"; }

ALL_SAMPLES="$OUTDIR/.all_samples.txt"
$BCFTOOLS query -l "$VCF" > "$ALL_SAMPLES"
N_TOTAL=$(wc -l < "$ALL_SAMPLES")
N_BATCHES=$(( (N_TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))
log "Total samples: $N_TOTAL. Batch size: $BATCH_SIZE. Number of batches: $N_BATCHES."

split -l "$BATCH_SIZE" -d -a 3 "$ALL_SAMPLES" "$OUTDIR/batch_sample_lists/batch_"
# rename split's numeric suffixes (batch_000, batch_001, ...) to a
# consistent 1-indexed, .samples.txt-suffixed naming.
i=1
for f in "$OUTDIR"/batch_sample_lists/batch_*; do
    [[ "$f" == *.samples.txt ]] && continue
    padded=$(printf "%03d" "$i")
    mv "$f" "$OUTDIR/batch_sample_lists/batch_${padded}.samples.txt"
    i=$((i + 1))
done

i=0
for sample_list in "$OUTDIR"/batch_sample_lists/batch_*.samples.txt; do
    i=$((i + 1))
    padded=$(printf "%03d" "$i")
    OUT_VCF="$OUTDIR/batch_${padded}.vcf.gz"
    N_SAMPLES=$(wc -l < "$sample_list")
    if [[ -s "$OUT_VCF" ]]; then
        log "[$i/$N_BATCHES] reuse $OUT_VCF ($N_SAMPLES samples)"
        continue
    fi
    log "[$i/$N_BATCHES] $sample_list ($N_SAMPLES samples) -> $OUT_VCF"
    $BCFTOOLS view -S "$sample_list" -Oz -o "$OUT_VCF" "$VCF"
    $BCFTOOLS index -t "$OUT_VCF"
done

log "[DONE] $N_BATCHES batch VCFs written to $OUTDIR/batch_*.vcf.gz"
