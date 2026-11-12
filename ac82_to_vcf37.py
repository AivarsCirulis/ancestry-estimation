#!/usr/bin/env python
"""
Convert the raw UKB AC82 allele-frequency tables (GRCh37, 50 gzipped tab-
separated parts) into a single sites-only GRCh37 VCF, so that Picard
LiftoverVcf can move it to GRCh38 the same way the AADR pipeline did.

Confirmed layout of each AC82_AF_split_partN.tab.gz (verified against real
data on 2026-08-27):
    col 1: CHROM   (no 'chr' prefix in the raw file, e.g. "1" -- but this
                    script writes it out as "chr1" etc., see below)
    col 2: POS
    col 3: ID      (rsID, or "." if none)
    col 4: REF
    col 5: ALT
    col 6: QUAL    (always ".")
    col 7: FILTER  (always ".")
    col 8: INFO    (internal coverage/count string, NOT used here)
    col 9..90: 82 population alt-allele-frequency columns (float)

Only part 1 has a header row; parts 2-50 do not.

Output: a plain-text (not bgzipped) sites-only VCF on stdout, with each of
the 82 populations packed as its own INFO tag (POP_<name>=<freq>). REF/ALT
are carried through unchanged -- they are needed for Picard's liftover
validation and swap-detection, even though aeon's own allele-freq file
format doesn't include REF/ALT columns.

Note: this does NOT require a local GRCh37 fasta/dict. The ##contig header
lines are cosmetic (Picard's LiftoverVcf validates the *target* (GRCh38)
reference, not the source), so standard, well-known GRCh37 primary-assembly
contig lengths are hardcoded below rather than depending on a multi-GB file
you may not actually have.

IMPORTANT: the UCSC hg19ToHg38.over.chain.gz chain file (the one already
used successfully for the AADR liftover) references chromosomes as "chr1",
"chr2", ..., "chrX", "chrY", "chrM" -- NOT the bare "1", "2", ... "MT" that
AC82's raw files use. Picard matches chain blocks by exact contig name, so
this script rewrites every CHROM to "chr"-prefixed form (and "MT" -> "chrM",
matching the same rename your AADR pipeline's own step [6] already does)
before writing the VCF. Without this, every single variant fails to
liftover with 0% success and no error -- which is exactly what happened the
first time this ran without the rename.

Usage:
    python ac82_to_vcf37.py <AC82_DIR> > ac82.hg19.vcf
"""
import sys
import gzip
import os

FIXED_N = 8  # CHROM POS ID REF ALT QUAL FILTER INFO

# Standard GRCh37/b37 primary-assembly contig lengths, already in the
# "chr"-prefixed form the UCSC chain file expects (note chrM, not chrMT).
# Cosmetic header only -- liftover correctness does not depend on lengths
# being exact, but the NAMES have to match the chain file exactly.
GRCH37_CONTIGS = [
    ("chr1", 249250621), ("chr2", 243199373), ("chr3", 198022430), ("chr4", 191154276),
    ("chr5", 180915260), ("chr6", 171115067), ("chr7", 159138663), ("chr8", 146364022),
    ("chr9", 141213431), ("chr10", 135534747), ("chr11", 135006516), ("chr12", 133851895),
    ("chr13", 115169878), ("chr14", 107349540), ("chr15", 102531392), ("chr16", 90354753),
    ("chr17", 81195210), ("chr18", 78077248), ("chr19", 59128983), ("chr20", 63025520),
    ("chr21", 48129895), ("chr22", 51304566), ("chrX", 155270560), ("chrY", 59373566),
    ("chrM", 16569),
]


def to_chr_name(raw_chrom):
    """AC82 raw chromosome ("1", "X", "MT", ...) -> UCSC chain-file naming."""
    if raw_chrom in ("MT", "M", "mt", "m"):
        return "chrM"
    if raw_chrom.startswith("chr"):
        return raw_chrom
    return f"chr{raw_chrom}"


def iter_ac82_rows(ac82_dir, n_parts=50):
    populations = None
    for i in range(1, n_parts + 1):
        path = os.path.join(ac82_dir, f"AC82_AF_split_part{i}.tab.gz")
        if not os.path.exists(path):
            raise SystemExit(f"FATAL: missing expected part file: {path}")
        with gzip.open(path, "rt") as fh:
            first_line = fh.readline()
            if i == 1:
                header = first_line.rstrip("\n").split("\t")
                if len(header) < FIXED_N + 1:
                    raise SystemExit(
                        f"FATAL: header in part1 only has {len(header)} columns; "
                        f"expected at least {FIXED_N + 1}. Layout assumption is wrong -- stop and check."
                    )
                populations = header[FIXED_N:]
                sys.stderr.write(f"[INFO] Detected {len(populations)} population columns from part1 header.\n")
            else:
                # parts 2..50 have NO header -- first_line is real data, don't discard it
                fields = first_line.rstrip("\n").split("\t")
                if len(fields) != FIXED_N + len(populations):
                    raise SystemExit(
                        f"FATAL: part{i} line 1 has {len(fields)} fields, "
                        f"expected {FIXED_N + len(populations)}. Aborting -- schema mismatch."
                    )
                yield fields
            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != FIXED_N + len(populations):
                    raise SystemExit(
                        f"FATAL: part{i} has a line with {len(fields)} fields, "
                        f"expected {FIXED_N + len(populations)}. Aborting -- schema mismatch mid-file."
                    )
                yield fields
    return populations


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    ac82_dir = sys.argv[1]

    # First pass: just grab the population names from part1's header without
    # reading the whole 25M-row dataset twice -- cheap, header line only.
    part1 = os.path.join(ac82_dir, "AC82_AF_split_part1.tab.gz")
    with gzip.open(part1, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    populations = header[FIXED_N:]
    n_pop = len(populations)
    sys.stderr.write(f"[INFO] {n_pop} populations detected.\n")

    out = sys.stdout
    out.write("##fileformat=VCFv4.2\n")
    out.write('##source=ac82_to_vcf37.py (UKB AC82, Hu et al. 2023, GRCh37)\n')
    for name, length in GRCH37_CONTIGS:
        out.write(f"##contig=<ID={name},length={length}>\n")
    for pop in populations:
        out.write(f'##INFO=<ID=POP_{pop},Number=1,Type=Float,Description="Alt allele freq in population {pop} (UKB AC82)">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    n_rows = 0
    for fields in iter_ac82_rows(ac82_dir):
        chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        chrom = to_chr_name(chrom)
        if vid in ("", "."):
            vid = f"{chrom}_{pos}_{ref}_{alt}"
        af_values = fields[FIXED_N:]
        info = ";".join(f"POP_{pop}={val}" for pop, val in zip(populations, af_values))
        out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n")
        n_rows += 1
        if n_rows % 2_000_000 == 0:
            sys.stderr.write(f"[INFO] {n_rows:,} rows written...\n")

    sys.stderr.write(f"[DONE] {n_rows:,} total rows written.\n")


if __name__ == "__main__":
    main()
