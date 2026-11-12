#!/usr/bin/env python
"""
Pivot per-population AADR allele-frequency tables (one CHROM/POS/REF/ALT/AF
TSV per population, produced by aadr_populations_to_af.sh via
`bcftools +fill-tags -t AF` + `bcftools query`) into one wide table in
aeon's -a format: CHROM START STOP VAR_ID POP_<name1> POP_<name2> ...

Enforces intersection-only: a locus is kept ONLY if every single population
has a real (non-missing, non-multiallelic) AF value there. A population can
have a missing AF at a locus if none of its samples had a called genotype
there (bcftools then reports AF as "."); rather than imputing that gap,
the locus is dropped entirely -- same policy already applied to combining
AC82 + 1000G + AADR at the top level, applied here first at the
AADR-internal level since some kept AADR populations are tiny (down to
n=1), which can make per-locus missingness large enough to matter.

Usage:
    python pivot_aadr_af.py <dir_with_pop_*.tsv files> <output_prefix>

Each input file must be named pop_<POPULATION_NAME>.tsv and contain columns
(no header): CHROM  POS  REF  ALT  AF
"""
import sys
import os
import glob
import gzip


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    in_dir, out_prefix = sys.argv[1], sys.argv[2]

    pop_files = sorted(glob.glob(os.path.join(in_dir, "pop_*.tsv")))
    if not pop_files:
        raise SystemExit(f"FATAL: no pop_*.tsv files found in {in_dir}")

    populations = []
    # locus_key -> {pop_name: af_float}
    locus_data = {}
    # locus_key -> (chrom, pos, ref, alt)  (kept once, from the first file
    # that mentions it -- REF/ALT are the same across all per-population
    # files since they all come from the same source VCF's own records)
    locus_fields = {}
    n_multiallelic_skipped = 0
    n_missing_af_events = 0

    for path in pop_files:
        base = os.path.basename(path)
        pop_name = base[len("pop_") : -len(".tsv")]
        populations.append(pop_name)
        n_lines = 0
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                n_lines += 1
                fields = line.split("\t")
                if len(fields) != 5:
                    sys.stderr.write(
                        f"[WARN] {path} line {n_lines}: expected 5 fields, got "
                        f"{len(fields)} -- skipping line: {line!r}\n"
                    )
                    continue
                chrom, pos, ref, alt, af = fields
                if "," in alt or "," in ref:
                    n_multiallelic_skipped += 1
                    continue
                key = f"{chrom}\t{pos}\t{ref}\t{alt}"
                if key not in locus_fields:
                    locus_fields[key] = (chrom, pos, ref, alt)
                    locus_data[key] = {}
                if af in (".", "", "nan", "NaN"):
                    n_missing_af_events += 1
                    continue  # leave this population absent for this locus -> dropped later
                try:
                    locus_data[key][pop_name] = float(af)
                except ValueError:
                    n_missing_af_events += 1
                    continue
        sys.stderr.write(f"[INFO] read {n_lines:,} lines from {base} (population={pop_name})\n")

    n_pops = len(populations)
    sys.stderr.write(f"[INFO] {n_pops} populations, {len(locus_data):,} candidate loci before intersection.\n")

    out_path = f"{out_prefix}.tsv.gz"
    n_written = 0
    n_dropped_incomplete = 0
    with gzip.open(out_path, "wt") as out:
        header = ["CHROM", "START", "STOP", "VAR_ID"] + populations
        out.write("\t".join(header) + "\n")
        for key, af_by_pop in locus_data.items():
            if len(af_by_pop) != n_pops:
                n_dropped_incomplete += 1
                continue  # at least one population missing data here -- drop per intersection-only policy
            chrom, pos, ref, alt = locus_fields[key]
            var_id = f"{chrom}_{pos}_{ref}_{alt}"
            start = str(int(pos) - 1)
            values = [str(af_by_pop[p]) for p in populations]
            out.write("\t".join([chrom, start, pos, var_id] + values) + "\n")
            n_written += 1

    sys.stderr.write(
        f"[DONE] {n_written:,} loci written (every population has real data), "
        f"{n_dropped_incomplete:,} dropped (incomplete across populations), "
        f"{n_multiallelic_skipped:,} multiallelic records skipped, "
        f"{n_missing_af_events:,} missing-AF events encountered.\n"
    )
    sys.stderr.write(f"[DONE] Output: {out_path}\n")


if __name__ == "__main__":
    main()