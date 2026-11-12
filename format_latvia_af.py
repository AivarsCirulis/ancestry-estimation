#!/usr/bin/env python
"""
Turn the raw CHROM/POS/REF/ALT/AF table (produced by extract_latvia_af.sh
via `bcftools +fill-tags -t AF` + `bcftools query`) into the
CHROM/START/STOP/VAR_ID/Latvian wide format used everywhere else in this
pipeline -- ready to inner-join onto
merged_ac82_1000g_aadr_final_allele_freqs.tsv.gz in add_latvian_population.py.

Applies the same correctness rules already used for AC82/AADR:
  - drops multiallelic records (comma in REF or ALT) -- aeon's VAR_ID
    convention is one ALT allele per row.
  - restricts to the 25 canonical chromosomes aeon's genotype-extraction
    code can parse (chr1-22, chrX, chrY, chrM) -- unplaced/alt contigs have
    underscores in their own names that break aeon's "split VAR_ID on _"
    logic.
  - normalizes chromosome naming to the "chrN" convention used by the rest
    of the merged reference (this VCF's own CHROM field might already say
    "chr1" or might say bare "1" depending on how it was built -- both are
    handled; a mismatch here would silently produce zero overlap at the
    final merge step instead of an obvious error, so this is worth getting
    right rather than assuming).
  - drops missing/non-numeric AF values (shouldn't occur given this VCF's
    upstream mis2ref correction, but checked defensively anyway).

Usage:
    python format_latvia_af.py <raw_chrom_pos_ref_alt_af.tsv> <out_prefix>

Output: <out_prefix>.tsv.gz with columns CHROM START STOP VAR_ID Latvian
"""
import sys
import gzip

ALLOWED_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
# Bare (no "chr" prefix) equivalents, for normalization.
ALLOWED_BARE = {c[3:] for c in ALLOWED_CHROMS}


def normalize_chrom(chrom):
    if chrom in ALLOWED_CHROMS:
        return chrom
    if chrom in ALLOWED_BARE:
        return f"chr{chrom}"
    # A couple of common alternate spellings for the mito/sex contigs.
    if chrom in ("MT", "chrMT"):
        return "chrM"
    return None  # unrecognized / non-canonical contig


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    in_path, out_prefix = sys.argv[1], sys.argv[2]

    n_total = 0
    n_written = 0
    n_multiallelic = 0
    n_bad_chrom = 0
    n_bad_af = 0
    chroms_seen_raw = set()
    normalized_any = False

    out_path = f"{out_prefix}.tsv.gz"
    with open(in_path) as fin, gzip.open(out_path, "wt") as fout:
        fout.write("CHROM\tSTART\tSTOP\tVAR_ID\tLatvian\n")
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            n_total += 1
            fields = line.split("\t")
            if len(fields) != 5:
                sys.stderr.write(f"[WARN] line {n_total}: expected 5 fields, got {len(fields)} -- skipping: {line!r}\n")
                continue
            chrom, pos, ref, alt, af = fields
            chroms_seen_raw.add(chrom)

            if "," in ref or "," in alt:
                n_multiallelic += 1
                continue

            norm_chrom = normalize_chrom(chrom)
            if norm_chrom is None:
                n_bad_chrom += 1
                continue
            if norm_chrom != chrom:
                normalized_any = True

            if af in (".", "", "nan", "NaN"):
                n_bad_af += 1
                continue
            try:
                af_val = float(af)
            except ValueError:
                n_bad_af += 1
                continue

            start = str(int(pos) - 1)
            var_id = f"{norm_chrom}_{pos}_{ref}_{alt}"
            fout.write(f"{norm_chrom}\t{start}\t{pos}\t{var_id}\t{af_val}\n")
            n_written += 1

            if n_total % 5_000_000 == 0:
                sys.stderr.write(f"[INFO] {n_total:,} lines read, {n_written:,} written so far...\n")

    if normalized_any:
        sys.stderr.write(
            f"[INFO] chromosome names were normalized (input didn't already use the 'chrN' convention). "
            f"Raw CHROM values seen: {sorted(chroms_seen_raw)[:5]}{'...' if len(chroms_seen_raw) > 5 else ''}\n"
        )

    sys.stderr.write(
        f"[DONE] {n_total:,} input lines, {n_written:,} written, "
        f"{n_multiallelic:,} multiallelic skipped, {n_bad_chrom:,} non-canonical-contig skipped, "
        f"{n_bad_af:,} missing/invalid AF skipped.\n"
    )
    if n_written == 0:
        raise SystemExit(
            "FATAL: 0 usable records written. Check that the input file actually has data and that "
            "chromosome naming wasn't rejected entirely (see 'Raw CHROM values seen' above)."
        )
    sys.stderr.write(f"[DONE] Output: {out_path}\n")


if __name__ == "__main__":
    main()
