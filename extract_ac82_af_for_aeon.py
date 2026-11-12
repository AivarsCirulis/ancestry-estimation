#!/usr/bin/env python
"""
Turn the GRCh38-lifted AC82 VCF (produced by Picard LiftoverVcf from
ac82_to_vcf37.py's output) into the two files aeon needs:

  1. ac82_allele_freqs_hg38.tsv.gz
       header: CHROM START STOP VAR_ID POP_<name1> POP_<name2> ...
       (tab-delimited, START=POS-1 / STOP=POS, exactly what aeon's -a flag wants)

  2. ac82_population_labels.DRAFT.tsv
       header: Superpopulation Population
       A best-effort grouping of the 82 populations into superpopulations,
       using known HGDP/1000G-style population-name conventions. THIS IS A
       DRAFT -- rows the script couldn't confidently classify are labeled
       "REVIEW_NEEDED" and every row should be sanity-checked before use,
       since this grouping affects how aeon aggregates/reports ancestry.

Critical correctness step: Picard tags any variant it had to allele-swap
during liftover with INFO/SwappedAlleles. Because aeon's AF file format has
no REF/ALT columns, an un-corrected swapped locus would silently hand aeon
an inverted allele frequency (1 - true_AF) for all 82 populations at that
site. This script detects SwappedAlleles and flips every POP_* value
(newval = 1 - oldval) before writing it out. Loci that Picard could only
place via ReverseComplementedAlleles (pure strand flip, no swap) need no
numeric correction -- the frequency still refers to the same (now
complemented) ALT allele.

Usage:
    python extract_ac82_af_for_aeon.py <lifted_hg38.vcf(.gz)> <out_prefix>
"""
import sys
import gzip
import re


# Keyword -> superpopulation, checked against the population name (prefix
# stripped) case-insensitively. First match wins. Extend/correct freely --
# this is a draft, not ground truth.
KEYWORD_SUPERPOP = [
    (["mandenka", "lwk", "amaxhosa", "yoruba", "amhara", "anuak"], "Africa"),
    (["bedouin", "syrian", "turkishe", "egypt", "morocco", "libya"], "Middle_East_North_Africa"),
    (["balochi", "india", "srilanka", "pathan", "sindhi", "gih", "bangladesh"], "South_Asia"),
    (["dai", "cambodian", "japanese", "chd", "malaysia", "philippines"], "East_Southeast_Asia"),
    (["pima", "maya"], "Americas"),
]


def guess_superpop(pop_name):
    if pop_name.startswith("G_"):
        return "British_Isles"
    if pop_name.startswith("E_"):
        return "Europe"
    if pop_name.startswith("W_"):
        lname = pop_name[2:].lower()
        for keywords, superpop in KEYWORD_SUPERPOP:
            if any(kw in lname for kw in keywords):
                return superpop
        return "REVIEW_NEEDED"
    return "REVIEW_NEEDED"


def openvcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    vcf_path, out_prefix = sys.argv[1], sys.argv[2]

    populations = None
    n_total = 0
    n_swapped = 0
    n_written = 0

    af_out_path = f"{out_prefix}.tsv.gz"
    with openvcf(vcf_path) as vin, gzip.open(af_out_path, "wt") as vout:
        for line in vin:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue  # we build our own header once we know population order
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            n_total += 1

            info_dict = {}
            flags = set()
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_dict[k] = v
                else:
                    flags.add(kv)

            if populations is None:
                # Recover population order from the first record's INFO keys,
                # in the order Picard preserved them.
                populations = [k[len("POP_"):] for k in info_dict if k.startswith("POP_")]
                header = ["CHROM", "START", "STOP", "VAR_ID"] + populations
                vout.write("\t".join(header) + "\n")
                sys.stderr.write(f"[INFO] Recovered {len(populations)} populations from first record.\n")

            swapped = "SwappedAlleles" in flags
            if swapped:
                n_swapped += 1

            try:
                values = []
                for pop in populations:
                    raw = float(info_dict[f"POP_{pop}"])
                    values.append(str(1.0 - raw) if swapped else str(raw))
            except (KeyError, ValueError) as e:
                sys.stderr.write(f"[WARN] skipping {chrom}:{pos} ({vid}) -- {e}\n")
                continue

            start = str(int(pos) - 1)
            stop = pos
            vout.write("\t".join([chrom, start, stop, vid] + values) + "\n")
            n_written += 1

            if n_total % 2_000_000 == 0:
                sys.stderr.write(f"[INFO] {n_total:,} records processed, {n_written:,} written...\n")

    if populations is None:
        raise SystemExit(
            f"FATAL: 0 usable records found in {vcf_path} (n_total={n_total} lines read). "
            "This means the input VCF has no data rows -- almost certainly a failed liftover "
            "upstream (check the LiftoverVcf log for '% of variants were not successfully "
            "lifted over'), not a bug in this script. Nothing to extract; aborting instead of "
            "writing empty/garbage output."
        )

    labels_path = f"{out_prefix.rsplit('/', 1)[0] if '/' in out_prefix else '.'}/ac82_population_labels.DRAFT.tsv"
    with open(labels_path, "w") as lf:
        lf.write("Superpopulation\tPopulation\n")
        n_review = 0
        for pop in populations:
            sp = guess_superpop(pop)
            if sp == "REVIEW_NEEDED":
                n_review += 1
            lf.write(f"{sp}\t{pop}\n")

    sys.stderr.write(f"[DONE] {n_total:,} records read, {n_written:,} written, "
                      f"{n_swapped:,} had SwappedAlleles (corrected).\n")
    sys.stderr.write(f"[DONE] Draft population labels written to {labels_path} "
                      f"({n_review} of {len(populations)} population(s) need manual review).\n")


if __name__ == "__main__":
    main()
