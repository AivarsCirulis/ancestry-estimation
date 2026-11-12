#!/usr/bin/env python
"""
Inner-join the new "Latvian" population (Latvia's ~2100-sample national
reference VCF, computed by extract_latvia_af.sh + format_latvia_af.py) onto
the already-resolved 231-population reference
(merged_ac82_1000g_aadr_final_allele_freqs.tsv.gz / _population_labels.tsv),
producing the actual final reference for aeon.

This is a national reference population (mixed ethnic composition, as
Latvia's population itself is), the same kind of "national aggregate" as
the AC82-derived "British" population already in this reference -- so it's
added as a single new column/label row, filed under Superpopulation
"Europe", same pattern as everything else.

Like every earlier merge step in this pipeline, this is an INNER join on
VAR_ID: a locus only survives if every population -- including the new
Latvian one -- has a real value there. This will shrink the locus count
from the current 36,196 down to however many of those loci are also
covered by Latvia's reference VCF.

Usage:
    python add_latvian_population.py \\
        --af merged_ac82_1000g_aadr_final_allele_freqs.tsv.gz \\
        --labels merged_ac82_1000g_aadr_final_population_labels.tsv \\
        --latvia-af latvia_allele_freqs_hg38.tsv.gz \\
        --out-prefix merged_ac82_1000g_aadr_latvia
"""
import argparse
import sys
import pandas as pd

NON_POP_COLS = {"CHROM", "START", "STOP", "VAR_ID"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--af", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--latvia-af", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    af = pd.read_table(args.af)
    labels = pd.read_table(args.labels)
    latvia = pd.read_table(args.latvia_af)

    if set(labels.columns) != {"Superpopulation", "Population"}:
        raise SystemExit(f"FATAL: {args.labels} columns are {list(labels.columns)}, expected Superpopulation, Population")
    if list(latvia.columns) != ["CHROM", "START", "STOP", "VAR_ID", "Latvian"]:
        raise SystemExit(f"FATAL: {args.latvia_af} columns are {list(latvia.columns)}, expected CHROM START STOP VAR_ID Latvian")
    if "Latvian" in af.columns:
        raise SystemExit(f"FATAL: {args.af} already has a 'Latvian' column -- refusing to overwrite, check you're not re-running on an already-merged file.")
    if "Latvian" in set(labels["Population"]):
        raise SystemExit(f"FATAL: {args.labels} already has a 'Latvian' population row.")

    n_before = len(af)
    n_pops_before = len([c for c in af.columns if c not in NON_POP_COLS])

    merged = af.merge(latvia, on=["CHROM", "START", "STOP", "VAR_ID"], how="inner")

    n_after = len(merged)
    sys.stderr.write(
        f"[INFO] Loci before adding Latvian: {n_before:,}. Latvia's reference VCF contributed "
        f"{len(latvia):,} usable loci. After intersection: {n_after:,} "
        f"({n_after / n_before * 100:.1f}% of the pre-Latvian loci retained).\n"
    )
    if n_after == 0:
        raise SystemExit(
            "FATAL: zero loci survived intersecting the Latvian reference with the existing merged "
            "reference. This almost certainly means a chromosome-naming or build mismatch (e.g. "
            "GRCh37 vs GRCh38, or 'chr1' vs '1') -- check format_latvia_af.py's normalization step "
            "and the VAR_ID samples from both files before re-running."
        )

    labels = pd.concat(
        [labels, pd.DataFrame([{"Superpopulation": "Europe", "Population": "Latvian"}])],
        ignore_index=True,
    )

    af_out_path = f"{args.out_prefix}_allele_freqs.tsv.gz"
    labels_out_path = f"{args.out_prefix}_population_labels.tsv"

    pop_cols_final = [c for c in merged.columns if c not in NON_POP_COLS]
    merged = merged[["CHROM", "START", "STOP", "VAR_ID"] + pop_cols_final]

    af_pop_set_final = set(pop_cols_final)
    label_pop_set_final = set(labels["Population"])
    if af_pop_set_final != label_pop_set_final:
        raise SystemExit(
            f"FATAL: post-merge AF table / labels mismatch. "
            f"In AF but not labels: {af_pop_set_final - label_pop_set_final}. "
            f"In labels but not AF: {label_pop_set_final - af_pop_set_final}."
        )

    merged.to_csv(af_out_path, sep="\t", index=False, compression="gzip")
    labels.to_csv(labels_out_path, sep="\t", index=False)

    sys.stderr.write(
        f"[DONE] Final reference: {n_pops_before + 1} populations "
        f"({n_pops_before} previous + Latvian), {n_after:,} loci.\n"
    )
    sys.stderr.write(f"[DONE] Wrote {af_out_path}\n")
    sys.stderr.write(f"[DONE] Wrote {labels_out_path}\n")


if __name__ == "__main__":
    main()
