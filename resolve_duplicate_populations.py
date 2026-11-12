#!/usr/bin/env python
"""
Resolve the 14 cross-source duplicate populations found after the initial
3-way merge (merged_ac82_1000g_aadr_allele_freqs.tsv.gz /
_population_labels.tsv) -- populations named differently across 1000G,
AC82, and AADR that actually refer to the same real-world group (e.g.
"E_polish" (AC82) and "Polish" (AADR) are both Polish).

Decision rule (per user instruction, 2026-09-04): PICK ONE WINNER per
group and drop the rest -- no averaging. The winner is chosen by:

  1. Direct evidence from the 2-sample sanity test
     (merged_test2_ancestry_ae.csv), when the evidence is decisive
     (one source's version carries real signal while the other is ~0 for
     both samples). This overrides the general default below -- concretely,
     this is what happens for Polish: AC82's "E_polish" fits the real test
     samples (0.20 / 0.17) while AADR's "Polish" does not (0.0 / 0.0), so
     AC82 wins there despite AADR being the "default" winner in general.

  2. Where the test gives no signal at all (both/all candidates are 0.0
     for both samples -- true for every group except Polish and, weakly,
     Finnish, since the 2 test samples carry only Baltic/Lithuanian-Polish
     ancestry), fall back to a general scientific-strength default:
     1000G (whole-genome sequenced, most rigorously validated) > AADR
     (well-curated meta-compilation, but often small per-population N) >
     AC82 (largest raw N via UK Biobank, but no documented per-population
     sample size and not an independently published reference panel).

Resulting per-group decision:
    Polish    -> AC82's E_polish wins (decisive test evidence overrides default)
    Finnish   -> 1000G's FIN wins (weak lean in test data; also matches default)
    Hungarian, Bulgarian, Norwegian, French, Turkish, Pathan, Sindhi
              -> AADR wins over AC82 (no test signal; default order)
    Japanese, Yoruba, Luhya, Tuscan, Gujarati
              -> 1000G wins over AC82 and AADR (no test signal; default order)

Usage:
    python resolve_duplicate_populations.py \\
        --af merged_ac82_1000g_aadr_allele_freqs.tsv.gz \\
        --labels merged_ac82_1000g_aadr_population_labels.tsv \\
        --out-prefix merged_ac82_1000g_aadr_final
"""
import argparse
import sys
import pandas as pd

NON_POP_COLS = {"CHROM", "START", "STOP", "VAR_ID"}

# final output name -> (column to KEEP, [columns to DROP])
# Column names exactly as they appear in the current merged AF table.
RESOLUTIONS = {
    "Polish":    ("E_polish", ["Polish"]),          # AC82 wins: decisive test evidence
    "Finnish":   ("FIN", ["E_Finland", "Finnish"]),  # 1000G wins: weak test lean + default
    "Hungarian": ("Hungarian", ["E_hungarian"]),     # AADR wins: default (no signal)
    "Bulgarian": ("Bulgarian", ["E_bulgarian"]),     # AADR wins: default
    "Norwegian": ("Norwegian", ["E_norwegian"]),     # AADR wins: default
    "French":    ("French", ["E_French"]),           # AADR wins: default
    "Turkish":   ("Turkish", ["W_turkishe"]),        # AADR wins: default
    "Pathan":    ("Pathan", ["W_Pathan"]),           # AADR wins: default
    "Sindhi":    ("Sindhi_Pakistan", ["W_Sindhi"]),  # AADR wins: default
    "Japanese":  ("JPT", ["W_Japanese", "Japanese"]),  # 1000G wins: default
    "Yoruba":    ("YRI", ["W_Yoruba", "Yoruba"]),      # 1000G wins: default
    "Luhya":     ("LWK", ["W_LWK", "Luhya"]),          # 1000G wins: default
    "Tuscan":    ("TSI", ["E_TSI", "Tuscan_1"]),       # 1000G wins: default
    "Gujarati":  ("GIH", ["W_GIH_2"]),                 # 1000G wins: default
}

# Only Polish needs a rename (its winning column, "E_polish", carries an
# AC82-style prefix that no longer makes sense once it's the sole survivor
# for this population). Everything else's winning column name is already
# clean (a bare 1000G code, or AADR's own bare population name).
RENAME_WINNER = {
    "Polish": "Polish",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--af", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    af = pd.read_table(args.af)
    labels = pd.read_table(args.labels)
    if set(labels.columns) != {"Superpopulation", "Population"}:
        raise SystemExit(f"FATAL: {args.labels} columns are {list(labels.columns)}, expected Superpopulation, Population")

    af_pop_cols = set(c for c in af.columns if c not in NON_POP_COLS)
    label_pop_set = set(labels["Population"])
    if af_pop_cols != label_pop_set:
        raise SystemExit(
            f"FATAL: input AF table / labels mismatch before resolution even starts. "
            f"In AF but not labels: {af_pop_cols - label_pop_set}. "
            f"In labels but not AF: {label_pop_set - af_pop_cols}."
        )

    n_resolved = 0
    for group_name, (keep, drop) in RESOLUTIONS.items():
        missing = [c for c in [keep] + drop if c not in af.columns]
        if missing:
            sys.stderr.write(
                f"[WARN] group '{group_name}': column(s) {missing} not found in {args.af} "
                f"(already resolved/renamed?) -- skipping this group.\n"
            )
            continue
        af = af.drop(columns=drop)
        labels = labels[~labels["Population"].isin(drop)]

        new_name = RENAME_WINNER.get(group_name, keep)
        if new_name != keep:
            af = af.rename(columns={keep: new_name})
            labels.loc[labels["Population"] == keep, "Population"] = new_name

        sys.stderr.write(f"[RESOLVE] {group_name}: kept '{new_name}' (was '{keep}'), dropped {drop}.\n")
        n_resolved += 1

    af_out_path = f"{args.out_prefix}_allele_freqs.tsv.gz"
    labels_out_path = f"{args.out_prefix}_population_labels.tsv"

    pop_cols_final = [c for c in af.columns if c not in NON_POP_COLS]
    af = af[["CHROM", "START", "STOP", "VAR_ID"] + pop_cols_final]

    af_pop_set_final = set(pop_cols_final)
    label_pop_set_final = set(labels["Population"])
    if af_pop_set_final != label_pop_set_final:
        raise SystemExit(
            f"FATAL: post-resolution AF table / labels mismatch. "
            f"In AF but not labels: {af_pop_set_final - label_pop_set_final}. "
            f"In labels but not AF: {label_pop_set_final - af_pop_set_final}."
        )
    if labels["Population"].duplicated().any():
        dupes = labels.loc[labels["Population"].duplicated(), "Population"].tolist()
        raise SystemExit(f"FATAL: duplicate population name(s) remain after resolution: {dupes}")

    af.to_csv(af_out_path, sep="\t", index=False, compression="gzip")
    labels.to_csv(labels_out_path, sep="\t", index=False)

    sys.stderr.write(f"\n[DONE] {n_resolved} of {len(RESOLUTIONS)} duplicate group(s) resolved.\n")
    sys.stderr.write(f"[DONE] Final population count: {len(pop_cols_final)} (was {len(af_pop_cols)} before resolution).\n")
    sys.stderr.write(f"[DONE] Wrote {af_out_path}\n")
    sys.stderr.write(f"[DONE] Wrote {labels_out_path}\n")


if __name__ == "__main__":
    main()