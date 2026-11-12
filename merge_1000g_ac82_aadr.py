#!/usr/bin/env python
"""
Merge aeon's built-in 1000G reference, the AC82 reference, and the AADR
reference into one allele-freq table + population_labels file for aeon.

Applies, in order, every merge decision made for this project
(2026-09-01):

  1. AC82's 23 British/Irish "G_*" sub-populations are averaged (simple
     unweighted mean -- no per-population sample size is available for
     AC82) into one new "British" population, filed under the "Europe"
     superpopulation. The 23 originals are dropped.
  2. AC82's competing entries for 8 populations that AADR also has its own
     (better-curated) version of are dropped in favor of AADR's version:
     E_Russian, W_Balochi, W_Bedouin_2, W_Cambodians, W_Dai, W_Mandenka,
     W_Maya, W_Pima.
  3. aeon's built-in GBR (1000 Genomes "British in England and Scotland")
     is dropped entirely -- superseded by AC82's merged "British" entry,
     per instruction to drop "the Brits from the other [sources]".
  4. 1000G's remaining 25 population codes keep their own allele
     frequencies untouched, but their superpopulation labels are remapped
     from the 1000G 5-code scheme (AFR/AMR/EAS/EUR/SAS) onto the same
     7-category scheme AC82 already uses (Africa/Americas/
     East_Southeast_Asia/Europe/Middle_East_North_Africa/South_Asia),
     plus "Oceania" (used only by AADR populations), so that aeon's
     superpopulation SUMMARY rollup is coherent across all three merged
     sources instead of mixing naming conventions.
  5. AADR's populations (already fully resolved by
     aadr_classify_populations.py + aadr_populations_to_af.sh +
     pivot_aadr_af.py) are taken as-is.
  6. Final loci = INTERSECTION of all three sources' VAR_IDs (inner join).
     No imputation for populations/loci that don't overlap -- see the
     project's "intersection only" decision. This is required by aeon's
     own code anyway (aeon.py loads the whole -a file into one dense
     float32 tensor -- there is no missing-value support).

Usage:
    python merge_1000g_ac82_aadr.py \\
        --g1k-af /path/to/aeon_ancestry/refs/g1k_allele_freqs.txt \\
        --g1k-labels /path/to/aeon_ancestry/refs/pop2super.txt \\
        --ac82-af ac82_allele_freqs_hg38.tsv.gz \\
        --ac82-labels ac82_population_labels.DRAFT.tsv \\
        --aadr-af aadr_allele_freqs_hg38.tsv.gz \\
        --aadr-classification aadr_population_classification.DRAFT.csv \\
        --out-prefix merged_ac82_1000g_aadr
"""
import argparse
import sys
import pandas as pd

AC82_BRITISH_COLUMNS = [
    "G_Anglia", "G_C_England", "G_Cornwall", "G_Cornwall_Tip", "G_Cumbria",
    "G_Devon", "G_Ireland", "G_Lincolnshire", "G_Merseyside", "G_NC_England",
    "G_NE_Scot", "G_NW_Scot", "G_NW_Wales", "G_N_Yorkshire", "G_Noi",
    "G_Northumberland", "G_Orkney", "G_SC_England", "G_SE_England",
    "G_SE_Wales", "G_SW_Wales", "G_S_England", "G_S_Yorkshire",
]

# AC82 columns dropped because AADR's own version of the same population
# wins (see module docstring, decision 2). Exact spellings confirmed
# against the real ac82_population_labels.DRAFT.tsv on 2026-09-01.
AC82_COLLISION_DROPS = [
    "E_Russian", "W_Balochi", "W_Bedouin_2", "W_Cambodians", "W_Dai",
    "W_Mandenka", "W_Maya", "W_Pima",
]

# 1000G's own 5-code superpopulation scheme remapped onto AC82's 7-category
# scheme (decision 4). GBR is deliberately absent -- it's dropped entirely
# (decision 3), not remapped.
G1K_SUPERPOP_REMAP = {
    "AFR": "Africa",
    "AMR": "Americas",
    "EAS": "East_Southeast_Asia",
    "EUR": "Europe",
    "SAS": "South_Asia",
}

NON_POP_COLS = {"CHROM", "START", "STOP", "VAR_ID"}


def load_af_table(path):
    df = pd.read_table(path)
    missing = NON_POP_COLS - set(df.columns)
    if missing:
        raise SystemExit(f"FATAL: {path} is missing expected columns {missing}")
    return df


def build_ac82_block(ac82_af_path, ac82_labels_path):
    af = load_af_table(ac82_af_path)
    labels = pd.read_table(ac82_labels_path)
    if set(labels.columns) != {"Superpopulation", "Population"}:
        raise SystemExit(
            f"FATAL: {ac82_labels_path} columns are {list(labels.columns)}, "
            "expected exactly Superpopulation, Population"
        )

    missing_british = [c for c in AC82_BRITISH_COLUMNS if c not in af.columns]
    if missing_british:
        raise SystemExit(f"FATAL: expected British sub-populations not found in AC82 AF table: {missing_british}")
    missing_collisions = [c for c in AC82_COLLISION_DROPS if c not in af.columns]
    if missing_collisions:
        raise SystemExit(f"FATAL: expected collision columns not found in AC82 AF table: {missing_collisions}")

    af = af.copy()
    af["British"] = af[AC82_BRITISH_COLUMNS].mean(axis=1)
    drop_cols = AC82_BRITISH_COLUMNS + AC82_COLLISION_DROPS
    af = af.drop(columns=drop_cols)

    labels = labels[~labels["Population"].isin(drop_cols)].copy()
    labels = pd.concat(
        [labels, pd.DataFrame([{"Superpopulation": "Europe", "Population": "British"}])],
        ignore_index=True,
    )

    kept_pop_cols = [c for c in af.columns if c not in NON_POP_COLS]
    n_expected = 82 - len(AC82_BRITISH_COLUMNS) - len(AC82_COLLISION_DROPS) + 1
    if len(kept_pop_cols) != n_expected:
        sys.stderr.write(
            f"[WARN] AC82 block has {len(kept_pop_cols)} population columns, "
            f"expected {n_expected} (82 - 23 British + 1 merged - 8 collisions). "
            "Check for unexpected duplicate/renamed columns.\n"
        )
    sys.stderr.write(f"[INFO] AC82 block: {len(kept_pop_cols)} populations, {len(af):,} loci.\n")
    return af, labels


def build_g1k_block(g1k_af_path, g1k_labels_path):
    af = load_af_table(g1k_af_path)
    labels = pd.read_table(g1k_labels_path)
    if "GBR" not in af.columns:
        raise SystemExit(f"FATAL: expected GBR column not found in {g1k_af_path}")
    af = af.drop(columns=["GBR"])
    labels = labels[labels["Population"] != "GBR"].copy()

    unmapped = set(labels["Superpopulation"]) - set(G1K_SUPERPOP_REMAP)
    if unmapped:
        raise SystemExit(f"FATAL: unexpected 1000G superpopulation code(s) not in remap table: {unmapped}")
    labels["Superpopulation"] = labels["Superpopulation"].map(G1K_SUPERPOP_REMAP)

    kept_pop_cols = [c for c in af.columns if c not in NON_POP_COLS]
    sys.stderr.write(f"[INFO] 1000G block: {len(kept_pop_cols)} populations (GBR dropped), {len(af):,} loci.\n")
    return af, labels


def build_aadr_block(aadr_af_path, aadr_classification_path):
    af = load_af_table(aadr_af_path)
    cls = pd.read_csv(aadr_classification_path)
    kept = cls[cls["decision"] == "KEEP"][["final_population", "superpopulation"]].copy()
    kept = kept.rename(columns={"final_population": "Population", "superpopulation": "Superpopulation"})

    review = kept[kept["Superpopulation"] == "REVIEW_NEEDED"]
    if len(review):
        raise SystemExit(
            f"FATAL: {len(review)} AADR population(s) still have superpopulation=REVIEW_NEEDED "
            f"in {aadr_classification_path}: {list(review['Population'])}. "
            "Resolve these (edit SUPERPOP_MAP in aadr_classify_populations.py and rerun it, "
            "or hand-edit the DRAFT csv) before merging."
        )

    af_pop_cols = set(c for c in af.columns if c not in NON_POP_COLS)
    label_pop_names = set(kept["Population"])
    only_in_af = af_pop_cols - label_pop_names
    only_in_labels = label_pop_names - af_pop_cols
    if only_in_af:
        raise SystemExit(f"FATAL: AADR AF table has population columns with no label: {only_in_af}")
    if only_in_labels:
        sys.stderr.write(
            f"[WARN] {len(only_in_labels)} AADR population(s) in the classification table never "
            f"produced a column in the AF table (likely zero surviving loci for that tiny "
            f"population after intersection) -- dropping their label rows: {only_in_labels}\n"
        )
        kept = kept[kept["Population"].isin(af_pop_cols)]

    sys.stderr.write(f"[INFO] AADR block: {len(af_pop_cols)} populations, {len(af):,} loci.\n")
    return af, kept


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--g1k-af", required=True)
    p.add_argument("--g1k-labels", required=True)
    p.add_argument("--ac82-af", required=True)
    p.add_argument("--ac82-labels", required=True)
    p.add_argument("--aadr-af", required=True)
    p.add_argument("--aadr-classification", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    g1k_af, g1k_labels = build_g1k_block(args.g1k_af, args.g1k_labels)
    ac82_af, ac82_labels = build_ac82_block(args.ac82_af, args.ac82_labels)
    aadr_af, aadr_labels = build_aadr_block(args.aadr_af, args.aadr_classification)

    all_pop_names = (
        [c for c in g1k_af.columns if c not in NON_POP_COLS]
        + [c for c in ac82_af.columns if c not in NON_POP_COLS]
        + [c for c in aadr_af.columns if c not in NON_POP_COLS]
    )
    dupes = {name for name in all_pop_names if all_pop_names.count(name) > 1}
    if dupes:
        raise SystemExit(
            f"FATAL: population name collision across sources after all merge rules applied: {dupes}. "
            "This means two sources ended up with the exact same final population name -- "
            "add it to the collision-handling rules instead of letting it clash silently."
        )

    sys.stderr.write("[INFO] Merging on VAR_ID (inner join = intersection across all three sources)...\n")
    merged = g1k_af.merge(ac82_af, on=["CHROM", "START", "STOP", "VAR_ID"], how="inner")
    merged = merged.merge(aadr_af, on=["CHROM", "START", "STOP", "VAR_ID"], how="inner")

    n_final_pops = len(all_pop_names)
    sys.stderr.write(
        f"[DONE] Merged reference: {n_final_pops} populations "
        f"({len([c for c in g1k_af.columns if c not in NON_POP_COLS])} 1000G + "
        f"{len([c for c in ac82_af.columns if c not in NON_POP_COLS])} AC82 + "
        f"{len([c for c in aadr_af.columns if c not in NON_POP_COLS])} AADR), "
        f"{len(merged):,} loci after full three-way intersection.\n"
    )
    if len(merged) == 0:
        raise SystemExit(
            "FATAL: zero loci survived the three-way intersection. This almost certainly means a "
            "build/coordinate mismatch between one of the sources -- do not proceed with an empty "
            "reference. Check each source's VAR_ID sample manually before re-running."
        )

    af_out_path = f"{args.out_prefix}_allele_freqs.tsv.gz"
    merged.to_csv(af_out_path, sep="\t", index=False, compression="gzip")

    labels_out = pd.concat([g1k_labels, ac82_labels, aadr_labels], ignore_index=True)
    labels_out = labels_out[["Superpopulation", "Population"]]
    label_pop_set = set(labels_out["Population"])
    af_pop_set = set(all_pop_names)
    if label_pop_set != af_pop_set:
        raise SystemExit(
            f"FATAL: population_labels/AF-table mismatch after merge. "
            f"In labels but not AF: {label_pop_set - af_pop_set}. "
            f"In AF but not labels: {af_pop_set - label_pop_set}. "
            "aeon does an inner-merge on population name between these two files -- any "
            "mismatch here means a population would silently vanish from aeon's output."
        )
    labels_out_path = f"{args.out_prefix}_population_labels.tsv"
    labels_out.to_csv(labels_out_path, sep="\t", index=False)

    sys.stderr.write(f"[DONE] Wrote {af_out_path}\n")
    sys.stderr.write(f"[DONE] Wrote {labels_out_path}\n")
    sys.stderr.write(
        f"[DONE] Superpopulation counts:\n{labels_out['Superpopulation'].value_counts().to_string()}\n"
    )


if __name__ == "__main__":
    main()