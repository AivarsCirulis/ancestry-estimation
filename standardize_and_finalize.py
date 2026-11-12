#!/usr/bin/env python
"""
Final polish pass on the merged reference
(merged_ac82_1000g_aadr_latvia_allele_freqs.tsv.gz / _population_labels.tsv):

  1. Drop 2 more cross-source duplicates found after a closer look at the
     final population list (missed in the earlier systematic review):
       - E_CzechRepublic (AC82) / Czech (AADR) -- same population.
         Test evidence was marginal (0.02 in one sample vs 0.0), not
         decisive like Polish was, so the established default applies:
         AADR wins, drop AC82's version.
       - E_North_Italian (AC82) / Italian_North (AADR) -- same population,
         words in reversed order (missed because the original scan checked
         for exact string matches, not reordered ones). No test signal
         either way -> default applies -> AADR wins.

  2. Drop 6 AC82 labels whose original meaning is ambiguous/undocumented
     (E_DutchlikeGerman, E_PortlikeSpain, E_Italy_3, E_Spain_3,
     E_Yugoslavia_1, W_India_2) rather than keep a guessed placeholder name
     on a permanent reference file -- per user instruction 2026-09-08:
     keep well-defined specific regional populations (Italian_North,
     Italian_Sardinian, Bergamo, Tuscan, etc.) as their own entries
     alongside the generic ones (Italian, Spanish), but drop the
     ambiguous/undocumented numbered sub-groups entirely.

  3. Rename every surviving population and superpopulation to a clean,
     human-readable name: aeon's own ae.csv output echoes these labels
     directly, so standardizing them here (rather than post-processing
     every future output) fixes it at the source. Covers: 1000G's 26
     cryptic 3-letter codes -> their standard published names, AC82's
     E_/W_/G_-prefixed CamelCase/underscore labels -> spaced Title Case,
     and superpopulation category names -> readable form (e.g.
     East_Southeast_Asia -> "East & Southeast Asia").

Usage:
    python standardize_and_finalize.py \\
        --af merged_ac82_1000g_aadr_latvia_allele_freqs.tsv.gz \\
        --labels merged_ac82_1000g_aadr_latvia_population_labels.tsv \\
        --out-prefix aeon_reference_final
"""
import argparse
import sys
import pandas as pd

NON_POP_COLS = {"CHROM", "START", "STOP", "VAR_ID"}

DROP = {
    "E_CzechRepublic",   # duplicate of AADR's Czech; AADR wins (default, marginal evidence)
    "E_North_Italian",   # duplicate of AADR's Italian_North (reversed word order); AADR wins (default)
    "E_DutchlikeGerman", # ambiguous/undocumented AC82 label
    "E_PortlikeSpain",   # ambiguous/undocumented AC82 label
    "E_Italy_3",         # ambiguous/undocumented AC82 numbering
    "E_Spain_3",         # ambiguous/undocumented AC82 numbering
    "E_Yugoslavia_1",    # ambiguous/undocumented AC82 numbering
    "W_India_2",         # ambiguous/undocumented AC82 numbering
}

# 1000 Genomes Phase 3 population codes -> standard published names.
G1K_NAMES = {
    "ACB": "African Caribbean (Barbados)", "ASW": "African Ancestry (SW USA)",
    "BEB": "Bengali (Bangladesh)", "CDX": "Chinese Dai (Xishuangbanna)",
    "CEU": "Northern & Western European (Utah)", "CHB": "Han Chinese (Beijing)",
    "CHS": "Southern Han Chinese", "CLM": "Colombian (Medellin)",
    "ESN": "Esan (Nigeria)", "FIN": "Finnish", "GIH": "Gujarati Indian (Houston)",
    "GWD": "Gambian (Western Division)", "IBS": "Iberian (Spain)",
    "ITU": "Indian Telugu (UK)", "JPT": "Japanese (Tokyo)",
    "KHV": "Kinh (Ho Chi Minh City)", "LWK": "Luhya (Webuye, Kenya)",
    "MSL": "Mende (Sierra Leone)", "MXL": "Mexican Ancestry (LA)",
    "PEL": "Peruvian (Lima)", "PJL": "Punjabi (Lahore)",
    "PUR": "Puerto Rican", "STU": "Sri Lankan Tamil (UK)",
    "TSI": "Tuscan (Italy)", "YRI": "Yoruba (Ibadan, Nigeria)",
}

SUPERPOP_RENAME = {
    "Africa": "Africa",
    "Americas": "Americas",
    "East_Southeast_Asia": "East & Southeast Asia",
    "Europe": "Europe",
    "Middle_East_North_Africa": "Middle East & North Africa",
    "South_Asia": "South Asia",
    "Oceania": "Oceania",
}

# Explicit overrides confirmed with the user (2026-09-08).
OVERRIDES = {
    "E_French_Basque": "French Basque",
    "E_Switzerland_1": "Swiss (group 1)",
    "E_sicilianw": "Sicilian",
    "E_croatian": "Croatian",
    "E_cypriot": "Cypriot",
    "E_german": "German",
    "E_greek": "Greek",
    "E_italians": "Italian",
    "E_ukrainian": "Ukrainian",
    "W_amaXhosa": "Xhosa",
    "W_syrian": "Syrian",
    "W_SriLanka": "Sri Lankan",
    "W_CHD_2": "Chinese (Denver, historical HapMap)",
    "Bedouin_A": "Bedouin A",
    "Bedouin_B": "Bedouin B",
    "Sindhi_Pakistan": "Sindhi (Pakistan)",
    "Italian_North": "Italian (North)",
    "Italian_Sardinian": "Italian (Sardinian)",
    "IBS_CanaryIslands": "Spanish (Canary Islands)",
    "Ju_hoan_North": "Ju/'hoan North",
    "Khomani_San": "Khomani San",
    "Russia_NorthOssetian": "North Ossetian (Russia)",
    "Kyrgyz_Kyrgyzstan": "Kyrgyz",
    "Mexico_Zapotec": "Zapotec (Mexico)",
    "Kenya_Somali": "Somali (Kenya)",
    "Greek_1": "Greek (1)",
    "Greek_2": "Greek (2)",
    "Relli_1": "Relli (1)",
    "Relli_2": "Relli (2)",
    "Mixtec_1": "Mixtec (1)",
    "Mixtec_2": "Mixtec (2)",
    "Cree1": "Cree (1)",
    "Cree2": "Cree (2)",
    "BantuSA": "Bantu (South Africa)",
    "BantuSA_Herero": "Bantu Herero (South Africa)",
    "BantuSA_Ovambo": "Bantu Ovambo (South Africa)",
    "BantuHerero": "Bantu Herero",
    "BantuTswana": "Bantu Tswana",
    "BantuKenya": "Bantu (Kenya)",
    "China_Lahu": "Lahu (China)",
    "Eskimo_Naukan": "Naukan Eskimo",
    "Eskimo_Sireniki": "Sireniki Eskimo",
    "Eskimo_Chaplin": "Chaplin Eskimo",
}


def clean_generic(name):
    n = name
    for prefix in ("E_", "W_", "G_"):
        if n.startswith(prefix):
            n = n[len(prefix):]
    n = n.replace("_", " ")
    out = []
    for i, ch in enumerate(n):
        if i > 0 and ch.isupper() and n[i - 1].islower():
            out.append(" ")
        out.append(ch)
    n = "".join(out)
    return " ".join(w.capitalize() if not w.isupper() else w for w in n.split())


def resolve_new_name(old_name):
    if old_name in G1K_NAMES:
        return G1K_NAMES[old_name]
    if old_name in OVERRIDES:
        return OVERRIDES[old_name]
    return clean_generic(old_name)


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
            f"FATAL: input AF table / labels mismatch before standardization even starts. "
            f"In AF but not labels: {af_pop_cols - label_pop_set}. "
            f"In labels but not AF: {label_pop_set - af_pop_cols}."
        )

    missing_drops = DROP - af_pop_cols
    if missing_drops:
        sys.stderr.write(f"[WARN] expected to drop {missing_drops} but they're not in the AF table (already removed?) -- continuing.\n")
    actual_drops = DROP & af_pop_cols

    af = af.drop(columns=list(actual_drops))
    labels = labels[~labels["Population"].isin(actual_drops)]
    sys.stderr.write(f"[DROP] removed {len(actual_drops)} population(s): {sorted(actual_drops)}\n")

    old_pop_cols = [c for c in af.columns if c not in NON_POP_COLS]
    rename_map = {old: resolve_new_name(old) for old in old_pop_cols}

    # Sanity check: no two old names should map to the same new name.
    new_name_counts = {}
    for old, new in rename_map.items():
        new_name_counts.setdefault(new, []).append(old)
    collisions = {new: olds for new, olds in new_name_counts.items() if len(olds) > 1}
    if collisions:
        raise SystemExit(f"FATAL: rename would create duplicate population names: {collisions}")

    af = af.rename(columns=rename_map)
    labels["Population"] = labels["Population"].map(lambda p: rename_map.get(p, p))

    unmapped_superpops = set(labels["Superpopulation"]) - set(SUPERPOP_RENAME)
    if unmapped_superpops:
        raise SystemExit(f"FATAL: unmapped superpopulation(s) found: {unmapped_superpops}")
    labels["Superpopulation"] = labels["Superpopulation"].map(SUPERPOP_RENAME)

    af_out_path = f"{args.out_prefix}_allele_freqs.tsv.gz"
    labels_out_path = f"{args.out_prefix}_population_labels.tsv"

    pop_cols_final = [c for c in af.columns if c not in NON_POP_COLS]
    af = af[["CHROM", "START", "STOP", "VAR_ID"] + pop_cols_final]

    af_pop_set_final = set(pop_cols_final)
    label_pop_set_final = set(labels["Population"])
    if af_pop_set_final != label_pop_set_final:
        raise SystemExit(
            f"FATAL: post-standardization AF table / labels mismatch. "
            f"In AF but not labels: {af_pop_set_final - label_pop_set_final}. "
            f"In labels but not AF: {label_pop_set_final - af_pop_set_final}."
        )
    if labels["Population"].duplicated().any():
        dupes = labels.loc[labels["Population"].duplicated(), "Population"].tolist()
        raise SystemExit(f"FATAL: duplicate population name(s) after rename: {dupes}")

    af.to_csv(af_out_path, sep="\t", index=False, compression="gzip")
    labels.to_csv(labels_out_path, sep="\t", index=False)

    sys.stderr.write(f"\n[DONE] Final population count: {len(pop_cols_final)} (dropped {len(actual_drops)} from the {len(af_pop_cols) + 0} before this step).\n")
    sys.stderr.write(f"[DONE] Superpopulation counts:\n{labels['Superpopulation'].value_counts().to_string()}\n")
    sys.stderr.write(f"[DONE] Wrote {af_out_path}\n")
    sys.stderr.write(f"[DONE] Wrote {labels_out_path}\n")


if __name__ == "__main__":
    main()
