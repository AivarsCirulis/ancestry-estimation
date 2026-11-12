#!/usr/bin/env python
"""
Classify every raw AADR population label (from the .ind-derived sample sheet)
into KEEP (a legitimate present-day population usable as an aeon ancestry
reference) or EXCLUDE, and assign a final population name + superpopulation
to every kept one.

Why this exists: AADR has no machine-readable "ancient vs present-day" field
in what's available on this cluster (no .anno file with a Date-BP column --
see README note below). The only structural signal available is AADR's own
suffix convention: ".DG" ("diploid genotype") is reserved almost exclusively
for high-coverage genomes, which in practice means present-day reference
individuals folded into AADR as modern comparanda; ancient/archaeological
samples are overwhelmingly ".AG"/".SG"/".IM" (low-coverage, pseudo-haploid).
This is NOT a perfect rule -- a manual pass over the actual .DG label list
found real exceptions (archaic hominins, and several genuinely ancient
samples that still carry a .DG suffix, e.g. explicit "_EBA"/"_LBK"/"_N"
culture-period tags, explicit "_BP" dates, well-known ancient individuals
like Ust'-Ishim). Those exceptions are hardcoded below as explicit
exclusions found by manual review of the real AADR label list on
2026-09-01, rather than caught by a generic pattern.

DRAFT status: this is exactly the same kind of draft-needs-review artifact
as ac82_population_labels.DRAFT.tsv was -- review the output CSV
(aadr_population_classification.DRAFT.csv) before trusting it, especially
any row with decision=REVIEW_NEEDED or an unfamiliar population name this
script had to guess a superpopulation for.

Usage:
    python aadr_classify_populations.py <aadr_samples.ind.tsv> <out_prefix>

Input format expected (tab-delimited, with header): sample  sex  population
(this is projected/aadr_samples.ind.tsv, or equivalently derived directly
from the eigenstrat .ind file's 3 whitespace-delimited columns).

Outputs:
  <out_prefix>_classification.DRAFT.csv
      One row per RAW population label (not per sample), columns:
      raw_label, n_samples, decision (KEEP/EXCLUDE), final_population,
      superpopulation, reason
  <out_prefix>_sample_map.DRAFT.tsv
      One row per KEPT sample, columns: sample, final_population
      (this is what the AF-extraction script actually consumes)
"""
import sys
import re
import csv
from collections import defaultdict, Counter

# --- The 26 1000 Genomes population codes aeon already ships natively in
# refs/g1k_allele_freqs.txt. AADR also carries these exact codes (its own
# copy of the 1000G individuals) with a .DG suffix -- e.g. "GBR.DG". Those
# are excluded here in favor of aeon's own built-in curated entry, not
# recomputed from AADR's copy (same underlying population, no reason to
# duplicate/potentially-degrade it).
G1K_CODES = {
    "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN",
    "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL",
    "PEL", "PJL", "PUR", "STU", "TSI", "YRI",
}

# Archaic hominins -- never appropriate as a "modern human ancestry" proxy,
# regardless of suffix.
ARCHAIC_SUBSTRINGS = ("neanderthal", "denisova")

# Manually identified during review of the real AADR .DG label list
# (2026-09-01): these carry a .DG suffix (so they'd otherwise pass the
# "DG = modern" heuristic) but are genuinely ancient/archaeological --
# identified by explicit period/culture/site names or BP dates. Exact
# label match, case-sensitive, as they appear in the AADR .ind file.
KNOWN_ANCIENT_DESPITE_DG = {
    "Bahamas_SouthAndros_Ceramic.DG",
    "Bulgaria_Dzhulyunitsa_EBA.DG",
    "Cameroon_ShumLaka_SMA.DG",
    "Ethiopia_4500BP.DG",
    "Germany_EN_LBK.DG",
    "Hungary_EN_Koros_oMesolithic.DG",
    "Hungary_MN_LBK.DG",
    "Iran_TepeHissar_C.DG",
    "Luxembourg_Mesolithic.DG",
    "Russia_Afanasievo.DG",
    "Russia_Afanasievo_son.I3388_son.I3950_brother.I3949.DG",
    "Russia_Afanasievo_son.I3388_son.I3950_brother.I6714.DG",
    "Russia_Irkutsk_N.DG",
    "Russia_Tyumen_HG.DG",
    "Russia_UstIshim_IUP.DG",
    "Russia_UstIshim_IUP_snpAD.DG",
    "Russia_WestSiberia_N.DG",
    "Scotland_LBA.DG",
    "Spain_IA_Celt.DG",
    "Turkey_Marmara_Barcin_N.DG",
    "Turkey_Southeast_Byzantine.DG",
    "Ukraine_Mesolithic_DG.DG",
    "USA_AK_Ancient_Athabaskan_1100BP.DG",
}

# QC-flagged outlier/discovery/relative subsets -- these are individuals
# flagged by AADR's curators as deviating from their nominal population
# (excess admixture, relatedness, being part of the original "discovery"
# panel later superseded, etc). Excluded rather than treated as their own
# population, and excluded rather than folded into their parent population
# (they were split out specifically *because* they're not representative).
OUTLIER_PATTERN = re.compile(
    r"_o[A-Z0-9]|_o$|\(relative\)|_relative|\(discovery\)|_discovery",
)
# A few more that don't match the regex cleanly but are the same idea,
# found by manual inspection.
KNOWN_OUTLIER_EXTRA = {
    "PEL_SomeAfrican.DG",
}

# Explicit renames for names that collide between AC82 and AADR (per
# project decision 2026-09-01: AADR's curated version wins in every case;
# AC82's corresponding entry is dropped when building the merged reference
# -- see merge_1000g_ac82_aadr.py). Mapping: raw AADR label -> clean final
# population name to use in the merged reference.
COLLISION_RENAMES = {
    "Russian.DG": "Russian",
    "BedouinA.DG": "Bedouin_A",
    "BedouinB.DG": "Bedouin_B",
    "Cambodian.DG": "Cambodian",
    "Dai.DG": "Dai",
    "Mandenka.DG": "Mandenka",
    "Mayan.DG": "Mayan",
    "Pima.DG": "Pima",
    "Balochi.DG": "Balochi",
}

# Superpopulation scheme: standardized on AC82's 7-category scheme (it's
# more granular than 1000G's 5 continental codes -- e.g. it splits
# Middle_East_North_Africa out of Africa/South_Asia), PLUS a new "Oceania"
# category for Papuan/Aboriginal-Australian/Melanesian populations that
# none of AC82's or 1000G's original categories cover well. 1000G's own
# 5-code scheme (AFR/AMR/EAS/EUR/SAS) gets remapped onto this same set in
# the merge script so aeon's superpopulation summary is coherent across
# all three merged sources rather than mixing naming conventions.
#
# This dict only needs to cover the labels that survive to KEEP status
# (see the exclusion rules above) -- anything not listed here that reaches
# the KEEP branch gets superpopulation=REVIEW_NEEDED so it surfaces for
# manual sign-off rather than silently guessing.
SUPERPOP_MAP = {
    # --- Europe ---
    "Albanian.DG": "Europe", "Basque.DG": "Europe", "Bergamo.DG": "Europe",
    "Bulgarian.DG": "Europe", "Cretan.DG": "Europe", "Czech.DG": "Europe",
    "Estonian.DG": "Europe", "Finnish.DG": "Europe", "French.DG": "Europe",
    "Greek_1.DG": "Europe", "Greek_2.DG": "Europe", "Hungarian.DG": "Europe",
    "Icelandic.DG": "Europe", "Italian_North.DG": "Europe",
    "Italian_Sardinian.DG": "Europe", "Norwegian.DG": "Europe",
    "Orcadian.DG": "Europe", "Polish.DG": "Europe", "Sardinian.DG": "Europe",
    "Saami.DG": "Europe", "Spanish.DG": "Europe", "Tuscan_1.DG": "Europe",
    "Turkish.DG": "Europe", "Russian.DG": "Europe",
    "English.DG": "Europe",  # only relevant if included in the British pool
    "IBS_CanaryIslands.DG": "Europe",  # small (n=2) but genuine modern Spanish subgroup, not an outlier
    # --- Middle East / North Africa ---
    "Abkhasian.DG": "Middle_East_North_Africa",
    "Adygei.DG": "Middle_East_North_Africa",
    "Armenian.DG": "Middle_East_North_Africa",
    "BedouinA.DG": "Middle_East_North_Africa",
    "BedouinB.DG": "Middle_East_North_Africa",
    "Chechen.DG": "Middle_East_North_Africa",
    "Druze.DG": "Middle_East_North_Africa",
    "Georgian.DG": "Middle_East_North_Africa",
    "Iranian.DG": "Middle_East_North_Africa",
    "Jew_Iraqi.DG": "Middle_East_North_Africa",
    "Jew_Yemenite.DG": "Middle_East_North_Africa",
    "Jordanian.DG": "Middle_East_North_Africa",
    "Lezgin.DG": "Middle_East_North_Africa",
    "Mozabite.DG": "Middle_East_North_Africa",
    "Palestinian.DG": "Middle_East_North_Africa",
    "Russia_NorthOssetian.DG": "Middle_East_North_Africa",
    "Saharawi.DG": "Middle_East_North_Africa",
    "Samaritan.DG": "Middle_East_North_Africa",
    "Tajik.DG": "Middle_East_North_Africa",
    # --- South Asia ---
    "Balochi.DG": "South_Asia", "Bengali.DG": "South_Asia",
    "Brahmin.DG": "South_Asia", "Brahui.DG": "South_Asia",
    "Burusho.DG": "South_Asia", "Hazara.DG": "South_Asia",
    "Irula.DG": "South_Asia", "Kalash.DG": "South_Asia",
    "Kapu.DG": "South_Asia", "Khonda_Dora.DG": "South_Asia",
    "Madiga.DG": "South_Asia", "Makrani.DG": "South_Asia",
    "Mala.DG": "South_Asia", "Pathan.DG": "South_Asia",
    "Punjabi.DG": "South_Asia", "Relli_1.DG": "South_Asia",
    "Relli_2.DG": "South_Asia", "Sindhi_Pakistan.DG": "South_Asia",
    "Yadava.DG": "South_Asia",
    # --- East / Southeast Asia ---
    "Ami.DG": "East_Southeast_Asia", "Atayal.DG": "East_Southeast_Asia",
    "Burmese.DG": "East_Southeast_Asia", "Cambodian.DG": "East_Southeast_Asia",
    "China_Lahu.DG": "East_Southeast_Asia", "Dai.DG": "East_Southeast_Asia",
    "Daur.DG": "East_Southeast_Asia", "Dusun.DG": "East_Southeast_Asia",
    "Han.DG": "East_Southeast_Asia", "Hezhen.DG": "East_Southeast_Asia",
    "Igorot.DG": "East_Southeast_Asia",
    "Indonesia_Alor.DG": "East_Southeast_Asia",
    "Indonesia_Bali_Gadon.DG": "East_Southeast_Asia",
    "Indonesia_Borneo_Maanyan.DG": "East_Southeast_Asia",
    "Indonesia_Flores_Bena.DG": "East_Southeast_Asia",
    "Indonesia_Flores_Bere.DG": "East_Southeast_Asia",
    "Indonesia_Flores_Cibol.DG": "East_Southeast_Asia",
    "Indonesia_Flores_Rampasasa.DG": "East_Southeast_Asia",
    "Indonesia_Java_Dieng.DG": "East_Southeast_Asia",
    "Indonesia_Kei_Faan.DG": "East_Southeast_Asia",
    "Indonesia_Kei_Ohoidertutu.DG": "East_Southeast_Asia",
    "Indonesia_Kei_Waur.DG": "East_Southeast_Asia",
    "Indonesia_Lembata_TimurKadakewa.DG": "East_Southeast_Asia",
    "Indonesia_Lembata_Waipukang.DG": "East_Southeast_Asia",
    "Indonesia_Mentawai.DG": "East_Southeast_Asia",
    "Indonesia_Nias_Gomo.DG": "East_Southeast_Asia",
    "Indonesia_Nias_Hilitobara.DG": "East_Southeast_Asia",
    "Indonesia_Sulawesi_Kajang.DG": "East_Southeast_Asia",
    "Indonesia_Sulawesi_Mandar.DG": "East_Southeast_Asia",
    "Indonesia_Sumatra_Toba.DG": "East_Southeast_Asia",
    "Indonesia_Tanimbar_AlusiKrawain.DG": "East_Southeast_Asia",
    "Indonesia_Tanimbar_Fordata.DG": "East_Southeast_Asia",
    "Indonesia_Tanimbar_Makatian.DG": "East_Southeast_Asia",
    "Indonesia_Tanimbar_SangliatDol.DG": "East_Southeast_Asia",
    "Indonesia_Tanimbar_Tumbur.DG": "East_Southeast_Asia",
    "Japanese.DG": "East_Southeast_Asia", "Kinh.DG": "East_Southeast_Asia",
    "Korean.DG": "East_Southeast_Asia", "Kusunda.DG": "East_Southeast_Asia",
    "Kyrgyz_Kyrgyzstan.DG": "East_Southeast_Asia",
    "Miao.DG": "East_Southeast_Asia", "Mongola.DG": "East_Southeast_Asia",
    "Naxi.DG": "East_Southeast_Asia", "Oroqen.DG": "East_Southeast_Asia",
    "She.DG": "East_Southeast_Asia", "Thai.DG": "East_Southeast_Asia",
    "Tu.DG": "East_Southeast_Asia", "Tubalar.DG": "East_Southeast_Asia",
    "Tujia.DG": "East_Southeast_Asia", "Uyghur.DG": "East_Southeast_Asia",
    "Xibo.DG": "East_Southeast_Asia", "Yakut.DG": "East_Southeast_Asia",
    "Yi.DG": "East_Southeast_Asia", "Altaian.DG": "East_Southeast_Asia",
    "Mansi.DG": "East_Southeast_Asia", "Itelmen.DG": "East_Southeast_Asia",
    "Chukchi.DG": "East_Southeast_Asia", "Even.DG": "East_Southeast_Asia",
    "Ulchi.DG": "East_Southeast_Asia",
    "Eskimo_Chaplin.DG": "East_Southeast_Asia",
    "Eskimo_Naukan.DG": "East_Southeast_Asia",
    "Eskimo_Sireniki.DG": "East_Southeast_Asia",
    "Aleut.DG": "East_Southeast_Asia",
    # --- Africa ---
    "BantuHerero.DG": "Africa", "BantuKenya.DG": "Africa",
    "BantuSA.DG": "Africa", "BantuSA_Herero.DG": "Africa",
    "BantuSA_Ovambo.DG": "Africa", "BantuTswana.DG": "Africa",
    "Biaka.DG": "Africa", "Dinka.DG": "Africa", "Gambian.DG": "Africa",
    "Igbo.DG": "Africa", "Ju_hoan_North.DG": "Africa",
    "Khomani_San.DG": "Africa", "Kongo.DG": "Africa",
    "Lemande.DG": "Africa", "Luhya.DG": "Africa", "Luo.DG": "Africa",
    "Mandenka.DG": "Africa", "Masai.DG": "Africa", "Mbuti.DG": "Africa",
    "Mende.DG": "Africa", "Yoruba.DG": "Africa",
    "Kenya_Somali.DG": "Middle_East_North_Africa",  # Somali -> MENA/Horn, not sub-Saharan Bantu belt
    # --- Americas ---
    "Chane.DG": "Americas", "Chipewyan.DG": "Americas",
    "Cree1.DG": "Americas", "Cree2.DG": "Americas",
    "Karitiana.DG": "Americas", "Mayan.DG": "Americas",
    "Mexico_Zapotec.DG": "Americas", "Mixe.DG": "Americas",
    "Mixtec_1.DG": "Americas", "Mixtec_2.DG": "Americas",
    "Nahua.DG": "Americas", "Piapoco.DG": "Americas", "Pima.DG": "Americas",
    "Quechua.DG": "Americas", "Surui.DG": "Americas",
    # --- Oceania (new category -- not in the original AC82/1000G schemes) ---
    "Australian.DG": "Oceania", "Bougainville.DG": "Oceania",
    "Hawaiian.DG": "Oceania", "Maori.DG": "Oceania",
    "Nasioi.DG": "Oceania", "PapuaNewGuinea.DG": "Oceania",
    "PapuaNewGuinea_NewBritain_Baining.DG": "Oceania",
    "PapuaNewGuinea_NewGuinea_MappiKorowai.DG": "Oceania",
    "Papuan.DG": "Oceania",
}


def classify(raw_label):
    """Return (decision, final_population, superpopulation, reason)."""
    if not raw_label.endswith(".DG"):
        return "EXCLUDE", None, None, "non-.DG suffix (default: ancient/archaeological)"

    if raw_label.startswith("Ignore_"):
        return "EXCLUDE", None, None, "AADR curator 'Ignore_' flag"

    bare = raw_label[: -len(".DG")]  # strip suffix for code comparisons

    if bare in G1K_CODES:
        return "EXCLUDE", None, None, "duplicate of aeon's built-in 1000G reference"

    low = raw_label.lower()
    if any(s in low for s in ARCHAIC_SUBSTRINGS):
        return "EXCLUDE", None, None, "archaic hominin, not modern human"

    if raw_label in KNOWN_ANCIENT_DESPITE_DG:
        return "EXCLUDE", None, None, "manually identified as ancient despite .DG suffix"

    # NOTE: match against `bare` (suffix stripped), not raw_label -- a bare
    # trailing "_o" immediately followed by ".DG" (e.g. "YRI_o.DG") does NOT
    # match "_o$" against the full raw_label, since the string actually ends
    # in ".DG". Found via real data: 13 "_o"-suffixed outlier labels slipped
    # through when this checked raw_label instead of bare.
    if raw_label in KNOWN_OUTLIER_EXTRA or OUTLIER_PATTERN.search(bare):
        return "EXCLUDE", None, None, "QC outlier/discovery/relative subset"

    final_name = COLLISION_RENAMES.get(raw_label, bare)
    superpop = SUPERPOP_MAP.get(raw_label, "REVIEW_NEEDED")
    reason = "kept: present-day reference population"
    if raw_label in COLLISION_RENAMES:
        reason += " (AADR version wins over AC82 duplicate)"
    return "KEEP", final_name, superpop, reason


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    ind_path, out_prefix = sys.argv[1], sys.argv[2]

    pop_counts = Counter()
    sample_pop = []  # (sample, raw_population)

    # Tolerate either the .tsv (tab-delimited, header "sample sex population")
    # or the raw eigenstrat .ind (whitespace-delimited, no header).
    with open(ind_path) as fh:
        first = fh.readline()
        is_tsv_with_header = first.strip().lower().replace(" ", "\t").startswith("sample")
        rows = []
        if not is_tsv_with_header:
            rows.append(first)
        rows.extend(fh.readlines())

    for line in rows:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) < 3:
            continue
        sample, sex, pop = parts[0], parts[1], parts[2]
        sample_pop.append((sample, pop))
        pop_counts[pop] += 1

    classification_rows = []
    final_sample_rows = []
    seen_raw = set()
    for raw_label, n in sorted(pop_counts.items(), key=lambda kv: -kv[1]):
        decision, final_name, superpop, reason = classify(raw_label)
        classification_rows.append(
            {
                "raw_label": raw_label,
                "n_samples": n,
                "decision": decision,
                "final_population": final_name or "",
                "superpopulation": superpop or "",
                "reason": reason,
            }
        )
        seen_raw.add(raw_label)

    keep_map = {
        row["raw_label"]: (row["final_population"], row["superpopulation"])
        for row in classification_rows
        if row["decision"] == "KEEP"
    }
    for sample, raw_label in sample_pop:
        if raw_label in keep_map:
            final_name, _ = keep_map[raw_label]
            final_sample_rows.append((sample, final_name))

    class_path = f"{out_prefix}_classification.DRAFT.csv"
    with open(class_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["raw_label", "n_samples", "decision", "final_population", "superpopulation", "reason"],
        )
        w.writeheader()
        w.writerows(classification_rows)

    map_path = f"{out_prefix}_sample_map.DRAFT.tsv"
    with open(map_path, "w") as f:
        f.write("sample\tfinal_population\n")
        for sample, final_name in final_sample_rows:
            f.write(f"{sample}\t{final_name}\n")

    n_keep = sum(1 for r in classification_rows if r["decision"] == "KEEP")
    n_exclude = len(classification_rows) - n_keep
    n_review = sum(1 for r in classification_rows if r["superpopulation"] == "REVIEW_NEEDED")
    final_pop_count = len(set(final_name for final_name, _ in keep_map.values()))
    sys.stderr.write(
        f"[DONE] {len(classification_rows):,} raw population labels seen "
        f"({len(sample_pop):,} total samples).\n"
    )
    sys.stderr.write(f"[DONE] {n_keep:,} raw labels KEPT, {n_exclude:,} EXCLUDED.\n")
    sys.stderr.write(
        f"[DONE] {n_keep:,} kept labels collapse to {final_pop_count:,} final "
        f"population names (collisions like Bedouin_A/B stay distinct; no "
        f"other merging applied here).\n"
    )
    if n_review:
        sys.stderr.write(
            f"[WARN] {n_review:,} kept population(s) have superpopulation=REVIEW_NEEDED "
            f"(no mapping in SUPERPOP_MAP yet) -- check {class_path} and either add them "
            f"to SUPERPOP_MAP and rerun, or hand-edit the DRAFT csv directly.\n"
        )
    sys.stderr.write(f"[DONE] Classification table: {class_path}\n")
    sys.stderr.write(f"[DONE] Kept-sample map:      {map_path}\n")


if __name__ == "__main__":
    main()