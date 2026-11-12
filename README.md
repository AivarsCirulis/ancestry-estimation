# ancestry-estimation
Ancestry estimation pipeline for a Latvian cohort using aeon against a custom 224-population reference (1000 Genomes + UK Biobank AC82 + AADR + a custom Latvian population), producing per-participant PDF reports in English or Latvian.

Estimates genetic ancestry for a cohort of Latvian participants using `aeon`,
against a custom 224-population reference merged from the 1000 Genomes
Project, UK Biobank's AC82 fine-scale panel, the Allen Ancient DNA Resource
(AADR), and a custom "Latvian" population built from this project's own
reference data. Produces one PDF ancestry report per participant — a pie
chart of their ancestry components with a color-swatch legend — in English
or Latvian.

> For the full narrative writeup (design rationale, verified script
> internals, example reports) see `Latvian_Ancestry_Pipeline.docx`. This file is the condensed, run-it version.

## Requirements

- An HPC cluster (or any machine) with PBS/Torque (`qsub`) and a conda/mamba
  environment named `aeon-ancestry` with `aeon` and `bcftools` installed.
- R with cairo support (`capabilities("cairo")` → `TRUE`; true on
  essentially every standard Linux/macOS/Windows R build) — required for
  Latvian-language reports, since diacritics (ā, č, ē, ģ, ī, ķ, ļ, ņ, š, ū, ž)
  don't render through R's base `pdf()` device at all, regardless of
  encoding.
- Python 3 with `pandas` (used throughout the merge/cleanup scripts).
- `bcftools`/`tabix`, for indexing VCFs.

## Repository layout

```
ac82_to_vcf37.py, run_ac82_for_aeon.pbs,          \
extract_ac82_af_for_aeon.py                        } Stage A — AC82 branch

aadr_classify_populations.py, patch_classifier.py, \
aadr_populations_to_af.sh, pivot_aadr_af.py         } Stage B — AADR branch

extract_latvia_af.sh, format_latvia_af.py           } Stage C — Latvian population

merge_1000g_ac82_aadr.py, resolve_duplicate_populations.py, \
add_latvian_population.py, standardize_and_finalize.py       } Stage D — merge & finalize reference

split_samples_into_batches.sh, run_batched_cohort_aeon.pbs, \
merge_batch_ancestry_results.py                              } Stage E — classify the cohort

ethnicity_name_fixes.csv, apply_population_renames.py        } Stage F — naming fix

generate_ancestry_pie_reports.R, population_names_lv.csv     } Stage G — PDF reports

run_new_sample_aeon.pbs                                       } Stage H — add one new sample
```

## Quickstart

### Part 1 — Build the reference (one-time)

Skip this section entirely if `aeon_reference_final_allele_freqs.tsv.gz` /
`aeon_reference_final_population_labels.tsv` already exist in your project
directory — that's the finalized, 224-population reference every
classification step below reads from, and rebuilding it is a multi-hour,
multi-source process you only do once (or after a source panel changes).

```bash
# Stage A — AC82 branch (VCF conversion, liftover, allele-freq extraction)
qsub run_ac82_for_aeon.pbs

# Stage B — AADR branch (see each script's own --help/header for exact flags)
python aadr_classify_populations.py ...
python patch_classifier.py ...
./aadr_populations_to_af.sh ...
python pivot_aadr_af.py ...

# Stage C — Latvian reference population
./extract_latvia_af.sh ...
python format_latvia_af.py ...

# Stage D — merge, deduplicate, add Latvian, finalize (run in this order)
python merge_1000g_ac82_aadr.py ...
python resolve_duplicate_populations.py \
    --af merged_ac82_1000g_aadr_allele_freqs.tsv.gz \
    --labels merged_ac82_1000g_aadr_population_labels.tsv ...

python add_latvian_population.py \
    --af merged_ac82_1000g_aadr_final_allele_freqs.tsv.gz \
    --labels merged_ac82_1000g_aadr_final_population_labels.tsv \
    --latvia-af latvia_allele_freqs_hg38.tsv.gz ...

python standardize_and_finalize.py \
    --af merged_ac82_1000g_aadr_latvia_allele_freqs.tsv.gz \
    --labels merged_ac82_1000g_aadr_latvia_population_labels.tsv \
    --out-prefix aeon_reference_final
```

End result: `aeon_reference_final_allele_freqs.tsv.gz` /
`aeon_reference_final_population_labels.tsv` — 224 populations across 7
superpopulation groups. These two files are the reference every step below
reads from and never modifies.

### Part 2 — Classify the full cohort and generate reports

```bash
# Stage E — classify the cohort (BATCHED — see Gotchas, the single-shot
# job OOM-kills on a full-size cohort VCF)
./split_samples_into_batches.sh <cohort.vcf.gz> 200 batches
qsub run_batched_cohort_aeon.pbs

# Stage F — fix bare-country-name labels (e.g. "Denmark" -> "Danish")
python apply_population_renames.py \
    --renames ethnicity_name_fixes.csv \
    --file latvia_full_cohort_ancestry_ae.csv \
    --format csv

# Stage G — generate one PDF report per participant
Rscript generate_ancestry_pie_reports.R \
    latvia_full_cohort_ancestry_ae.csv ancestry_reports 1.0        # English

Rscript generate_ancestry_pie_reports.R \
    latvia_full_cohort_ancestry_ae.csv ancestry_reports_lv 1.0 lv  # Latvian
```

`population_names_lv.csv` (the English→Latvian population-name lookup) must
sit next to `generate_ancestry_pie_reports.R`, or its path can be passed as
a 5th argument. End result: one `<participant>_ancestry_report.pdf` per
cohort member in `ancestry_reports/` (and/or `ancestry_reports_lv/`).

### Part 3 — Add a new sample later (no cohort/reference rebuild needed)

For a participant recruited after the cohort above was finalized. The
reference is never touched — this sample is just classified against it,
the same way the cohort was.

```bash
# the new sample's VCF needs an index
bcftools index -t /path/to/the/new/sample.vcf.gz

# edit VCF=/OUT_PREFIX= in run_new_sample_aeon.pbs for this sample, then:
qsub run_new_sample_aeon.pbs

Rscript generate_ancestry_pie_reports.R \
    new_sample_ancestry_ae.csv new_sample_report 1.0 lv
```

## Gotchas

- **The single-shot cohort job OOM-kills.** `run_full_cohort_aeon.pbs` /
  `run_full_cohort_aeon_1.pbs` run `aeon` once against the entire cohort VCF
  in one process; aeon's genotype-extraction step builds an unvectorized,
  per-sample-per-locus structure that doesn't scale to a cohort this size
  even at 32GB, and the job gets killed by the OOM killer partway through
  "Extracting genotypes" (confirmed on this project's own job log). Use the
  batched path (`split_samples_into_batches.sh` → `run_batched_cohort_aeon.pbs`
  → `merge_batch_ancestry_results.py`) instead — that's what actually
  classified the real cohort here, in 11 batches of 200 (last batch 191).
- **Latvian output needs `encoding="UTF-8"`, not `fileEncoding="UTF-8"`,**
  when reading `population_names_lv.csv` — already handled inside
  `generate_ancestry_pie_reports.R`, but worth knowing if you ever touch
  that code. Under a non-UTF-8 system locale (common on HPC clusters),
  `fileEncoding=` silently re-encodes the file into the system's native
  (ASCII) encoding while reading, corrupting every diacritic and
  potentially misparsing the whole table.
- **Latvian output needs `cairo_pdf()`, not `pdf()`** — R's base `pdf()`
  device has no glyphs for Latvian's diacritics at all, regardless of
  string encoding. Already handled inside the script, which checks
  `capabilities("cairo")` up front and stops with a clear error if it's
  somehow unavailable rather than silently producing garbled text.
- **A new sample's VCF needs a tabix/CSI index** before Stage H's `aeon`
  call — it fails with `"Index file required for VCF.GZ or BCF file"`
  otherwise. `bcftools index -t <file>.vcf.gz` fixes it; if that itself
  errors, the file is likely plain `gzip`-compressed rather than
  BGZF-compressed (check with `htsfile`) — re-compress with `bgzip` first.

## Citations

This pipeline builds on several external tools and reference panels. If you
use it, please cite:

- **AEon** (the ancestry-estimation tool itself) — *"AEon: A global genetic
  ancestry estimation tool"* (2024), bioRxiv,
  https://www.biorxiv.org/content/10.1101/2024.06.18.599246. Software:
  [`aeon-ancestry`](https://pypi.org/project/aeon-ancestry/) on PyPI.
- **1000 Genomes Project** — The 1000 Genomes Project Consortium, *"A global
  reference for human genetic variation,"* Nature 526, 68–74 (2015).
  doi:[10.1038/nature15393](https://doi.org/10.1038/nature15393).
- **AADR (Allen Ancient DNA Resource)** — Mallick et al., *"The Allen
  Ancient DNA Resource (AADR): a curated compendium of ancient human
  genomes,"* Scientific Data (2024). Check the
  [Reich Lab's AADR page](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data)
  for their currently preferred citation wording before publishing.
- **AC82** (the 82-fine-scale-population UK Biobank panel) — Lawson, Hu &
  Myers, *"Frequency of variants in the UK Biobank for 82 fine-scale
  populations,"* University of Bristol (2023).
  doi:[10.5523/bris.3g5oatl682kz82as80jakjrq91](https://doi.org/10.5523/bris.3g5oatl682kz82as80jakjrq91).
  **This dataset's own terms-of-use page states "all rights reserved" with
  no confirmed open-redistribution license** — before publishing or sharing
  anything derived from it (including this pipeline's merged
  allele-frequency tables), check that dataset's actual terms directly
  rather than assuming they permit redistribution.

The Latvian reference population and the cohort itself are this project's
own participant data, governed by its own ethics approval/consent terms —
not included in this repository, and not covered by any of the citations
above.

## Citing this repository

If you use this pipeline itself (as opposed to, or in addition to, the
external tools and reference panels cited above), please cite:

> Cīrulis, A. (2026). *Ancestry Estimation Pipeline for a Latvian Cohort*
> [Computer software]. NIRI. https://github.com/AivarsCirulis/ancestry-estimation

A machine-readable [`CITATION.cff`](CITATION.cff) is included at the repo
root — GitHub picks this up automatically and shows a "Cite this
repository" button in the sidebar with ready-made APA/BibTeX citations, no
further setup needed. If you'd like a permanent, version-specific DOI
instead of a repo-URL citation (recommended once this is cited in a paper),
link the repo to [Zenodo](https://zenodo.org) — it's a free GitHub login,
and every tagged release then gets its own archived DOI automatically.

## License

The scripts in this repository are released under the [MIT License](LICENSE)
(see `LICENSE`). This covers the pipeline code only — it does not grant any
rights to the external reference panels or tools listed under
[Citations](#citations) above, each of which has its own terms.
