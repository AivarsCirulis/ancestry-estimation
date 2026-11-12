#!/usr/bin/env Rscript
###############################################################################
# generate_ancestry_pie_reports.R
#
# Reads aeon's ancestry-fraction output (wide format: one row per reference
# population, one column per participant, values 0-1 summing to ~1 per
# column -- exactly the shape of latvia_full_cohort_ancestry_ae.csv) and
# writes ONE PDF report PER PARTICIPANT: a pie chart of their ancestry
# components (by %) with a color-swatch legend, no per-slice leader lines
# and no detailed table/disclaimer -- just the chart and bigger legend text.
#
# Deliberately writes one file per participant rather than one combined
# multi-page PDF, since a report "given to participants" should show each
# person only their own results, not the whole cohort's in one file.
#
# Populations below --threshold (default 1%) are folded into a single grey
# "Other" slice rather than cluttering the chart with dozens of near-zero
# wedges -- aeon's output is naturally sparse (most populations come out
# at exactly 0 for any given person), so most reports will show a small,
# readable number of slices.
#
# LANGUAGE: pass "lv" as the 4th argument to render in Latvian -- population
# names and the report's own text (title, "Participant ID", the "Other"
# slice label) are translated via population_names_lv.csv (a reviewed
# English -> Latvian lookup sitting next to this script; anything not found
# there falls back to the English name, with a one-time warning listing
# what's missing so a gap is visible rather than silently inconsistent).
# Default is "en" (unchanged English behavior).
#
# This always renders through cairo_pdf() rather than base pdf(): base
# pdf()'s built-in Type1 fonts don't have glyphs for Latvian's macron/
# cedilla letters (ā, č, ē, ģ, ī, ķ, ļ, ņ, š, ū, ž) -- confirmed by testing,
# they render as missing-glyph boxes and even throw off text centering.
# cairo_pdf() renders through the system's actual font (Pango/Cairo), which
# has full Unicode coverage, and produces visually equivalent output for
# plain English text too, so this isn't a language-conditional branch.
#
# Separately: population_names_lv.csv is read with encoding="UTF-8", not
# fileEncoding="UTF-8" -- on a non-UTF-8 system locale (seen on this
# project's own cluster), fileEncoding= re-encodes into the system's
# native (here, ASCII) encoding while reading, which silently corrupts
# every diacritic and can even misparse the whole table. encoding="UTF-8"
# just tags the file's already-correct UTF-8 bytes as UTF-8, with no lossy
# re-encode, and works regardless of the system locale.
# Requires capabilities("cairo") to be TRUE on the R installation in use
# (true on virtually every Linux/macOS/Windows R build; the script checks
# and fails with a clear message if not).
#
# Usage:
#   Rscript generate_ancestry_pie_reports.R <ae.csv> <output_dir> [threshold_pct] [lang] [lv_lookup_csv]
#
# Examples:
#   Rscript generate_ancestry_pie_reports.R latvia_full_cohort_ancestry_ae.csv ancestry_reports 1.0
#   Rscript generate_ancestry_pie_reports.R latvia_full_cohort_ancestry_ae.csv ancestry_reports_lv 1.0 lv
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript generate_ancestry_pie_reports.R <ae.csv> <output_dir> [threshold_pct] [lang: en|lv] [lv_lookup_csv]")
}
ae_path <- args[1]
out_dir <- args[2]
threshold_pct <- if (length(args) >= 3) as.numeric(args[3]) else 1.0
lang <- if (length(args) >= 4) args[4] else "en"
lv_lookup_path <- if (length(args) >= 5) args[5] else "population_names_lv.csv"

if (!capabilities("cairo")) {
  stop("This R installation was not built with cairo support (capabilities(\"cairo\") is FALSE). ",
       "cairo_pdf() is required for correct text rendering (see the comment block at the top of this ",
       "script) -- install a cairo-enabled R build, e.g. via your distro's r-base-core package or conda-forge's r-base.")
}
if (!lang %in% c("en", "lv")) stop("lang must be 'en' or 'lv', got: ", lang)

if (!file.exists(ae_path)) stop("Input file not found: ", ae_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ae <- read.csv(ae_path, row.names = 1, check.names = FALSE, encoding = "UTF-8")
if (!"Superpopulation" %in% colnames(ae)) {
  stop("Expected a 'Superpopulation' column in the input CSV -- got: ",
       paste(colnames(ae), collapse = ", "))
}

population_names <- rownames(ae)
sample_cols <- setdiff(colnames(ae), "Superpopulation")
if (length(sample_cols) == 0) stop("No participant columns found (only Superpopulation present).")

# --- UI text, by language ---------------------------------------------------
ui_text <- if (lang == "lv") {
  list(title = "Tautību izcelsmes pārskats",
       participant_fmt = "Dalībnieka ID: %s",
       other_label = "citi (zem robežvērtības)")
} else {
  list(title = "Ancestry Estimate Report",
       participant_fmt = "Participant ID: %s",
       other_label = "Other (below threshold)")
}
for (nm in names(ui_text)) Encoding(ui_text[[nm]]) <- "UTF-8"

# --- Population name translation (lang == "lv" only) -------------------------
translate_population <- function(name) name
if (lang == "lv") {
  if (!file.exists(lv_lookup_path)) {
    stop("lang=lv but lookup file not found: ", lv_lookup_path,
         " (expected next to this script, or pass its path as the 5th argument).")
  }
  # encoding="UTF-8" (NOT fileEncoding="UTF-8") is required here: under a
  # non-UTF-8 system locale (common on HPC clusters -- confirmed on this
  # project's cluster), fileEncoding= re-encodes the file's bytes into the
  # system's native encoding while reading, which is lossy/silently wrong
  # when that native encoding is ASCII -- it corrupted every diacritic and
  # truncated the table to a single row in testing. encoding="UTF-8" instead
  # just tags the already-correct UTF-8 bytes as UTF-8 without re-encoding,
  # which parses correctly regardless of the system locale.
  lv_lookup_tbl <- read.csv(lv_lookup_path, stringsAsFactors = FALSE, encoding = "UTF-8")
  if (!all(c("english_name", "latvian_name") %in% colnames(lv_lookup_tbl))) {
    stop("lv_lookup_csv must have columns english_name,latvian_name -- got: ",
         paste(colnames(lv_lookup_tbl), collapse = ", "))
  }
  lv_map <- lv_lookup_tbl$latvian_name
  names(lv_map) <- lv_lookup_tbl$english_name
  Encoding(lv_map) <- "UTF-8"

  missing_from_lookup <- setdiff(population_names, names(lv_map))
  if (length(missing_from_lookup) > 0) {
    cat(sprintf(
      "[WARN] %d population name(s) have no Latvian translation in %s -- falling back to English for: %s\n",
      length(missing_from_lookup), lv_lookup_path, paste(missing_from_lookup, collapse = ", ")
    ))
  }

  translate_population <- function(name) {
    if (name %in% names(lv_map)) unname(lv_map[name]) else name
  }
}

cat(sprintf(
  "Found %d reference populations, %d participant(s). 'Other' threshold: %.1f%%. Language: %s\n",
  length(population_names), length(sample_cols), threshold_pct, lang
))

# A curated, print-friendly, mutually distinct palette (Tableau-10-derived)
# for up to 12 slices; falls back to an HCL rainbow for anything larger.
get_palette <- function(n) {
  base <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#86BCB6", "#D37295"
  )
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::hcl.colors(n, palette = "Dark 3")
}
other_color <- "#B0B0B0"

title_text      <- ui_text$title
legend_cex_max  <- 2.8   # ceiling on legend text size -- auto-shrinks per participant
                         # (a long population list needs smaller text to fit than a short one)
pie_radius      <- 1.0
max_margin_frac <- 0.55  # cap: legend's margin can use at most this fraction of panel
                         # width, so the pie itself can never be squeezed to a dot

# Draws one participant's title + pie + legend into whatever device/panel is
# currently active. Verified across a narrow device, an extreme 420x300
# stress case, and this exact 8.5x11 PDF shape in explore_ancestry_pie.R.
draw_report <- function(sample, major, labels, colors) {
  layout(matrix(c(1, 2), nrow = 2), heights = c(lcm(2.2), 1))

  # Title block -- fixed-height panel regardless of device size, so
  # Participant ID always has clear room below the title. A small top
  # margin plus a lower y for the title text keeps the bold title clear
  # of the page's top edge instead of crowding right up against it.
  par(mar = c(0, 1, 0.6, 1))
  plot.new()
  text(0.5, 0.55, title_text, cex = 2.1, font = 2)
  text(0.5, 0.12, sprintf(ui_text$participant_fmt, sample), cex = 0.95)

  # Pie + legend. Rather than guess formulas for how wide/tall legend()
  # will render its box at a given cex, ask it directly: legend(...,
  # plot = FALSE) returns the exact box size without drawing anything.
  # That's used to find the largest cex that (a) fits inside a right
  # margin capped at max_margin_frac of the panel width -- so the pie can
  # never be squeezed down to a dot -- and (b) fits inside the panel's
  # height, so a long population list doesn't run off the bottom either.
  #
  # A placeholder margin is used just to get a plot region open (legend
  # measurements need one); the real margin, sized to the final chosen
  # cex, is applied afterward and the pie re-drawn with it.
  edge <- max(1, pie_radius)
  panel_w_in <- par("fin")[1]
  in_per_line <- par("csi")
  max_margin_in <- panel_w_in * max_margin_frac

  par(mar = c(0.3, 1, 0.3, 10), xpd = NA)
  pie(major, labels = NA, col = colors, main = "",
      radius = pie_radius, clockwise = TRUE, init.angle = 90)
  scale_in_per_unit <- par("pin")[1] / diff(par("usr")[1:2])
  avail_h_in <- par("pin")[2] * 0.97

  measure <- function(cex) {
    lg <- legend(x = edge + 0.08, y = 0, yjust = 0.5, legend = labels, fill = colors,
                 bty = "n", cex = cex, y.intersp = 1.15, xjust = 0, plot = FALSE)
    c(w_in = lg$rect$w * scale_in_per_unit, h_in = lg$rect$h * scale_in_per_unit)
  }

  cex <- legend_cex_max
  m <- measure(cex)
  if (m["w_in"] > max_margin_in) cex <- unname(cex * max_margin_in / m["w_in"])
  m <- measure(cex)
  if (m["h_in"] > avail_h_in) cex <- unname(cex * avail_h_in / m["h_in"])
  final_cex <- max(0.5, cex)
  m <- measure(final_cex)

  side_margin_eff <- max(3, ceiling(min(max_margin_in, unname(m["w_in"]) + 0.05) / in_per_line))

  # new = TRUE is essential here: pie() calls plot.new() internally, and
  # without this we'd be issuing a THIRD plot.new()-equivalent inside a
  # layout() that only reserved 2 cells (title + body) -- the first two
  # were consumed by the title panel and the measurement pie() draw above.
  # A 3rd call wraps the layout back to cell 1 and clears the device. But
  # new = TRUE only stops the CELL from advancing -- it does NOT erase what
  # the measurement pie() already drew, so blank the whole figure region
  # (mar = 0 makes the plot region exactly equal to the figure region)
  # with a white rect first, THEN draw the real pie + legend into it.
  par(mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  rect(0, 0, 1, 1, col = "white", border = NA)

  par(mar = c(0.3, 1, 0.3, side_margin_eff), xpd = NA, new = TRUE)
  pie(major, labels = NA, col = colors, main = "",
      radius = pie_radius, clockwise = TRUE, init.angle = 90)

  legend(x = edge + 0.08, y = 0, yjust = 0.5, legend = labels, fill = colors,
         bty = "n", cex = final_cex, y.intersp = 1.15, xjust = 0)
}

n_written <- 0
n_skipped <- 0

for (sample in sample_cols) {
  fractions <- suppressWarnings(as.numeric(ae[[sample]]))
  names(fractions) <- population_names

  # Drop exact zeros/NAs -- populations aeon assigned literally no weight
  # to for this person, not "small" ones.
  fractions <- fractions[!is.na(fractions) & fractions > 0]
  if (length(fractions) == 0) {
    cat(sprintf("[WARN] %s: no nonzero ancestry fractions -- skipping.\n", sample))
    n_skipped <- n_skipped + 1
    next
  }

  # Renormalize to exactly 100% (aeon's own fractions can sum to slightly
  # under 1 due to rounding/imputation loss) so the pie chart's slices are
  # internally consistent and always add up to 100%.
  pct <- fractions / sum(fractions) * 100
  pct <- sort(pct, decreasing = TRUE)

  keep <- pct >= threshold_pct
  major <- pct[keep]
  minor_sum <- sum(pct[!keep])
  if (minor_sum > 0) {
    major <- c(major, minor_sum)
    names(major)[length(major)] <- "__OTHER__"
  }

  display_names <- vapply(names(major), function(n) {
    if (identical(n, "__OTHER__")) ui_text$other_label else translate_population(n)
  }, character(1))
  # Latvian uses a comma as the decimal separator (e.g. "65,3%"); R's sprintf
  # always uses a period regardless of locale, so swap it in for lang == "lv"
  # after formatting (names stay untouched -- only the numeric part can
  # contain a '.', since '.' never otherwise appears in a translated name).
  pct_str <- sprintf("%.1f", major)
  if (lang == "lv") pct_str <- gsub(".", ",", pct_str, fixed = TRUE)
  labels <- sprintf("%s (%s%%)", display_names, pct_str)
  colors <- get_palette(length(major))
  if ("__OTHER__" %in% names(major)) {
    colors[length(colors)] <- other_color
  }

  safe_name <- gsub("[^A-Za-z0-9._-]", "_", sample)
  pdf_path <- file.path(out_dir, sprintf("%s_ancestry_report.pdf", safe_name))

  cairo_pdf(pdf_path, width = 8.5, height = 11)
  draw_report(sample, major, labels, colors)
  dev.off()

  n_written <- n_written + 1
  if (n_written %% 100 == 0) cat(sprintf("  ... %d report(s) written so far\n", n_written))
}

cat(sprintf(
  "\nDone. %d report(s) written to %s (%d participant(s) skipped: no nonzero data).\n",
  n_written, out_dir, n_skipped
))
