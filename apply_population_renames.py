#!/usr/bin/env python
"""
Patch population-name labels in an ALREADY-BUILT reference or ALREADY-
COMPUTED result file, without rerunning aeon or standardize_and_finalize.py.

Why this exists: standardize_and_finalize.py's ETHNICITY_NAME_FIXES (added
2026-09-10) fixes the naming at the SOURCE for anything built or rebuilt
from now on. But it can't retroactively fix files that were already
generated before that fix existed:
  - aeon_reference_final_population_labels.tsv / _allele_freqs.tsv.gz --
    for these, just re-run standardize_and_finalize.py from your existing
    merged_ac82_1000g_aadr_latvia_* files; it's a fast pandas pass over
    already-merged data, not a pipeline rerun, and will regenerate both
    files with the corrected names baked in.
  - The already-computed per-participant aeon output for the full 2191-
    person cohort (e.g. *_ae.csv, and the matching sheet in an exported
    .xlsx) -- these were produced by an aeon run against the OLD reference
    labels, so re-running standardize_and_finalize.py doesn't touch them.
    THIS script is for exactly that case: a pure string relabel of the
    population axis, since the underlying ancestry fractions haven't
    changed at all -- only what we call each population.

*_ae_stats.csv does NOT need this -- it has no population column (it's
per-sample QC metrics: FractionLociImputed, FractionLociNonCalled, etc.).

Usage:
    # CSV (e.g. a PREFIX_ae.csv, population names in the first column)
    python apply_population_renames.py \\
        --renames ethnicity_name_fixes.csv \\
        --file latvia_full_cohort_ancestry_ae.csv --format csv

    # A population-labels TSV (Population column, not necessarily first)
    python apply_population_renames.py \\
        --renames ethnicity_name_fixes.csv \\
        --file aeon_reference_final_population_labels.tsv --format csv \\
        --sep '\\t' --column Population

    # An .xlsx sheet (writes a NEW file so you keep the original)
    python apply_population_renames.py \\
        --renames ethnicity_name_fixes.csv \\
        --file cohort_results.xlsx --format xlsx \\
        --sheet latvia_full_cohort_ancestry_ae --column 0 \\
        --out cohort_results.renamed.xlsx

By default this edits CSV/TSV files in place (after printing what it's
about to do) and refuses to touch an .xlsx without an explicit --out, so
you never silently lose the original of a file that's harder to diff.
"""
import argparse
import csv
import sys


def load_renames(path):
    renames = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        cols = r.fieldnames or []
        if "old_population" not in cols or "new_population" not in cols:
            raise SystemExit(
                f"FATAL: {path} must have columns old_population,new_population -- got {cols}"
            )
        for row in r:
            old, new = row["old_population"], row["new_population"]
            if old in renames and renames[old] != new:
                raise SystemExit(f"FATAL: conflicting renames for {old!r} in {path}")
            renames[old] = new
    return renames


def apply_to_csv(path, renames, sep, column, out_path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter=sep))
    if not rows:
        raise SystemExit(f"FATAL: {path} is empty")

    header = rows[0]
    col_idx = 0
    if column is not None:
        if column.isdigit():
            col_idx = int(column)
        elif column in header:
            col_idx = header.index(column)
        else:
            raise SystemExit(f"FATAL: column {column!r} not found in header {header}")

    seen_old = set()
    changed = 0
    for row in rows[1:]:
        if col_idx >= len(row):
            continue
        old = row[col_idx]
        if old in renames:
            row[col_idx] = renames[old]
            seen_old.add(old)
            changed += 1

    missing = set(renames) - seen_old
    if missing:
        sys.stderr.write(
            f"[WARN] {len(missing)} rename(s) not found in {path} (already renamed, or wrong file?): "
            f"{sorted(missing)}\n"
        )

    new_values = [row[col_idx] for row in rows[1:] if col_idx < len(row)]
    dupes = {v for v in new_values if new_values.count(v) > 1}
    if dupes:
        raise SystemExit(
            f"FATAL: rename would create duplicate value(s) in column {col_idx} of {path}: {sorted(dupes)}"
        )

    dest = out_path or path
    with open(dest, "w", newline="") as f:
        csv.writer(f, delimiter=sep).writerows(rows)
    sys.stderr.write(f"[DONE] {changed} row(s) relabeled in {dest}\n")


def apply_to_xlsx(path, renames, sheet, column, out_path):
    import openpyxl

    if not out_path:
        raise SystemExit("FATAL: --out is required for --format xlsx (won't overwrite an .xlsx in place)")

    wb = openpyxl.load_workbook(path)
    if sheet not in wb.sheetnames:
        raise SystemExit(f"FATAL: sheet {sheet!r} not found -- available: {wb.sheetnames}")
    ws = wb[sheet]

    col_idx = int(column) if str(column).isdigit() else None
    if col_idx is None:
        header = [c.value for c in ws[1]]
        if column not in header:
            raise SystemExit(f"FATAL: column {column!r} not found in header {header}")
        col_idx = header.index(column)
    excel_col = col_idx + 1  # openpyxl is 1-indexed

    seen_old = set()
    changed = 0
    new_values = []
    for row in ws.iter_rows(min_row=2):
        cell = row[col_idx] if col_idx < len(row) else None
        if cell is None:
            continue
        old = cell.value
        if old in renames:
            cell.value = renames[old]
            seen_old.add(old)
            changed += 1
        new_values.append(cell.value)

    missing = set(renames) - seen_old
    if missing:
        sys.stderr.write(
            f"[WARN] {len(missing)} rename(s) not found in sheet {sheet!r} of {path}: {sorted(missing)}\n"
        )
    dupes = {v for v in new_values if new_values.count(v) > 1}
    if dupes:
        raise SystemExit(f"FATAL: rename would create duplicate value(s): {sorted(dupes)}")

    wb.save(out_path)
    sys.stderr.write(f"[DONE] {changed} row(s) relabeled -> {out_path}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--renames", required=True, help="CSV with old_population,new_population columns")
    p.add_argument("--file", required=True)
    p.add_argument("--format", choices=["csv", "xlsx"], required=True)
    p.add_argument("--sep", default=",", help="CSV/TSV delimiter (default ',', use $'\\t' for TSV)")
    p.add_argument("--sheet", help="xlsx sheet name (required for --format xlsx)")
    p.add_argument("--column", default="0", help="Column name or 0-based index holding population names (default: first column)")
    p.add_argument("--out", help="Write to a new path instead of editing --file in place (required for xlsx)")
    args = p.parse_args()

    renames = load_renames(args.renames)
    sys.stderr.write(f"[INFO] loaded {len(renames)} rename(s) from {args.renames}\n")

    if args.format == "csv":
        apply_to_csv(args.file, renames, args.sep, args.column, args.out)
    else:
        if not args.sheet:
            raise SystemExit("FATAL: --sheet is required for --format xlsx")
        apply_to_xlsx(args.file, renames, args.sheet, args.column, args.out)


if __name__ == "__main__":
    main()
