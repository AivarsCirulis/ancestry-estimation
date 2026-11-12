#!/usr/bin/env python
"""
Merge per-batch aeon ancestry output (_ae.csv and _ae_stats.csv) back into
one combined result covering all samples -- the flip side of
split_samples_into_batches.sh.

Each batch's _ae.csv has the same Population/Superpopulation rows (since
every batch was run against the same reference) with different sample
columns, so batches merge column-wise. Each batch's _ae_stats.csv has one
row per sample, so those merge row-wise (concatenation).

Usage:
    python merge_batch_ancestry_results.py \\
        --ae-glob "batch_*_ancestry_ae.csv" \\
        --stats-glob "batch_*_ancestry_ae_stats.csv" \\
        --out-prefix latvia_full_cohort_ancestry
"""
import argparse
import glob
import sys
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ae-glob", required=True, help="Glob pattern matching every batch's _ae.csv")
    p.add_argument("--stats-glob", required=True, help="Glob pattern matching every batch's _ae_stats.csv")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    ae_files = sorted(glob.glob(args.ae_glob))
    stats_files = sorted(glob.glob(args.stats_glob))
    if not ae_files:
        raise SystemExit(f"FATAL: no files matched --ae-glob '{args.ae_glob}'")
    if not stats_files:
        raise SystemExit(f"FATAL: no files matched --stats-glob '{args.stats_glob}'")
    sys.stderr.write(f"[INFO] {len(ae_files)} ae.csv batch file(s), {len(stats_files)} stats batch file(s) found.\n")

    # --- Merge the ancestry-fraction tables column-wise. ---
    merged_ae = None
    seen_samples = set()
    key_cols = None
    for path in ae_files:
        df = pd.read_csv(path, index_col=0)
        this_key_cols = [c for c in df.columns if c == "Superpopulation"]
        sample_cols = [c for c in df.columns if c != "Superpopulation"]
        dupes = set(sample_cols) & seen_samples
        if dupes:
            raise SystemExit(f"FATAL: sample(s) {dupes} appear in more than one batch file -- check for overlapping batches.")
        seen_samples.update(sample_cols)

        if merged_ae is None:
            merged_ae = df
            key_cols = this_key_cols
        else:
            if not df.index.equals(merged_ae.index):
                raise SystemExit(
                    f"FATAL: {path} has a different population row set/order than earlier batches -- "
                    f"every batch must have been run against the exact same reference."
                )
            if not (df["Superpopulation"] == merged_ae["Superpopulation"]).all():
                raise SystemExit(f"FATAL: {path}'s Superpopulation column doesn't match earlier batches for the same population.")
            merged_ae = pd.concat([merged_ae, df[sample_cols]], axis=1)
        sys.stderr.write(f"[INFO] merged {path}: {len(sample_cols)} sample(s), running total {len(seen_samples)}.\n")

    ae_out_path = f"{args.out_prefix}_ae.csv"
    merged_ae.to_csv(ae_out_path)

    # --- Merge the per-sample stats tables row-wise. ---
    stats_frames = [pd.read_csv(path) for path in stats_files]
    merged_stats = pd.concat(stats_frames, ignore_index=True)
    if merged_stats["Sample"].duplicated().any():
        dupes = merged_stats.loc[merged_stats["Sample"].duplicated(), "Sample"].tolist()
        raise SystemExit(f"FATAL: duplicate sample(s) across stats batch files: {dupes}")

    stats_out_path = f"{args.out_prefix}_ae_stats.csv"
    merged_stats.to_csv(stats_out_path, index=False)

    sys.stderr.write(
        f"\n[DONE] Combined {len(seen_samples)} samples across {len(ae_files)} batch(es) into:\n"
        f"       {ae_out_path}\n"
        f"       {stats_out_path}\n"
    )


if __name__ == "__main__":
    main()
