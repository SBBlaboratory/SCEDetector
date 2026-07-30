#!/usr/bin/env python3
"""Detect SCE and recurrent translocation candidates from StrandPhaseR output."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

# SCE-valid class transitions (pattern matching table)
SCE_TRANSITIONS: dict[str, frozenset[str]] = {
    "CC": frozenset({"CW", "WC"}),
    "CW": frozenset({"CC", "WW"}),
    "WC": frozenset({"CC", "WW"}),
    "WW": frozenset({"CW", "WC"}),
}

INPUT_COLUMNS = ("chrom", "start", "end", "sample", "cell", "class")
OUTPUT_COLUMNS = ("Sample", "Cell_ID", "chr", "start", "Event", "Shared_cell_percent")


def is_sce_transition(prev_class: str, next_class: str) -> bool:
    """Return True if prev_class -> next_class matches an SCE pattern."""
    if prev_class == next_class:
        return False
    allowed = SCE_TRANSITIONS.get(prev_class)
    return allowed is not None and next_class in allowed


def _merge_consecutive_classes(
    classes: list[str], starts: list[int]
) -> list[tuple[str, int]]:
    """Collapse consecutive identical classes; keep start of each run."""
    if not classes:
        return []
    merged: list[tuple[str, int]] = [(classes[0], starts[0])]
    for cls, start in zip(classes[1:], starts[1:]):
        if cls != merged[-1][0]:
            merged.append((cls, start))
    return merged


def classify_recurrent_breakpoints(
    events: pd.DataFrame,
    source: pd.DataFrame,
    tolerance: int = 10_000,
    min_cell_fraction: float = 0.05,
) -> pd.DataFrame:
    """
    Relabel recurrent SCE-like breakpoints as translocation candidates.

    Within each sample and chromosome, an event is recurrent when breakpoints
    from at least 5% of all cells in that sample occur within +/- tolerance.
    At least two distinct cells are required to define a recurrent subclone.
    Translocation rows also carry the percentage of sample cells sharing the
    breakpoint; SCE rows leave that value empty.
    """
    if events.empty:
        events["Event"] = pd.Series(dtype="object")
        events["Shared_cell_percent"] = pd.Series(dtype="float")
        return events

    total_cells = source.groupby("sample")["cell"].nunique().to_dict()
    events = events.copy()
    events["Event"] = "SCE"
    events["Shared_cell_percent"] = pd.NA

    for sample, sample_events in events.groupby("Sample", sort=False):
        sample_cells = total_cells[sample]
        required_cells = max(2, math.ceil(sample_cells * min_cell_fraction))
        for _, chromosome_events in sample_events.groupby("chr", sort=False):
            starts = chromosome_events["start"]
            cells = chromosome_events["Cell_ID"]
            for index, breakpoint in starts.items():
                nearby_cells = cells[(starts - breakpoint).abs() <= tolerance]
                shared_cells = nearby_cells.nunique()
                if shared_cells >= required_cells:
                    events.at[index, "Event"] = "Translocation"
                    events.at[index, "Shared_cell_percent"] = round(
                        shared_cells / sample_cells * 100, 2
                    )

    return events


def detect_sce(
    df: pd.DataFrame,
    translocation_tolerance: int = 10_000,
    translocation_min_fraction: float = 0.05,
) -> pd.DataFrame:
    """
    Detect SCE breakpoints per (sample, cell, chrom).

    SCE requires exactly one state-class change along the chromosome:
    a valid SCE transition, after which the new class is maintained to
    the chromosome end (no further state changes).
    """
    missing = [c for c in INPUT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    work = df.copy()
    work["start"] = work["start"].astype(int)
    work["end"] = work["end"].astype(int)
    work["class"] = work["class"].astype(str).str.strip().str.upper()

    records: list[dict[str, object]] = []

    grouped = work.groupby(["sample", "cell", "chrom"], sort=False)
    for (sample, cell, chrom), group in grouped:
        segments = group.sort_values(["start", "end"], kind="mergesort")
        merged = _merge_consecutive_classes(
            segments["class"].tolist(),
            segments["start"].tolist(),
        )

        # Exactly one class change, then held to chr end => two runs only
        if len(merged) != 2:
            continue

        prev_class, _ = merged[0]
        next_class, breakpoint = merged[1]
        if is_sce_transition(prev_class, next_class):
            records.append(
                {
                    "Sample": sample,
                    "Cell_ID": cell,
                    "chr": chrom,
                    "start": breakpoint,
                }
            )

    result = pd.DataFrame.from_records(records, columns=list(OUTPUT_COLUMNS[:-2]))
    result = classify_recurrent_breakpoints(
        result,
        work,
        tolerance=translocation_tolerance,
        min_cell_fraction=translocation_min_fraction,
    )
    if result.empty:
        return result

    return result.sort_values(
        ["Sample", "Cell_ID", "chr", "start"], kind="mergesort"
    ).reset_index(drop=True)


def load_strandphaser_output(path: Path) -> pd.DataFrame:
    """Load StrandPhaseR-style TSV (chrom, start, end, sample, cell, class)."""
    return pd.read_csv(path, sep="\t")


def write_excel(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(path, index=False, engine="openpyxl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Detect SCE and recurrent translocation candidates from StrandPhaseR "
            "final output and write an Excel summary "
            "(Sample, Cell_ID, chr, start, Event)."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Input TSV (e.g. StrandPhaseR_final_output.txt)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("SCE_detected.xlsx"),
        help="Output Excel path (default: SCE_detected.xlsx)",
    )
    parser.add_argument(
        "--translocation-tolerance",
        type=int,
        default=10_000,
        help="Breakpoint tolerance in bp (default: 10000, i.e. +/-10 kb)",
    )
    parser.add_argument(
        "--translocation-min-fraction",
        type=float,
        default=0.05,
        help="Minimum fraction of sample cells sharing a breakpoint (default: 0.05)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.translocation_tolerance < 0:
        raise ValueError("--translocation-tolerance must be >= 0")
    if not 0 < args.translocation_min_fraction <= 1:
        raise ValueError("--translocation-min-fraction must be in (0, 1]")

    df = load_strandphaser_output(args.input)
    events = detect_sce(
        df,
        translocation_tolerance=args.translocation_tolerance,
        translocation_min_fraction=args.translocation_min_fraction,
    )
    write_excel(events, args.output)
    counts = events["Event"].value_counts()
    print(f"Detected {len(events)} breakpoint candidate(s)")
    print(f"SCE: {counts.get('SCE', 0)}")
    print(f"Translocation: {counts.get('Translocation', 0)}")
    print(f"Wrote: {args.output.resolve()}")


if __name__ == "__main__":
    main()