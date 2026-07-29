#!/usr/bin/env python3
"""SCEDetector: detect Sister Chromatid Exchange breakpoints from StrandPhaseR output."""

from __future__ import annotations

import argparse
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
OUTPUT_COLUMNS = ("Sample", "Cell_ID", "chr", "start")


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


def detect_sce(df: pd.DataFrame) -> pd.DataFrame:
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

    result = pd.DataFrame.from_records(records, columns=list(OUTPUT_COLUMNS))
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
            "Detect SCE breakpoints from StrandPhaseR final output and write "
            "an Excel summary (Sample, Cell_ID, chr, start)."
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = load_strandphaser_output(args.input)
    sce = detect_sce(df)
    write_excel(sce, args.output)
    print(f"Detected {len(sce)} SCE breakpoint(s)")
    print(f"Wrote: {args.output.resolve()}")


if __name__ == "__main__":
    main()