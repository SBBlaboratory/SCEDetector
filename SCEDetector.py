#!/usr/bin/env python3
"""Detect SCE and recurrent translocation candidates from 200 kb bin strand states."""

from __future__ import annotations

import argparse
import gzip
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd

# Valid single-switch SCE transitions on filtered 200 kb bins
SCE_TRANSITIONS: dict[str, frozenset[str]] = {
    "CC": frozenset({"WC"}),
    "WW": frozenset({"WC"}),
    "WC": frozenset({"WW", "CC"}),
}
# Male chrX is hemizygous: drop WC segments and allow only WW↔CC switches.
MALE_X_SCE_TRANSITIONS: dict[str, frozenset[str]] = {
    "CC": frozenset({"WW"}),
    "WW": frozenset({"CC"}),
}

VALID_CLASSES = frozenset(SCE_TRANSITIONS)
MALE_X_VALID_CLASSES = frozenset(MALE_X_SCE_TRANSITIONS)
OUTPUT_COLUMNS = (
    "Sample",
    "Cell_ID",
    "chr",
    "start",
    "end",
    "Event",
    "Shared_cell_percent",
)
DEFAULT_MASK = Path(__file__).resolve().parent / "HGSVC.200000.txt"
DEFAULT_HUMAN_ARM_POSITIONS = (
    Path(__file__).resolve().parent / "chromosome_arm_positions_grch38.txt"
)
# Species → default chromosome arm / centromere table (p-end .. q-start = centro).
SPECIES_ARM_POSITIONS: dict[str, Path] = {
    "human": DEFAULT_HUMAN_ARM_POSITIONS,
}
# Good islands shorter than this (in 200 kb bins) and flanked by None are
# treated as None (e.g. sparse mappability speckles inside chr9 heterochromatin).
DEFAULT_MAX_SPARSE_GOOD_BINS = 5
# Large None runs (>= this many bins) also absorb up to DEFAULT_MAX_SPARSE_GOOD_BINS
# flanking good bins so edge stubs next to heterochromatin are skipped.
DEFAULT_MIN_NONE_RUN_FOR_EDGE = 10
# Fallback when a chrom has no arm-table centro: contiguous None this long (bp).
DEFAULT_MIN_CENTROMERE_NONE_BP = 1_000_000
# Merge adjacent large-None runs separated by at most this gap into one barrier.
DEFAULT_CENTROMERE_MERGE_GAP_BP = 1_000_000
# Minimum flank length on each side of a centromere barrier for a
# cross-centromere SCE exception call.
DEFAULT_MIN_CENTRO_FLANK_BP = 5_000_000
# After SV skip, drop short edge stubs that abut an SV hole (e.g. idup)
# so flanking states can merge instead of forming false SCE sandwiches.
DEFAULT_MAX_SV_EDGE_STUB_BINS = 2
DEFAULT_BIN_BP = 200_000
# Drop WC as duplication-like when (vs local homozygous flank dominant depth
# ``dom``): max(c,w)/dom, min(c,w)/dom, and (c+w)/dom all clear thresholds.
# True WC ≈ both strands ~0.5×dom and tot~dom; unmarked dup often keeps one
# strand near full (max~dom) while adding the other (min≳0.35×dom) so tot≳1.35×dom.
DEFAULT_WC_DUP_MAX_STRAND_RATIO = 0.85
DEFAULT_WC_DUP_MIN_STRAND_RATIO = 0.35
DEFAULT_WC_DUP_TOT_RATIO = 1.35
DEFAULT_WC_DUP_FLANK_BINS = 25
# Tip slack: SV that overlaps a tip window (±1 Mb) is tip-reaching.
DEFAULT_TELOMERE_SV_SLACK_BP = 1_000_000
# Tip-직전: SV ending within this gap of qter (or starting within this of pter)
# is also tip-linked for inv/complex keep (covers e.g. complex stopping ~1–2 Mb
# short of chrom_end with a tiny uncovered tip stub).
DEFAULT_TELOMERE_PRE_TIP_GAP_BP = 2_000_000
ACROCENTRIC_CHROMS = frozenset(
    {f"chr{i}" for i in (13, 14, 15, 21, 22)}
)


def is_sce_transition(
    prev_class: str, next_class: str, *, male_chrx: bool = False
) -> bool:
    """Return True if prev_class -> next_class matches an SCE pattern."""
    if prev_class == next_class:
        return False
    table = MALE_X_SCE_TRANSITIONS if male_chrx else SCE_TRANSITIONS
    allowed = table.get(prev_class)
    return allowed is not None and next_class in allowed


def _is_male_chrx(chrom: str, sex: str) -> bool:
    return chrom == "chrX" and sex == "male"


def _normalize_chrom(chrom: object) -> str:
    text = str(chrom).strip()
    if text.startswith("chr"):
        return text
    return f"chr{text}"


def load_centromere_barriers_from_arms(
    path: Path,
) -> dict[str, list[tuple[int, int]]]:
    """
    Build per-chromosome centromere intervals from p/q arm coordinates.

    Centromere = ``[p_arm.End, q_arm.Start)`` (gap between p and q).
    """
    arms = pd.read_csv(path, sep="\t")
    required = {"Chrom", "Start", "End"}
    missing = required - set(arms.columns)
    if missing:
        raise ValueError(f"Arm positions file missing columns: {sorted(missing)}")

    by_chrom: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    for row in arms.itertuples(index=False):
        chrom = _normalize_chrom(getattr(row, "Chrom"))
        start = int(row.Start)
        end = int(row.End)
        idf = str(getattr(row, "Idf", "") or "")
        arm = idf[-1].lower() if idf else ""
        if arm not in {"p", "q"}:
            # Infer: arm starting at 0 is p, otherwise q.
            arm = "p" if start <= 0 else "q"
        by_chrom[chrom][arm] = (start, end)

    centros: dict[str, list[tuple[int, int]]] = {}
    for chrom, parts in by_chrom.items():
        if "p" not in parts or "q" not in parts:
            continue
        _, p_end = parts["p"]
        q_start, _ = parts["q"]
        if q_start <= p_end:
            continue
        centros[chrom] = [(p_end, q_start)]
    return centros


def resolve_species_arm_path(
    species: str, arm_positions: Path | None = None
) -> Path:
    """Return arm-position table path for ``species`` (or explicit override)."""
    if arm_positions is not None:
        return arm_positions
    key = species.strip().lower()
    if key not in SPECIES_ARM_POSITIONS:
        known = ", ".join(sorted(SPECIES_ARM_POSITIONS))
        raise ValueError(f"Unknown species {species!r}; known: {known}")
    return SPECIES_ARM_POSITIONS[key]


def _large_none_barriers(
    none_bins: set[tuple[str, int, int]],
    chrom: str,
    min_run_bp: int = DEFAULT_MIN_CENTROMERE_NONE_BP,
    merge_gap_bp: int = DEFAULT_CENTROMERE_MERGE_GAP_BP,
) -> list[tuple[int, int]]:
    """Fallback centro/arm barriers from large HGSVC None runs (e.g. chrX)."""
    runs: list[tuple[int, int]] = []
    for start, end in sorted(
        (s, e) for c, s, e in none_bins if c == chrom
    ):
        if runs and start <= runs[-1][1] + merge_gap_bp:
            runs[-1] = (runs[-1][0], max(runs[-1][1], end))
        else:
            runs.append((start, end))

    barriers: list[tuple[int, int]] = []
    for start, end in runs:
        length = end - start
        if length < min_run_bp:
            continue
        if (
            chrom not in ACROCENTRIC_CHROMS
            and start < 10_000_000
            and length < 20_000_000
        ):
            continue
        barriers.append((start, end))
    return barriers


def _centromere_barriers_for_chrom(
    chrom: str,
    arm_centros: dict[str, list[tuple[int, int]]],
    none_bins: set[tuple[str, int, int]],
) -> list[tuple[int, int]]:
    """Prefer arm-table centro; else fall back to large-None barriers."""
    barriers = arm_centros.get(chrom)
    if barriers:
        return barriers
    return _large_none_barriers(none_bins, chrom)


def _interval_overlaps_barrier(
    start: int, end: int, barriers: list[tuple[int, int]]
) -> bool:
    """Return True if [start, end) overlaps any barrier interval."""
    for b0, b1 in barriers:
        if not (end <= b0 or start >= b1):
            return True
    return False


def _split_bins_into_arms(
    classes: list[str],
    starts: list[int],
    ends: list[int],
    barriers: list[tuple[int, int]],
) -> list[tuple[list[str], list[int], list[int]]]:
    """
    Split kept bins into chromosome arms at large-None barriers.

    Small None/SV holes do not split arms; only gaps that overlap a large
    centromere-scale None barrier do.
    """
    if not classes:
        return []
    if not barriers:
        return [(classes, starts, ends)]

    arms: list[tuple[list[str], list[int], list[int]]] = []
    cur_c = [classes[0]]
    cur_s = [starts[0]]
    cur_e = [ends[0]]
    for cls, start, end, prev_end in zip(
        classes[1:], starts[1:], ends[1:], ends[:-1]
    ):
        if start > prev_end and _interval_overlaps_barrier(
            prev_end, start, barriers
        ):
            arms.append((cur_c, cur_s, cur_e))
            cur_c, cur_s, cur_e = [cls], [start], [end]
        else:
            cur_c.append(cls)
            cur_s.append(start)
            cur_e.append(end)
    arms.append((cur_c, cur_s, cur_e))
    return arms


def _call_events_on_arm(
    sample: object,
    cell: object,
    chrom: object,
    classes: list[str],
    starts: list[int],
    ends: list[int],
    single_records: list[dict[str, object]],
    raw_starts: list[int] | None = None,
    raw_ends: list[int] | None = None,
    raw_classes: list[str] | None = None,
    skip_sv: list[tuple[int, int]] | None = None,
    *,
    male_chrx: bool = False,
) -> None:
    """Call a single-switch SCE on one arm (after A-B-A cleanup)."""
    runs = _merge_runs(classes, starts, ends)
    if len(runs) != 2:
        return
    left_cls, _, left_end = runs[0]
    right_cls, right_start, _ = runs[1]
    if not is_sce_transition(left_cls, right_cls, male_chrx=male_chrx):
        return
    breakpoint = right_start
    if raw_starts is not None and raw_ends is not None and raw_classes is not None:
        breakpoint = _refine_sce_breakpoint_for_sv(
            breakpoint,
            left_cls,
            right_cls,
            left_end,
            right_start,
            raw_starts,
            raw_ends,
            raw_classes,
            skip_sv or [],
            male_chrx=male_chrx,
        )
    single_records.append(
        {
            "Sample": sample,
            "Cell_ID": cell,
            "chr": chrom,
            "start": breakpoint,
            "end": pd.NA,
            "Event": "SCE",
            "Shared_cell_percent": pd.NA,
        }
    )


def _sv_touches_interval(
    start: int, end: int, skip_sv: list[tuple[int, int]]
) -> bool:
    """True if [start, end) overlaps an SV or an endpoint lies on an SV edge."""
    if not skip_sv:
        return False
    if end < start:
        start, end = end, start
    if end == start:
        end = start + 1
    if _overlaps_any(start, end, skip_sv):
        return True
    for a, b in skip_sv:
        if a <= start <= b or a <= end <= b:
            return True
        if start == b or end == a:
            return True
    return False


def _refine_sce_breakpoint_for_sv(
    bp: int,
    left_cls: str,
    right_cls: str,
    left_end: int,
    right_start: int,
    raw_starts: list[int],
    raw_ends: list[int],
    raw_classes: list[str],
    skip_sv: list[tuple[int, int]],
    *,
    male_chrx: bool = False,
) -> int:
    """
    If the filtered SCE breakpoint (or the gap that created it) touches an SV
    interval, replace it with the raw left_cls→right_cls switch in that gap.
    Judgment that an SCE exists is unchanged; only the coordinate is corrected.
    """
    if not skip_sv:
        return bp
    gap_lo = min(left_end, right_start)
    gap_hi = max(left_end, right_start)
    if not (
        _sv_touches_interval(gap_lo, gap_hi, skip_sv)
        or _sv_touches_interval(bp, bp + 1, skip_sv)
    ):
        return bp
    # Peek raw without SV skip so the true switch inside the SV remains visible.
    raw_bp = _raw_sce_breakpoint_in_window(
        raw_starts,
        raw_ends,
        raw_classes,
        gap_lo,
        gap_hi if gap_hi > gap_lo else bp + 1,
        left_cls,
        right_cls,
        skip_sv=[],
        seed_prev_cls=left_cls,
        male_chrx=male_chrx,
    )
    return raw_bp if raw_bp is not None else bp


def _barrier_flank_runs(
    runs: list[tuple[str, int, int]],
    barrier: tuple[int, int],
) -> tuple[tuple[str, int, int] | None, tuple[str, int, int] | None]:
    """Return (left_run, right_run) flanking a centromere barrier."""
    b0, b1 = barrier
    left: tuple[str, int, int] | None = None
    for cls, start, end in runs:
        if end <= b0:
            left = (cls, start, end)
    right: tuple[str, int, int] | None = None
    for cls, start, end in runs:
        if start >= b1:
            right = (cls, start, end)
            break
    if right is None:
        # Unmasked island inside the barrier: use first kept run past b0.
        for cls, start, end in runs:
            if start > b0:
                right = (cls, start, end)
                break
    return left, right


def _raw_sce_breakpoint_in_window(
    starts: list[int],
    ends: list[int],
    classes: list[str],
    window_start: int,
    window_end: int,
    left_cls: str,
    right_cls: str,
    skip_sv: list[tuple[int, int]],
    seed_prev_cls: str | None = None,
    *,
    male_chrx: bool = False,
) -> int | None:
    """
    Peek raw (including None-masked) bins inside ``[window_start, window_end)``
    for a single valid SCE transition matching left_cls → right_cls.

    ``seed_prev_cls`` (usually ``left_cls``) provides the class immediately
    before the window when the left flank ends exactly at ``window_start``.
    """
    allowed = MALE_X_VALID_CLASSES if male_chrx else VALID_CLASSES
    prev_cls: str | None = seed_prev_cls
    hits: list[tuple[int, str, str]] = []
    for start, end, cls in zip(starts, ends, classes):
        if end <= window_start or start >= window_end:
            continue
        if skip_sv and _overlaps_any(start, end, skip_sv):
            continue
        if cls not in allowed:
            continue
        if prev_cls is not None and cls != prev_cls:
            if is_sce_transition(prev_cls, cls, male_chrx=male_chrx):
                hits.append((start, prev_cls, cls))
        prev_cls = cls

    if len(hits) != 1:
        return None
    bp, a, b = hits[0]
    if a == left_cls and b == right_cls:
        return bp
    return None


def _call_centromere_crossing_sce(
    sample: object,
    cell: object,
    chrom: object,
    kept_classes: list[str],
    kept_starts: list[int],
    kept_ends: list[int],
    raw_starts: list[int],
    raw_ends: list[int],
    raw_classes: list[str],
    barriers: list[tuple[int, int]],
    skip_sv: list[tuple[int, int]],
    single_records: list[dict[str, object]],
    min_flank_bp: int = DEFAULT_MIN_CENTRO_FLANK_BP,
    *,
    male_chrx: bool = False,
) -> None:
    """
    Exception: SCE across a centromere-scale None barrier.

    When kept flanks on either side of a barrier differ by a valid SCE
    transition and each flank is long enough:
    1. Prefer the raw breakpoint inside the barrier if there is exactly one
       matching switch (clean pericentromeric SCE).
    2. Otherwise fall back to the start of the right kept flank (None-masked
       path is already a single valid switch across the gap).
    """
    if not barriers or not kept_classes:
        return
    runs = _merge_runs(kept_classes, kept_starts, kept_ends)
    for barrier in barriers:
        b0, b1 = barrier
        left, right = _barrier_flank_runs(runs, barrier)
        if left is None or right is None:
            continue
        left_cls, left_start, left_end = left
        right_cls, right_start, right_end = right
        if left_cls == right_cls:
            continue
        if not is_sce_transition(left_cls, right_cls, male_chrx=male_chrx):
            continue
        if (left_end - left_start) < min_flank_bp:
            continue
        if (right_end - right_start) < min_flank_bp:
            continue
        bp = _raw_sce_breakpoint_in_window(
            raw_starts,
            raw_ends,
            raw_classes,
            b0,
            b1,
            left_cls,
            right_cls,
            skip_sv,
            male_chrx=male_chrx,
        )
        if bp is None:
            bp = right_start
        bp = _refine_sce_breakpoint_for_sv(
            bp,
            left_cls,
            right_cls,
            left_end,
            right_start,
            raw_starts,
            raw_ends,
            raw_classes,
            skip_sv,
            male_chrx=male_chrx,
        )
        single_records.append(
            {
                "Sample": sample,
                "Cell_ID": cell,
                "chr": chrom,
                "start": bp,
                "end": pd.NA,
                "Event": "SCE",
                "Shared_cell_percent": pd.NA,
            }
        )


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


def _merge_runs(
    classes: list[str], starts: list[int], ends: list[int]
) -> list[tuple[str, int, int]]:
    """Collapse consecutive identical classes into (class, start, end) runs."""
    if not classes:
        return []
    runs: list[tuple[str, int, int]] = [(classes[0], starts[0], ends[0])]
    for cls, start, end in zip(classes[1:], starts[1:], ends[1:]):
        prev_cls, prev_start, _ = runs[-1]
        if cls == prev_cls:
            runs[-1] = (prev_cls, prev_start, end)
        else:
            runs.append((cls, start, end))
    return runs


def _run_touches_sv(
    start: int, end: int, skip_sv: list[tuple[int, int]]
) -> bool:
    """True if [start, end) overlaps or abuts an SV skip interval."""
    for a, b in skip_sv:
        if not (end <= a or start >= b):
            return True
        if end == a or start == b:
            return True
    return False


def _drop_sv_adjacent_stubs(
    classes: list[str],
    starts: list[int],
    ends: list[int],
    skip_sv: list[tuple[int, int]],
    max_stub_bins: int = DEFAULT_MAX_SV_EDGE_STUB_BINS,
    bin_bp: int = DEFAULT_BIN_BP,
) -> tuple[list[str], list[int], list[int]]:
    """
    Remove short state runs that touch an SV skip hole.

    After SV bins are skipped, a 1–2 bin stub often remains at the hole edge
    and creates false A-B-A / SCE calls. Drop those stubs so flanking identical
    classes merge across the gap.
    """
    if not skip_sv or not classes or max_stub_bins <= 0:
        return classes, starts, ends

    max_stub_bp = max_stub_bins * bin_bp
    while True:
        runs = _merge_runs(classes, starts, ends)
        drop_start: int | None = None
        drop_end: int | None = None
        for _, start, end in runs:
            if end - start > max_stub_bp:
                continue
            if not _run_touches_sv(start, end, skip_sv):
                continue
            drop_start, drop_end = start, end
            break
        if drop_start is None or drop_end is None:
            break
        kept = [
            (cls, start, end)
            for cls, start, end in zip(classes, starts, ends)
            if end <= drop_start or start >= drop_end
        ]
        if len(kept) == len(classes):
            break
        classes = [item[0] for item in kept]
        starts = [item[1] for item in kept]
        ends = [item[2] for item in kept]

    return classes, starts, ends


def _cell_homozygous_depth_refs(
    work: pd.DataFrame,
    none_bins: set[tuple[str, int, int]],
) -> dict[tuple[str, str], float]:
    """
    Per-cell homozygous depth baseline for WC duplication filtering.

    For each ``(sample, cell)``, take the median dominant-strand depth over
    non-None ``CC`` (``c``) and ``WW`` (``w``) bins genome-wide (chrY excluded).
    """
    if "c" not in work.columns or "w" not in work.columns:
        return {}
    homo = work[
        work["class"].isin(["CC", "WW"]) & (work["chrom"].astype(str) != "chrY")
    ]
    if homo.empty:
        return {}

    refs: dict[tuple[str, str], float] = {}
    for (sample, cell), group in homo.groupby(["sample", "cell"], sort=False):
        vals: list[float] = []
        for chrom, start, end, cls, c_val, w_val in zip(
            group["chrom"].astype(str),
            group["start"].astype(int),
            group["end"].astype(int),
            group["class"].astype(str),
            group["c"].astype(float),
            group["w"].astype(float),
        ):
            if (chrom, start, end) in none_bins:
                continue
            dom = w_val if cls == "WW" else c_val
            if dom > 0:
                vals.append(dom)
        if vals:
            vals.sort()
            refs[(str(sample), str(cell))] = vals[len(vals) // 2]
    return refs


def _contiguous_raw_wc_depth(
    wc_start: int,
    wc_end: int,
    raw_starts: list[int],
    raw_ends: list[int],
    raw_classes: list[str],
    raw_c: list[float],
    raw_w: list[float],
) -> tuple[float, float] | None:
    """
    Mean ``c``/``w`` over the full contiguous raw WC run covering ``[wc_start, wc_end)``.

    Includes bins that may later be skipped by the None mask so low-mappability
    tips of a true WC segment still contribute to the duplication check.
    """
    n = len(raw_classes)
    if n == 0 or len(raw_c) != n or len(raw_w) != n:
        return None
    seed: int | None = None
    for idx, (start, end, cls) in enumerate(
        zip(raw_starts, raw_ends, raw_classes)
    ):
        if cls != "WC":
            continue
        if not (end <= wc_start or start >= wc_end):
            seed = idx
            break
    if seed is None:
        return None

    lo = seed
    while lo > 0 and raw_classes[lo - 1] == "WC" and raw_ends[lo - 1] == raw_starts[lo]:
        lo -= 1
    hi = seed
    while (
        hi + 1 < n
        and raw_classes[hi + 1] == "WC"
        and raw_ends[hi] == raw_starts[hi + 1]
    ):
        hi += 1

    count = hi - lo + 1
    mean_c = sum(raw_c[lo : hi + 1]) / count
    mean_w = sum(raw_w[lo : hi + 1]) / count
    return mean_c, mean_w


def _homozygous_flank_dom(
    side: str,
    wc_start: int,
    wc_end: int,
    raw_starts: list[int],
    raw_ends: list[int],
    raw_classes: list[str],
    raw_c: list[float],
    raw_w: list[float],
    *,
    max_bins: int,
) -> float | None:
    """
    Mean dominant-strand depth of the homozygous run abutting a WC interval.

    ``side`` ``"L"`` uses bins immediately before ``wc_start``; ``"R"`` uses
    bins immediately after ``wc_end``. WW → ``w``, CC → ``c``.
    """
    n = len(raw_classes)
    if n == 0 or max_bins <= 0:
        return None
    if side == "L":
        j = None
        for idx, end in enumerate(raw_ends):
            if end == wc_start:
                j = idx
                break
        if j is None:
            return None
        cls = raw_classes[j]
        if cls not in {"WW", "CC"}:
            return None
        vals: list[float] = []
        while j >= 0 and raw_classes[j] == cls and len(vals) < max_bins:
            vals.append(raw_w[j] if cls == "WW" else raw_c[j])
            if j == 0 or raw_ends[j - 1] != raw_starts[j]:
                break
            j -= 1
    else:
        j = None
        for idx, start in enumerate(raw_starts):
            if start == wc_end:
                j = idx
                break
        if j is None:
            return None
        cls = raw_classes[j]
        if cls not in {"WW", "CC"}:
            return None
        vals = []
        while j < n and raw_classes[j] == cls and len(vals) < max_bins:
            vals.append(raw_w[j] if cls == "WW" else raw_c[j])
            if j + 1 >= n or raw_ends[j] != raw_starts[j + 1]:
                break
            j += 1
    if not vals:
        return None
    return sum(vals) / len(vals)


def _is_duplication_like_wc_depth(
    mean_c: float,
    mean_w: float,
    dom: float,
    *,
    max_strand_ratio: float,
    min_strand_ratio: float,
    tot_ratio: float,
) -> bool:
    """True if WC depths match the keep-one-strand / add-other-strand dup pattern."""
    if dom <= 0:
        return False
    hi = max(mean_c, mean_w) / dom
    lo = min(mean_c, mean_w) / dom
    tot = (mean_c + mean_w) / dom
    return (
        hi >= max_strand_ratio
        and lo >= min_strand_ratio
        and tot >= tot_ratio
    )


def _drop_duplication_like_wc(
    classes: list[str],
    starts: list[int],
    ends: list[int],
    *,
    ref_depth: float,
    raw_starts: list[int],
    raw_ends: list[int],
    raw_classes: list[str],
    raw_c: list[float],
    raw_w: list[float],
    max_strand_ratio: float = DEFAULT_WC_DUP_MAX_STRAND_RATIO,
    min_strand_ratio: float = DEFAULT_WC_DUP_MIN_STRAND_RATIO,
    tot_ratio: float = DEFAULT_WC_DUP_TOT_RATIO,
    flank_bins: int = DEFAULT_WC_DUP_FLANK_BINS,
) -> tuple[list[str], list[int], list[int]]:
    """
    Drop WC runs that look like unmarked duplication vs local homozygous flanks.

    WC depth uses the full contiguous raw WC run (including None-masked bins).
    Dominant flank depth ``dom`` is the mean WW-``w`` / CC-``c`` of abutting
    homozygous runs (up to ``flank_bins`` each side; average if both sides),
    floored by the cell-wide homozygous median ``ref_depth`` so shallow
    peri-centromere / noisy flanks do not inflate ratios. If no homozygous
    flank exists, use ``ref_depth`` alone.

    Drop when all hold: ``max(c,w)/dom >= max_strand_ratio``,
    ``min(c,w)/dom >= min_strand_ratio``, ``(c+w)/dom >= tot_ratio``.
    """
    if (
        not classes
        or tot_ratio <= 0
        or max_strand_ratio <= 0
        or min_strand_ratio <= 0
        or flank_bins <= 0
    ):
        return classes, starts, ends

    drop_idx: set[int] = set()
    i = 0
    n = len(classes)
    while i < n:
        if classes[i] != "WC":
            i += 1
            continue
        j = i + 1
        while j < n and classes[j] == "WC":
            j += 1
        depth = _contiguous_raw_wc_depth(
            starts[i],
            ends[j - 1],
            raw_starts,
            raw_ends,
            raw_classes,
            raw_c,
            raw_w,
        )
        if depth is not None:
            mean_c, mean_w = depth
            flank_doms = [
                d
                for d in (
                    _homozygous_flank_dom(
                        "L",
                        starts[i],
                        ends[j - 1],
                        raw_starts,
                        raw_ends,
                        raw_classes,
                        raw_c,
                        raw_w,
                        max_bins=flank_bins,
                    ),
                    _homozygous_flank_dom(
                        "R",
                        starts[i],
                        ends[j - 1],
                        raw_starts,
                        raw_ends,
                        raw_classes,
                        raw_c,
                        raw_w,
                        max_bins=flank_bins,
                    ),
                )
                if d is not None and d > 0
            ]
            if flank_doms:
                flank_mean = sum(flank_doms) / len(flank_doms)
                # Floor by cell ref: shallow local flanks (centro/tip noise)
                # otherwise inflate max/min/tot and falsely drop true WC.
                dom = max(flank_mean, ref_depth) if ref_depth > 0 else flank_mean
            else:
                dom = ref_depth
            if _is_duplication_like_wc_depth(
                mean_c,
                mean_w,
                dom,
                max_strand_ratio=max_strand_ratio,
                min_strand_ratio=min_strand_ratio,
                tot_ratio=tot_ratio,
            ):
                drop_idx.update(range(i, j))
        i = j

    if not drop_idx:
        return classes, starts, ends
    kept = [
        (cls, start, end)
        for idx, (cls, start, end) in enumerate(zip(classes, starts, ends))
        if idx not in drop_idx
    ]
    return (
        [item[0] for item in kept],
        [item[1] for item in kept],
        [item[2] for item in kept],
    )


def _overlaps_any(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    """Return True if [start, end) overlaps any [a, b) interval."""
    for a, b in intervals:
        if not (end <= a or start >= b):
            return True
    return False


def _is_inversion_call(name: object) -> bool:
    return isinstance(name, str) and "inv" in name.lower()


def _is_complex_call(name: object) -> bool:
    return isinstance(name, str) and name.lower().startswith("complex")


def _is_deletion_call(name: object) -> bool:
    return isinstance(name, str) and name.lower().startswith("del")


def _is_duplication_call(name: object) -> bool:
    """True for dup_* / idup_* copy-number gain calls."""
    if not isinstance(name, str):
        return False
    lower = name.lower()
    return lower.startswith("dup") or lower.startswith("idup")


def _is_cn_skip_sv(name: object) -> bool:
    """Deletion / duplication intervals are always used as skip masks."""
    return _is_deletion_call(name) or _is_duplication_call(name)


def _is_inv_or_complex_call(name: object) -> bool:
    return _is_inversion_call(name) or _is_complex_call(name)


def _reaches_chromosome_telomere(
    start: int,
    end: int,
    chrom_end: int | None,
    slack_bp: int = DEFAULT_TELOMERE_SV_SLACK_BP,
) -> bool:
    """True if an SV interval includes/abuts a chromosome tip (±slack)."""
    if start <= slack_bp:
        return True
    if chrom_end is not None and end >= chrom_end - slack_bp:
        return True
    return False


def _is_tip_linked_anchor(
    start: int,
    end: int,
    chrom_end: int | None,
    *,
    tip_slack_bp: int = DEFAULT_TELOMERE_SV_SLACK_BP,
    pre_tip_gap_bp: int = DEFAULT_TELOMERE_PRE_TIP_GAP_BP,
) -> bool:
    """
    True if an SV is tip-reaching or tip-직전.

    Tip-reaching: overlaps pter/qter within ``tip_slack_bp``.
    Tip-직전: starts within ``pre_tip_gap_bp`` of pter, or ends within
    ``pre_tip_gap_bp`` of ``chrom_end`` (small uncovered tip stub after SV).
    """
    if _reaches_chromosome_telomere(start, end, chrom_end, tip_slack_bp):
        return True
    if start <= pre_tip_gap_bp:
        return True
    if chrom_end is not None and 0 <= chrom_end - end <= pre_tip_gap_bp:
        return True
    return False


def _intervals_touch(a0: int, a1: int, b0: int, b1: int) -> bool:
    """True if half-open intervals overlap or abut."""
    if a0 < b1 and b0 < a1:
        return True
    return a1 == b0 or b1 == a0


def _tip_linked_inv_complex(
    inv_complex: list[tuple[int, int]],
    cn_skips: list[tuple[int, int]],
    chrom_end: int | None,
) -> set[int]:
    """
    Indices of inv/complex intervals that stay for SCE judgment.

    Keep when the interval itself is tip-reaching or tip-직전, or when it
    abuts (transitively) a tip-linked del/dup/idup or another kept
    inv/complex.
    """
    keep = {
        i
        for i, (start, end) in enumerate(inv_complex)
        if _is_tip_linked_anchor(start, end, chrom_end)
    }
    tip_cn = [
        (start, end)
        for start, end in cn_skips
        if _is_tip_linked_anchor(start, end, chrom_end)
    ]
    changed = True
    while changed:
        changed = False
        for i, (start, end) in enumerate(inv_complex):
            if i in keep:
                continue
            if any(_intervals_touch(start, end, ts, te) for ts, te in tip_cn):
                keep.add(i)
                changed = True
                continue
            if any(
                j in keep
                and _intervals_touch(start, end, inv_complex[j][0], inv_complex[j][1])
                for j in range(len(inv_complex))
            ):
                keep.add(i)
                changed = True
    return keep


def _is_two_step_opposite(
    left_cls: str, mid_cls: str, right_cls: str, *, male_chrx: bool = False
) -> bool:
    """WW→WC→CC or CC→WC→WW (two SCE steps to the opposite homozygous state)."""
    if male_chrx:
        return False
    if mid_cls != "WC" or left_cls == right_cls:
        return False
    if {left_cls, right_cls} != {"WW", "CC"}:
        return False
    return is_sce_transition(left_cls, mid_cls) and is_sce_transition(
        mid_cls, right_cls
    )


def _extract_aba_sandwiches(
    classes: list[str],
    starts: list[int],
    ends: list[int],
    barriers: list[tuple[int, int]] | None = None,
    *,
    male_chrx: bool = False,
) -> tuple[list[str], list[int], list[int], list[tuple[int, int, bool, str]]]:
    """
    Extract A-B-A and two-step-opposite sandwiches.

    Returns duals as (left_bp, right_bp, spans_centromere, kind) where kind is
    ``"aba"`` or ``"twostep"``.

    - ``aba``: same-flank sandwich. ``spans_centromere`` is retained for
      diagnostics but resolution is by subclone recurrence only (see
      ``resolve_double_switch_patterns``).
    - ``twostep``: WW→WC→CC or CC→WC→WW. Always two SCEs (not used on male
      chrX, where WC is removed and only WW↔CC applies).

    Invalid flank-matched sandwiches and other A-B-C middles are dropped
    without emitting a dual.
    """
    barriers = barriers or []
    duals: list[tuple[int, int, bool, str]] = []

    while True:
        runs = _merge_runs(classes, starts, ends)
        if len(runs) < 3:
            break

        drop_start: int | None = None
        drop_end: int | None = None
        # For two-step, drop left+middle (keep final state run).
        drop_through_right_start: int | None = None

        for i in range(len(runs) - 2):
            left_cls, left_start, _ = runs[i]
            mid_cls, mid_start, mid_end = runs[i + 1]
            right_cls, right_start, right_end = runs[i + 2]

            spans = _interval_overlaps_barrier(mid_start, right_start, barriers)

            if (
                left_cls == right_cls
                and is_sce_transition(left_cls, mid_cls, male_chrx=male_chrx)
                and is_sce_transition(mid_cls, right_cls, male_chrx=male_chrx)
            ):
                duals.append((mid_start, right_start, spans, "aba"))
                drop_start, drop_end = mid_start, mid_end
                break

            if _is_two_step_opposite(
                left_cls, mid_cls, right_cls, male_chrx=male_chrx
            ):
                # Always double SCE; centro in the final run does not make inversion.
                duals.append((mid_start, right_start, False, "twostep"))
                # Remove left+mid; keep the final homozygous run.
                drop_through_right_start = right_start
                drop_start, drop_end = left_start, right_start
                break

            if left_cls == right_cls:
                # Flank-matched but not a valid SCE sandwich: clear path only.
                drop_start, drop_end = mid_start, mid_end
                break

        if drop_start is None:
            if len(runs) == 3 and runs[0][0] != runs[2][0]:
                _, drop_start, drop_end = runs[1]
            else:
                break

        if drop_start is None or drop_end is None:
            break

        if drop_through_right_start is not None:
            kept = [
                (cls, start, end)
                for cls, start, end in zip(classes, starts, ends)
                if start >= drop_through_right_start
            ]
        else:
            kept = [
                (cls, start, end)
                for cls, start, end in zip(classes, starts, ends)
                if end <= drop_start or start >= drop_end
            ]
        if len(kept) == len(classes):
            break
        classes = [item[0] for item in kept]
        starts = [item[1] for item in kept]
        ends = [item[2] for item in kept]

    return classes, starts, ends, duals


def load_none_mask(
    path: Path,
    max_sparse_good_bins: int = DEFAULT_MAX_SPARSE_GOOD_BINS,
    min_none_run_for_edge: int = DEFAULT_MIN_NONE_RUN_FOR_EDGE,
) -> tuple[set[tuple[str, int, int]], dict[str, int]]:
    """
    Load fixed None bins and per-chromosome ends from HGSVC-style mask.

    Without rewriting the mask file:
    1. Short ``good`` runs flanked on both sides by ``None`` are treated as None.
    2. Large ``None`` runs also absorb up to ``max_sparse_good_bins`` flanking
       ``good`` bins on each side, so short state stubs at heterochromatin
       edges (e.g. chr9 ~38 Mb) are skipped.
    """
    mask = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = {"chrom", "start", "end", "class"}
    missing = required - set(mask.columns)
    if missing:
        raise ValueError(f"Mask file missing columns: {sorted(missing)}")

    work = mask.copy()
    work["chrom"] = work["chrom"].astype(str)
    work["start"] = work["start"].astype(int)
    work["end"] = work["end"].astype(int)
    work["class"] = work["class"].astype(str).str.strip()

    if max_sparse_good_bins > 0:
        for _, chrom_bins in work.groupby("chrom", sort=False):
            chrom_bins = chrom_bins.sort_values(["start", "end"], kind="mergesort")
            classes = chrom_bins["class"].tolist()
            indices = chrom_bins.index.tolist()
            n = len(classes)

            # Fill short good islands inside None deserts.
            i = 0
            while i < n:
                if classes[i] != "good":
                    i += 1
                    continue
                j = i
                while j < n and classes[j] == "good":
                    j += 1
                left_none = i > 0 and classes[i - 1] == "None"
                right_none = j < n and classes[j] == "None"
                if (
                    left_none
                    and right_none
                    and (j - i) <= max_sparse_good_bins
                ):
                    for k in range(i, j):
                        classes[k] = "None"
                    work.loc[indices[i:j], "class"] = "None"
                i = j

            # Expand large None runs into flanking good edge bins.
            i = 0
            while i < n:
                if classes[i] != "None":
                    i += 1
                    continue
                j = i
                while j < n and classes[j] == "None":
                    j += 1
                if (j - i) >= min_none_run_for_edge:
                    left = i - 1
                    absorbed = 0
                    while (
                        left >= 0
                        and classes[left] == "good"
                        and absorbed < max_sparse_good_bins
                    ):
                        classes[left] = "None"
                        work.loc[indices[left], "class"] = "None"
                        absorbed += 1
                        left -= 1
                    right = j
                    absorbed = 0
                    while (
                        right < n
                        and classes[right] == "good"
                        and absorbed < max_sparse_good_bins
                    ):
                        classes[right] = "None"
                        work.loc[indices[right], "class"] = "None"
                        absorbed += 1
                        right += 1
                i = j

    none = work[work["class"] == "None"]
    none_bins = {
        (str(chrom), int(start), int(end))
        for chrom, start, end in zip(none["chrom"], none["start"], none["end"])
    }
    chrom_ends = work.groupby("chrom", sort=False)["end"].max().to_dict()
    return none_bins, chrom_ends


# chrY total (c+w) above this → male; QC HPNE cells are all well above.
DEFAULT_MIN_CHRY_DEPTH_MALE = 50.0


def classify_cells_sex(
    raw: pd.DataFrame,
    min_chrY_depth: float = DEFAULT_MIN_CHRY_DEPTH_MALE,
) -> dict[tuple[str, str], str]:
    """
    Infer male/female per cell from chrY mapping depth (sum of c+w).

    Cells with chrY depth >= ``min_chrY_depth`` are ``male``; otherwise
    ``female`` (including missing chrY or no c/w columns).
    """
    keys = (
        raw[["sample", "cell"]]
        .astype(str)
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    sex: dict[tuple[str, str], str] = {key: "female" for key in keys}
    if "c" not in raw.columns or "w" not in raw.columns:
        return sex

    y = raw[raw["chrom"].astype(str) == "chrY"].copy()
    if y.empty:
        return sex
    y["c"] = pd.to_numeric(y["c"], errors="coerce").fillna(0.0)
    y["w"] = pd.to_numeric(y["w"], errors="coerce").fillna(0.0)
    y["depth"] = y["c"] + y["w"]
    totals = y.groupby(["sample", "cell"], sort=False)["depth"].sum()
    for (sample, cell), depth in totals.items():
        key = (str(sample), str(cell))
        if depth >= min_chrY_depth:
            sex[key] = "male"
    return sex


def load_sv_intervals(
    path: Path,
    chrom_ends: dict[str, int] | None = None,
) -> dict[tuple[str, str, str], list[tuple[int, int, str]]]:
    """
    Load per-(sample, cell, chrom) SV intervals.

    Each interval is ``(start, end, sv_call_name)``. ``chrom_ends`` is accepted
    for call-site compatibility and is unused here; tip-linking of
    inversion / complex is decided later in ``_sv_skip_intervals``.
    """
    del chrom_ends  # retained for API compatibility
    sv = pd.read_csv(path, sep="\t")
    required = {"chrom", "start", "end", "sample", "cell"}
    missing = required - set(sv.columns)
    if missing:
        raise ValueError(f"SV file missing columns: {sorted(missing)}")

    has_sv_name = "sv_call_name" in sv.columns
    intervals: dict[tuple[str, str, str], list[tuple[int, int, str]]] = defaultdict(
        list
    )
    for row in sv.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.start)
        end = int(row.end)
        name = (
            str(getattr(row, "sv_call_name", "") or "")
            if has_sv_name
            else ""
        )
        key = (str(row.sample), str(row.cell), chrom)
        intervals[key].append((start, end, name))
    return intervals


def _sv_skip_intervals(
    intervals: list[tuple[int, int, str]],
    chrom: str,
    sex: str,
    chrom_end: int | None = None,
) -> list[tuple[int, int]]:
    """
    Build skip mask from SV intervals.

    Always skip deletion / duplication (``del_*``, ``dup_*``, ``idup_*``),
    except male chrX deletions (kept in the stitch path). Inversion and
    complex are skipped unless tip-linked: tip-reaching (±1 Mb), tip-직전
    (within 2 Mb of a tip), or abutting a tip-linked CN skip / tip-linked
    inv/complex. Empty call names are skipped conservatively.
    """
    cn_skips: list[tuple[int, int]] = []
    inv_complex: list[tuple[int, int]] = []
    for start, end, name in intervals:
        if chrom == "chrX" and sex == "male" and _is_deletion_call(name):
            continue
        if _is_inv_or_complex_call(name):
            inv_complex.append((start, end))
        elif not name or _is_cn_skip_sv(name):
            cn_skips.append((start, end))
        else:
            cn_skips.append((start, end))

    keep_idx = _tip_linked_inv_complex(inv_complex, cn_skips, chrom_end)
    skip = list(cn_skips)
    for i, (start, end) in enumerate(inv_complex):
        if i not in keep_idx:
            skip.append((start, end))
    return skip


def load_raw_bins(path: Path) -> pd.DataFrame:
    """Load 200 kb bin strand-state table (plain or gzip TSV)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        df = pd.read_csv(handle, sep="\t")

    required = {"chrom", "start", "end", "sample", "cell", "class"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Raw input missing columns: {sorted(missing)}")

    return df


def load_qc_cells(path: Path) -> set[tuple[str, str]]:
    """Load QC-passed (sample, cell) pairs from StrandPhaseR final output."""
    qc = pd.read_csv(path, sep="\t")
    required = {"sample", "cell"}
    missing = required - set(qc.columns)
    if missing:
        raise ValueError(f"QC file missing columns: {sorted(missing)}")

    return {
        (str(sample), str(cell))
        for sample, cell in zip(qc["sample"], qc["cell"])
    }


def filter_to_qc_cells(raw: pd.DataFrame, qc_cells: set[tuple[str, str]]) -> pd.DataFrame:
    """Keep only bins from QC-passed (sample, cell) pairs."""
    keys = pd.MultiIndex.from_arrays(
        [raw["sample"].astype(str), raw["cell"].astype(str)]
    )
    keep = keys.isin(qc_cells)
    filtered = raw.loc[keep].copy()
    if filtered.empty:
        raise ValueError("No raw bins remain after filtering to QC-passed cells")
    return filtered


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
    breakpoint; SCE rows leave that value empty. Inversion rows are unchanged.
    """
    if events.empty:
        events = events.copy()
        if "Event" not in events.columns:
            events["Event"] = pd.Series(dtype="object")
        if "Shared_cell_percent" not in events.columns:
            events["Shared_cell_percent"] = pd.Series(dtype="float")
        return events

    total_cells = source.groupby("sample")["cell"].nunique().to_dict()
    events = events.copy()
    sce_mask = events["Event"].isna() | (events["Event"] == "SCE")
    events.loc[sce_mask, "Event"] = "SCE"
    events.loc[sce_mask, "Shared_cell_percent"] = pd.NA

    sce_events = events[events["Event"] == "SCE"]
    for sample, sample_events in sce_events.groupby("Sample", sort=False):
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


def resolve_double_switch_patterns(
    duals: pd.DataFrame,
    source: pd.DataFrame,
    tolerance: int = 10_000,
    min_cell_fraction: float = 0.05,
) -> pd.DataFrame:
    """
    Resolve dual switches into Inversion or two SCE records.

    - ``twostep`` (WW→WC→CC / CC→WC→WW): always two SCEs.
    - ``aba``: ambiguous (true inversion vs two SCEs). Emit ``Inversion`` only
      when both breakpoints are shared by ≥5% of QC cells (subclone);
      otherwise two SCE rows. Centromere overlap does not force Inversion.
    """
    empty = pd.DataFrame(columns=list(OUTPUT_COLUMNS))
    if duals.empty:
        return empty

    total_cells = source.groupby("sample")["cell"].nunique().to_dict()
    records: list[dict[str, object]] = []

    for sample, sample_duals in duals.groupby("Sample", sort=False):
        sample_cells = total_cells[sample]
        required_cells = max(2, math.ceil(sample_cells * min_cell_fraction))
        for chrom, chrom_duals in sample_duals.groupby("chr", sort=False):
            starts = chrom_duals["start"]
            ends = chrom_duals["end"]
            cells = chrom_duals["Cell_ID"]
            for _, row in chrom_duals.iterrows():
                left = int(row["start"])
                right = int(row["end"])
                kind = str(row.get("kind", "aba"))
                if kind == "twostep":
                    for breakpoint in (left, right):
                        records.append(
                            {
                                "Sample": sample,
                                "Cell_ID": row["Cell_ID"],
                                "chr": row["chr"],
                                "start": breakpoint,
                                "end": pd.NA,
                                "Event": "SCE",
                                "Shared_cell_percent": pd.NA,
                            }
                        )
                    continue

                shared_mask = (starts - left).abs() <= tolerance
                shared_mask &= (ends - right).abs() <= tolerance
                shared_cells = cells[shared_mask].nunique()
                if shared_cells >= required_cells:
                    records.append(
                        {
                            "Sample": sample,
                            "Cell_ID": row["Cell_ID"],
                            "chr": row["chr"],
                            "start": left,
                            "end": right,
                            "Event": "Inversion",
                            "Shared_cell_percent": round(
                                shared_cells / sample_cells * 100, 2
                            ),
                        }
                    )
                else:
                    for breakpoint in (left, right):
                        records.append(
                            {
                                "Sample": sample,
                                "Cell_ID": row["Cell_ID"],
                                "chr": row["chr"],
                                "start": breakpoint,
                                "end": pd.NA,
                                "Event": "SCE",
                                "Shared_cell_percent": pd.NA,
                            }
                        )

    if not records:
        return empty
    return pd.DataFrame.from_records(records, columns=list(OUTPUT_COLUMNS))


def detect_sce(
    raw: pd.DataFrame,
    none_bins: set[tuple[str, int, int]],
    sv_intervals: dict[tuple[str, str, str], list[tuple[int, int, str]]],
    cell_sex: dict[tuple[str, str], str] | None = None,
    centromere_barriers: dict[str, list[tuple[int, int]]] | None = None,
    chrom_ends: dict[str, int] | None = None,
    translocation_tolerance: int = 10_000,
    translocation_min_fraction: float = 0.05,
    wc_dup_tot_ratio: float = DEFAULT_WC_DUP_TOT_RATIO,
    wc_dup_max_strand_ratio: float = DEFAULT_WC_DUP_MAX_STRAND_RATIO,
    wc_dup_min_strand_ratio: float = DEFAULT_WC_DUP_MIN_STRAND_RATIO,
    wc_dup_flank_bins: int = DEFAULT_WC_DUP_FLANK_BINS,
) -> pd.DataFrame:
    """
    Detect SCE, inversion, and recurrent translocation candidates.

    Skip HGSVC None bins and copy-number SV intervals (deletion /
    duplication / idup). Inversion and complex are skipped unless tip-linked
    (tip-reaching ±1 Mb, tip-직전 within 2 Mb of a tip, or abutting a
    tip-linked CN skip / tip-linked inv/complex). On male chrX, deletion SV
    intervals are not skipped.
    Centromere / arm barriers come from the species arm-position table when
    available (fallback: large HGSVC None runs). Non-SV WC runs that match
    the local-flank duplication depth pattern (one strand near-full, the
    other present, elevated total) are dropped. A-B-A sandwiches are
    resolved by subclone recurrence. Remaining single switches are scored
    per arm; a centromere-crossing exception recovers masked single switches.
    """
    work = raw.copy()
    work["start"] = work["start"].astype(int)
    work["end"] = work["end"].astype(int)
    work["class"] = work["class"].astype(str).str.strip()
    has_depth = "c" in work.columns and "w" in work.columns
    if has_depth:
        work["c"] = pd.to_numeric(work["c"], errors="coerce").fillna(0.0)
        work["w"] = pd.to_numeric(work["w"], errors="coerce").fillna(0.0)

    if cell_sex is None:
        cell_sex = classify_cells_sex(work)
    arm_centros = centromere_barriers or {}
    ends_by_chrom = chrom_ends or {}
    cell_depth_refs = (
        _cell_homozygous_depth_refs(work, none_bins) if has_depth else {}
    )

    single_records: list[dict[str, object]] = []
    dual_records: list[dict[str, object]] = []
    barriers_by_chrom: dict[str, list[tuple[int, int]]] = {}

    grouped = work.groupby(["sample", "cell", "chrom"], sort=False)
    for (sample, cell, chrom), group in grouped:
        chrom_key = str(chrom)
        if chrom_key == "chrY":
            continue
        segments = group.sort_values(["start", "end"], kind="mergesort")
        sex = cell_sex.get((str(sample), str(cell)), "female")
        sv_raw = sv_intervals.get((str(sample), str(cell), chrom_key), [])
        chrom_end = ends_by_chrom.get(chrom_key)
        if chrom_end is None and len(segments):
            chrom_end = int(segments["end"].max())
        skip_sv = _sv_skip_intervals(sv_raw, chrom_key, sex, chrom_end)
        if chrom_key not in barriers_by_chrom:
            barriers_by_chrom[chrom_key] = _centromere_barriers_for_chrom(
                chrom_key, arm_centros, none_bins
            )
        barriers = barriers_by_chrom[chrom_key]

        raw_starts = segments["start"].tolist()
        raw_ends = segments["end"].tolist()
        raw_classes = segments["class"].tolist()
        raw_c = segments["c"].tolist() if has_depth else None
        raw_w = segments["w"].tolist() if has_depth else None

        male_chrx = _is_male_chrx(chrom_key, sex)
        allowed_classes = MALE_X_VALID_CLASSES if male_chrx else VALID_CLASSES

        classes: list[str] = []
        starts: list[int] = []
        ends: list[int] = []
        for start, end, cls in zip(raw_starts, raw_ends, raw_classes):
            if (chrom_key, start, end) in none_bins:
                continue
            # Drop bins that sit inside the centromere barrier.
            if barriers and _interval_overlaps_barrier(start, end, barriers):
                continue
            if skip_sv and _overlaps_any(start, end, skip_sv):
                continue
            if cls not in allowed_classes:
                continue
            classes.append(cls)
            starts.append(start)
            ends.append(end)

        ref_depth = cell_depth_refs.get((str(sample), str(cell)), 0.0)
        if (
            has_depth
            and classes
            and ref_depth > 0
            and raw_c is not None
            and raw_w is not None
        ):
            classes, starts, ends = _drop_duplication_like_wc(
                classes,
                starts,
                ends,
                ref_depth=ref_depth,
                raw_starts=raw_starts,
                raw_ends=raw_ends,
                raw_classes=raw_classes,
                raw_c=[float(x) for x in raw_c],
                raw_w=[float(x) for x in raw_w],
                tot_ratio=wc_dup_tot_ratio,
                max_strand_ratio=wc_dup_max_strand_ratio,
                min_strand_ratio=wc_dup_min_strand_ratio,
                flank_bins=wc_dup_flank_bins,
            )

        classes, starts, ends = _drop_sv_adjacent_stubs(
            classes, starts, ends, skip_sv
        )

        # Detect A-B-A / two-step-opposite on the full chromosome.
        classes, starts, ends, aba_duals = _extract_aba_sandwiches(
            classes, starts, ends, barriers=barriers, male_chrx=male_chrx
        )
        for left_bp, right_bp, spans_centromere, kind in aba_duals:
            dual_records.append(
                {
                    "Sample": sample,
                    "Cell_ID": cell,
                    "chr": chrom,
                    "start": left_bp,
                    "end": right_bp,
                    "spans_centromere": spans_centromere,
                    "kind": kind,
                }
            )

        # Single-switch SCE only within one arm after cleanup.
        for arm_c, arm_s, arm_e in _split_bins_into_arms(
            classes, starts, ends, barriers
        ):
            _call_events_on_arm(
                sample,
                cell,
                chrom,
                arm_c,
                arm_s,
                arm_e,
                single_records,
                raw_starts=raw_starts,
                raw_ends=raw_ends,
                raw_classes=raw_classes,
                skip_sv=skip_sv,
                male_chrx=male_chrx,
            )

        # Exception: SCE whose breakpoint sits inside the centromere.
        _call_centromere_crossing_sce(
            sample,
            cell,
            chrom,
            classes,
            starts,
            ends,
            raw_starts,
            raw_ends,
            raw_classes,
            barriers,
            skip_sv,
            single_records,
            male_chrx=male_chrx,
        )

    singles = pd.DataFrame.from_records(
        single_records, columns=list(OUTPUT_COLUMNS)
    )
    duals = pd.DataFrame.from_records(
        dual_records,
        columns=[
            "Sample",
            "Cell_ID",
            "chr",
            "start",
            "end",
            "spans_centromere",
            "kind",
        ],
    )
    resolved_duals = resolve_double_switch_patterns(
        duals,
        work,
        tolerance=translocation_tolerance,
        min_cell_fraction=translocation_min_fraction,
    )

    combined = pd.concat([singles, resolved_duals], ignore_index=True)
    if combined.empty:
        return pd.DataFrame(columns=list(OUTPUT_COLUMNS))

    inversions = combined[combined["Event"] == "Inversion"].copy()
    sce_like = combined[combined["Event"] != "Inversion"].copy()
    sce_like = classify_recurrent_breakpoints(
        sce_like,
        work,
        tolerance=translocation_tolerance,
        min_cell_fraction=translocation_min_fraction,
    )
    result = pd.concat([sce_like, inversions], ignore_index=True)
    return result.sort_values(
        ["Sample", "Cell_ID", "chr", "start", "end"], kind="mergesort"
    ).reset_index(drop=True)


def write_excel(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(path, index=False, engine="openpyxl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Detect SCE and recurrent translocation candidates from 200 kb "
            "strand-state bins, skipping SV intervals and fixed None bins, "
            "then write an Excel summary (Sample, Cell_ID, chr, start, Event)."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Raw 200 kb bin TSV/TSV.GZ (chrom,start,end,sample,cell,...,class)",
    )
    parser.add_argument(
        "--sv",
        type=Path,
        required=True,
        help="SV interval TSV (e.g. lenient_filterFALSE.tsv)",
    )
    parser.add_argument(
        "--qc",
        type=Path,
        required=True,
        help=(
            "StrandPhaseR final output TSV listing QC-passed cells "
            "(e.g. StrandPhaseR_final_output.txt)"
        ),
    )
    parser.add_argument(
        "--species",
        type=str,
        default="human",
        help=(
            "Species for centromere/arm coordinates (default: human). "
            f"Known: {', '.join(sorted(SPECIES_ARM_POSITIONS))}"
        ),
    )
    parser.add_argument(
        "--arm-positions",
        type=Path,
        default=None,
        help=(
            "Optional override TSV of chromosome arm coordinates "
            "(Idf, Chrom, Start, End, ...). Default comes from --species "
            f"(human → {DEFAULT_HUMAN_ARM_POSITIONS.name})"
        ),
    )
    parser.add_argument(
        "--mask",
        type=Path,
        default=DEFAULT_MASK,
        help=f"Fixed None-bin mask (default: {DEFAULT_MASK.name})",
    )
    parser.add_argument(
        "--max-sparse-good-bins",
        type=int,
        default=DEFAULT_MAX_SPARSE_GOOD_BINS,
        help=(
            "Treat None-flanked good islands of at most this many 200 kb bins "
            f"as None (default: {DEFAULT_MAX_SPARSE_GOOD_BINS}; 0 disables)"
        ),
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
    parser.add_argument(
        "--wc-dup-tot-ratio",
        type=float,
        default=DEFAULT_WC_DUP_TOT_RATIO,
        help=(
            "Dup-like WC: require (mean c+w)/flank_dom >= this "
            f"(default: {DEFAULT_WC_DUP_TOT_RATIO}; 0 disables WC-dup filter)"
        ),
    )
    parser.add_argument(
        "--wc-dup-max-strand-ratio",
        type=float,
        default=DEFAULT_WC_DUP_MAX_STRAND_RATIO,
        help=(
            "Dup-like WC: require max(c,w)/flank_dom >= this "
            f"(default: {DEFAULT_WC_DUP_MAX_STRAND_RATIO})"
        ),
    )
    parser.add_argument(
        "--wc-dup-min-strand-ratio",
        type=float,
        default=DEFAULT_WC_DUP_MIN_STRAND_RATIO,
        help=(
            "Dup-like WC: require min(c,w)/flank_dom >= this "
            f"(default: {DEFAULT_WC_DUP_MIN_STRAND_RATIO})"
        ),
    )
    parser.add_argument(
        "--wc-dup-flank-bins",
        type=int,
        default=DEFAULT_WC_DUP_FLANK_BINS,
        help=(
            "Homozygous flank bins per side for WC dup depth "
            f"(default: {DEFAULT_WC_DUP_FLANK_BINS})"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.translocation_tolerance < 0:
        raise ValueError("--translocation-tolerance must be >= 0")
    if not 0 < args.translocation_min_fraction <= 1:
        raise ValueError("--translocation-min-fraction must be in (0, 1]")
    if args.max_sparse_good_bins < 0:
        raise ValueError("--max-sparse-good-bins must be >= 0")
    if not args.input.is_file():
        raise FileNotFoundError(f"Input not found: {args.input}")
    if not args.sv.is_file():
        raise FileNotFoundError(f"SV file not found: {args.sv}")
    if not args.qc.is_file():
        raise FileNotFoundError(f"QC file not found: {args.qc}")
    if not args.mask.is_file():
        raise FileNotFoundError(f"Mask file not found: {args.mask}")

    raw = load_raw_bins(args.input)
    qc_cells = load_qc_cells(args.qc)
    raw = filter_to_qc_cells(raw, qc_cells)
    none_bins, chrom_ends = load_none_mask(
        args.mask, max_sparse_good_bins=args.max_sparse_good_bins
    )
    arm_path = resolve_species_arm_path(args.species, args.arm_positions)
    if not arm_path.is_file():
        raise FileNotFoundError(f"Arm positions file not found: {arm_path}")
    centromere_barriers = load_centromere_barriers_from_arms(arm_path)
    sv_intervals = load_sv_intervals(args.sv, chrom_ends)
    cell_sex = classify_cells_sex(raw)
    n_male = sum(1 for s in cell_sex.values() if s == "male")
    n_female = sum(1 for s in cell_sex.values() if s == "female")
    events = detect_sce(
        raw,
        none_bins=none_bins,
        sv_intervals=sv_intervals,
        cell_sex=cell_sex,
        centromere_barriers=centromere_barriers,
        chrom_ends=chrom_ends,
        translocation_tolerance=args.translocation_tolerance,
        translocation_min_fraction=args.translocation_min_fraction,
        wc_dup_tot_ratio=args.wc_dup_tot_ratio,
        wc_dup_max_strand_ratio=args.wc_dup_max_strand_ratio,
        wc_dup_min_strand_ratio=args.wc_dup_min_strand_ratio,
        wc_dup_flank_bins=args.wc_dup_flank_bins,
    )
    write_excel(events, args.output)
    counts = events["Event"].value_counts()
    print(f"QC-passed cells: {len(qc_cells)}")
    print(f"Species: {args.species} (arms: {arm_path.name})")
    print(f"Centromeres loaded: {len(centromere_barriers)} chroms")
    print(f"Sex (chrY depth): male={n_male}, female={n_female}")
    print(f"Detected {len(events)} breakpoint candidate(s)")
    print(f"SCE: {counts.get('SCE', 0)}")
    print(f"Translocation: {counts.get('Translocation', 0)}")
    print(f"Inversion: {counts.get('Inversion', 0)}")
    print(f"Wrote: {args.output.resolve()}")


if __name__ == "__main__":
    main()
