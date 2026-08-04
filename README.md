# SCEDetector

A **lightweight SCE and recurrent translocation candidate inference** tool for
200 kb bin strand-state output.

Results are heuristic candidates intended for quick screening rather than
definitive SCE or structural-variant calling. Double crossovers are not modeled.
The translocation label is based on recurrent breakpoints across cells and does
not identify translocation partners.

## Overview

Given:

1. a raw 200 kb strand-state table (e.g. `*.txt.raw.gz`)
2. a per-cell SV interval table (e.g. `lenient_filterFALSE.tsv`)
3. a StrandPhaseR final output listing QC-passed cells (e.g. `StrandPhaseR_final_output.txt`)
4. a fixed low-mappability / None-bin mask (`HGSVC.200000.txt`, shipped with the tool)
5. a species chromosome arm table for centromeres (default human:
   `chromosome_arm_positions_grch38.txt`)

Only QC-passed cells from the StrandPhaseR final output are analyzed. SCEDetector
skips SV, None, and centromere-gap bins (copy-number SVs always; inversion /
complex only when they are not tip-linked), cleans nested
inversion-like sandwiches, and calls single-switch SCE-like candidates.
Whole-chromosome double switches (`A-B-A`) become `Inversion` when shared by
≥5% of cells, otherwise two `SCE` records. Recurrent single breakpoints shared
by ≥5% of QC-passed cells are labeled as translocation candidates.

### Scope and limitations

| Included | Not considered |
|----------|----------------|
| Single valid class switch per chromosome after skipping SV/None bins | Double (or multiple) crossovers |
| Recurrent-breakpoint translocation heuristic | Translocation partner identification |
| Pattern matching on `CC` / `WC` / `WW` | Full structural-variant or haplotype-aware validation |

Chromosomes with multiple state changes (after filtering) are skipped by design.

## Requirements

- Python 3.9+
- `pandas`
- `openpyxl`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Input formats

### 1. Raw 200 kb bins (`-i`)

Tab-separated file (plain or `.gz`) with at least:

| Column  | Description |
|---------|-------------|
| `chrom` | Chromosome (e.g. `chr1`) |
| `start` | Bin start |
| `end`   | Bin end |
| `sample`| Sample name |
| `cell`  | Cell ID |
| `class` | Strand state: `CC`, `WC`, `WW`, or `None` |
| `c` / `w` | Crick / Watson read counts (used to drop duplication-like WC) |

### 2. SV intervals (`--sv`)

Tab-separated file with SV spans already expressed on the same 200 kb grid
(or overlapping it), e.g. `lenient_filterFALSE.tsv`:

| Column  | Description |
|---------|-------------|
| `chrom` | Chromosome |
| `start` | SV interval start |
| `end`   | SV interval end |
| `sample`| Sample name |
| `cell`  | Cell ID |

Bins overlapping these intervals for the same `(sample, cell, chrom)` are
skipped as follows:

- **deletion / duplication** (`del_*`, `dup_*`, `idup_*`): always skipped,
  except **male chrX deletions** (kept when stitching; sex from chrY `c+w`).
- **inversion / complex**: skipped **unless tip-linked** — the interval
  reaches a chromosome tip (±1 Mb), sits tip-직전 (starts within 2 Mb of
  pter or ends within 2 Mb of qter), or abuts (transitively) a tip-linked
  CN skip or another tip-linked inv/complex. Tip-linked inv/complex strand
  states stay available for SCE (e.g. WC→WW at a complex start that abuts a
  tip deletion, or a tip-직전 complex with a short uncovered tip stub).

**Male chrX:** WC segments are removed; SCE / A-B-A (inversion vs double SCE)
use only `WW↔CC` transitions.

### 3. QC-passed cells (`--qc`)

StrandPhaseR final output TSV, e.g. `StrandPhaseR_final_output.txt`. Only
`(sample, cell)` pairs present in this file are retained for SCE calling and
for the translocation shared-cell denominator.

| Column  | Description |
|---------|-------------|
| `sample`| Sample name |
| `cell`  | Cell ID |

Other columns (`chrom`, `start`, `end`, `class`, …) may be present and are ignored
for QC filtering.

### 4. Fixed None mask (`--mask`, default `HGSVC.200000.txt`)

Shipped reference file of low-mappability / centromere-like bins:

| Column  | Description |
|---------|-------------|
| `chrom` | Chromosome |
| `start` | Bin start |
| `end`   | Bin end |
| `class` | `None` or `good` |

Only `class == None` bins are skipped by default. In addition, without modifying
this file:

1. Short `good` islands flanked on both sides by `None` (≤5 bins by default)
   are treated as `None`.
2. Large `None` runs (≥10 bins) also absorb up to 5 flanking `good` bins on
   each side, so short state stubs at heterochromatin edges (e.g. chr9 ~38 Mb)
   are skipped.

Use `--max-sparse-good-bins` to change or disable (`0`) that fill-in.

### 5. Chromosome arm / centromere table (`--species` / `--arm-positions`)

Default for `--species human` is `chromosome_arm_positions_grch38.txt`
(GRCh38 p/q arms for chr1–22). Centromere barrier per chrom is the gap
between p-arm end and q-arm start: `[p.End, q.Start)`.

Chromosomes missing from the arm table (e.g. chrX) fall back to large
HGSVC `None` runs (≥1 Mb) as barriers.

## SCE candidate rule

For each (`sample`, `cell`, `chrom`) — **`chrY` is skipped** (still used only
to infer male/female via mapping depth):

1. Sort 200 kb bins by genomic start.
2. Skip bins that are `None` in the fixed mask, overlap the centromere
   barrier, overlap an SV skip interval (del/dup/idup always; inv/complex
   only when not tip-linked; male chrX deletions kept), or lack a valid
   strand class (`CC` / `WC` / `WW`; male chrX keeps only `CC` / `WW`). Then
   drop **duplication-like WC** runs using the local-flank asymmetry rule (see
   below). Remaining true WC keeps the original SCE logic. Then drop short
   (≤2 bin) state stubs that abut an SV skip hole so flanking states can
   merge (guards idup/dup edge artifacts).
3. Stitch remaining bins on the full chromosome and find `A-B-A` sandwiches
   and two-step opposite paths (`WW→WC→CC` / `CC→WC→WW`).
4. Valid `A-B-A` (ambiguous inversion vs two SCEs): both breakpoints shared
   by ≥5% of QC cells (subclone) → `Inversion`; otherwise two `SCE` rows.
   Centromere overlap does not force `Inversion`. Middle run is removed after
   extraction.
5. Two-step opposite (`WW→WC→CC` / `CC→WC→WW`): always two `SCE` rows (even if
   the final homozygous run includes the centromere — not `Inversion`).
6. Matching-flank sandwiches with invalid SCE transitions (e.g. `CC-WW-CC`) and
   other non-sandwich `A-B-C` middles are dropped without emitting a dual.
7. After cleanup, split at the centromere barrier and call a remaining
   **exactly one** valid class change held to the **arm end** as SCE.
8. **Centromere exception:** if kept flanks on either side of the centromere
   barrier differ by a valid SCE transition and each flank is ≥5 Mb, emit one
   `SCE`: prefer the raw breakpoint inside the barrier when there is exactly
   one matching switch; otherwise use the start of the right kept flank.
9. If a called SCE breakpoint (or the filtered gap that produced it) touches an
   SV skip interval, replace the coordinate with the raw `left→right` switch in
   that gap (judgment unchanged; only the reported position is corrected).

### Valid SCE transitions

| From | To |
|------|----|
| `WC` | `WW`, `CC` |
| `WW` | `WC` |
| `CC` | `WC` |

**Male chrX only:** WC bins are dropped; allowed switches are `WW↔CC` only.

### Double-switch rules

**`A-B-A`** (return to the original flank class), including nested sandwiches —
ambiguous between a true inversion and two SCEs:

1. If both breakpoints are shared within the same `sample` and `chrom`
   (±10 kb) by ≥5% of QC cells (subclone), emit `Inversion`.
2. Else emit two `SCE` rows.

Centromere overlap is not used to force `Inversion`.

**Two-step opposite** (`WW→WC→CC` or `CC→WC→WW`):

Always emit two `SCE` rows. Centromere inside the final homozygous state does
**not** make this an `Inversion` (contrast with `A-B-A`).

### Translocation candidate rule

After SCE-like breakpoints are collected (including non-shared double-switch
SCEs):

1. Events are compared only within the same `sample` and `chrom`.
2. For each breakpoint, the detector counts distinct cells with a breakpoint
   within ±10 kb.
3. If those cells represent at least 5% of all QC-passed cells in that sample,
   their events are labeled `Translocation` instead of `SCE`, and the shared
   percentage is reported in `Shared_cell_percent`.
4. At least two distinct cells are required, even when one cell alone would
   exceed 5% in a small sample.
5. `Inversion` rows are not relabeled as translocation.

The 10 kb tolerance and 5% threshold can be changed through command-line
arguments. Recurrence-based labels indicate a possible subclone; they are not
proof of a translocation or inversion.

## Usage

```bash
python3 SCEDetector.py \
  -i fastq0022_HPNE_M.txt.raw.gz \
  --sv lenient_filterFALSE.tsv \
  --qc StrandPhaseR_final_output.txt \
  -o SCE_detected.xlsx
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Raw 200 kb bin TSV / TSV.GZ (required) |
| `--sv` | SV interval TSV (required) |
| `--qc` | StrandPhaseR final output with QC-passed cells (required) |
| `--mask` | Fixed None-bin mask (default: `HGSVC.200000.txt` next to the script) |
| `--species` | Species for arm/centromere table (default: `human`) |
| `--arm-positions` | Optional override TSV of p/q arm coordinates |
| `-o`, `--output` | Output Excel path (default: `SCE_detected.xlsx`) |
| `--translocation-tolerance` | Breakpoint tolerance in bp (default: `10000`) |
| `--translocation-min-fraction` | Minimum shared-cell fraction (default: `0.05`) |
| `--wc-dup-tot-ratio` | Dup filter: require `(c+w)/dom ≥` this (default: `1.35`; `0` disables the filter) |
| `--wc-dup-max-strand-ratio` | Dup filter: require `max(c,w)/dom ≥` this (default: `0.85`) |
| `--wc-dup-min-strand-ratio` | Dup filter: require `min(c,w)/dom ≥` this (default: `0.35`) |
| `--wc-dup-flank-bins` | Homozygous flank bins per side used for `dom` (default: `25`) |

### WC duplication filter (local-flank asymmetry)

Goal: remove unmarked duplication that looks like WC (both colors present) but
is not a true sister-chromatid WC state.

**Measurements**
- WC depth: mean `c` and `w` on the **full contiguous raw WC run**, including
  bins that the None mask would otherwise skip.
- Flank depth `dom`: mean dominant strand of abutting homozygous run(s)
  (WW→`w`, CC→`c`), up to `--wc-dup-flank-bins` bins per side; average if both
  sides exist, then **floored by** the cell-wide homozygous median so shallow
  peri-centromere / noisy flanks do not inflate ratios. If no homozygous
  flank, use the cell-wide median alone.

**Drop (all three must hold)**
1. `max(c,w) / dom ≥ 0.85` — at least one strand still near full homozygous
2. `min(c,w) / dom ≥ 0.35` — the other strand is clearly present
3. `(c+w) / dom ≥ 1.35` — total depth is elevated vs the flank

**Why this matches Strand-seq biology**
- True WC: each strand ≈ half of `dom`, so `max` and `min` both ≈ 0.5 and
  `tot/dom ≈ 1` → keep.
- Tip duplication (keep original strand + add the opposite): `max ≈ 1`,
  `min ≈ 0.5`, `tot/dom ≈ 1.5` → drop (e.g. G61 chr17 q-tip).

## Output format

Excel file with columns:

| Column    | Description                         |
|-----------|-------------------------------------|
| `Sample`  | Sample name                         |
| `Cell_ID` | Cell ID                             |
| `chr`     | Chromosome                          |
| `start`   | Breakpoint start (left end of an inversion) |
| `end`     | Right breakpoint for `Inversion`; empty for `SCE` / `Translocation` |
| `Event`   | `SCE`, `Translocation`, or `Inversion` |
| `Shared_cell_percent` | Shared-cell percentage for `Translocation` / `Inversion`; empty for `SCE` |

Example:

| Sample | Cell_ID | chr | start | Event | Shared_cell_percent |
|--------|---------|-----|-------|-------|---------------------|
| fastq0022_HPNE_M | ...G02 | chr7 | 62200000 | SCE | |
| fastq0022_HPNE_M | ...G03 | chr7 | 17600000 | Translocation | 8.33 |

## License

Use freely for research and analysis.

## Contributor

* Hyunjin Cho

[hyun-jin891](https://github.com/hyun-jin891)
