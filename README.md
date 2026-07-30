# SCEDetector

A **lightweight SCE and recurrent translocation candidate inference** tool for
StrandPhaseR strand-state output.

Results are heuristic candidates intended for quick screening rather than
definitive SCE or structural-variant calling. Double crossovers are not modeled.
The translocation label is based on recurrent breakpoints across cells and does
not identify translocation partners.

## Overview

Given a tab-separated StrandPhaseR final output file, SCEDetector first finds
single-switch SCE-like candidates. It then compares their breakpoint positions
across cells in the same sample and chromosome. Recurrent positions are labeled
as translocation candidates; all remaining events are labeled as SCE candidates.

### Scope and limitations

| Included | Not considered |
|----------|----------------|
| Single valid class switch per chromosome | Double (or multiple) crossovers |
| Recurrent-breakpoint translocation heuristic | Translocation partner identification |
| Pattern matching on `CC` / `CW` / `WC` / `WW` | Full structural-variant or haplotype-aware validation |

Chromosomes with multiple state changes are skipped by design.

## Requirements

- Python 3.9+
- `pandas`
- `openpyxl`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Input format

Tab-separated file with the following columns (as in `StrandPhaseR_final_output.txt`):

| Column  | Description                          |
|---------|--------------------------------------|
| `chrom` | Chromosome (e.g. `chr1`)             |
| `start` | Segment start coordinate             |
| `end`   | Segment end coordinate               |
| `sample`| Sample name                          |
| `cell`  | Cell ID                              |
| `class` | Strand state: `CC`, `CW`, `WC`, `WW` |

Example:

```text
chrom	start	end	sample	cell	class
chr1	0	22800000	fastq0024_HPNE_T3	...A40	CW
chr1	22800000	248956422	fastq0024_HPNE_T3	...A40	WW
```

## SCE candidate rule

For each (`sample`, `cell`, `chrom`):

1. Sort segments by genomic start.
2. Merge consecutive segments with the same class.
3. Call an **SCE candidate** only if there is **exactly one** class change along the chromosome, and that change matches a valid SCE transition.
4. After the change, the new class must continue to the chromosome end (no further state changes).

### Valid SCE transitions

| From | To        |
|------|-----------|
| `CC` | `CW`, `WC` |
| `CW` | `CC`, `WW` |
| `WC` | `CC`, `WW` |
| `WW` | `CW`, `WC` |

Transitions such as `CC ↔ WW` or `CW ↔ WC`, and chromosomes with multiple class changes, are **not** counted as SCE candidates.

The candidate start breakpoint is the genomic `start` of the segment where the new class begins.

## Translocation candidate rule

After SCE-like breakpoints are detected:

1. Events are compared only within the same `sample` and `chrom`.
2. For each breakpoint, the detector counts distinct cells with a breakpoint
   within ±10 kb.
3. If those cells represent at least 5% of all cells in that sample, their
   events are labeled `Translocation` instead of `SCE`, and the shared
   percentage is reported in `Shared_cell_percent`.
4. At least two distinct cells are required, even when one cell alone would
   exceed 5% in a small sample.

The 10 kb tolerance and 5% threshold can be changed through command-line
arguments. This recurrence-based label indicates a possible subclone; it is not
proof of a translocation.

## Usage

```bash
python3 SCEDetector.py -i StrandPhaseR_final_output.txt -o SCE_detected.xlsx
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Input TSV path (required) |
| `-o`, `--output` | Output Excel path (default: `SCE_detected.xlsx`) |
| `--translocation-tolerance` | Breakpoint tolerance in bp (default: `10000`) |
| `--translocation-min-fraction` | Minimum shared-cell fraction (default: `0.05`) |

## Output format

Excel file with columns:

| Column    | Description                         |
|-----------|-------------------------------------|
| `Sample`  | Sample name                         |
| `Cell_ID` | Cell ID                             |
| `chr`     | Chromosome                          |
| `start`   | SCE candidate start breakpoint      |
| `Event`   | `SCE` or `Translocation`            |
| `Shared_cell_percent` | Percentage of sample cells sharing the breakpoint; empty for `SCE` rows |

Example:

| Sample | Cell_ID | chr | start | Event | Shared_cell_percent |
|--------|---------|-----|-------|-------|---------------------|
| fastq0024_HPNE_T3 | ...A10 | chr15 | 22600000 | Translocation | 26.47 |
| fastq0024_HPNE_T3 | ...A02 | chr10 | 50600000 | SCE | |

## License

Use freely for research and analysis.

## Contributor

* Hyunjin Cho

[hyun-jin891](https://github.com/hyun-jin891)
