# SCEDetector

A **lightweight SCE candidate inference** tool for StrandPhaseR strand-state output.

It does **not** model complex rearrangements such as **double crossovers** or **translocations**. Results are heuristic SCE candidates based on a simple single-switch strand-state pattern, intended for quick screening rather than definitive SCE calling.

## Overview

Given a tab-separated StrandPhaseR final output file, SCEDetector groups segments by sample, cell, and chromosome, then looks for a single valid strand-state class change that is maintained to the end of the chromosome. Matching events are written to an Excel file as SCE candidates.

### Scope and limitations

| Included | Not considered |
|----------|----------------|
| Single valid class switch per chromosome | Double (or multiple) crossovers |
| Switch held to chromosome end | Translocations / inter-chromosomal events |
| Pattern matching on `CC` / `CW` / `WC` / `WW` | Full structural-variant or haplotype-aware SCE validation |

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

## Usage

```bash
python3 SCEDetector.py -i StrandPhaseR_final_output.txt -o SCE_detected.xlsx
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Input TSV path (required) |
| `-o`, `--output` | Output Excel path (default: `SCE_detected.xlsx`) |

## Output format

Excel file with columns:

| Column    | Description                         |
|-----------|-------------------------------------|
| `Sample`  | Sample name                         |
| `Cell_ID` | Cell ID                             |
| `chr`     | Chromosome                          |
| `start`   | SCE candidate start breakpoint      |

## License

Use freely for research and analysis.

## Contributor

* Hyunjin Cho

[hyun-jin891](https://github.com/hyun-jin891)
