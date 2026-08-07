# SCEDetector

A lightweight tool that finds **SCE**, **Inversion**, and recurrent **Translocation** candidates from 200 kb bin strand states (`CC` / `WC` / `WW`).

Calls are **screening candidates**, not definitive events. The translocation label only reflects breakpoint recurrence across cells; it does not identify translocation partners.

---

## Pipeline at a glance

For each cell × chromosome (except `chrY`, used only for sex inference):

```
[1] Bin masking
      Drop None / centromere / SV / invalid strand classes
        ↓
[2] Depth-based artifact removal
      (a) Drop duplication-like WC
      (b) Drop short stubs at SV hole edges
      (c) Drop lopsided A-WC-A WC islands
      (d) Drop short depth-dropped WC-A-WC homozygous islands
        ↓
[3] Double-switch extraction
      A-B-A sandwich → Inversion or two SCEs
      WW→WC→CC / CC→WC→WW → always two SCEs
        ↓
[4] Single-switch calling
      Exactly one state switch per chromosome arm → SCE
      (+ centromere-crossing exception)
        ↓
[5] Breakpoint coordinate refinement
      Move stitched breakpoints to the true switch inside an SV hole
        ↓
[6] Recurrence relabeling
      Shared breakpoints across cells → Translocation
```

---

## Requirements

- Python 3.9+
- `pandas`, `openpyxl`

```bash
pip install -r requirements.txt
```

## Usage

```bash
python3 SCEDetector.py \
  -i fastq0022_HPNE_M.txt.raw.gz \
  --sv lenient_filterFALSE.tsv \
  --qc StrandPhaseR_final_output.txt \
  -o SCE_detected.xlsx
```

---

## Inputs

| Argument | File | Role |
|----------|------|------|
| `-i` | `*.txt.raw.gz` | 200 kb bin strand states (`chrom, start, end, sample, cell, class, c, w`) |
| `--sv` | `lenient_filterFALSE.tsv` | Per-cell SV intervals |
| `--qc` | `StrandPhaseR_final_output.txt` | QC-passed cells (only these are analyzed) |
| `--mask` | `HGSVC.200000.txt` (default) | Low-mappability / None-bin mask |
| `--species` | `chromosome_arm_positions_grch38.txt` (via `human`) | Centromere barriers (`[p.End, q.Start)`) |

**Sex:** inferred from chrY read depth. On male chrX (hemizygous), `WC` is removed and only `WW↔CC` is allowed.

---

## Algorithm (step by step)

### Step 1 — Bin masking

Each 200 kb bin is checked in order and **dropped** if any of the following applies:

| Condition | Description |
|-----------|-------------|
| None mask | `None` bins from `HGSVC.200000.txt`. Short `good` islands (≤5 bins) inside None are treated as None. Large None runs (≥10 bins) also absorb up to 5 flanking `good` bins on each side |
| Centromere | `[p.End, q.Start)` from the arm table; if missing, fall back to large None runs (≥1 Mb) |
| SV interval | See SV rules below |
| Invalid class | Keep `CC` / `WC` / `WW` normally; male chrX keeps only `CC` / `WW` |

**When are SVs removed?**

- **Deletion / duplication** (`del_*`, `dup_*`, `idup_*`): **always** removed  
  (exception: male chrX deletions are kept)
- **Inversion / complex**: removed by default, **kept only if tip-linked**
  - reaches a chromosome tip (±1 Mb), or
  - starts within 2 Mb of pter / ends within 2 Mb of qter, or
  - abuts another tip-linked SV

Surviving bins are stitched into **state runs**.  
Example: `WW WWW WW WC WC CC CC` → `WW | WC | CC`

---

### Step 2 — Depth-based artifact removal

Some regions look like SCE from strand state alone but are duplications or noise by **read depth**. Key measurements:

| Symbol | Meaning |
|--------|---------|
| `c`, `w` | Median Crick / Watson depth on the **full contiguous raw WC run** (including bins that the None mask would skip) |
| `dom` | Dominant-strand depth of abutting homozygous run(s) (WW→`w`, CC→`c`). Average of both sides when present; floored by the cell-wide homozygous median when flanks are missing or too shallow |
| `balance` | `min(c,w) / max(c,w)` — closer to 1 means more even strands |
| Expected coverage | Per-bin **across-cell median of relative depth** × this cell’s genome-wide median. Prevents mappability-hot bands (deep in every cell) from being mistaken for gains |

#### 2a. Drop duplication-like WC

An unmarked duplication often keeps one strand near-full and adds the opposite strand. True WC has both strands near half.

**Flank asymmetry — drop only if all three hold**

1. `max(c,w) / dom ≥ 0.85` — one strand still near homozygous depth
2. `min(c,w) / dom ≥ 0.35` — the other strand is clearly present
3. `(c+w) / dom ≥ 1.35` — total depth is elevated vs the flank

**Coverage ceiling — enough to drop on its own**

- `(c+w) > expected coverage × 1.50`

| | True WC (keep) | Tip duplication (drop) |
|--|----------------|------------------------|
| Pattern | both ≈ 0.5×dom, total ≈ dom | one ≈ dom, other ≈ 0.5×dom, total ≈ 1.5×dom |

#### 2b. Drop short stubs at SV hole edges

After SV removal, 1–2 bin state stubs at hole edges create false sandwiches.  
Any run of **≤2 bins** that touches an SV hole is dropped so flanking states can merge.

#### 2c. Drop lopsided A-WC-A WC islands

A WC island between two identical homozygous runs (`WW-WC-WW` / `CC-WC-CC`) often resolves to a spurious pair of SCEs. A real WC island should be strand-balanced and copy-neutral.  
**Keep if any rule below holds**; otherwise drop the island so the flanks merge.  
(Single WC transitions are never touched.)

| # | Keep when | Intuition |
|---|-----------|-----------|
| 1 | `balance ≥ 0.75` | Strands are sufficiently even |
| 2 | length ≥ 10 Mb **and** `balance ≥ 0.50` | Longer islands get a relaxed balance bar |
| 3 | `(c+w) < 0.55 × dom` | Depth-depleted islands are a different artifact — leave them alone here |
| 4 | Flank-dominant strand inside the island is `≤ 0.47 ×` the deeper flank **and** `balance ≥ 0.35` | Clean halving of the dominant strand (copy-neutral SCE) |
| 5 | Both flanks exist, length ≤ 5 Mb, `balance ≥ 0.55`, island total depth lies **between** the two flanks, dominant strand still `> 0.55 ×` the deeper flank, minor strand `≥ 0.45 ×` the shallower flank | Uneven flanks (shallower ≈ island, deeper flank larger) |

#### 2d. Drop short WC-A-WC homozygous islands

Mirror case: a short `WW`/`CC` island inside WC (`WC-CC-WC` / `WC-WW-WC`).  
Two SCEs only a few Mb apart with a depth drop are implausible.

**Drop only when both hold**

1. Island length `< 5 Mb`
2. Island `(c+w) < deeper flanking WC × 0.80`

Short islands that stay depth-neutral are kept.

---

### Step 3 — Double-switch extraction and resolution

Double switches are pulled out of the filtered state runs and resolved separately.

#### A-B-A sandwich (e.g. `WW-WC-WW`, `CC-WC-CC`)

Return to the original flank class → could be a **true inversion** or **two SCEs**.

| Priority | Condition | Result |
|----------|-----------|--------|
| 1 | Both breakpoints shared within ±10 kb by ≥5% of QC cells in the same sample·chrom | `Inversion` (+ shared %) |
| 2 | Exactly **one** inversion SV call lies inside the sandwich and covers ≥40% of it | `Inversion` (shared % empty — marks SV-backed) |
| 3 | Otherwise | Two `SCE` rows |

> Example for rule 2: the SV caller places `inv_h2` at 28.2–31.8 Mb while the raw WC block runs 28.2–34.6 Mb.  
> A clipped boundary leaves a tail that would otherwise become two SCEs.  
> Long sandwiches that merely contain **several** scattered inversion calls stay as two SCEs.

Crossing the centromere alone does **not** force `Inversion`.

#### Two-step opposite (`WW→WC→CC` / `CC→WC→WW`)

Always **two SCEs**. A centromere inside the final homozygous run does not make this an inversion.

#### Other patterns

- Matching flanks but not a valid SCE sandwich (e.g. `CC-WW-CC`): clear the middle, emit nothing
- Other `A-B-C` middles: clear the middle

---

### Step 4 — Single-switch SCE

After double switches are removed, split at the centromere into arms. For each arm:

- exactly **two** state runs, and
- a valid SCE transition  
→ one breakpoint = `SCE`

**Centromere exception:** if the kept flanks on either side of the barrier form a valid SCE transition and each flank is ≥5 Mb, emit one cross-centromere SCE. Prefer the raw switch inside the barrier when there is exactly one; otherwise use the start of the right flank.

**Valid SCE transitions**

| From | To |
|------|----|
| `WC` | `WW`, `CC` |
| `WW` | `WC` |
| `CC` | `WC` |

Male chrX: `WW↔CC` only.

---

### Step 5 — Breakpoint refinement across SV holes

After stitching, the default breakpoint is the start of the first surviving bin **after** the removed gap.  
If the gap (or the breakpoint) touches an SV, peek the raw bins inside the gap **without** the SV mask:

- exactly one `left → right` switch → move the coordinate there
- none or more than one → keep the default (far / right edge)

Sandwich breakpoints use the same correction, but **only when SV intervals tile the gap completely**.  
(Sandwich runs are merged over the whole chromosome; without that guard the peek can wander into a centromere or None stretch.)

> Without correction, a stitched breakpoint snaps to the **far edge** of the removed SV.  
> That can collide with a real breakpoint in another cell and fake a recurrent subclone.  
> Example: G63 chr10 — after removing complex 127.8–131.0 Mb the call sat at 131.0, matching G04’s true 131.0, so both were labeled Translocation. After correction: 128.0 (true WC→CC).

---

### Step 6 — Recurrence → Translocation

Among collected SCE-like breakpoints:

1. Compare only within the same `sample` + `chrom`
2. Count distinct cells with a breakpoint within ±10 kb
3. If that count is **≥5%** of QC cells (at least 2 cells) → `Translocation` + shared %
4. `Inversion` rows are never relabeled as translocation

---

## Output

| Column | Meaning |
|--------|---------|
| `Sample` | Sample name |
| `Cell_ID` | Cell ID |
| `chr` | Chromosome |
| `start` | Breakpoint (left end for inversions) |
| `end` | Right breakpoint for `Inversion`; empty for SCE / Translocation |
| `Event` | `SCE` / `Translocation` / `Inversion` |
| `Shared_cell_percent` | Shared-cell % for Translocation and recurrent Inversion; empty for SCE and SV-backed Inversion |

If `Shared_cell_percent` is empty on an `Inversion`, it was called from an **SV-backed** rule; if filled, it was called from **cross-cell recurrence**.

---

## Command-line arguments

Only file / species paths are exposed. Filter thresholds are fixed in the code to
the values used for `SCE_detected.xlsx` (see Algorithm above).

| Argument | Default | Description |
|----------|---------|-------------|
| `-i` / `--input` | (required) | Raw 200 kb bins |
| `--sv` | (required) | SV intervals |
| `--qc` | (required) | QC-passed cells |
| `--mask` | `HGSVC.200000.txt` | None-bin mask |
| `--species` | `human` | Arm / centromere table |
| `-o` / `--output` | `SCE_detected.xlsx` | Output path |

---

## Scope and limitations

| Included | Not included |
|----------|--------------|
| Single valid switch per arm after filtering | Full model of double / multiple crossovers |
| Breakpoint recurrence → translocation **candidates** | Translocation partner identification |
| `CC` / `WC` / `WW` patterns + depth | Full SV or haplotype-aware validation |

Arms that still have multiple state changes after filtering are not called as a single SCE (double switches are already handled in Step 3).

---

## License

Use freely for research and analysis.

## Contributor

* Hyunjin Cho

[hyun-jin891](https://github.com/hyun-jin891)
