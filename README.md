# SCEDetector

**SCEDetector** screens single-cell Strand-seq data for candidate **sister chromatid exchanges (SCEs)**, **inversions**, and recurrent **translocation-like** breakpoints. It works from fixed **200 kb** genomic bins labeled with Watson/Crick strand states (`CC`, `WC`, `WW`) and uses read depth plus structural-variant (SV) calls to remove common artifacts.

> **Important.** Every call is a **screening candidate**, not a definitive biological event. The label `Translocation` only means that the same breakpoint coordinate recurs across cells; the tool does **not** identify translocation partners or prove that a rearrangement occurred.

---

## Biological background (read this first)

### What Strand-seq measures

In Strand-seq, DNA is labeled so that, after replication, each sister chromatid can be read as mostly **Watson (W)** or mostly **Crick (C)** template. For every cell and every 200 kb bin the pipeline reports:

| Class | Meaning in a diploid autosome |
|-------|-------------------------------|
| `WW` | Both homologs / both sisters read as Watson-dominated (homozygous Watson) |
| `CC` | Both read as Crick-dominated (homozygous Crick) |
| `WC` | One Watson-like and one Crick-like contribution (heterozygous / mixed) |

The raw file also carries per-bin depths `c` and `w` (Crick and Watson read counts). Those depths are essential: the same letter pattern can be a true exchange, a duplication, a deletion, or noise next to the centromere.

### What an SCE looks like

A **sister chromatid exchange** swaps template strands between sisters. In Strand-seq this appears as a **single, clean change of state** along a chromosome arm, for example:

```text
WW WWW WWW | WC WC WC WC     →  breakpoint = start of the WC run
CC CC | WC WC | WW WW WW     →  not a single SCE (two switches; handled separately)
```

Valid single-switch SCE transitions used here:

| From | To |
|------|----|
| `WC` | `WW` or `CC` |
| `WW` | `WC` |
| `CC` | `WC` |

Direct `WW ↔ CC` without a `WC` middle is **not** treated as a normal autosomal SCE (except on **male chrX**, which is hemizygous: only `WW ↔ CC` is allowed, and `WC` is discarded).

### Why inversions and “two SCEs” look similar

A heterozygous **inversion** can flip the strand state of one interval and then restore the original flank state, producing an **A–B–A sandwich**:

```text
WW ── WC ── WW     or     CC ── WC ── CC
```

The same pattern can also be **two genuine SCEs** a few megabases apart. SCEDetector cannot always tell them apart from one cell alone, so it uses:

1. **Recurrence across cells** at both breakpoints, and/or  
2. A **single inversion SV call** that covers most of the sandwich  

before labeling `Inversion`; otherwise it reports two `SCE` rows.

### Why many filters exist

Strand-seq chromosomes are full of regions that **look** like state switches but are not SCEs:

| Artifact | Why it fools a naïve switch caller |
|----------|-------------------------------------|
| Centromeres / heterochromatin | Sparse or `None` bins; state stubs at edges |
| Deletions | Homozygous-looking dips inside WC |
| Duplications | Extra WC-like signal with elevated total depth |
| Inversions / complex SVs | Large masked holes; tiny state scraps on either side |
| Acrocentric p-arms | Masking leaves a short homozygous stub before a long WC |

The pipeline therefore **masks unreliable bins**, **drops depth/SV-supported artifacts**, then calls switches on what remains.

---

## Pipeline at a glance

For each QC-passed cell × chromosome (`chrY` is used only to infer sex):

```text
[1] Bin masking
      Drop None / centromere / SV / invalid classes
        ↓
[2] Artifact removal (depth + SV geometry)
      (a)  Duplication-like WC
      (b)  Short stubs at SV hole edges
      (b′) Sparse remnants after large SV holes
      (c)  Lopsided A–WC–A WC islands
      (c′) Short A–WC–A WC scraps beside large SV-tiled gaps
      (d)  Short / deletion-backed / long-shallow WC–A–WC homozygous islands
      (e)  Shallow homozygous tips (deletion-like)
      (f)  Short leading homozygous stubs on acrocentrics
        ↓
[3] Double-switch extraction
      A–B–A → Inversion or two SCEs
      WW→WC→CC / CC→WC→WW → always two SCEs
        ↓
[4] Single-switch calling
      Exactly one valid switch per arm → SCE
      (+ centromere-crossing exception)
        ↓
[5] Breakpoint refinement
      Move stitched breakpoints to the true switch inside an SV hole
        ↓
[6] Recurrence relabeling
      Shared SCE-like breakpoints across cells → Translocation
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
| `-i` | `*.txt.raw.gz` | 200 kb bins: `chrom, start, end, sample, cell, class, c, w` |
| `--sv` | `lenient_filterFALSE.tsv` | Per-cell SV intervals (deletion, duplication, inversion, complex, …) |
| `--qc` | `StrandPhaseR_final_output.txt` | QC-passed cells only these are analyzed |
| `--mask` | `HGSVC.200000.txt` (default) | Low-mappability / `None`-bin mask |
| `--species` | `human` → `chromosome_arm_positions_grch38.txt` | Centromere barriers `[p.End, q.Start)` |

**Sex.** Inferred from chrY read depth. On **male chrX** (hemizygous), `WC` bins are removed and only `WW ↔ CC` transitions are allowed. Deletion SVs on male chrX are **not** masked (hemizygous deletions would otherwise erase the only informative state).

---

## Algorithm (step by step)

### Step 1 — Bin masking

Each 200 kb bin is dropped if any of the following applies.

#### 1.1 None / low-mappability mask (`HGSVC.200000.txt`)

| Rule | Detail |
|------|--------|
| Explicit `None` | Always dropped |
| Short `good` islands | ≤5 bins of `good` fully flanked by `None` → treated as `None` (mappability speckles) |
| Edge absorption | A large `None` run (≥10 bins) also absorbs up to 5 flanking `good` bins on each side so heterochromatin edges do not leave tiny state stubs |

#### 1.2 Centromere barrier

Bins overlapping `[p.End, q.Start)` from the arm table are dropped.  
If a chromosome has no arm-table entry, contiguous `None` stretches ≥1 Mb (merged across gaps ≤1 Mb) are used as a fallback barrier.

#### 1.3 SV intervals

| SV class | Masked? |
|----------|---------|
| Deletion / duplication / inverted duplication (`del_*`, `dup_*`, `idup_*`) | **Always** (except male-chrX deletions — kept) |
| Inversion / complex | Masked by default; **kept** (not masked) only if **tip-linked** |

**Tip-linked** means any of:

- the SV reaches a chromosome tip within ±1 Mb, or  
- it starts within 2 Mb of pter / ends within 2 Mb of qter, or  
- it abuts another tip-linked SV (including tip-linked CN events).

**Why tip-linked inversions are kept.** A tip inversion often *is* the biological event that should remain visible as a state change. Masking it would erase the signal. Interior inversions / complex events are usually masked so their messy interiors do not invent extra SCE breakpoints; their effect is reconsidered later when resolving A–B–A sandwiches.

#### 1.4 Invalid strand class

Keep only `CC` / `WC` / `WW` (male chrX: only `CC` / `WW`).

Surviving bins are merged into **state runs**. Example:

```text
WW WW WW WC WC CC CC  →  WW | WC | CC
```

---

### Step 2 — Artifact removal

After masking, some runs still look like SCE switches but are better explained by copy-number change or mask geometry. Depth terms used below:

| Symbol | Meaning |
|--------|---------|
| `c`, `w` | Median Crick / Watson depth on the **full contiguous raw WC run** (including bins the None mask would skip) |
| `dom` | Dominant-strand depth of flanking homozygous run(s): WW → `w`, CC → `c`. Average of both sides when present; floored by the cell-wide homozygous median when flanks are missing or too shallow. Flanks may be taken **across** an SV/None gap when nothing abuts the WC edge. |
| `balance` | `min(c,w) / max(c,w)` — near 1 means strands are even |
| Expected coverage | Per-bin across-cell median of relative depth × this cell’s genome-wide median. Stops mappability-hot bands (deep in every cell) from looking like gains |

---

#### 2a. Drop duplication-like WC

**Biology.** An unmarked (or incompletely called) duplication often keeps one strand near full homozygous depth and adds signal on the other strand. True copy-neutral WC should sit near **half** depth on each strand and **~1×** the homozygous flank in total.

**Drop if either path fires.**

**Path A — flank asymmetry (all three required):**

1. `max(c,w) / dom ≥ 0.85` — one strand still near homozygous depth  
2. `min(c,w) / dom ≥ 0.35` — the other strand is clearly present  
3. `(c+w) / dom ≥ 1.35` — total depth is elevated vs the flank  

**Path B — coverage ceiling (alone enough):**

- `(c+w) > expected coverage × 1.50`

| | True WC (keep) | Duplication-like WC (drop) |
|--|----------------|----------------------------|
| Strands | both ≈ 0.5 × `dom` | one ≈ `dom`, other ≈ 0.5 × `dom` |
| Total | ≈ `dom` | ≈ 1.5 × `dom` |

---

#### 2b. Drop short stubs at SV hole edges

**Biology.** After a large SV is removed, 1–2 bins of a foreign state often cling to the hole edge and create a false A–B–A sandwich.

**Rule.** Any run of **≤2 bins (≤400 kb)** that **touches** an SV skip interval is dropped so the true flanks can merge.

---

#### 2b′. Drop sparse remnants after large SV holes

**Biology.** `_merge_runs` reports genomic `start`/`end` **across** skipped bins. A scrap on either side of a del/dup/inv block can therefore look megabases long even though almost no kept sequence remains — enough to seed a false SCE.

**Drop a WC or homozygous run when all hold:**

1. The run **touches** an SV skip interval  
2. Kept bin length ≤ **2 Mb** (WC) or ≤ **0.8 Mb** (WW/CC)  
3. Genomic hole (`span − kept`) is **strictly larger** than kept length  

Solid short flanks with `kept ≈ span` are retained.

**How 2b / 2b′ / 2c′ differ** (easy to confuse — they all touch SV holes):

| | 2b | 2b′ | 2c′ |
|--|----|-----|-----|
| Target | Any class | WC or WW/CC run | **A–WC–A only** (`WW-WC-WW` / `CC-WC-CC`) |
| Length idea | ≤ **2 bins (0.4 Mb)** — shortest of the three | Small **kept** vs large **hole** (WC ≤2 Mb / homo ≤0.8 Mb) | WC ≤ **1 Mb** (can be longer than a 2b stub) |
| Geometry | Touches SV | Touches SV and `hole > kept` | Flank gap ≥ **10 Mb** almost fully SV-tiled |
| Role | Dust on the hole edge | Run that looks long only because it spans a hole | Fake dual-SCE from a short WC island after a huge SV desert |

So “short” in 2c′ means short **relative to a real ABA WC island / dual-SCE spacing**, not shorter than 2b. A 0.6–1.0 Mb WC between identical flanks across an ~18 Mb SV gap is a 2c′ case that 2b misses (`>2` bins).

---

#### 2c. Drop lopsided A–WC–A WC islands

**Biology.** A WC island between two identical homozygous runs (`WW–WC–WW` / `CC–WC–CC`) is the classic “two SCE or one inversion” pattern. If the island is really a homozygous stretch with background on the minor strand (lopsided, near-flank total depth), treating it as two exchanges is wrong.

**Keep the island if any escape rule holds; otherwise drop it** (flanks merge). Single WC transitions (not sandwiched) are never touched here.

| # | Keep when | Intuition |
|---|-----------|-----------|
| 1 | `balance ≥ 0.75` | Strands are even enough for real WC |
| 2 | length ≥ 10 Mb **and** `balance ≥ 0.50` | Long islands get a relaxed balance bar |
| 3 | `(c+w) < 0.55 × dom` | Depth-depleted islands are a different artifact — handled elsewhere |
| 4 | Flank-dominant strand inside the island is `≤ 0.47 ×` the deeper flank **and** `balance ≥ 0.35` | Clean halving (copy-neutral exchange) |
| 5 | Both flanks exist, length ≤ 5 Mb, `balance ≥ 0.55`, island total lies **between** the two flanks, dominant strand still `> 0.55 ×` deeper flank, minor strand `≥ 0.45 ×` shallower flank | True exchange next to uneven flanks |

---

#### 2c′. Drop short A–WC–A WC scraps beside large SV-tiled gaps

**Biology.** After a multi-megabase SV/complex mask, a short WC remnant (≤1 Mb) can remain between identical homozygous flanks and invent a nonsense dual SCE a few hundred kb apart. Depth-depleted islands escape rule 2c (keep when `(c+w) < 0.55 × dom`), so a separate **ABA + gap + SV-tiling** rule is needed. This is **not** “shorter than 2b”: 2b already removes ≤2-bin edge stubs; 2c′ targets short **WC islands** next to a **huge SV-tiled hole** (including 3–5 bin scraps that 2b leaves).

**Drop `WW–WC–WW` / `CC–WC–CC` when all hold:**

1. WC length ≤ **1 Mb**  
2. At least one flank gap ≥ **10 Mb**  
3. That gap is almost fully tiled by SV skip (**uncovered ≤ 1 Mb**)  

Depth-carved holes **without** SV tiling (for example inversion sandwiches where WC was trimmed by other filters) are **kept**, so true ABA inversion calls are not destroyed. Overlap with 2b is possible when the WC scrap is also ≤2 bins; 2b runs earlier and usually clears those first.

---

#### 2d. Drop short / deletion-backed WC–A–WC homozygous islands

**Biology.** Mirror of 2c: a short `WW`/`CC` island inside WC (`WC–CC–WC` / `WC–WW–WC`) looks like two SCEs. Two exchanges only a few Mb apart are rare; a real island would stay roughly depth-neutral. A shallow island overlapping a deletion call is usually lost sequence, not two exchanges.

**Short island — drop when all hold:**

1. Island length `< 5 Mb`  
2. Island `(c+w) < deeper flanking WC × 0.80`  
3. Island `(c+w) < shallower flanking WC` (depleted vs **both** sides)

**Deletion-backed island (no length cap) — drop when either holds:**

1. Deletion SV covers ≥ **20%** of the island **and** island `(c+w) < deeper WC × 0.80`  
2. Deletion SV covers ≥ **15%** **and** island `(c+w) < deeper WC × 0.40` (already deeply depleted)

**Long shallow island without deletion support** — drop when all hold:

1. Island length ≥ **5 Mb**  
2. Deletion SV covers **< 15%** of the island (weak/no del call; otherwise the deletion-backed path applies)  
3. Island `(c+w) < deeper WC × 0.40`  
4. Island `(c+w) < shallower WC`

This closes the gap where a long unmarked deletion-like island is too long for the short path and has no `del_*` call for the deletion-backed path. The depth bar matches the soft deletion-backed tier so only deeply depleted islands are removed.

Depth-neutral islands, or islands deeper than one WC flank, are kept.

---

#### 2e. Drop shallow homozygous tips (deletion-like)

**Biology.** A terminal `WW`/`CC` tip next to inward `WC` can be a real tip SCE, or a deletion/noise stub that seeds a false `WW→WC→CC` two-step. Male chrX is skipped.

**A. Deletion-backed (exactly two runs on the chromosome) — drop when all hold:**

1. Tip **abuts** the WC flank (no masked gap)  
2. Tip `(c+w) < WC flank × 0.70`  
3. Deletion SV calls cover ≥ **10%** of the tip in total  
4. At least one deletion lies within **2 Mb** of the WC↔homo junction  

**B. Short shallow tip (any run count) — drop when both hold:**

1. Tip length ≤ **15 Mb**  
2. Tip `(c+w) ≤ WC flank × 0.70`  

Path B catches compact terminal homozygous blocks even when no deletion call sits on the tip.

---

#### 2f. Drop short leading homozygous stubs on acrocentrics

**Biology.** On acrocentric chromosomes (chr13/14/15/21/22), p-arm and centromere masking often leaves a few bins of `WW`/`CC` abutting a long `WC` on the q-arm. That stub looks like a tip SCE but is almost always mask-edge noise.

**Rule.** Leading `WW`/`CC` ≤ **1 Mb** that **abuts** a following `WC` → drop the stub.  
Stubs separated from WC by a masked gap are kept.

---

### Step 3 — Double-switch extraction and resolution

Double switches are removed from the filtered run list and resolved on their own.

#### 3.1 A–B–A sandwich (`WW–WC–WW`, `CC–WC–CC`, …)

Return to the original flank class → either one **inversion** or **two SCEs**.

| Priority | Condition | Result |
|----------|-----------|--------|
| 1 | Both breakpoints shared within ±10 kb by ≥5% of QC cells (same sample·chrom; ≥2 cells) | `Inversion` (+ shared %) |
| 2 | Exactly **one** inversion SV call lies inside the sandwich and covers ≥ **40%** of it | `Inversion` (shared % empty = SV-backed) |
| 3 | Otherwise | Two `SCE` rows |

Crossing the centromere alone does **not** force `Inversion`.

Long sandwiches that merely contain **several** scattered inversion calls stay as two SCEs (rule 2 requires exactly one contained call).

#### 3.2 Two-step opposite (`WW→WC→CC` / `CC→WC→WW`)

Always **two SCEs**. A centromere inside the final homozygous run does not convert this into an inversion. (Not used on male chrX, where WC is absent.)

#### 3.3 Other patterns

- Matching flanks but not a valid SCE sandwich (e.g. `CC–WW–CC`): clear the middle, emit nothing  
- Other `A–B–C` middles: clear the middle  

---

### Step 4 — Single-switch SCE

After double switches are stripped, each chromosome is split at the centromere into arms. For each arm:

- exactly **two** state runs, and  
- a valid SCE transition  

→ one breakpoint = `SCE` (the start of the second run, after refinement in Step 5).

**Centromere-crossing exception.** If the kept flanks on either side of the barrier form a valid SCE transition and **each flank is ≥ 5 Mb**, emit one cross-centromere SCE. Prefer the raw switch inside the barrier when there is exactly one; otherwise use the start of the right flank.

Arms that still have multiple state changes after filtering are left without a single-SCE call (their double switches were already handled in Step 3).

---

### Step 5 — Breakpoint refinement across SV holes

After stitching, the default breakpoint is the start of the first surviving bin **after** a removed gap. If that gap (or the breakpoint) touches an SV, the tool peeks at raw bins **inside** the gap without the SV mask:

| Peek result | Action |
|-------------|--------|
| Exactly one `left → right` switch | Move the coordinate there |
| None or more than one | Keep the default (far / right edge) |

Sandwich breakpoints use the same correction, but **only when SV intervals tile the gap completely**. Without that guard, a chromosome-wide sandwich peek can wander into the centromere or a None desert and land tens of Mb away.

**Why this matters.** Without correction, a stitched breakpoint snaps to the **far edge** of a removed SV. That coordinate can collide with a real breakpoint in another cell and fake a recurrent subclone (`Translocation`).

---

### Step 6 — Recurrence → Translocation

Among collected SCE-like breakpoints (not already resolved as SV-backed-only logic for inversions):

1. Compare only within the same `sample` + `chrom`  
2. Count distinct cells with a breakpoint within **±10 kb**  
3. If that count is **≥ 5%** of QC cells (**and at least 2 cells**) → relabel as `Translocation` and store the shared-cell percent  
4. Rows already labeled `Inversion` are never relabeled as translocation  

Recurrence is evidence of a **shared breakpoint** (clonal SV, recurrent fragile site, or mapping artifact). It is **not** proof of a translocation partner.

---

## Output

Excel table with one row per candidate breakpoint (two rows for dual SCEs; one row with `start`/`end` for inversions).

| Column | Meaning |
|--------|---------|
| `Sample` | Sample name |
| `Cell_ID` | Cell barcode / ID |
| `chr` | Chromosome |
| `start` | Breakpoint (left end for inversions) |
| `end` | Right breakpoint for `Inversion`; empty for SCE / Translocation |
| `Event` | `SCE` / `Translocation` / `Inversion` |
| `Shared_cell_percent` | Shared-cell % for Translocation and recurrent Inversion; empty for ordinary SCE and for **SV-backed** Inversion |

**How to read empty `Shared_cell_percent` on an Inversion**

| `Shared_cell_percent` | Meaning |
|-----------------------|---------|
| Filled | Called from **cross-cell recurrence** of both sandwich ends |
| Empty | Called from the **single inversion SV** covering ≥40% of the sandwich |

---

## Command-line arguments

Only file / species paths are exposed on the CLI. Filter thresholds are fixed in `SCEDetector.py` to the values documented above (tuned against reference outputs such as `SCE_detected.xlsx`).

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
| Single valid switch per arm after filtering | Full generative model of multiple crossovers |
| A–B–A / two-step double-switch resolution | Haplotype-resolved validation of every call |
| Depth- and SV-aware artifact filters | Partner chromosome for translocation candidates |
| Breakpoint recurrence → translocation **candidates** | Proof that a recurrent site is a true rearrangement |

Practical caveats:

- Calls remain **candidates**; always inspect depth and SV context for important events.  
- Large duplications can erase a true nearby SCE if the whole WC block is dropped as duplication-like.  
- True dual SCEs with WC ≤1 Mb beside a large SV-tiled hole are intentionally suppressed (rule 2c′).  
- Thresholds were tuned on HPNE Strand-seq cohorts; new chemistries or bin sizes may need retuning.

---

## License

Use freely for research and analysis.

## Contributor

* Hyunjin Cho

[hyun-jin891](https://github.com/hyun-jin891)
