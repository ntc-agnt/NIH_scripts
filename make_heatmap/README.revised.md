# make_heatmap

`make_heatmap` summarizes a genome-wide **signal** (e.g. PRO-seq / ChIP-seq
coverage, or a set of mapped intervals) over a list of **features** (genes, TSSs,
peaks, …), producing a **feature × bin matrix**: one row per feature, one column
per positional bin, each cell holding the aggregated signal in that bin for that
feature. That matrix is what metagene profiles and heatmaps are built from.

```
            bin -2   bin -1   bin 0(anchor)  bin +1   bin +2
 gene A       3.0      5.5        12.1         4.0      1.0
 gene B       0.0      1.2         8.7         6.3      2.1
 gene C       ...
```

Columns are placed **relative to an anchor point** in each feature (e.g. the TSS),
so every feature is measured on the same ruler and the rows can be averaged
(metagene) or stacked (heatmap).

---

## Terminology

This tool predates some now-standard names. The equivalents:

| This document | In the original help text | What it is |
|---|---|---|
| **signal track** / signal file | "hit file" | the data being summarized — a bedGraph of coverage, or a BED/cppmatch set of intervals |
| **feature list** | "gene list file" | the regions that become matrix **rows** (genes, TSSs, peaks …) |
| **anchor** | anchor | the reference point *inside each feature* that column offsets are measured from |
| **bin** / column | bin | one output column: a window at a fixed offset from the anchor |
| **value** | bin value | the number put in a cell (sum / average / density / coverage of the signal in that bin) |

---

## Run modes

How you invoke `make_heatmap` depends first on **how the signal is supplied**.
There are two modes; the first has three sub-forms (by how the bins are defined).

### Mode 1 — a single signal file

The signal file is the first positional argument. Pick **one** of
three types of trailing argument depending on how you want to define the bins

| Sub-mode | Bins | Command |
|---|---|---|
| **1a — uniform** (`-b c`, default) | `count` equal bins of width `size`, the first starting `start` from the anchor | `make_heatmap [opts] -b c -- SIGNAL FEATURES OUT start size count` |
| **1b — scaled** (`-b v`) | each feature's own span split into `count` bins (features normalized to the same number of columns) | `make_heatmap [opts] -b v SIGNAL FEATURES OUT count` |
| **1c — explicit** (`-b f`) | bins read from a file of `start<TAB>end` offsets (gaps, varying widths, negatives allowed) | `make_heatmap [opts] -b f SIGNAL FEATURES OUT BINFILE` |

- `SIGNAL` = signal track, `FEATURES` = feature list, `OUT` = output matrix.
- **The `--` separates options from positional arguments.** It is needed only when
  a positional value can begin with `-` — i.e. `-b c` with a **negative** `start`
  (e.g. `-2000`), which otherwise gets parsed as an option. `-b v` and `-b f` never
  need it (their first positional is a filename, and `count` is positive).
- `-b c` and `-b f` place bins at **fixed offsets from an anchor** (good for
  "signal ±2 kb around the TSS"). `-b v` instead **divides each feature** (good for
  "signal across the whole gene body, 5′→3′"); the anchor is ignored in `-b v`.

### Mode 2 — strand-split signal (`-p` / `-m`)

The signal is supplied as **two files, one per strand**: `-p PLUS -m MINUS`. The
feature file is then the first positional argument (there is no positional signal
file). The three bin sub-forms are identical; only the leading arguments change:

```bash
# fixed size bins (end coord = binSize X numBins; e.g. 25nt bins from -150/+150 bp relative to TSS)
# '--' is required here: StartCoord can be negative (e.g. -150), which would otherwise be read as an option
make_heatmap [opts] -p PLUS -m MINUS -b c -- FEATURES OUT StartCoord BinSize NumberOfBins

# Divide features into a fixed number of bins regardless of size (e.g. 100 bins per gene)
# no '--': NumberOfBins is positive and FEATURES is a filename
make_heatmap [opts] -p PLUS -m MINUS -b v    FEATURES OUT NumberOfBins

# count across bins specified in file
# no '--': the trailing argument is a filename
make_heatmap [opts] -p PLUS -m MINUS -b f    FEATURES OUT BINFILE
```

**When you need Mode 2.** bedGraph and basic-BED signal carry no strand column, so
from a single such file the tool can only do non-strand-specific, physical-position
matching. To do **strand-specific matching** (`-s s`/`-s o`) or to measure a signal
interval by its **genetic** (5′/3′) end (`-l s`/`-l e`) with strandless data, split
the signal into a plus-strand file and a minus-strand file and pass them with
`-p`/`-m` — the flag assigns the strand. (Extended BED and cppmatch already carry a
strand column, so they can do strand-specific work from a single Mode-1 file.)
`-p`/`-m` is the usual form for PRO-seq/GRO-seq, whose coverage is already in
separate `+`/`−` bedGraphs.

---

## Choosing what gets measured

These options apply in both modes. Defaults in **bold**.

### Anchor — where columns are measured from (`-a`)
The point in each feature that bin offset `0` refers to:
`s` genetic start (**default**, the feature's 5′ end), `e` genetic end (3′ end),
`p` physical start (leftmost), `d` physical end (rightmost), `u` user-defined
(an integer in **column 2** of the feature file). Genetic anchors (`s`/`e`) need a
strand column on the features.

### Measured point — how a signal interval is reduced to one position (`-l`)
Each signal interval is collapsed to a single point before being binned (except in
coverage mode, below): `p` physical start (**default**), `d` physical end,
`s` genetic start (5′), `e` genetic end (3′), `c` center. Genetic points (`s`/`e`)
need strand on the signal — i.e. Mode 2, or a stranded signal format.

### Bin direction — genetic vs physical (`-d`)
`g` genetic (**default**): offsets run 5′→3′, so `+100` means 100 nt **downstream
in the direction of transcription** for both strands (for `−` features the tool
subtracts from the anchor). `p` physical: offsets always run left→right along the
chromosome. Use `g` to align `+`/`−` features consistently; `g` needs feature
strand.

### Value — what goes in each cell (`-v`)
`t` total (**default**): sum of signal values in the bin (for BED input every
interval counts as 1, so this is a **count**). `a` average: total ÷ number of
intervals in the bin. `d` density: total ÷ bin width. `c` coverage: treats each
signal interval as spanning its full length and adds `value × (overlapping nt)` to
every bin it crosses, then divides by bin width — i.e. mean per-nt coverage; `-l`
is ignored in this mode.

### Strand matching (`-s`) and `--nostrand`
`b` both (**default**): every signal interval is matched to features on either
strand. `s` same-strand: `+` signal → `+` features only, `−` → `−`. `o` opposite.
`s`/`o` require strand on **both** the signal (Mode 2, or a stranded format) and the
features. `--nostrand` is a shortcut for fully strand-agnostic operation
(`-s b -l p -a p -d p`) and reads features without a strand column.

### Other
`-t N` worker threads (default 1). `--nohead` omit the descriptive header block.
`--help` usage.

---

## Input formats

### Signal track (`-h` selects the format)
All use UCSC coordinates: **0-based start, 1-based (exclusive) end**, except where
noted.

| `-h` | Format | Columns | Value used | Strand |
|---|---|---|---|---|
| `G` | standard bedGraph (**default**) | `chr  start  end  value` | the `value` column | none |
| `g` | bedGraph, **1-based** start | `chr  start  end  value` | the `value` column | none |
| `b` | basic BED | `chr  start  end` | 1 per interval | none |
| `e` | extended BED | `chr  start  end  name  score  [strand]` | 1 per interval (`score` ignored) | from col 6 |
| `c` | cppmatch | `id  value  chr  start  end  [strand]` | the `value` column | from col 6 |

- `G` vs `g` differ only in the start coordinate convention (the historical
  `bowtie2bedgraph`/`extract_fragments` tools emit 1-based starts → use `g`).
- BED formats contribute a count (value = 1), so `-v t` over BED counts intervals.
- cppmatch strand and extended-BED strand accept `+`, `-`, `plus`, or `minus`.
- A signal file may be a **comma-separated list** (see *Multiple signal files*).

### Feature list
Tab-separated, one feature per row — these six columns become the first six output
columns:

```
ID    Description/anchor    Chromosome    PhysStart    PhysEnd    [Strand]
```

- **ID** must be unique (duplicates are skipped with a warning).
- **Column 2** is free text *unless* `-a u`, in which case it must be an integer
  anchor position inside `[PhysStart, PhysEnd]`.
- **Strand** (`+`/`-`/`plus`/`minus`) is required for any genetic option
  (`-a s|e`, `-l s|e`, `-d g`, `-s s|o`); omit it only with physical/`--nostrand`
  operation.
- `PhysStart ≤ PhysEnd`; features may overlap.

### Bin file (`-b f`)
Two tab-separated columns, `start<TAB>end`, one bin per row, as **offsets from the
anchor**:

```
-100   -1
0      99
100    299
400    500
```

Gaps, differing widths, and negative coordinates are allowed; bins may **not**
overlap, and each `end ≥ start`. **Do not leave a trailing newline** — the parser
treats a blank final line as a formatting error. (This is a quirk of the original;
the bundled `*_chr_size.txt` files follow it.)

---

## Output

A header block (suppressed by `--nohead`) echoing the chosen options and input
filenames, then a column header row, then one row per feature:

```
Gene ID   Description   Chromosome   Gene Start   Gene End   Strand   <bin 1> <bin 2> ...
```

- For `-b c`/`-b f`, each bin column is labeled by its offset range `start:end`;
  for `-b v`, by its index `0 … count-1`.
- With genetic bin direction (`-d g`), the bins of `−`-strand features are written
  in reversed order so every row reads 5′→3′ left-to-right.

---

## Multiple signal files (advanced)

The signal argument (or each of `-p`/`-m`) may be a **comma-separated list** of
files. Two behaviors:

- **One output file** (default): the matrix gains one block of bin columns **per
  input file**, side by side — a quick way to tabulate several tracks against the
  same features.
- **N output files** (comma-separated `OUT` list, count matching the inputs): each
  input file is written to its own matrix.

---

## Examples

### Single signal file (Mode 1)

ChIP-seq coverage (one bedGraph) over whole gene bodies scaled to 50 bins:
```
make_heatmap -h G -b v -s b -a s -d g sample.bedGraph genes.txt out.matrix 50
```

Count mapped BED intervals in custom bins around peak centers, strand-agnostic:
```
make_heatmap -h b -b f --nostrand reads.bed peaks.txt out.matrix bins.txt
```

Mean per-nucleotide coverage (`-v c`) in 25 nt bins around the TSS:
```bash
# '--' required: the -500 StartCoord is negative
make_heatmap -h G -b c -v c -a s -d g -- sample.bedGraph genes.txt out.matrix -500 25 40
```

### Strand-split signal (Mode 2, `-p` / `-m`)

One example per bin sub-form — PRO-seq coverage in separate `+`/`−` bedGraphs,
same-strand matching, 5′ end of each signal interval (`-l s`):

`-b c` — fixed 100 nt bins, 2 kb each side of the TSS:
```bash
# '--' required: the -2000 StartCoord is negative
make_heatmap -p fwd.bedGraph -m rev.bedGraph -b c -s s -a s -l s -d g -- \
    genes.txt out.matrix -2000 100 40
```

`-b v` — each gene body scaled to 100 bins (anchor ignored; 5′→3′ via `-d g`):
```
make_heatmap -p fwd.bedGraph -m rev.bedGraph -b v -s s -l s -d g genes.txt out.matrix 100
```

`-b f` — custom bins from a file (offsets relative to the anchor):
```
make_heatmap -p fwd.bedGraph -m rev.bedGraph -b f -s s -a s -l s -d g genes.txt out.matrix bins.txt
```

---

## Defaults at a glance

`-h G` (standard bedGraph) · `-b c` (uniform bins) · `-s b` (both strands) ·
`-a s` (genetic start / TSS anchor) · `-l p` (physical start of each signal
interval) · `-d g` (genetic bin direction) · `-v t` (total) · `-t 1`.

Note the defaults assume **stranded features** (because `-a s` and `-d g` are
genetic). For an unstranded feature list, add `--nostrand` or use `-a p -d p`.

---

## Build

```
g++ -O3 -o make_heatmap make_heatmap.cpp -lpthread          # Linux
g++ -O3 -o make_heatmap make_heatmap_OSX.cpp -stdlib=libstdc++   # macOS
g++ -O3 -DSINGLE -o make_heatmap make_heatmap.cpp           # no pthreads (single-threaded)
```

A faster, byte-identical drop-in (`make_heatmap_optimized.cpp`) is available in
this directory — see `README_optimized.md`.

---

## Caveats

- **Strandless signal + single file** limits you to `-s b` and physical hit
  locations; use `-p`/`-m` (Mode 2) or a stranded format for anything else.
- **Bin files must not end in a trailing newline** (see above).
- A signal file consisting of a very large number of **consecutive malformed
  lines** can exhaust the stack (each bad line recurses); well-formed input is
  unaffected.
- Under `-t N > 1`, floating-point bin totals may differ in the last digits from a
  `-t 1` run (summation order depends on thread interleaving). Use `-t 1` when you
  need bit-reproducible output.
