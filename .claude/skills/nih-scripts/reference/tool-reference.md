# NIH_scripts — per-tool CLI reference

Seven standalone tools, run directly (no pipeline wrapper). All coordinates in standard
bedGraph output are 0-based start / 1-based end. Chromosome-name rewriting is on by default
(disable with `-D`).

---

## 1. trim_and_filter_SE.pl / trim_and_filter_PE.pl
Crop reads to a position window and drop reads below an average base-quality score.

**SE:**
```
trim_and_filter_SE.pl [-h] -i <input_fastq> -a <Keep_Pos1> -b <Keep_PosN> \
  -m <MinimumAvgBQS> -q <illumina|sanger> -o <OutputNameRoot>
```
**PE** (separate crop window per mate):
```
trim_and_filter_PE.pl [-h] -1 <Mate1_fastq> -2 <Mate2_fastq> \
  -a <M1_Pos1> -b <M1_PosN> -c <M2_Pos1> -d <M2_PosN> \
  -m <MinimumAvgBQS> -q <illumina|sanger> -o <OutputNameRoot>
```
| Opt | Meaning |
|--|--|
| `-i` / `-1`,`-2` | input FASTQ (SE / PE mates); `.gz`/`.bz2` read via zcat/bzcat |
| `-a`,`-b` (`-c`,`-d` PE mate2) | first/last read base to **retain** (1-based inclusive) |
| `-m` | minimum average base quality score |
| `-q` | quality encoding: `illumina` (+64, old 1.3–1.5) or `sanger` (+33, modern 1.8+) |
| `-o` | output root |

- **In:** FASTQ. **Out:** `<root>.trim_<a>_<b>.minQS_<m>.fastq` (PE: `.1.`/`.2.` per mate)
  + a `.FilterStats.txt` counting how many reads/pairs failed and why.
- Gotcha: wrong `-q` silently mis-filters; empty input can divide-by-zero.
- `trim_and_filter_SE.pl` was added later (2024) than the long-standing PE script.

---

## 2. bowtie2stdBedGraph.pl
Bowtie-v1 text (or name-sorted SAM) alignment → **5'-end stranded** standard bedGraph.
Uses Perl `threads`.
```
bowtie2stdBedGraph.pl [options] <input file> <output file prefix>
```
| Opt | Meaning |
|--|--|
| `-p` | input is paired-end alignment (reports End1 5' ends by default) |
| `-a` | input is **SAM** (must be sorted by query name) instead of Bowtie text |
| `-o [nbmv]` | extra outputs: `n`=normalized, `b`=binned, `m`=merged, `v`=binned+merged |
| `-b <int≥1>` | bin size (needs `-o b` or `-o v`; default 25) |
| `-s <int≥0>` | bp to shift ChIP-seq hits before merging (needs `-o m`/`v`; default 75) |
| `-l <file>` | chromosome name/length list (prevents shifting past chrom end) |
| `-r <genome>` | fetch chrom lengths from UCSC MySQL instead of `-l` (needs DBI/DBD::mysql) |
| `-t <int≥1>` | nt trimmed prior to alignment |
| `-x` | swap strand outputs (forward↔reverse); for PE, report End2 5' ends |
| `-D` | **disable** `chr`-prefixing and mitochondrial renaming |
| `-M <str>` | mitochondrial-name match string (default `mit`) |
| `-u <int≥1>` | de-duplicate, retaining up to N copies per position (`-u 1` = full dedup) |
| `-i` | strand-independent dedup (requires `-p`) |

- **In:** Bowtie-v1 output / name-sorted SAM. **Out:** `<prefix>_forward.bedGraph`,
  `<prefix>_reverse.bedGraph` (+ `_norm1m_<strand>`, `_binned`, `_merged` variants when
  requested). Output is 0-based start / 1-based end, `chr`-prefixed.
- For Start-seq/startRNA PE: the 5prRNA file (End1) is default; `-x` + `-p` gives the 3prRNA
  (End2) file with swapped strands.

---

## 3. extract_fragments.pl
PE alignment → fragment BED and/or binned (fragment-center) bedGraph, length-filtered.
```
perl extract_fragments.pl [options] <input> <output prefix> <miseq|hiseq>
```
| Opt | Meaning (default) |
|--|--|
| `-a` | input is SAM (sorted by query name) |
| `-o [b\|g\|a]` | output: `b` BED only (default), `g` bedGraph only, `a` both |
| `-min` / `-max` | min/max fragment length to keep (100 / 200) |
| `-b` | bin size for bedGraph (25) |
| `-l <file>` / `-r <genome>` | chrom lengths (file) or UCSC fetch — clamp final bin to chrom end |
| `-D` | disable chr-name fixing |
| `-M <str>` | mitochondrial match string (default `mit`) |
| `-u <int>` | de-duplicate, retain up to N per position |

- **Out:** `<prefix>.bed` (one line per fragment) and/or `<prefix>.bedgraph` (counts of
  fragment centers per bin). bedGraph is standard 0-based start / 1-based end.

---

## 4. bedgraphs2stdBedGraph
Merge/sum all bedGraphs over **identical intervals** found in the working directory into one
standard bedGraph; format inferred from file name. Also converts a single 1-based `.bedgraph`
→ standard `.bedGraph`.
```
bedgraphs2stdBedGraph <output prefix>
```
- Operates on files in the cwd (no explicit input args). Writes the merged standard bedGraph
  + a log listing the inputs used.
- Name convention it relies on: `.bedgraph` = 1-based start (shifted −1 to standardize),
  `.bedGraph` = already standard 0-based.

---

## 5. normalize_bedGraph
Scale the value column (col 4) of one or more bedGraphs by a constant; factor is embedded in
the output filename for provenance.
```
perl normalize_bedGraph <factor> <multiply|divide> <bg1> [bg2 ... bgN]
```
- **Out:** `<name>_normM<factor>.bedGraph` (multiply) or `<name>_normD<factor>.bedGraph`
  (divide). e.g. `normalize_bedGraph 0.85867 multiply S2_F.bedGraph` → `S2_F_normM0.86.bedGraph`.

---

## 6. make_heatmap (C++)
Gene-list × bin metagene/heatmap matrix: for each gene-list feature and each relative bin,
compute total/average/density/coverage of hit-file values intersecting that bin.

**Four command-line forms:**
```
make_heatmap [opts] <hit file> <gene list> <output> <bin file>          # explicit bins
make_heatmap [opts] -b c ... <output> <bin start> <bin size> <bin count># fixed-size bins
make_heatmap [opts] -b v ... <output> <bin count>                       # variable bins/feature
make_heatmap [opts] -p <plus hits> -m <minus hits> <gene list> ...      # split-strand hits
```
Negative bin starts require a `--` before positional args:
```
make_heatmap -b c -p F.bedgraph -m R.bedgraph -- gene_list.txt out.matrix -100 1 201
```

### Option matrix
| Opt | Arg (default) | Choices |
|--|--|--|
| `-p` / `-m` | file | plus- / minus-strand hit file (for strand-less hit formats) |
| `-t` | int (1) | threads |
| `-b --bintype` | (`c`) | `c` cmd-line fixed size · `v` cmd-line variable · `f` bin file |
| `-h --hittype` | (`G`) | `G` std bedGraph (0-based) · `g` bedgraph (1-based) · `b` basic bed · `e` extended bed · `c` cppmatch |
| `-l --hitloc` | (`p`) | `p` phys start · `d` phys end · `s` genetic start* · `e` genetic end* · `c` center |
| `-a --anchor` | (`s`) | `s` genetic start* · `e` genetic end* · `p` phys start · `d` phys end · `u` user col-2 |
| `-v --binvalue` | (`t`) | `t` total · `a` average · `d` density (total/binsize) · `c` mean per-nt coverage |
| `-d --binloc` | (`g`) | `g` genetic distance* · `p` physical distance |
| `-s --strands` | (`b`) | `b` both · `s` same-strand* · `o` opposite-strand* |
| `--nostrand` | — | strand-agnostic shortcut: applies `-s b -l p -a p -d p` |
| `--nohead` | — | suppress output header |

\* Genetic/same/opposite-strand options require strand info — i.e. `-p`/`-m`, or hit type
`c` (cppmatch) or `e` (extended bed) — and a stranded gene list.

### make_heatmap file formats
- **Gene list** (TSV, 5 required + 1 optional col):
  `[Gene ID] [Description or user anchor] [Chromosome] [Physical Start] [Physical End] [Strand]`.
  Gene ID must be unique; with `-a u` col 2 must be an integer anchor inside start..end;
  strand ∈ `plus|minus|+|-`.
- **Hit file** types: standard bedGraph `G` (0-based start), bedgraph `g` (1-based start, for
  compat with bowtie2bedgraph/extract_fragments output), basic bed `b` (chrom,start,end),
  extended bed `e` (≥5 fields; only strand is used beyond that), cppmatch `c`
  (`ID value chrom physStart physEnd strand`). Chromosome IDs must match the gene list exactly.
- **Bin file** (`-b f`): two tab cols `start<TAB>end`, may be negative/gapped/variable-size,
  no overlaps, end ≥ start.
- **Output:** header (options + inputs) then the matrix — the 6 gene-list columns + one
  column per bin.

### Compiling
```
g++ -O3 -o make_heatmap make_heatmap.cpp -lpthread            # Linux
g++ -O3 -o make_heatmap make_heatmap_OSX.cpp -stdlib=libstdc++ # macOS (make_heatmap_OSX.cpp)
g++ -O3 -o make_heatmap make_heatmap.cpp -DSINGLE             # no-pthreads single-threaded
```

---

## Bundled chromosome-size files
Two-column `name<TAB>length`. Shipped for **hg19, hg38, mm9, mm10, dm3**, duplicated in both
`bowtie2stdbedgraph/` and `extract_fragments/` (e.g. `mm10_chr_size.txt`, `dm3_chr_size.txt`).
Pass with `-l` to keep bins/shifts from running off a chromosome end; prefer these over
`-r <genome>` UCSC MySQL fetches on compute nodes (no network, no DBD::mysql dependency).

---

## Evolution, rationale & gotchas
- **Provenance & stability:** written at NIEHS (Burkholder/Grimm/Fargo) for early
  Adelman-lab PRO-seq/Start-seq/ChIP-seq papers (Nechaev 2010 → Henriques 2018). ~56 commits,
  mostly stable 2010–2016 "boring old code"; primary later maintenance BenjaminMartin with
  K. Adelman / G. Nelson / sgoldman101. `trim_and_filter_SE.pl` is the newest addition
  (2024-05) — the PE variant predates it. Treat the core Perl/C++ as a frozen, trusted
  dependency rather than something to refactor.
- **Committed make_heatmap binary is old:** built for Red Hat Enterprise Linux 5.11 with
  gcc 4.1.2 (the shipped ELF targets GNU/Linux 2.6.18). It is checked in for convenience but
  can fail to load against newer glibc/libstdc++ — recompile from `make_heatmap.cpp` when it
  does, rather than assuming the tool is broken.
- **Chr-name rewriting is a deliberate UCSC-compat hack:** bowtie2stdBedGraph /
  extract_fragments prepend `chr` and rename `mit`→`chrM` so alignments done against bare-
  numbered indexes (e.g. the default fly genome) load in the UCSC browser. It is silent and
  on by default — the frequent source of chromosome-mismatch errors between NIH_scripts output
  and other tools. Use `-D` to turn it off; `-M` to change the mito match string.
- **`.bedgraph` (1-based) vs `.bedGraph` (0-based) is a real coordinate convention, encoded in
  the file extension's case.** The legacy 1-based `.bedgraph` exists only to stay compatible
  with older bowtie2bedgraph/extract_fragments output; `make_heatmap` accepts it via `-h g`,
  and `bedgraphs2stdBedGraph` standardizes it to 0-based. Prose in some READMEs says
  `.bedgraph` where the code actually writes `.bedGraph` — trust the extension semantics and
  the code, not the prose, and keep the case correct when handing files to downstream tools.
