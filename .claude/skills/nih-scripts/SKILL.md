---
name: nih-scripts
description: >
  Use when working with the NTC low-level NGS toolbox — converting Bowtie/SAM
  alignments to 5'-end stranded bedGraphs, merging/summing bedGraphs on identical
  intervals, normalizing a bedGraph by a scale factor, trimming+quality-filtering
  FASTQ, extracting PE fragments, or building gene-by-bin metagene/heatmap matrices
  with make_heatmap — or when a task in a related NTC repo calls one of these tools
  (often as a hardcoded /n/data1 binary): GeneAnnotationScripts, pipeline2, DOBBE_F,
  metagene_heatmap_tool, pauseIndex, proTSScall, prepareBrowser, peakMorph_toolkit all
  consume them; alignments come from RNAseqMapping/pipeline2. Covers the 7 standalone
  Perl/C++ scripts, their CLIs, the bundled chrom-size files, and the bedGraph
  coordinate conventions shared across the ecosystem. Aliases: make_heatmap,
  bowtie2stdBedGraph, bedgraphs2stdBedGraph, normalize_bedGraph, trim_and_filter,
  extract_fragments, NIH_scripts.
---

# NIH_scripts (NTC shared toolbox)

## Purpose
NIH_scripts is the **shared low-level utility toolbox** of the NTC ecosystem: a
collection of standalone Perl scripts plus one C++ program (`make_heatmap`) for the
common steps of PRO-seq/Start-seq/ChIP-seq analysis — FASTQ trimming, alignment→bedGraph
conversion, bedGraph merging/normalization, and metagene/heatmap matrix building. Written
at NIEHS (Adam Burkholder, Sara Grimm, David Fargo) for Adelman-lab studies. Repo:
`NIH_scripts` (branch `master`); upstream `github.com/AdelmanLab/NIH_scripts`. Deployed
on O2 and frequently referenced by **hardcoded absolute `/n/data1` paths** from other NTC
repos (not a conda package).

## Architecture
**Not a pipeline — not smartSlurm.** There is no `runAsPipeline`, no `ssbatch`, no flag
folder, no DAG. It is **7 independent command-line tools**, each run directly:

- `trim_and_filter_SE.pl`, `trim_and_filter_PE.pl` — crop reads to a position window +
  average-base-quality filter (FASTQ → FASTQ + FilterStats).
- `bowtie2stdBedGraph.pl` — Bowtie-v1 text or name-sorted SAM → 5'-end **stranded**
  standard bedGraph (`_forward`/`_reverse`); optional normalized / binned / merged /
  ChIP-shifted outputs; threaded Perl.
- `extract_fragments.pl` — PE alignment → fragment BED and/or binned bedGraph, length
  filtered (`-min`/`-max`).
- `bedgraphs2stdBedGraph` — merge/sum all bedGraphs over identical intervals in the cwd
  into one standard bedGraph (also converts 1-based→0-based).
- `normalize_bedGraph` — multiply/divide col-4 values by a factor (→ `_normM`/`_normD`).
- `make_heatmap` (C++) — gene-list × bin metagene/heatmap matrix from hit files.

Full per-tool CLI, option matrix, and I/O: `reference/tool-reference.md`.

## I/O contract
**Inputs:** FASTQ (phred+33 *sanger* or +64 *illumina* — the `-q` scale **must** match the
data); Bowtie-v1 text output or name-sorted SAM (`-a`); standard (`.bedGraph`, 0-based) or
1-based (`.bedgraph`) bedGraphs; `make_heatmap` gene-list TSV
(`ID, anchor/desc, chr, start, end, strand`); chromosome-size files (two-column
`name<TAB>length`).

**Outputs:** trimmed FASTQ + `.FilterStats.txt`; stranded standard bedGraphs (0-based
start, 1-based end, `chr`-prefixed, `_forward`/`_reverse`, optional `_norm1m_`,
`_binned`, `_merged`); fragment `.bed` / `.bedgraph`; merged / summed bedGraph;
normalized bedGraph; heatmap matrix TSV (6 gene-list columns + one column per bin).

Typical invocations:
```
trim_and_filter_PE.pl -1 R1.fq -2 R2.fq -a 1 -b 36 -c 1 -d 36 -m 20 -q sanger -o SAMPLE
bowtie2stdBedGraph.pl -a -o nm -b 25 -l mm10_chr_size.txt alignment.sam SAMPLE
bedgraphs2stdBedGraph merged_prefix          # operates on cwd
perl normalize_bedGraph 0.859 multiply S2_F.bedGraph S2_R.bedGraph
make_heatmap -p F.bedGraph -m R.bedGraph genelist.txt out.matrix -b c -- -500 5 201
```

## Dependencies & environment
- **No conda env, no modules required.** Perl 5 (with `threads` for
  `bowtie2stdBedGraph.pl`); `zcat`/`bzcat` for gzipped/bzipped FASTQ.
- `DBI` + `DBD::mysql` **only** if you use `-r <UCSC genome>` to fetch chromosome lengths
  from `genome-mysql.cse.ucsc.edu`. On compute nodes prefer `-l <local chrom-size file>`
  to avoid network calls.
- `make_heatmap` is C++: a **committed prebuilt Linux binary** is shipped; recompile with
  `g++ -O3 -o make_heatmap make_heatmap.cpp -lpthread` if it won't run.
- Reference data: bundled chrom-size files for **hg19, hg38, mm9, mm10, dm3** (in both
  `bowtie2stdbedgraph/` and `extract_fragments/`). Full list in `reference/tool-reference.md`.

## Modifying it safely (gotchas & assumptions)
- **No rerun/resume/flag folder** — these are one-shot scripts; wrap them yourself for
  reproducibility.
- **Chromosome names are silently rewritten** by default: `chr` is prepended to bare names
  and anything matching `mit` → `chrM` (match string `-M`, override). Disable with **`-D`**.
  This is the usual cause of "my chr names don't match" bugs between tools.
- **`.bedgraph` vs `.bedGraph` is a load-bearing distinction**: lowercase `.bedgraph` =
  **1-based** start (legacy bowtie2bedgraph/extract_fragments); `.bedGraph` = **standard
  0-based** start. `bedgraphs2stdBedGraph` and `make_heatmap -h g/-h G` infer/choose by
  this. README of `bedgraphs2stdBedGraph` says `.bedgraph` but the bowtie script writes
  `.bedGraph` — trust the extension semantics, not prose.
- **Quality scale must match the data** (`-q`): modern Illumina 1.8+ is **sanger** (+33);
  only old Illumina 1.3–1.5 is `illumina` (+64). Wrong scale silently corrupts filtering.
- `make_heatmap` negative bin starts need a `--` separator before positional args (see
  reference).
- Trim scripts can divide-by-zero on empty input; committed `make_heatmap` binary was built
  for RHEL5 and may fail on newer libc (recompile).

## Integration points
- **Upstream:** alignments from `RNAseqMapping` / `pipeline2` (Bowtie/SAM); trimmed FASTQ
  from the `trim_and_filter` scripts feeds the aligner.
- **Downstream / consumers:** `GeneAnnotationScripts` (calls `make_heatmap` for metagene QC,
  `bedgraphs2stdBedGraph` for 5' bedGraph inputs), `pipeline2`, `DOBBE_F`,
  `metagene_heatmap_tool`, `pauseIndex`, `proTSScall`, `prepareBrowser`, `peakMorph_toolkit`
  — several reference the binaries by hardcoded `/n/data1` path or vendor copies.
- **Shared formats:** the stranded standard-bedGraph convention and the `make_heatmap`
  gene-list format (same shape as GGA `formakeheatmap`) are used across the ecosystem.

## Reference files
- `reference/tool-reference.md` — per-tool CLI (all 7), key options & I/O, bundled
  chrom-size files, the full `make_heatmap` option matrix, and evolution/rationale/gotchas.
- `.claude/reference/runAsPipeline-grammar.md` — smartSlurm grammar (NIH_scripts is *not*
  rAP, but its many consumer pipelines are; useful context).
