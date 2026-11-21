# NIH_scripts — revision changelog (for the NTC ecosystem manager)

<!-- THIS FILE IS OWNED BY THIS REPO. The ecosystem manager (ntc-ecosystem-context) READS it to
track pipeline revisions and keep its overseer context current; it does NOT overwrite it.
Append a new entry at the TOP whenever you make a notable change (a revision / version bump /
results-affecting fix), so the manager has visibility. Keep entries to one or two lines. -->

Entry format (newest first):
```
## <YYYY-MM-DD> · <version/tag or claude/branch> · <NTC-### if any>
- <one-line summary>; affects results: <yes|no>
```

<!-- entries below, newest first -->

## 2026-11-21 · claude/make_heatmap-parallel
- Added make_heatmap_optimized_parallel.cpp (additive; original + _optimized untouched): multi-threaded, BYTE-IDENTICAL to the original at EVERY thread count (not just -t1). Parallelizes parse+match+output while replaying each cell's bin updates in strict input-file order, so float summation order == -t1 (unlike the original/_optimized -tN, whose single global accumulation lock makes -tN float sums thread-interleaving-dependent AND — now that _optimized's windowed lookup made match cheap — actually SLOWER than -t1). Motivation grounded in actual usage: pipeline2/pipeline.sh and pauseIndex both call make_heatmap with `-t 6` on `-c 6` sbatch jobs. Measured (12M hits × 20k TSS, 4 cores): original -t6 43s → _optimized -t1 4.6s (9.4×) but _optimized -t6 6.8s (lock regression) → _parallel -t6 3.3s (13× vs original -t6, ~2× vs _optimized -t6), as a drop-in needing no pipeline edits. make_heatmap is memory-bandwidth-bound so scaling is ~1.4–2× over _optimized serial (sub-linear); a functional/per-thread-accumulate variant was prototyped and rejected (tied on speed, 2× RAM, lost byte-identity). + README_optimized_parallel.md, tests/run_parallel.sh (byte-identity at -t1/2/4 across the full option matrix incl. multi-index), NUMERICS.md updated. affects results: no (byte-identical at every -t)

## 2026-11-21 · claude/revert-summation-variant
- Removed make_heatmap_optimized_summation.cpp + its harness (run_summation.sh, compare_matrix.py). Benchmarking (hyperfine) showed it was ~20% SLOWER single-threaded than _optimized and only ~1.3–1.5× on ≤6 cores (make_heatmap is not accumulation-bound — the costs are match ~47% and output formatting ~40%, with read ~1%), and it altered output bits. Strictly dominated by _optimized for the lab's ≤6-core use. Added tests/benchmark.sh (hyperfine matrix: value-type, depth, gene-count, binning, threads). Docs updated to two builds. affects results: no (removes an unused functional-identity variant; _optimized unchanged)

## 2026-11-21 · claude/make_heatmap-summation-variant
- Added make_heatmap_optimized_summation.cpp (3rd variant; relaxes byte-identity to FUNCTIONAL identity to allow a better summation strategy): lock-free per-thread accumulation + Neumaier compensated summation + batched reader. Results are thread-count-independent & run-to-run reproducible (-t1=-t2=-t8 byte-identical), more numerically robust (handles signed/cancellation), and functionally identical to the original (rel ~1e-13 internally; <=1 unit in the 6-sigfig printout; integer/count matrices exact). Multicore speedup is MODEST (~1.5x on 4 cores; capped by serial input reading — the big win is still the shared Phase 1-3 algorithmic work). Per-thread memory ~24B x features x bins x threads. + NUMERICS.md (lab doc: summation-order reproducibility, the 3 builds' guarantees) and tests/run_summation.sh + compare_matrix.py (functional-tolerance + thread-independence harness; 25/25). affects results: functionally identical, NOT byte-identical — use make_heatmap_optimized for bit-match to historical output

## 2026-11-20 · claude/make_heatmap-modernize-cpp
- make_heatmap_optimized.cpp: replaced the deprecated GNU `<tr1/unordered_map>` (and `using tr1::`) with the standard C++11 `<unordered_map>`/`std::unordered_map` — the same thing make_heatmap_OSX.cpp already did, so one source now compiles on Linux and macOS/clang (no OSX variant, no `-stdlib=libstdc++`). Updated README_optimized.md build flags (`-std=c++11 -pthread`; don't use `-ffast-math`/`-march=native`). Validated byte-identical: golden-master harness 108/108 (unordered_map iteration order never reaches output). affects results: no

## 2026-11-20 · claude/add-make_heatmap_optimized
- Added make_heatmap_optimized.cpp (additive; original untouched): byte-identical drop-in at -t 1 with three output-preserving optimizations — (1) windowed binary-search feature lookup replacing the O(hits×genes) per-hit linear scan (the big win), (2) a cached per-gene table pointer removing a std::map<string> lookup from the hot loop, (3) manual strtol/strtod line parsing replicating operator>> exactly (no per-line istringstream). Measured 2.1× (200 genes) → 6.3× (2k) → 13.1× (20k/GGA-scale), 3M hits, -t1; no small-list regression. + README_optimized.md and a golden-master harness tests/ (108 checks: full option matrix + adversarial parser fuzz; builds reference from make_heatmap.cpp and byte-diffs). Threading/global-lock left unchanged (so -tN float sums are, as before, not bit-equal to -t1). Harness also surfaced two PRE-EXISTING original bugs (documented, not changed): bin files must not end in a trailing newline; a file of very many consecutive bad lines can overflow the stack via update() recursion. affects results: no (byte-identical at -t1)

## 2026-11-20 · claude/add-extract_fragments_optimized
- Added extract_fragments_optimized.pl (additive; original untouched): BED sort rewritten as a Guttman-Rösler transform (pack+plain-sort instead of arrayref-Schwartzian) — ~1.57x faster + ~1.57x lower peak RSS on the BED-heavy path, and it makes the BED output DETERMINISTIC by adding `end` as the final sort key (the original emitted equal-start fragments in randomized hash order). bedGraph output byte-identical; BED is the same line-set, now reproducible. + README_optimized.md, OPTIMIZATION_PLAN.md, tests/ (54/54 pass). affects results: BED line *order* only (same set of lines; now stable run-to-run) — bedGraph unchanged

## 2026-11-20 · claude/bedgraphs2stdBedGraph_optimized_v2
- Added bedgraphs2stdBedGraph_optimized (streaming k-way merge + atomic write-temp-then-rename): byte-identical .bedGraph output to the original with far lower peak RAM on large merges, plus a guard for unsorted input. Additive (original untouched) + byte-diff harness. Re-validated: 25/25 pass (incl. 1.2M-record input). affects results: no

## 2026-11-18 · claude/add-trim_and_filter_optimized
- Added trim_and_filter_PE_optimized.pl + trim_and_filter_SE_optimized.pl (byte-identical drop-in additions; originals unchanged): quality-average hot loop rewritten with unpack("%32C*") instead of per-base substr/ord, ~8.5x faster; + README_optimized.md and tests/ (byte-diff + throughput incl. cutadapt reference); affects results: no
