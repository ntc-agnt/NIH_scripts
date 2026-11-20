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

## 2026-11-20 · claude/add-make_heatmap_optimized
- Added make_heatmap_optimized.cpp (additive; original untouched): byte-identical drop-in at -t 1 with three output-preserving optimizations — (1) windowed binary-search feature lookup replacing the O(hits×genes) per-hit linear scan (the big win), (2) a cached per-gene table pointer removing a std::map<string> lookup from the hot loop, (3) manual strtol/strtod line parsing replicating operator>> exactly (no per-line istringstream). Measured 2.1× (200 genes) → 6.3× (2k) → 13.1× (20k/GGA-scale), 3M hits, -t1; no small-list regression. + README_optimized.md and a golden-master harness tests/ (108 checks: full option matrix + adversarial parser fuzz; builds reference from make_heatmap.cpp and byte-diffs). Threading/global-lock left unchanged (so -tN float sums are, as before, not bit-equal to -t1). Harness also surfaced two PRE-EXISTING original bugs (documented, not changed): bin files must not end in a trailing newline; a file of very many consecutive bad lines can overflow the stack via update() recursion. affects results: no (byte-identical at -t1)

## 2026-11-20 · claude/add-extract_fragments_optimized
- Added extract_fragments_optimized.pl (additive; original untouched): BED sort rewritten as a Guttman-Rösler transform (pack+plain-sort instead of arrayref-Schwartzian) — ~1.57x faster + ~1.57x lower peak RSS on the BED-heavy path, and it makes the BED output DETERMINISTIC by adding `end` as the final sort key (the original emitted equal-start fragments in randomized hash order). bedGraph output byte-identical; BED is the same line-set, now reproducible. + README_optimized.md, OPTIMIZATION_PLAN.md, tests/ (54/54 pass). affects results: BED line *order* only (same set of lines; now stable run-to-run) — bedGraph unchanged

## 2026-11-20 · claude/bedgraphs2stdBedGraph_optimized_v2
- Added bedgraphs2stdBedGraph_optimized (streaming k-way merge + atomic write-temp-then-rename): byte-identical .bedGraph output to the original with far lower peak RAM on large merges, plus a guard for unsorted input. Additive (original untouched) + byte-diff harness. Re-validated: 25/25 pass (incl. 1.2M-record input). affects results: no

## 2026-11-18 · claude/add-trim_and_filter_optimized
- Added trim_and_filter_PE_optimized.pl + trim_and_filter_SE_optimized.pl (byte-identical drop-in additions; originals unchanged): quality-average hot loop rewritten with unpack("%32C*") instead of per-base substr/ord, ~8.5x faster; + README_optimized.md and tests/ (byte-diff + throughput incl. cutadapt reference); affects results: no
