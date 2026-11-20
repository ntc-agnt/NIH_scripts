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

## 2026-11-20 · claude/bedgraphs2stdBedGraph_optimized_v2
- Added bedgraphs2stdBedGraph_optimized (streaming k-way merge + atomic write-temp-then-rename): byte-identical .bedGraph output to the original with far lower peak RAM on large merges, plus a guard for unsorted input. Additive (original untouched) + byte-diff harness. Re-validated: 25/25 pass (incl. 1.2M-record input). affects results: no

## 2026-11-18 · claude/add-trim_and_filter_optimized
- Added trim_and_filter_PE_optimized.pl + trim_and_filter_SE_optimized.pl (byte-identical drop-in additions; originals unchanged): quality-average hot loop rewritten with unpack("%32C*") instead of per-base substr/ord, ~8.5x faster; + README_optimized.md and tests/ (byte-diff + throughput incl. cutadapt reference); affects results: no
