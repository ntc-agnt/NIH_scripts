## bedgraphs2stdBedGraph_streaming

A memory-lean, streaming alternative to `bedgraphs2stdBedGraph` (same directory).
It is an **addition, not a replacement**: the original script is unchanged, and
this one is a drop-in with the **same command line and byte-identical output**
on well-formed input.

```
bedgraphs2stdBedGraph_streaming [output prefix]
```

Like the original, it finds all `*.bedGraph` (standard 0-based) and `*.bedgraph`
(1-based start) files in the working directory, merges/sums them over identical
intervals, and writes `<prefix>.bedGraph` plus a `<prefix>.log`.

### Why

The original reads **every** input fully into one in-memory hash
(`%table{chr}{end} = [start, sum]`) and only sorts and writes at the end, so its
peak memory is the size of the **union of all intervals** — for 1-nt
whole-genome PRO-seq/Start-seq bedGraphs that is tens of millions of entries and
multiple GB, all resident at once, and it grows as you merge more coverage.

This version exploits the fact that the inputs are already sorted by
`(chr in cmp order, start ascending)` — the order `bowtie2stdBedGraph` emits and
`normalize_bedGraph` preserves — and does a **k-way merge** across the files,
streaming output as it goes. It holds only one "front" record per input file, so
peak memory is `O(number of files)` instead of `O(union of intervals)`.

### Measured (1.2M-record union, 8 files, ~137 MB input)

| implementation            | peak RSS | wall time |
|---------------------------|---------:|----------:|
| `bedgraphs2stdBedGraph` (original, hash) | 638 MB | 21.8 s |
| `bedgraphs2stdBedGraph_streaming`        | 5.8 MB | 16.1 s |

≈ **110× lower peak memory** and modestly faster, with byte-identical output.
(The speedup comes from avoiding the giant hash build + final sort; the hot loop
is fully inlined with no per-record subroutine calls.)

### Guarantees

- **Byte-identical `.bedGraph`** to the original on well-formed input, including
  the float case: values are summed in the same file order and printed with
  Perl's same default stringification, so normalized (floating-point) bedGraphs
  reproduce bit-for-bit. Verified by the harness in `tests/` across integer,
  float, mixed 0/1-based, single-file (1→0 conversion), 50-file, sparse-union,
  track/blank-line, multi-base-width, and chromosome-ordering scenarios.
- The `.log` matches the original except for its `Initiated <localtime>` line
  (a wall-clock timestamp, never reproducible run to run).

### Two deliberate differences from the original

1. **Requires sorted input.** Each file must be sorted by `(chr cmp, start)`.
   Pipeline-produced inputs already are. If a file is not, the streaming merge
   would mis-merge, so instead it **dies with a clear message** rather than
   produce a silently wrong result. The hash-based original tolerates arbitrary
   input order; if you have unsorted inputs, use the original (or sort first).
2. **Atomic output.** Output is written to a temporary file and `rename`d onto
   `<prefix>.bedGraph` only after a clean finish. A mid-run error (e.g. a
   misaligned interval) therefore **never leaves a partial `.bedGraph`** and
   never clobbers a good file from a previous run — matching the original, which
   only creates its output after all input has been read successfully.

For genuinely misaligned input (an interval whose `end` recurs with a different
`start`), both scripts abort with a non-zero exit, exactly as the original does.

### Tests

`tests/` contains a self-contained byte-diff harness and a synthetic bedGraph
generator (no cluster or real data needed):

```
bash tests/run_harness.sh
```

It generates synthetic inputs, runs both the original and the streaming script,
`cmp`s their `.bedGraph` output byte-for-byte across all scenarios, checks the
misaligned/ unsorted/ atomic error paths, and prints a peak-RSS + wall-time
comparison on a large input.
