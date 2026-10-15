## bowtie2stdBedGraph_optimized.pl

A performance-optimized, **byte-identical** drop-in for `bowtie2stdBedGraph.pl`
(same directory). It is an **addition, not a replacement**: the original is
unchanged, and this variant takes the **same command line, options, and inputs**
and produces **byte-for-byte identical output** on every mode.

```
perl bowtie2stdBedGraph_optimized.pl [options] [input file] [output prefix]
```

(See `bowtie2stdBedGraph_readme.txt` / the original `README.md` for the full
option list — `-o nbmv`, `-b`, `-s`, `-l`/`-r`, `-t`, `-u`, `-i`, `-p`, `-x`,
`-D`, `-M`. They behave identically here.)

### What differs (and why it's faster / leaner)

Only two things change, neither of which affects output:

1. **No ithreads.** The original spawns a thread per output table and hands it
   the count hash. Perl ithreads **deep-copy** each hash into the new thread, so
   the process holds ~2x the table memory (more under `-o nbmv`) while sorting.
   This version sorts and writes the tables **sequentially in the main process**,
   so only one copy of the data is ever resident.
2. **Lean sort.** The original sorts keys with a Perl comparator over arrayrefs
   built by a regex (`sort { cmp || <=> } map { [/regex/] } keys`), which builds
   an arrayref per key (large transient scratch) and runs a Perl callback on
   every comparison. This version uses a **Guttman-Rösler transform**: each key
   is packed as `chr\t<zero-padded position>` and sorted with a plain, C-level
   string sort (no callback), reproducing `(chr cmp, pos <=>)` exactly with far
   less scratch and less CPU.

The parser and the paired-end accounting are left exactly as in the original, so
byte-identity is preserved in all cases (including `-u` dedup and `-i`).

### Measured (synthetic deep coverage; see `tests/`)

| output set | original | optimized | change |
|---|---|---|---|
| `-o n`  (12M reads / 4M sites) | 5450 MB, 115 s | **1528 MB, 81 s** | **0.28× RAM, 0.70× time** |
| `-o nbmv` (5M / 1.5M)          | 11887 MB, 153 s | **1717 MB, 78 s** | **0.14× RAM, 0.51× time** |
| `-o nbmv` (12M / 4M)           | **OOM-killed (>13.6 GB)** | 4.3 GB, completes | stays memory-bounded |

Roughly **3.6–7× lower peak memory and 1.3–2× faster**, byte-identical — and it
stays memory-bounded under `-o nbmv`, where the ithreads original can exhaust RAM.

### Notes
- The position pack uses a 10-digit zero pad, valid for any position `< 10^10`
  (human max ~2.5×10^8 = 9 digits).
- Same dependency footprint as the original (`Getopt::Std`, `threads` is no
  longer used, `DBI` only for the `-r` UCSC lookup).

### Tests

`tests/` has a self-contained byte-diff + benchmark harness and a synthetic
bowtie-format generator (no cluster or real data needed):

```
bash tests/run_benchmark.sh
```

It generates synthetic alignments, runs both scripts across all option modes,
`cmp`s the `.bedGraph` outputs byte-for-byte, and prints a peak-RSS + wall-time
comparison on a larger input.
