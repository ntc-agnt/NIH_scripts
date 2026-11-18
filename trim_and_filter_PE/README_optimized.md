## trim_and_filter_PE_optimized.pl

A performance-optimized, **byte-identical** drop-in for `trim_and_filter_PE.pl`
(same directory). It is an **addition, not a replacement**: the original is
unchanged, and this variant takes the **same command line, options, and inputs**
and produces **byte-for-byte identical** trimmed FASTQ and `FilterStats.txt` on
every mode (both `sanger` and `illumina` scales, plain or `.gz`).

```
perl trim_and_filter_PE_optimized.pl -1 <m1.fq> -2 <m2.fq> \
     -a <m1_pos1> -b <m1_posN> -c <m2_pos1> -d <m2_posN> \
     -m <minAvgBQS> -q <illumina|sanger> -o <outRoot>
```

### What differs (and why it's faster)

Only the per-read quality-averaging loop changes. The original summed the
trimmed quality window one base at a time (`substr`+`ord` in a Perl loop, with
the `-q` scale branch re-tested on every base). This version computes the same
sum in a single C-level operation:

```perl
$off = ($qscale eq "ilmn") ? 64 : 33;            # resolved once, not per base
...
$avg = (unpack("%32C*", $trimQ) - $off*length($trimQ)) / length($trimQ);
```

`unpack("%32C*", …)` returns the sum of the byte values of the trimmed quality
string; subtracting `offset × length` yields exactly the original
`sum(ord(b) − offset)`. Same integer sum, same division → identical averages,
pass/fail decisions, output, and stats. Everything else (parsing, trimming,
I/O, stats) is unchanged.

### Measured (synthetic, 4M read pairs, readlen 50, single thread; see `tests/`)

| tool | wall | reads kept |
|---|---:|---|
| `trim_and_filter_PE.pl` (original) | 50.2 s | 1,675,770 |
| **`trim_and_filter_PE_optimized.pl`** | **5.8 s** | 1,675,770 (identical) |
| cutadapt `-j1` (best-matching flags) | 13.2 s | 1,414,376 |
| cutadapt `-j4` | 4.3 s | — |

**~8.7× faster than the original**, byte-identical. For reference, the
single-threaded optimized script also beats cutadapt at 1–2 threads on this
fixed-window-trim + quality-filter task; cutadapt only overtakes it using all
cores. (cutadapt's nearest filter, `--max-average-error-rate`, is *not* the same
metric as this tool's mean-Phred filter — hence the different kept-read count —
so that row is a throughput reference only, not a drop-in substitute.)

### Notes
- Same dependency footprint as the original (core Perl; `zcat`/`bzcat` for
  compressed input).
- Pre-existing behavior is preserved exactly, including the division-by-zero if a
  read is shorter than `-a` (empty trimmed window) — not introduced here.

### Tests
`tests/run_benchmark.sh` (self-contained; needs only perl + coreutils, and
cutadapt if present) generates synthetic FASTQ, verifies byte-identity of both
the PE and SE optimized scripts across all modes, and prints the throughput
comparison. See `tests/gen_fastq.pl` for the generator.
