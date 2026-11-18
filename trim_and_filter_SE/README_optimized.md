## trim_and_filter_SE_optimized.pl

A performance-optimized, **byte-identical** drop-in for `trim_and_filter_SE.pl`
(same directory). It is an **addition, not a replacement**: the original is
unchanged, and this variant takes the **same command line, options, and inputs**
and produces **byte-for-byte identical** trimmed FASTQ and `FilterStats.txt` on
every mode (both `sanger` and `illumina` scales, plain or `.gz`).

```
perl trim_and_filter_SE_optimized.pl -i <in.fq> -a <pos1> -b <posN> \
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
pass/fail decisions, output, and stats. Everything else is unchanged.

### Performance

On the same synthetic benchmark as the PE tool (`../trim_and_filter_PE/tests/`),
the optimized script is several-fold faster than the original (the per-base
`substr`/`ord` loop was the dominant per-read cost), byte-identical. See that
harness for the PE and SE comparison and the cutadapt reference.

### Notes
- Same dependency footprint as the original (core Perl; `zcat` for `.gz` input).
- Pre-existing behavior is preserved exactly, including the division-by-zero if a
  read is shorter than `-a` (empty trimmed window) — not introduced here.
