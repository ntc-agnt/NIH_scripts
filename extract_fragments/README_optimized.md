## extract_fragments_optimized

A faster, lighter, and **deterministic** alternative to `extract_fragments.pl`
(same directory). It is an **addition, not a replacement**: the original script
is unchanged, and this one is a drop-in with the **same command line, options,
inputs, output filenames, track headers, and dedup/chr-renaming/length-list
semantics**.

```
extract_fragments_optimized.pl [options] [input file] [output prefix]
```

See the original `README.md` for the full option list (`-a -o -min -max -b -l
-r -D -M -u`); they are identical here.

### What changed

Only the **BED-writing step**. The original sorts distinct fragments with the
arrayref-Schwartzian idiom

```perl
sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } map { [/regex/] } keys(%entries)
```

which builds a 3-element array reference for **every** distinct fragment (a large
transient that coexists with the fragment hash) and runs a Perl comparator on
every comparison. This version uses a **Guttman-Rösler transform**: it packs each
key as `chr\t<zero-padded start>\t<zero-padded end>`, sorts the packed strings
with a plain C-level `sort` (no comparator, no arrayref-per-key), then recovers
the fields. Positions are zero-padded to 10 digits (ample for any genome — a
wider coordinate would break the fixed-width ordering, so the script dies with a
clear message if one ever appears).

The bedGraph path is **unchanged** (it was already lean and deterministic — bins
are unique integers, already fully sorted).

DBI is now loaded lazily (only on the `-r` UCSC-fetch path) so the script runs on
hosts without DBI installed unless `-r` is used. Output is unaffected.

### Two output guarantees

1. **bedGraph output is byte-identical** to the original, across `-o g`/`-o a`,
   all bin sizes, `-l` length lists, `-D`/`-M` chr handling, and bowtie vs SAM
   input.
2. **BED output is the same *set* of lines as the original, and is now
   deterministic.** The original tiebreaks only on `(chr, start)`, so fragments
   that share a start were emitted in randomized Perl hash-key order — the **same
   lines in a different order on every run** (verified). Adding `end` as the
   final sort key makes the BED output reproducible byte-for-byte run to run. The
   *set* of lines (and every count) is identical to the original; only the order
   of equal-start fragments changes, and is now stable.

> **The one deliberate behavior change.** For fragments that share `(chr,
> start)`, the BED line *order* relative to any single original run changes
> (because the original's order was random). No BED consumer depends on the order
> of equal-start lines, and reproducible output is strictly better — but it is a
> change, so it is called out here. If you need bit-for-bit reproduction of a
> *specific historical* `extract_fragments.pl` BED file, that file was not
> reproducible to begin with; the line set is preserved.

### Measured (synthetic, 400k distinct in-range fragments, `-o a`)

| implementation            | peak RSS | wall time |
|---------------------------|---------:|----------:|
| `extract_fragments.pl` (original) | 322 MB | 5.8 s |
| `extract_fragments_optimized.pl`  | 206 MB | 3.7 s |

≈ **1.57× lower peak memory** and **1.57× faster** on the BED-heavy path, with a
byte-identical bedGraph and a set-equal, deterministic BED. The win grows with
the number of distinct in-range fragments (the memory floor is still
`O(distinct fragments)` — see "Scope" below).

### Scope / non-goals

- The parser, dedup/accounting (`-u`), fragment-length and bin-center math, and
  the leftmost-mate assumption are **untouched** — this is a sort/write change
  only.
- Peak memory is still `O(distinct fragments)` (the `%entries` hash), just with a
  smaller constant. A deep, complex library with tens of millions of distinct
  fragments still needs many GB on the BED path. Making memory **bounded** would
  require streaming fragments to `LC_ALL=C sort` instead of the in-RAM hash — a
  real rewrite, proposed separately in `OPTIMIZATION_PLAN.md §6`, not done here.
- No multithreading.

### Tests

`tests/` is a self-contained harness (perl + coreutils; no cluster, no real
data, no DBI):

```
bash tests/run_benchmark.sh
```

It generates synthetic bowtie and SAM paired-end inputs (`gen_fragments.pl`) and,
across a mode matrix (`-o b|g|a`; varied `-min/-max/-b`; `-u 0/1/N` dedup;
`-l` length list; `-D`/`-M` chr handling; bowtie vs SAM; a same-start/different-
end tie fixture; empty/all-filtered and single-chromosome cases):

- `cmp`s the bedGraph output byte-for-byte (optimized == original),
- checks the BED **line set** equals the original's (`sort | cmp`) and that the
  optimized BED is **byte-stable across two runs** (determinism),
- reports a peak-RSS + wall-time comparison on a large input.

Last run: **54/54 pass**. (`-r` is not exercised — it needs network access to
UCSC; the parse path it feeds is covered by the `-l` length-list cases.)
