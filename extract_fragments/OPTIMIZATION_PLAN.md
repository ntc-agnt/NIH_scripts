# extract_fragments.pl — optimization plan

Status: **proposed (not yet implemented)**. Author: Claude Code session, 2026-11-18.
This is a design/intent document only; no behavior has changed yet.

## 1. Motivation

`extract_fragments.pl` reads a paired-end bowtie/SAM alignment file, computes each
fragment's `(chr, start, end)` + length, and for fragments within `[min,max]`
emits a BED file (one line per fragment) and/or a binned-fragment-center
bedGraph. Two issues, both measured on synthetic data:

1. **Cost.** The final BED sort uses the arrayref-Schwartzian idiom
   (`sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } map { [/regex/] }
   keys(%entries)`), which builds a 3-element arrayref per distinct fragment (large
   transient scratch that coexists with the fragment hash) and runs a Perl
   comparator on every comparison. Peak memory is dominated by `%entries`
   (O(distinct fragments) — high cardinality, since fragments are keyed on *both*
   endpoints) plus that sort scratch.
2. **Non-determinism (latent bug).** The BED sort tiebreaks only on
   `(chr, start)`, not on `end`. Fragments that share a start are therefore
   emitted in Perl hash-key order, which is randomized per process. Verified: the
   same input produces the **same set of BED lines in a different order** on
   repeated runs. The bedGraph path is unaffected (bins are unique integers,
   fully sorted → deterministic).

### Measured baseline (synthetic, worst case: all fragments distinct + in range)

| scale | wall | peak RSS |
|---|---:|---:|
| 2 M pairs, `-o a` | 37 s | 1.58 GB |
| 4 M pairs, `-o a` | 82 s | 3.12 GB |

→ ~19 s and ~0.8 GB **per million distinct in-range fragments**, single-threaded.
Scales linearly; a deep complex library (tens of millions of distinct fragments)
reaches tens of GB on the BED path. (`-o g` *without* `-u` never builds
`%entries` and already runs in well under 1 GB — not a target.)

### Prototype result (GRT lean sort, measured at 4 M pairs, `-o a`)

| variant | wall | peak RSS |
|---|---:|---:|
| original | 82.5 s | 3121 MB |
| GRT-optimized | 53.9 s | 1800 MB |
| | **~35 % faster** | **~42 % lower peak** |

Scaled estimate for a full dataset:

| pairs | time 19→~12 s/M | peak 0.8→~0.46 GB/M distinct |
|---|---:|---:|
| 20 M | ~6 → ~4 min | ~16 → ~9 GB |
| 50 M | ~16 → ~11 min | ~40 → ~23 GB |
| 100 M | ~32 → ~21 min | ~80 → ~46 GB |

## 2. Goals / non-goals

**Goals**
- Ship a faster, lighter variant as an **addition** (`extract_fragments_optimized.pl`),
  leaving `extract_fragments.pl` untouched (same pattern as the other NIH_scripts
  optimized additions).
- **bedGraph output: byte-identical** to the original.
- **BED output: same set of lines as the original, and deterministic**
  (reproducible byte-for-byte run-to-run) by adding `end` as the final sort key.
- Same CLI, options, inputs, filenames, track headers, dedup (`-u`) semantics,
  chr-renaming (`-D`/`-M`), `-l`/`-r` length handling, and SAM (`-a`) vs bowtie
  parsing.

**Non-goals (this round)**
- No change to the parser, the dedup/accounting logic, the fragment-length or
  bin-center math, or the leftmost-mate assumption in the bowtie path.
- No external-sort rewrite (see §6, listed as a separate future option).
- No multithreading.

## 3. The determinism decision (needs sign-off before coding)

Making the BED output deterministic changes the **order** of fragments that share
`(chr, start)` relative to any given original run (the original's order was
random). The *set* of lines is identical. This is a (benign) behavior change:
reproducible output is strictly better, and no sane downstream BED consumer
depends on the order of equal-start fragments — but it is a change, so it must be
approved. Options if full determinism is not wanted:
- **(b)** Optimize **only** the bedGraph path (strictly byte-identical), leave BED
  untouched. Smaller win, zero behavior change.
- **(c)** Do nothing (gains are smaller than the trim/bowtie cases).

Recommendation: **(a)** — optimize + make BED deterministic. The determinism fix
is the main reason this is worth doing at all.

## 4. Design

New file `extract_fragments/extract_fragments_optimized.pl`, copied from the
original with exactly these changes:

1. **BED sort → Guttman-Rösler transform.** Replace
   `sort{cmp||<=>} map{[/regex/]} keys(%entries)` with: pack each key as
   `chr\t<zero-padded start>\t<zero-padded end>` (10-digit pad, valid for any
   position < 10^10), sort with a plain C-level `sort @packed` (no comparator,
   no arrayref-per-key), then recover `chr`/`start`/`end` from the packed string
   and look up the count. This gives a fully-ordered `(chr, start, end)` sort
   (deterministic) with far less scratch.
2. **bedGraph path unchanged** (already lean + deterministic).
3. Everything else (option parsing, parse, dedup, bin math, stats, headers,
   filenames) copied verbatim.

Sketch (illustrative, not final):
```perl
my @packed = map { my($c,$s,$e)=/^([^\t]+)\t(\d+)\t(\d+)/;
                   sprintf("%s\t%010d\t%010d",$c,$s,$e) } keys(%entries);
@packed = sort @packed;
foreach (@packed) {
    my ($c,$sp,$ep)=/^([^\t]+)\t(\d+)\t(\d+)/; my $s=$sp+0; my $e=$ep+0;
    my $cnt=$entries{"$c\t$s\t$e"}; my $s0=$s-1;
    print BEDFILE "$c\t$s0\t$e\n" for 1..$cnt;
}
```

## 5. Validation (test harness, to live in `extract_fragments/tests/`)

A synthetic bowtie/SAM PE generator + `run_benchmark.sh` that:
- **bedGraph**: strict byte-diff (optimized == original) across modes.
- **BED**: assert (i) optimized line-set == original line-set (`sort | cmp`), and
  (ii) optimized output is byte-stable across two runs (determinism).
- **Mode matrix**: `-o b|g|a`; varied `-min/-max/-b`; `-u` dedup (0 / 1 / N);
  `-l` length list and (if reachable) `-r`; `-D`/`-M`; SAM (`-a`) vs bowtie input;
  same-start/different-end tie fixtures; empty-output and all-filtered cases.
- **Throughput + peak RSS** comparison at a large size (original vs optimized),
  reusing the `/proc` VmRSS sampler from the other NIH_scripts harnesses.
- Self-contained (perl + coreutils); `DBI` only needed for `-r` (harness uses a
  lazy-require copy so it runs without DBI installed).

## 6. Future option (out of scope here): external-sort for bounded memory

The in-language GRT refactor trims ~⅓ off time/memory but does **not** remove the
O(distinct-fragments) memory floor (tens of GB at high depth). Streaming fragments
to `LC_ALL=C sort -k1,1 -k2,2n -k3,3n` instead of the in-RAM `%entries` hash would
make memory **bounded** (a sort buffer of a few hundred MB regardless of depth →
sub-GB), is deterministic by construction, and can use multiple cores via
`sort --parallel`. Cost: a real rewrite of the accumulate/dedup/output logic and a
dependency on coreutils `sort`. Worth it only if the tens-of-GB ceiling is a real
operational problem; propose as a separate effort.

## 7. Risks & mitigations
- *Byte-identity vs determinism confusion* → harness encodes the exact contract
  (byte-identical bedGraph; set-equal + stable BED).
- *Zero-pad width* → 10 digits covers any genome; assert positions fit.
- *Dedup / `-u` interaction with the new writer* → covered by the `-u 0/1/N` matrix.
- *SAM vs bowtie parse parity* → both exercised in the matrix (parser is unchanged,
  so this is a guard, not a new risk).

## 8. Deliverables (when approved)
- `extract_fragments/extract_fragments_optimized.pl` (addition; original untouched)
- `extract_fragments/README_optimized.md`
- `extract_fragments/tests/` (generator + `run_benchmark.sh`)
- `.claude/ntc-changelog.md` entry (affects results: BED line *order* only, see §3)
- On its own `claude/add-extract_fragments_optimized` branch, like the others.

## 9. Open questions for the maintainer
1. Approve the BED determinism change (§3 option a), or restrict to bedGraph-only
   (option b)?
2. Is the tens-of-GB BED memory at full depth an operational pain worth the
   external-sort rewrite (§6), or is the ~40 % GRT reduction enough for now?
3. Should the optimized variant eventually *replace* the original (once trusted),
   or remain a permanent sibling?
