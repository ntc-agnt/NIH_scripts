# make_heatmap — numerical reproducibility notes

_For the NTC lab. Explains what "identical output" means for make_heatmap, where
floating-point summation order does and doesn't matter, and which of the three
builds to use when bit-for-bit reproducibility matters._

## Two different "orders" — don't conflate them

1. **Row order** — which feature is on which output line.
2. **Summation order** — the sequence in which individual signal values are
   *added together* to compute one matrix cell.

These are independent. Everything below is about **summation order** (#2), which
affects only the last-bit rounding of a cell's *value* — it never moves a feature
to a different row.

**Row order is fixed and safe.** make_heatmap emits rows **sorted lexicographically
by the feature-ID column** (an internal ordered map), regardless of input order.
This is identical in both builds (original and `_optimized`).
Note two consequences, true of every build:
- make_heatmap does **not** preserve *input* order — "input row N ↔ output row N"
  holds only if the feature list is already ID-sorted. A consumer that joins by the
  ID column is always safe; one that joins by line position must feed an ID-sorted
  list.
- The sort is lexicographic (string), so `gene10` sorts before `gene2`.

## Why summation order matters at all

Floating-point addition is **not associative**: `(a+b)+c` can differ from
`a+(b+c)` in the last bit, because each `+` rounds to the nearest representable
`double`. So summing the same values in a different order can change a cell by a
few units in the last place (ULPs). This is what breaks *byte*-identity even when
the result is numerically the same.

## How big is the difference? (measured)

Summing the same set of values in ascending / descending / shuffled / reversed
order and taking the spread:

| Values in a cell | Relative spread across orderings | Verdict |
|---|---|---|
| **Non-negative coverage** (100k values, ~U(0,100)) | **2.5e-14** | last ~2 ULPs — numerically identical |
| **Integer counts** (BED `value=1`, integer coverage) | **exactly 0** | order-*independent*, bit-identical in any order |
| **Mixed-sign with cancellation** (true sum ≪ the terms being added) | **4.6e-5**, rising as the true sum → 0 | **significant divergence** |
| **Extreme dynamic range** (1e16 + many 1.0s in one cell) | 1e-11 relative, 1e5 **absolute** | small terms silently absorbed |

Rule of thumb: for naive sequential summation the orderings can differ by about
`n · ε · Σ|xᵢ|`, where `n` is the number of values in the cell and
`ε ≈ 1.1e-16` (double epsilon). For non-negative coverage this stays at the
`~1e-14` level; it only blows up when large terms cancel (so `Σ|xᵢ|` ≫ the true
result) or when the dynamic range approaches `1/ε ≈ 1e16`.

## What this means for real make_heatmap data

- **Raw / normalized coverage, counts:** all values are non-negative and modest in
  range, and the number per cell is in the hundreds–thousands. Summation order
  changes a cell by ~1e-13–1e-14 relative — effectively identical. **Integer-valued
  matrices (BED, integer coverage) are exactly identical in any order.**
- **The one case that genuinely diverges: signed / differential tracks.** If a
  bedGraph fed to make_heatmap holds *signed* values whose per-bin total is a small
  residual of large opposing terms — e.g. a WT−KO subtraction track or a log2
  fold-change track — then summation order (and therefore any parallel or
  multi-threaded run) can move results at the 1e-4 level or worse. Treat such tracks
  specially (see builds, below).
- **Downstream amplification.** Even the benign 1e-14 can be magnified when a later
  step divides or subtracts near-zero quantities — e.g. a pauseIndex ratio with a
  near-silent gene body in the denominator. That quantity is numerically unstable
  regardless of summation order (it's noise-dominated), so this is *exposed*
  pre-existing instability, not new error. When you set an acceptance tolerance,
  make it **relative** (≈1e-4), not a fixed number of decimal places — a fixed "3
  dp" is meaningless across cells that range from 0.0003 to 300000 — and validate it
  against the *downstream* analysis, not just the matrix.

## Does the result depend on thread count?

**In the original and the `_optimized` build: yes, for float values.**
Those builds accumulate each cell under a single global lock, in whatever order the
threads happen to acquire it:
- `-t 1` is **deterministic** (values are summed in input-file order).
- `-t 2`, `-t 6`, … each interleave differently, so a float cell can differ from
  the `-t 1` result **and run-to-run** — at the ~1e-14 level for coverage (Case A).
- **Integer/count matrices are identical at every thread count** (exact, Case B).

So on the same float data, those builds give trivially different numbers at `-t 1`,
`-t 2`, `-t 6`. For bit-for-bit reproducibility from them, use `-t 1`.

**In the `_optimized_parallel` build: no — it is byte-identical at every thread count.**
It preserves each cell's input-file summation order even when threaded (workers emit
file-ordered deltas that a serial pass replays in order), so `-t 1 == -t 2 == -t 6`
bit-for-bit, for float data too. (It is also the only threaded build whose `-t N` is
actually *faster* than its own `-t 1`; see the note below on the lock.)

## The builds and their guarantees

| Build | Speed | Summation | Reproducibility |
|---|---|---|---|
| `make_heatmap` (original) | baseline | input-file order at `-t 1` | bit-identical only at `-t 1`; `-t N` varies in last bits (float) |
| `make_heatmap_optimized` | ~2–13× (scales with gene count) | **byte-identical to the original at `-t 1`** | same as original (`-t 1` deterministic; `-t N` varies) |
| `make_heatmap_optimized_parallel` | `_optimized` + thread scaling (~1.4–2× over its `-t 1` on ≤6 cores) | **byte-identical to the original at *every* `-t`** | `-t 1 == -t N` bit-for-bit, float data included |

- Use **`_optimized`** — the faster, byte-identical serial drop-in and the safe default:
  for standard non-negative coverage/count data it reproduces historical `make_heatmap`
  output exactly, and it's 2–13× faster.
- Use **`_optimized_parallel`** when the job reserves cores and calls `-t N > 1` (the
  `pipeline2` / `pauseIndex` case, both hard-coded to `-t 6`): it is byte-identical at
  every thread count — so a safe drop-in even for signed/differential tracks — and it is
  the only threaded build whose `-t N` beats its own `-t 1`.
- For any **signed/differential** input track with `_optimized` (not `_parallel`), pin to
  `-t 1` for bit-reproducibility and validate numerically.

> **The lock, and why threading the lock-based builds backfires.** The original and
> `_optimized` accumulate every bin under one global mutex, so only the *match* runs in
> parallel. That paid off when the match was the original's O(genes) linear scan, but
> `_optimized` made the match cheap (windowed lookup), so the lock now dominates: on a
> 12M-hit × 20k-TSS run, `_optimized -t 6` (6.8 s) is *slower* than `_optimized -t 1`
> (4.6 s). `_optimized_parallel` removes the lock (file-ordered delta replay, parallel
> output) and brings that `-t 6` to 3.3 s — ~2× over `_optimized -t 6`, ~13× over the
> original `-t 6` (43 s) — while staying byte-identical.
>
> A compensated-summation parallel variant (`_optimized_summation`) and a
> per-thread-accumulate variant were both prototyped and dropped: they tied the
> byte-identical build on speed (make_heatmap is memory-bandwidth-bound, not
> accumulation-bound) while giving up byte-identity and, for the per-thread one, doubling
> RAM. The byte-identical `_optimized_parallel` dominates them.

## Verifying on your data

`tests/run_harness.sh` builds `_optimized` and the reference from source and byte-diffs
them across the full option matrix; `tests/run_parallel.sh` does the same for
`_optimized_parallel` at `-t 1/2/4` (asserting byte-identity at every thread count);
`tests/benchmark.sh` reports absolute wall times across depth / gene-count / binning /
value-type. Re-run them on O2 with your toolchain before trusting a swapped `/n/data1`
binary.
