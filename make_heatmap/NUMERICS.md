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

**Yes, in the original and the byte-identical `_optimized` build — for float values.**
Those builds accumulate each cell under a single global lock, in whatever order the
threads happen to acquire it:
- `-t 1` is **deterministic** (values are summed in input-file order).
- `-t 2`, `-t 6`, … each interleave differently, so a float cell can differ from
  the `-t 1` result **and run-to-run** — at the ~1e-14 level for coverage (Case A).
- **Integer/count matrices are identical at every thread count** (exact, Case B).

So yes: on the same float data, `-t 1`, `-t 2`, `-t 6` can give trivially different
numbers. If you need bit-for-bit reproducibility from those builds, use `-t 1`.

## The builds and their guarantees

| Build | Speed | Summation | Reproducibility |
|---|---|---|---|
| `make_heatmap` (original) | baseline | input-file order at `-t 1` | bit-identical only at `-t 1`; `-t N` varies in last bits (float) |
| `make_heatmap_optimized` | ~2–13× (scales with gene count) | **byte-identical to the original at `-t 1`** | same as original (`-t 1` deterministic; `-t N` varies) |

- Use **`_optimized`** — it is the faster, byte-identical drop-in and the safe default:
  for standard non-negative coverage/count data it reproduces historical
  `make_heatmap` output exactly, and it's 2–13× faster.
- For any **signed/differential** input track (the one case summation order can matter,
  above), pin to `-t 1` for bit-reproducibility and validate numerically.

> A compensated-summation parallel variant (`_optimized_summation`) was prototyped and
> then dropped: it was ~20% *slower* single-threaded and only ~1.3–1.5× on ≤6 cores
> (make_heatmap is not accumulation-bound), and it altered the output bits. The
> remaining speedup lever is a **byte-identical** parallelization of the match + output
> stages; see the benchmark/plan work.

## Verifying on your data

`tests/run_harness.sh` builds `_optimized` and the reference from source and byte-diffs
them across the full option matrix; `tests/benchmark.sh` reports absolute wall times
across depth / gene-count / binning / value-type. Re-run both on O2 with your toolchain
before trusting a swapped `/n/data1` binary.
