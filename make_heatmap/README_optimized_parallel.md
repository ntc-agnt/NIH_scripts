## make_heatmap_optimized_parallel

A **multi-threaded, byte-identical** build of `make_heatmap`. It layers real thread
scaling on top of `make_heatmap_optimized` (same windowed feature lookup, cached table
pointers, and manual parsing) **without changing a single output bit** — the output is
byte-for-byte identical to the original `make_heatmap` at *every* thread count, not just
`-t 1`.

It is an **addition, not a replacement**: `make_heatmap.cpp` and
`make_heatmap_optimized.cpp` are untouched.

```bash
# Linux and macOS, one source (std::unordered_map; needs C++11 + pthreads):
g++ -O3 -std=c++11 -pthread -o make_heatmap_optimized_parallel make_heatmap_optimized_parallel.cpp

# self-contained for /n/data1 deployment (no libstdc++ module coupling):
g++ -O3 -std=c++11 -pthread -static-libstdc++ -static-libgcc \
    -o make_heatmap_optimized_parallel make_heatmap_optimized_parallel.cpp

# single-threaded build (no pthreads; behaves exactly like make_heatmap_optimized):
g++ -O3 -std=c++11 -DSINGLE -o make_heatmap_optimized_parallel make_heatmap_optimized_parallel.cpp
```

Build flags, toolchain and FP guidance (no `-ffast-math`, avoid `-march=native`/FMA)
are identical to `README_optimized.md` — see there.

### Why this exists: the original's `-t N` doesn't scale (it can regress)

The original `make_heatmap` threads under **two global mutexes** — one around reading
each input line, one around every per-bin accumulation. Only the *match* runs in
parallel. That was a reasonable trade when the match was the original's **O(genes)
linear scan per hit** (expensive, so worth parallelizing). But `make_heatmap_optimized`
replaced that scan with an `O(log G + window)` lookup, which makes the match cheap — so
now the two locks dominate, and adding threads mostly adds lock contention. On a real
lab workload (12M hits × 20k TSS) `make_heatmap_optimized` at `-t 6` is **slower** than
the same binary at `-t 1`.

This matters because the pipelines hard-code `-t 6`: `pipeline2/pipeline.sh` and
`pauseIndex` both invoke make_heatmap with `-t 6` on a `-c 6` sbatch reservation. This
build makes that `-t 6` actually pay off, as a drop-in, with no change to the callers.

### How it stays byte-identical while scaling

Three stages are parallelized without touching any value or any row's position:

1. **Parse + match** run in worker threads over a block of input lines.
2. **Accumulation stays in strict input-file order.** Workers emit each intersection as
   a `(cell, delta)` into per-thread, file-order lists; a short serial pass replays those
   deltas into the bin totals in exactly the order `-t 1` would, so every cell's
   floating-point summation order is identical to the serial run. (This is the opposite
   of a per-thread-accumulate-then-reduce design, which would reorder the sums and only
   be *functionally* identical — that variant was prototyped and rejected: same speed,
   double the RAM, weaker guarantee.)
3. **Output formatting** is parallelized — each thread formats a contiguous range of the
   (sorted) output rows into its own buffer, concatenated in row order.

The read of raw lines and the ordered delta-apply are the serial floor; everything that
can move off the critical path without reordering a sum has. `strand_map` lookups were
made non-inserting so parallel match has no shared-state write at all.

### The identity contract

| invocation | guarantee |
|---|---|
| **`-t 1`** | byte-for-byte identical to `make_heatmap` (serial path). |
| **`-t N`, any N** | byte-for-byte identical to `make_heatmap` at `-t 1` — for *all* data (float bedGraph coverage included), because the per-cell summation order is preserved. |

This is strictly stronger than `make_heatmap_optimized`, whose `-t N` float sums (like
the original's) vary in the last bits with thread interleaving. Here `-t N` is both
faster *and* bit-reproducible.

### Performance (one call, 12M hits × 20k TSS; 4-core box)

Measured against the original and `make_heatmap_optimized` on the pauseIndex-style
workload (`-b v`, stranded `-p/-m`):

| binary @ setting | wall | vs the original at `-t 6` |
|---|--:|--:|
| `make_heatmap` (original) `-t 1` | 145.5 s | — |
| `make_heatmap` (original) `-t 6` *(what the pipelines run today)* | 43.2 s | 1.0× |
| `make_heatmap_optimized` `-t 1` | 4.6 s | 9.4× |
| `make_heatmap_optimized` `-t 6` | 6.8 s | 6.4× *(lock regression: slower than its own `-t 1`)* |
| **`make_heatmap_optimized_parallel` `-t 6`** | **3.3 s** | **13×** |

The bulk of the win is the windowed lookup already in `make_heatmap_optimized` (the
9.4× at `-t 1`). This build adds thread scaling on top: as a **drop-in for the
pipelines' hard-coded `-t 6`** it is ~2× faster than `make_heatmap_optimized` at `-t 6`
and ~1.4× faster than `make_heatmap_optimized` at `-t 1`, with no caller changes.

The workload is memory-bandwidth-bound once the lock is gone, so scaling is sub-linear
(~1.4–2× over the single-thread `_optimized`, not N×) and flattens past a handful of
cores — do not expect large gains beyond ~4–6 threads. Peak memory is modest: the
per-thread delta lists are sized to the number of intersections in a block, not to the
whole matrix (543 MB vs 404 MB single-thread on the run above).

### When to use which build

- **`make_heatmap_optimized`** — the default serial build. 2–13× faster than the
  original, byte-identical at `-t 1`. Use it (or this build at `-t 1`) when you are not
  giving it multiple cores.
- **`make_heatmap_optimized_parallel`** — when the job already reserves cores and calls
  `-t N > 1` (the pipelines' case). Byte-identical at every thread count, so it is a safe
  drop-in for the existing `-t 6` call sites and also fixes their lock regression.

### Tests

`tests/run_parallel.sh` builds the reference (`make_heatmap.cpp`) and this candidate and
asserts the candidate at `-t 1/2/4` is `cmp`-identical to the reference across a
path-complete option matrix, including the multi-index (comma-separated hit files)
output path:

```bash
bash tests/run_parallel.sh                      # -t {1,2,4} all byte-identical to ref
MH_THREADS="1 2 3 4 6" bash tests/run_parallel.sh
```

Rebuild both binaries on O2 and re-run there before trusting a swapped `/n/data1`
binary — it rebuilds the reference from source, so it validates byte-identity under
whatever compiler/libstdc++ the lab uses.
