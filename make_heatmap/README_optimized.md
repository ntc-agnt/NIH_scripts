## make_heatmap_optimized

A faster, **byte-identical** drop-in for `make_heatmap` (same directory). It is an
**addition, not a replacement**: `make_heatmap.cpp` is unchanged, and
`make_heatmap_optimized.cpp` takes the **same command line, options, input/output
formats, and (at `-t 1`) produces byte-for-byte identical output** across every
supported mode.

```bash
# Linux and macOS, one source (uses std::unordered_map; needs C++11):
g++ -O3 -std=c++11 -pthread -o make_heatmap_optimized make_heatmap_optimized.cpp

# single-threaded build (no pthreads):
g++ -O3 -std=c++11 -DSINGLE -o make_heatmap_optimized make_heatmap_optimized.cpp
```

Build notes:
- **`std::unordered_map`, not `<tr1/...>`.** Unlike `make_heatmap.cpp` (which uses
  the deprecated GNU `tr1` header and needs a separate `make_heatmap_OSX.cpp` plus
  `-stdlib=libstdc++` on macOS), this source uses the standard C++11 container, so
  **one file compiles on both Linux and macOS/clang** — no OSX variant, no
  `-stdlib=libstdc++`.
- **`-std=c++11`** is only needed for compilers that still default to C++03
  (modern g++ defaults to a later standard and builds it without the flag); it is
  harmless to pass always.
- **`-pthread`** (not just `-lpthread`) sets the threading flags for both the
  compile and link steps.
- **Do not add `-ffast-math`** — it reorders floating-point operations and would
  break bit-for-bit reproducibility (and the byte-identity guarantee). `-O3` and
  `-O2` are both fine.
- **Avoid `-march=native` / `-mfma`.** Besides cluster portability (a binary built
  on a newer node may fault on older ones), enabling FMA is the one thing that can
  shift coverage-mode (`-v c`) results by a last bit, by fusing `value × length`
  into the running total. Baseline x86-64 has no FMA instruction, so the default
  build is contraction-stable; if you ever do tune the arch, add
  `-ffp-contract=off` and prefer an explicit baseline (e.g. `-march=x86-64-v2`)
  over `native`.

### Recommended toolchain
- **Any GCC from 4.8.1 up is fine** (4.8.1 is the first complete C++11; pass
  `-std=c++11`). Newer is good, not a problem: a recent GCC — e.g. **`gcc/14` on
  O2** — gives a stronger optimizer and faster libstdc++ at no correctness cost.
  Tested here on GCC 13.3. The `-std` level is **semantic only** for this code
  (C++11/14/17/20 produce byte-identical output at the same speed), so stay at
  `-std=c++11`.
- **The caveat is the *runtime* C++ library, not the compiler version.** A binary
  built with `gcc/14` needs that module's `libstdc++.so` at launch. Under lmod that
  library lives in the module, so a plain dynamic build works **only if the
  environment that runs make_heatmap also has the same `gcc/NN` module loaded**;
  otherwise you get `GLIBCXX_3.4.xx not found`. Because this binary is a hardcoded
  `/n/data1/...` dependency invoked by many pipelines (which won't all load the
  same module), build it **self-contained**:
  ```bash
  g++ -O3 -std=c++11 -pthread -static-libstdc++ -static-libgcc \
      -o make_heatmap_optimized make_heatmap_optimized.cpp
  ```
  That removes the libstdc++/libgcc runtime dependency entirely (larger binary,
  no module coupling). Avoid full `-static` (it can break NSS lookups).
- Output is byte-identical *for the same flags* regardless of GCC version — IEEE
  `double` arithmetic is deterministic and FMA is left off (above). Confirm on the
  target toolchain by rebuilding both binaries and running `tests/run_harness.sh`
  there.

See the original `README.md` for the full option reference — nothing there
changes.

### What changed (three layers, all output-preserving)

1. **Windowed feature lookup (the big win).** The original matches each hit by a
   **linear scan over every gene on the chromosome** (`for i in 0..genes`,
   bounding-box test each). This version keeps, per chromosome, the genes **sorted
   by span-start** plus the widest span, and answers each hit with a **binary
   search over a bounded window** instead of a full scan — `O(log G + window)` in
   place of `O(G)`. The *set* of matched genes is identical to the full scan; only
   the visit order within a single hit changes, and that is irrelevant because each
   gene is an independent output-table cell (so no per-cell sum changes). Both the
   point path (`-v t/a/d`) and the coverage path (`-v c`) are indexed.
2. **Direct table pointers.** The hot loop updated the output table with a
   `std::map<std::string>` lookup **per intersection** (`table[gene_id]…`). Each
   gene now caches a direct pointer to its table entry (resolved once after
   parsing), removing that lookup from the inner loop.
3. **Manual line parsing.** Each input line was parsed by constructing a
   `std::istringstream` per line. This version extracts fields with
   `strtol`/`strtod` plus a hand-written decimal scanner that **replicates
   `operator>>`'s semantics exactly** — including where it sets failbit (no digits,
   incomplete exponent like `1e`), and its rejection of hex floats (`0x1F`→`0`) and
   `inf`/`nan`. Bad-line detection, `track`/`\r` skipping, 0-vs-1-based handling,
   and strand normalization are unchanged.

### Guarantees

- **Byte-identical to `make_heatmap` at `-t 1`** across the full option matrix —
  bin types `c`/`v`/`f`, hit types `G`/`g`/`b`/`e`/`c` (± stranded), all
  `-l`/`-a`/`-v`/`-d`/`-s` combinations, `-p`/`-m` split inputs, multi-file and
  multi-output, `--nostrand`, `--nohead`, and the `-DSINGLE` build. Verified by
  `tests/run_harness.sh` (108 checks: a path-complete matrix + an adversarial
  parser fuzz) which builds the **reference from `make_heatmap.cpp`** and the
  candidate from source and `cmp`s their output.
- **Threaded output (`-t N`) is unchanged in kind, not guaranteed bit-equal to
  `-t 1`.** This is pre-existing: the original accumulates bin totals in
  thread-interleaving order under one global lock, so `-t N` float sums already
  differ from `-t 1` in the last ULPs run-to-run. This version does not change that
  (the threading/locking is untouched); the stable contract is `-t 1`.

### Performance (synthetic, 3M hits, `-t 1`; speedup grows with gene count)

| gene list | `make_heatmap` | optimized | speedup |
|----------:|---------------:|----------:|--------:|
| 10 genes   |  1.22 s | 0.50 s |  2.4× |
| 200 genes  |  1.58 s | 0.76 s |  2.1× |
| 2,000 genes |  7.42 s | 1.17 s |  6.3× |
| 20,000 genes (GGA scale) | 70.6 s | 5.39 s | **13.1×** |

The linear scan is `O(hits × genes)`, so the win scales with the gene list; a
typical metagene run over a ~20k-gene GGA dominant-TSS list sees ~**13×**. Peak
memory is essentially unchanged (the index adds `O(genes)` small arrays per
chromosome). Worst case — one gene whose span covers most of a chromosome —
degrades gracefully toward the original scan cost, never worse.

### Scope / non-goals

- The **global table lock** that serializes `-t N` accumulation is **not** changed
  (changing it would alter threaded float sums, which are already non-reproducible;
  out of the byte-identical mandate). Thread *scaling* is therefore still limited by
  that lock — the single-thread path is where these gains land.
- No change to the bin math, dedup/accounting, strand logic, or output formatting.

### Two pre-existing issues surfaced by the harness (NOT changed here)

Both are properties of the **original** `make_heatmap`; this port preserves them to
stay byte-identical, but they are worth a maintainer's attention:

1. **Bin definition files must not end with a trailing newline.** The original
   `bin_parser` treats a trailing newline as an empty final record and exits with
   "Bin file contains formatting errors". The bundled `*_chr_size.txt` files and the
   test bin file follow this. (A one-line fix is possible but would be a behavior
   change, so it is deliberately not made here.)
2. **A file of very many consecutive malformed lines can overflow the stack.**
   `update()` recurses once per bad line (`return update(he)`), so e.g. tens of
   thousands of consecutive unparseable lines can exhaust the stack (the original
   crashes; the optimized version's tail-recursive parse happens to survive, but the
   underlying unbounded recursion remains). A future hardening would convert that
   recursion to a loop.

### Tests

`tests/` is a self-contained golden-master harness (`g++` + `perl` + coreutils; no
cluster, no real data):

```
bash tests/run_harness.sh          # builds ref+candidate, runs the full matrix + fuzz
CXX=g++ CXXFLAGS=-O3 bash tests/run_harness.sh   # override compiler/flags
```

- `gen_data.pl` — synthetic genes (stranded + unstranded, overlapping spans, in-range
  user anchors) and the same hit set in every format (G/g/b/e/c + stranded + plus/minus
  split), a bin file, and a messy/edge-case file.
- `gen_fuzz.pl` — adversarial line generator (wrong field counts, non-numeric values,
  floats in integer slots, signs, scientific notation, hex-looking tokens, mixed
  whitespace, CRLF, `track`, blanks) for the parser-equivalence check.
- `run_harness.sh` — builds reference (`make_heatmap.cpp`) + candidate, byte-diffs
  output at `-t 1` across the matrix and the fuzz, characterizes `-t 1` vs `-t N`,
  and reports a wall-time + peak-RSS comparison. Last run: **108/108 pass**.

> Before deploying on O2, rebuild **both** binaries with the cluster's toolchain and
> run `tests/run_harness.sh` there — it rebuilds the reference from source, so it
> validates byte-identity under whatever compiler/libstdc++ the lab uses.
