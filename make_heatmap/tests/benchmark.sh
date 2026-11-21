#!/usr/bin/env bash
# make_heatmap performance benchmark (hyperfine-based). Characterizes wall time across
# the dimensions that actually move it, so you can decide on real data whether a given
# build/flag is worth it. Measures ABSOLUTE seconds (not just ratios) -- if a scenario
# is already a couple of seconds, optimizing it further is not worth much.
#
# Dimensions covered (representative points, not a full cartesian blow-up):
#   * value type ...... integer counts vs normalized (float) values
#   * hit mode ........ single-position (-v t) vs per-nt coverage (-v c)
#   * depth ........... low vs high hit-line count
#   * gene-list size .. small / medium / large
#   * binning ......... fixed-size (-b c) vs scaled-per-feature (-b v), on short+long features
#   * threads ......... a sweep, for any build that honors -t (summation/parallel variants)
#
# Requires: hyperfine, g++, perl, coreutils.
# Env:
#   BENCH_HIGH   high-depth hit-line count (default 10000000)
#   BENCH_LOW    low-depth  hit-line count (default 2000000)
#   BENCH_THREADS space-separated thread counts for the sweep (default "1 2 4 6")
#   BENCH_RUNS   hyperfine runs per command (default 3)
#   BENCH_EXTRA_BIN  path to another binary to include (e.g. a parallel build); honors -t
#   BENCH_OUT    markdown output dir (default ./bench_results)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$HERE/.."
HIGH="${BENCH_HIGH:-10000000}"; LOW="${BENCH_LOW:-2000000}"
THREADS="${BENCH_THREADS:-1 2 4 6}"; RUNS="${BENCH_RUNS:-3}"
OUT="${BENCH_OUT:-$PWD/bench_results}"; mkdir -p "$OUT"
command -v hyperfine >/dev/null || { echo "need hyperfine (apt install hyperfine or grab a release binary)"; exit 1; }
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT

echo "=== build binaries ==="
g++ -O3 -o "$W/orig" "$ROOT/make_heatmap.cpp" -lpthread && echo "  orig"
g++ -O3 -std=c++11 -pthread -o "$W/opt" "$ROOT/make_heatmap_optimized.cpp" && echo "  optimized"
EXTRA="${BENCH_EXTRA_BIN:-}"; [ -n "$EXTRA" ] && echo "  extra: $EXTRA"

echo "=== generate synthetic fixtures (high=$HIGH low=$LOW lines) ==="
gen(){ perl "$HERE/gen_data.pl" "$@" >/dev/null; }
# large + small/medium gene lists, high + low depth, short + long features
gen --outdir "$W/hi" --seed 1 --ngenes 20000 --nhits "$HIGH" --coord-max 60000000
gen --outdir "$W/lo" --seed 1 --ngenes 20000 --nhits "$LOW"  --coord-max 60000000
head -5000 "$W/hi/genes_stranded.txt" > "$W/hi/genes_5k.txt"
head -1000 "$W/hi/genes_stranded.txt" > "$W/hi/genes_1k.txt"
# long-feature gene list (for scaled-bin cost on big features): widen start..end
awk -F'\t' 'BEGIN{OFS="\t"}{ $5=$4+50000; print }' "$W/hi/genes_stranded.txt" > "$W/hi/genes_long.txt"
# integer- and float-valued copies of the high-depth hit file (same positions)
awk -F'\t' 'BEGIN{OFS="\t"}{v=($4<0?-$4:$4); print $1,$2,$3,(int(v*100)%97+1)}' "$W/hi/hits.G" > "$W/hi/hits_int.G"
awk -F'\t' 'BEGIN{OFS="\t"}{v=($4<0?-$4:$4); printf "%s\t%s\t%s\t%.4f\n",$1,$2,$3,(int(v*100)%9700+1)/100.0}' "$W/hi/hits.G" > "$W/hi/hits_flt.G"
awk -F'\t' 'BEGIN{OFS="\t"}{v=($4<0?-$4:$4); printf "%s\t%s\t%s\t%.4f\n",$1,$2,$3,(int(v*100)%9700+1)/100.0}' "$W/lo/hits.G" > "$W/lo/hits_flt.G"

echo "=== DIMENSION SWEEP (optimized, single build; absolute seconds) ==="
G20="$W/hi/genes_stranded.txt"
HF=()                                   # hyperfine -n/command argument array (no eval)
add(){ local tag="$1"; shift; HF+=(-n "$tag" "$*"); }   # paths carry no spaces -> $* is safe
add "value:int   point 20k 400 HIGH"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/hi/hits_int.G $G20 $W/o -2000 10 400"
add "value:float point 20k 400 HIGH"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 400"
add "mode:coverage float 20k 400 HIGH" "$W/opt -h G -b c -v c -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 400"
add "depth:LOW  float point 20k 400"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/lo/hits_flt.G $G20 $W/o -2000 10 400"
add "genes:1k  float point 400 HIGH"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $W/hi/genes_1k.txt $W/o -2000 10 400"
add "genes:5k  float point 400 HIGH"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $W/hi/genes_5k.txt $W/o -2000 10 400"
add "bins:100  float point 20k HIGH"   "$W/opt -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 100"
add "bins:scaled 20k shortfeat HIGH"   "$W/opt -h G -b v -s b -a s -d g $W/hi/hits_flt.G $G20 $W/o 400"
add "bins:scaled 20k LONGfeat HIGH"    "$W/opt -h G -b v -s b -a s -d g $W/hi/hits_flt.G $W/hi/genes_long.txt $W/o 400"
hyperfine --warmup 1 --runs "$RUNS" --export-markdown "$OUT/dimensions.md" "${HF[@]}" 2>&1 | tail -3
echo "  -> $OUT/dimensions.md"

echo "=== optimized vs original (absolute win), high depth ==="
hyperfine --warmup 1 --runs "$RUNS" --export-markdown "$OUT/opt_vs_orig.md" \
  -n 'original  20k point HIGH' "$W/orig -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 400" \
  -n 'optimized 20k point HIGH' "$W/opt  -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 400" \
  2>&1 | tail -3
echo "  -> $OUT/opt_vs_orig.md"

if [ -n "$EXTRA" ]; then
  echo "=== THREAD SWEEP ($EXTRA; high depth, 20k genes, point) ==="
  TH=()
  for t in $THREADS; do TH+=(-n "extra -t$t" "$EXTRA -t $t -h G -b c -s b -l p -a s -d g -- $W/hi/hits_flt.G $G20 $W/o -2000 10 400"); done
  hyperfine --warmup 1 --runs "$RUNS" --export-markdown "$OUT/threads.md" "${TH[@]}" 2>&1 | tail -3
  echo "  -> $OUT/threads.md"
fi

echo
echo "=== results in $OUT/ (absolute seconds; decide per your real depth/gene-list) ==="
cat "$OUT"/dimensions.md 2>/dev/null
