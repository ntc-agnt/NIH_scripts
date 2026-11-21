#!/usr/bin/env bash
# Validation harness for make_heatmap_optimized_summation.cpp (the parallel /
# compensated-summation variant). Unlike run_harness.sh (which proves BYTE identity
# for the _optimized build), this checks the summation variant's weaker-but-defined
# contract:
#   1. FUNCTIONAL identity to the original at -t 1 (relative tolerance; the printed
#      output is only ~6 significant figures, so a handful of cells may differ by one
#      unit in the last digit).
#   2. THREAD-COUNT INDEPENDENCE: -t 1 == -t 2 == -t 4 byte-for-byte (compensated
#      summation makes the result independent of how hits are distributed to threads).
#   3. a wall-time comparison across thread counts.
# Self-contained: g++ + perl + coreutils. Env: CXX (g++), CXXFLAGS (-O3), TOL (2e-5).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC_REF="$HERE/../make_heatmap.cpp"
SRC_SUM="$HERE/../make_heatmap_optimized_summation.cpp"
GEN="$HERE/gen_data.pl"; CMP="$HERE/compare_matrix.py"
CXX="${CXX:-g++}"; CXXFLAGS="${CXXFLAGS:--O3}"; TOL="${TOL:-2e-5}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
PASS=0; FAIL=0
pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

echo "=== build reference + summation variant ($CXX $CXXFLAGS) ==="
$CXX $CXXFLAGS -o "$WORK/ref" "$SRC_REF" -lpthread 2>"$WORK/rb" || { echo "ref build failed"; cat "$WORK/rb"; exit 1; }
$CXX $CXXFLAGS -std=c++11 -pthread -o "$WORK/sum" "$SRC_SUM" 2>"$WORK/sb" || { echo "summation build failed"; cat "$WORK/sb"; exit 1; }
echo "  built"
FX="$WORK/fx"; perl "$GEN" --outdir "$FX" --seed 7 --ngenes 1500 --nhits 300000 >/dev/null
REF="$WORK/ref"; SUM="$WORK/sum"

# $1 name; rest = argv with OUT token where the output prefix goes
check(){
  local name="$1"; shift
  local aR=() a1=() a2=() a4=() tok
  for tok in "$@"; do
    case "$tok" in
      OUT) aR+=("$WORK/$name.r"); a1+=("$WORK/$name.1"); a2+=("$WORK/$name.2"); a4+=("$WORK/$name.4");;
      *@*) local f="${tok//@/$FX/}"; aR+=("$f"); a1+=("$f"); a2+=("$f"); a4+=("$f");;
      *)   aR+=("$tok"); a1+=("$tok"); a2+=("$tok"); a4+=("$tok");;
    esac
  done
  "$REF"           "${aR[@]}" >/dev/null 2>&1; local rcR=$?
  "$SUM" -t 1      "${a1[@]}" >/dev/null 2>&1; local rc1=$?
  "$SUM" -t 2      "${a2[@]}" >/dev/null 2>&1; local rc2=$?
  "$SUM" -t 4      "${a4[@]}" >/dev/null 2>&1; local rc4=$?
  if [ $rcR -ne 0 ] || [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ] || [ $rc4 -ne 0 ]; then
    fail "$name  nonzero exit (ref=$rcR t1=$rc1 t2=$rc2 t4=$rc4)"; return; fi
  # thread-count independence
  if cmp -s "$WORK/$name.1" "$WORK/$name.2" && cmp -s "$WORK/$name.1" "$WORK/$name.4"; then
    pass "$name  thread-independent (-t1 = -t2 = -t4)"
  else
    fail "$name  thread-DEPENDENT output"
  fi
  # functional identity vs original
  if python3 "$CMP" "$WORK/$name.r" "$WORK/$name.1" "$TOL" >"$WORK/$name.cmp" 2>&1; then
    pass "$name  functional vs original ($(cat "$WORK/$name.cmp"))"
  else
    fail "$name  exceeds tolerance ($(cat "$WORK/$name.cmp"))"
  fi
}

echo "=== functional identity + thread-count independence across modes ==="
check total        -h G -b c -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -2000 10 400
check average      -h G -b c -v a -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -2000 10 400
check density      -h G -b c -v d -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -2000 10 400
check coverage     -h G -b c -v c -a s -d g -- @hits.G @genes_stranded.txt OUT -2000 10 400
check variable     -h G -b v -s b -a s -d g @hits.G @genes_stranded.txt OUT 100
check binfile      -h G -b f -s b -l p -a s -d g @hits.G @genes_stranded.txt OUT @bins.txt
check strand_same  -p @hits_plus.G -m @hits_minus.G -h G -s s -a s -l s -d g -b c -- @genes_stranded.txt OUT -2000 10 400
check strand_opp   -p @hits_plus.G -m @hits_minus.G -h G -s o -a s -l s -d g -b c -- @genes_stranded.txt OUT -2000 10 400
check ext_bed      -h e -s s -l s -a s -d g -b c -- @hits.e @genes_stranded.txt OUT -300 30 21
check bed_counts   -h b -b c -s b -l p -a p -d p -- @hits.b @genes_stranded.txt OUT 0 100 40
check nostrand     --nostrand -h G -b c @hits.G @genes_nostrand.txt OUT 0 100 20
check multi_hits   -h G -b c -s b -l p -a p -d p -- "@hits.G,@hits.G" @genes_stranded.txt OUT 0 100 20

echo
echo "=== wall-time across thread counts (large, compute-heavy input) ==="
BIG="$WORK/big"; perl "$GEN" --outdir "$BIG" --seed 42 --ngenes 8000 --nhits 3000000 --coord-max 20000000 >/dev/null
BA=(-h G -b c -s b -l p -a s -d g -- "$BIG/hits.G" "$BIG/genes_stranded.txt")
wt(){ local t0 t1; t0=$(date +%s.%N); "$@" >/dev/null 2>&1; t1=$(date +%s.%N); awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}'; }
w1=$(wt "$SUM" -t 1 "${BA[@]}" "$BIG/o1" -2000 10 400)
w2=$(wt "$SUM" -t 2 "${BA[@]}" "$BIG/o2" -2000 10 400)
w4=$(wt "$SUM" -t 4 "${BA[@]}" "$BIG/o4" -2000 10 400)
cmp -s "$BIG/o1" "$BIG/o2" && cmp -s "$BIG/o1" "$BIG/o4" && pass "big input thread-independent" || fail "big input thread-dependent"
printf "  -t1=%ss  -t2=%ss  -t4=%ss  (cores here: %s)\n" "$w1" "$w2" "$w4" "$(nproc 2>/dev/null || echo ?)"
awk -v a="$w1" -v b="$w4" 'BEGIN{ if(b>0) printf "  speedup -t1->-t4: %.2fx (bounded by serial input reading; scales further on more cores)\n", a/b }'

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
