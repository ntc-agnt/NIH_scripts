#!/usr/bin/env bash
# Golden-master harness for make_heatmap_optimized.cpp.
#
# Builds the reference (make_heatmap.cpp) and candidate (make_heatmap_optimized.cpp)
# from source with the SAME compiler/flags, then asserts the candidate's output is
# BYTE-IDENTICAL to the reference at -t 1 across a path-complete option matrix and a
# set of edge-case fixtures. Also characterizes -t 1 vs -t N reproducibility (the
# tool's accumulation order is thread-interleaving dependent, so threaded float
# output can differ in the last ULPs -- that is pre-existing, not introduced here)
# and reports a wall-time + peak-RSS comparison on a large input.
#
# The contract this gates: optimized == original, byte-for-byte, at -t 1, for every
# supported mode. Self-contained: g++ + perl + coreutils.
#
# Env: CXX (default g++), CXXFLAGS (default "-O3"), MH_KEEP=1 to keep the workdir.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC_REF="$HERE/../make_heatmap.cpp"
SRC_CAND="$HERE/../make_heatmap_optimized.cpp"
GEN="$HERE/gen_data.pl"
CXX="${CXX:-g++}"; CXXFLAGS="${CXXFLAGS:--O3}"
WORK="$(mktemp -d)"; if [ "${MH_KEEP:-0}" = 1 ]; then trap 'echo "kept: $WORK"' EXIT; else trap 'rm -rf "$WORK"' EXIT; fi

PASS=0; FAIL=0
pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

echo "=== build reference + candidate ($CXX $CXXFLAGS) ==="
$CXX $CXXFLAGS -o "$WORK/mh_ref"  "$SRC_REF"  -lpthread 2>"$WORK/ref.build"  || { echo "ref build FAILED"; cat "$WORK/ref.build"; exit 1; }
$CXX $CXXFLAGS -o "$WORK/mh_cand" "$SRC_CAND" -lpthread 2>"$WORK/cand.build" || { echo "candidate build FAILED"; cat "$WORK/cand.build"; exit 1; }
echo "  built mh_ref + mh_cand"
[ -s "$WORK/cand.build" ] && { echo "  candidate build warnings:"; sed 's/^/    /' "$WORK/cand.build" | head -20; }
REF="$WORK/mh_ref"; CAND="$WORK/mh_cand"

echo "=== generate fixtures ==="
FX="$WORK/fx"
perl "$GEN" --outdir "$FX" --seed 7 --ngenes 80 --nhits 40000 >/dev/null || { echo "gen failed"; exit 1; }
echo "  fixtures in $FX"

# run both binaries with identical args (OUT substituted per arm) and byte-diff.
# $1 name; rest = argv after the binary, using the token OUT for the output path
# and file paths relative to $FX.
check(){
  local name="$1"; shift
  local args=("$@")
  local aref=() acand=()
  local tok sub
  for tok in "${args[@]}"; do
    if [ "$tok" = OUT ]; then
      aref+=("$WORK/$name.ref"); acand+=("$WORK/$name.cand")
    elif [[ "$tok" == *@* ]]; then           # one or more @fixture refs (comma lists ok)
      sub="${tok//@/$FX/}"; aref+=("$sub"); acand+=("$sub")
    else
      aref+=("$tok"); acand+=("$tok")
    fi
  done
  ( cd "$WORK" && "$REF"  "${aref[@]}"  >"$name.ref.log"  2>&1 ); local rr=$?
  ( cd "$WORK" && "$CAND" "${acand[@]}" >"$name.cand.log" 2>&1 ); local rc=$?
  if [ $rr -ne 0 ] || [ $rc -ne 0 ]; then fail "$name (exit ref=$rr cand=$rc)"; return; fi
  if [ ! -f "$WORK/$name.ref" ]; then fail "$name (no reference output produced)"; return; fi
  if cmp -s "$WORK/$name.ref" "$WORK/$name.cand"; then
    pass "$name  ($(grep -vc -e '^Match' -e '^Anchor' -e '^Hit' -e '^Bin' -e '^Gene' -e '^$' "$WORK/$name.ref" 2>/dev/null) data rows)"
  else
    fail "$name  (output differs)"; diff "$WORK/$name.ref" "$WORK/$name.cand" | head -6 | sed 's/^/      /'
  fi
}

echo
echo "=== single unstranded file (G/g/b): -s b, -l p/d/c only ==="
# bin form -b c needs: [hits] [genes] [out] [start] [size] [count]
check G_total_phys   -h G -b c -s b -l p -a p -d p -- @hits.G @genes_stranded.txt OUT 0 100 20
check G_total_gen    -h G -b c -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -500 50 25
check G_center       -h G -b c -s b -l c -a p -d p -- @hits.G @genes_stranded.txt OUT 0 100 20
check G_physend      -h G -b c -s b -l d -a e -d g -- @hits.G @genes_stranded.txt OUT -1000 100 21
check g_oneBased     -h g -b c -s b -l p -a p -d p -- @hits.g @genes_stranded.txt OUT 0 100 20
check b_bed3         -h b -b c -s b -l p -a p -d p -- @hits.b @genes_stranded.txt OUT 0 100 20
check G_avg          -v a -h G -b c -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -500 50 25
check G_density      -v d -h G -b c -s b -l p -a s -d g -- @hits.G @genes_stranded.txt OUT -500 50 25
check G_coverage     -v c -h G -b c -s b         -a s -d g -- @hits.G @genes_stranded.txt OUT -500 50 25
check G_anchor_user  -h G -b c -s b -l p -a u -d p -- @hits.G @genes_stranded.txt OUT -100 20 15

echo
echo "=== stranded single file: bed6 (e) and cppmatch-stranded (cs): -s s/o, -l s/e ==="
check e_samestrand   -h e -s s -b c -l s -a s -d g -- @hits.e  @genes_stranded.txt OUT -300 30 21
check e_oppstrand    -h e -s o -b c -l e -a e -d g -- @hits.e  @genes_stranded.txt OUT -300 30 21
check e_both         -h e -s b -b c -l p -a p -d p -- @hits.e  @genes_stranded.txt OUT 0 100 20
check cs_samestrand  -h c -s s -b c -l s -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21
check cs_oppstrand   -h c -s o -b c -l e -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21
check c_unstranded   -h c -s b -b c -l p -a p -d p -- @hits.c  @genes_stranded.txt OUT 0 100 20
check cs_coverage    -v c -h c -s s -b c -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21

echo
echo "=== -p/-m split (plus/minus, G format): strand-specific without stranded single file ==="
# -p/-m + -b c needs: [genes] [out] [start] [size] [count]
check pm_same        -p @hits_plus.G -m @hits_minus.G -h G -s s -b c -l s -a s -d g -- @genes_stranded.txt OUT -300 30 21
check pm_opp         -p @hits_plus.G -m @hits_minus.G -h G -s o -b c -l e -a e -d g -- @genes_stranded.txt OUT -300 30 21
check pm_both        -p @hits_plus.G -m @hits_minus.G -h G -s b -b c -l p -a p -d p -- @genes_stranded.txt OUT 0 100 20
check pm_coverage    -p @hits_plus.G -m @hits_minus.G -h G -s s -v c -b c -a s -d g -- @genes_stranded.txt OUT -300 30 21
check p_only         -p @hits_plus.G                  -h G -s s -b c -l s -a s -d g -- @genes_stranded.txt OUT -300 30 21

echo
echo "=== bin forms: file (-b f) and variable (-b v) ==="
check binfile        -h G -b f -s b -l p -a p -d p    @hits.G @genes_stranded.txt OUT @bins.txt
check binfile_gen    -h G -b f -s b -l p -a s -d g    @hits.G @genes_stranded.txt OUT @bins.txt
check var_bins       -h G -b v -s b -l p -a p -d p    @hits.G @genes_stranded.txt OUT 50
check var_bins_gen   -h G -b v -s b -l p -a s -d g    @hits.G @genes_stranded.txt OUT 50
check var_bins_frac  -h G -b v -s b -l p -a p -d p    @hits.G @genes_stranded.txt OUT 5000   # count > span -> fractional bins
check var_bins_pm    -p @hits_plus.G -m @hits_minus.G -h G -s s -b v -l s -a s -d g @genes_stranded.txt OUT 40

echo
echo "=== nostrand gene list + --nostrand / physical strand-independent branch ==="
check nostrand_flag  --nostrand -h G -b c    @hits.G @genes_nostrand.txt OUT 0 100 20
check phys_indep     -h G -b c -s b -l p -a p -d p -- @hits.G @genes_nostrand.txt OUT 0 100 20
check nostrand_var   --nostrand -h G -b v    @hits.G @genes_nostrand.txt OUT 50

echo
echo "=== header suppression + edge cases ==="
check nohead         --nohead -h G -b c -s b -l p -a p -d p -- @hits.G @genes_stranded.txt OUT 0 100 20
check messy_input    -h G -b c -s b -l p -a p -d p -- @hits_messy.G @genes_stranded.txt OUT -300 50 20
check empty_result   -h G -b c -s b -l p -a p -d p -- @hits.G @genes_stranded.txt OUT 9000000 100 5   # bins far from any hit

echo
echo "=== multi-file (comma-separated) hit inputs, single output ==="
check multi_hits     -h G -b c -s b -l p -a p -d p -- "@hits.G,@hits.G" @genes_stranded.txt OUT 0 100 20

echo
echo "=== -t 1 vs -t N reproducibility characterization (candidate) ==="
# integer-valued bed totals: addition is exact -> should be -tN stable.
"$CAND" -h b -b c -s b -l p -a p -d p -t 1 -- "$FX/hits.b" "$FX/genes_stranded.txt" "$WORK/tn_int_t1" 0 100 20 >/dev/null 2>&1
"$CAND" -h b -b c -s b -l p -a p -d p -t 4 -- "$FX/hits.b" "$FX/genes_stranded.txt" "$WORK/tn_int_t4" 0 100 20 >/dev/null 2>&1
if cmp -s "$WORK/tn_int_t1" "$WORK/tn_int_t4"; then echo "  integer-value bed: -t1 == -t4 (exact, expected)"; else echo "  NOTE integer-value bed: -t1 != -t4 (unexpected -- investigate)"; fi
# float bedGraph totals: summation order is interleaving-dependent -> may differ.
"$CAND" -h G -b c -s b -l p -a p -d p -t 1 -- "$FX/hits.G" "$FX/genes_stranded.txt" "$WORK/tn_flt_t1" 0 100 20 >/dev/null 2>&1
"$CAND" -h G -b c -s b -l p -a p -d p -t 4 -- "$FX/hits.G" "$FX/genes_stranded.txt" "$WORK/tn_flt_t4" 0 100 20 >/dev/null 2>&1
if cmp -s "$WORK/tn_flt_t1" "$WORK/tn_flt_t4"; then echo "  float bedGraph:    -t1 == -t4 on this input"; else echo "  float bedGraph:    -t1 != -t4 (pre-existing thread accumulation-order noise; contract is -t1)"; fi

echo
echo "=== wall-time + peak-RSS: reference vs candidate on a large input (-t 1) ==="
BIG="$WORK/big"; perl "$GEN" --outdir "$BIG" --seed 42 --ngenes 2000 --nhits 3000000 --coord-max 20000000 >/dev/null
sample_rss(){
  local t0 t1 pid max hwm rc
  t0=$(date +%s.%N); "$@" >/dev/null 2>&1 & pid=$!; max=0
  while kill -0 "$pid" 2>/dev/null; do
    hwm=$(awk '/VmHWM/{print $2}' "/proc/$pid/status" 2>/dev/null)
    [ -n "$hwm" ] && [ "$hwm" -gt "$max" ] && max=$hwm; sleep 0.02
  done
  wait "$pid"; rc=$?; t1=$(date +%s.%N)
  awk -v m="$max" -v a="$t0" -v b="$t1" 'BEGIN{printf "%d %.2f %d", m, b-a, '"$rc"'}'
}
BARGS=(-h G -b c -s b -l p -a s -d g -- "$BIG/hits.G" "$BIG/genes_stranded.txt")
read RO_RSS RO_WALL RO_RC < <(sample_rss "$REF"  "${BARGS[@]}" "$WORK/big.ref"  -2000 10 400)
read RP_RSS RP_WALL RP_RC < <(sample_rss "$CAND" "${BARGS[@]}" "$WORK/big.cand" -2000 10 400)
if [ "$RO_RC" -eq 0 ] && [ "$RP_RC" -eq 0 ] && cmp -s "$WORK/big.ref" "$WORK/big.cand"; then
  pass "big input: byte-identical ($(grep -c '	' "$WORK/big.ref") lines)"
else
  fail "big input: differ or nonzero (ref=$RO_RC cand=$RP_RC)"
fi
printf "  %-10s %12s %10s\n" "impl" "maxRSS(MB)" "wall(s)"
awk -v r="$RO_RSS" -v w="$RO_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","reference", r/1024, w}'
awk -v r="$RP_RSS" -v w="$RP_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","candidate", r/1024, w}'
awk -v o="$RO_WALL" -v p="$RP_WALL" 'BEGIN{ if(p>0) printf "  wall-time speedup:  %.2fx\n", o/p }'
awk -v o="$RO_RSS" -v p="$RP_RSS" 'BEGIN{ if(p>0) printf "  peak-RSS ratio:     %.2fx\n", o/p }'

echo
echo "=== parser-equivalence fuzz (manual parse vs operator>>) ==="
# Feed the SAME adversarial file to reference and candidate under each parser and
# byte-diff the output. Any divergence in which lines are accepted/rejected, or in
# a parsed value, changes the matrix -> caught. Run several seeds.
fuzz_check(){
  local name="$1"; shift   # rest = parser flags (hit file + genes + OUT + bins appended)
  local args=("$@")
  local aref=() acand=() tok sub
  for tok in "${args[@]}"; do
    if [ "$tok" = OUT ]; then aref+=("$WORK/$name.ref"); acand+=("$WORK/$name.cand")
    elif [[ "$tok" == *@* ]]; then sub="${tok//@/$FZ/}"; aref+=("$sub"); acand+=("$sub")
    else aref+=("$tok"); acand+=("$tok"); fi
  done
  "$REF"  "${aref[@]}"  >/dev/null 2>&1; local rr=$?
  "$CAND" "${acand[@]}" >/dev/null 2>&1; local rc=$?
  # both may legitimately exit nonzero on some inputs; what matters is identical output
  if [ -f "$WORK/$name.ref" ] && cmp -s "$WORK/$name.ref" "$WORK/$name.cand"; then
    pass "fuzz $name  (ref==cand; $(grep -vc -e '^Match' -e '^Anchor' -e '^Hit' -e '^Bin' -e '^Gene' "$WORK/$name.ref" 2>/dev/null) rows)"
  elif [ ! -f "$WORK/$name.ref" ] && [ ! -f "$WORK/$name.cand" ]; then
    pass "fuzz $name  (both produced no output, identical)"
  else
    fail "fuzz $name  (output differs; ref_rc=$rr cand_rc=$rc)"
    diff "$WORK/$name.ref" "$WORK/$name.cand" 2>/dev/null | head -6 | sed 's/^/      /'
  fi
}
# NOTE: both scripts share a pre-existing fragility -- update() recurses once per
# bad line (return update(he)), so a file of tens of thousands of CONSECUTIVE bad
# lines can overflow the stack (the original crashes; the candidate's manual parse
# tail-optimizes and survives). That is orthogonal to parser-value equivalence, so
# the fuzz uses files small enough never to approach that depth, across many seeds.
for seed in 1 2 3 4 5 6 7 8; do
  FZ="$WORK/fuzz$seed"; mkdir -p "$FZ"
  perl "$HERE/gen_fuzz.pl" --out "$FZ/f" --lines 1200 --seed "$seed" --chroms chr1,chr2,chr3 >/dev/null
  cp "$FX/genes_stranded.txt" "$FZ/g"
  fuzz_check "G_s$seed"   -h G -b c -s b -l p -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "g_s$seed"   -h g -b c -s b -l p -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "b_s$seed"   -h b -b c -s b -l p -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "e_s$seed"   -h e -b c -s b -l p -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "es_s$seed"  -h e -b c -s s -l s -a s -d g -- @f @g OUT -200 20 20
  fuzz_check "c_s$seed"   -h c -b c -s b -l p -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "cs_s$seed"  -h c -b c -s s -l s -a s -d g -- @f @g OUT -200 20 20
  fuzz_check "Gcov_s$seed" -h G -b c -s b -v c -a p -d p -- @f @g OUT 0 100 20
  fuzz_check "ccov_s$seed" -h c -b c -s b -v c -a p -d p -- @f @g OUT 0 100 20
done

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
