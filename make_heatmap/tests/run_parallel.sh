#!/usr/bin/env bash
# Byte-identity harness for make_heatmap_optimized_parallel.cpp.
#
# CONTRACT: this build is BYTE-FOR-BYTE identical to the reference (make_heatmap.cpp)
# at EVERY thread count -- not just -t 1. It reaches that by reading hit lines serially
# in file order, parsing+matching them in parallel, replaying each cell's updates in
# strict file order (so the per-bin float summation order is exactly the -t 1 order),
# and formatting the finished output rows in parallel into per-thread buffers that are
# concatenated in sorted row order. No value and no row position changes, so -t N output
# is identical to the serial reference bit-for-bit (unlike the original/_optimized -t N,
# whose single global accumulation lock makes float sums thread-interleaving dependent).
#
# Builds the reference from make_heatmap.cpp and the candidate from
# make_heatmap_optimized_parallel.cpp with the same compiler/flags, then asserts the
# candidate at -t 1, -t 2 and -t 4 is `cmp`-identical to the reference (-t 1) across a
# path-complete option matrix (mirrors run_harness.sh), including the multi-index path
# (comma-separated hit files -> one output with per-file column blocks), which exercises
# every branch of the parallel row formatter.
#
# Self-contained: g++ + perl + coreutils. Env: CXX (default g++), CXXFLAGS
# (default "-O3 -std=c++11"), MH_THREADS (default "1 2 4"), MH_KEEP=1 to keep the workdir.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC_REF="$HERE/../make_heatmap.cpp"
SRC_CAND="$HERE/../make_heatmap_optimized_parallel.cpp"
GEN="$HERE/gen_data.pl"
CXX="${CXX:-g++}"; CXXFLAGS="${CXXFLAGS:--O3 -std=c++11}"
WORK="$(mktemp -d)"; if [ "${MH_KEEP:-0}" = 1 ]; then trap 'echo "kept: $WORK"' EXIT; else trap 'rm -rf "$WORK"' EXIT; fi
THREADS="${MH_THREADS:-1 2 4}"

PASS=0; FAIL=0
pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

echo "=== build reference + parallel candidate ($CXX $CXXFLAGS) ==="
$CXX $CXXFLAGS -o "$WORK/mh_ref"  "$SRC_REF"  -lpthread 2>"$WORK/ref.build"  || { echo "ref build FAILED"; cat "$WORK/ref.build"; exit 1; }
$CXX $CXXFLAGS -pthread -o "$WORK/mh_par" "$SRC_CAND" 2>"$WORK/cand.build" || { echo "candidate build FAILED"; cat "$WORK/cand.build"; exit 1; }
echo "  built mh_ref + mh_par"
[ -s "$WORK/cand.build" ] && { echo "  candidate build warnings:"; sed 's/^/    /' "$WORK/cand.build" | head -20; }
REF="$WORK/mh_ref"; CAND="$WORK/mh_par"

echo "=== generate fixtures ==="
FX="$WORK/fx"
perl "$GEN" --outdir "$FX" --seed 7 --ngenes 80 --nhits 40000 >/dev/null || { echo "gen failed"; exit 1; }
echo "  fixtures in $FX"

# check NAME ARGS... : reference at -t1 is the gold; candidate at each THREADS count must
# be byte-identical to it. Token OUT -> per-run output path; token with @ -> fixture path(s).
check(){
  local name="$1"; shift; local args=("$@")
  local base=() tok sub
  for tok in "${args[@]}"; do
    if [ "$tok" = OUT ]; then base+=("OUTSENTINEL")
    elif [[ "$tok" == *@* ]]; then sub="${tok//@/$FX/}"; base+=("$sub")
    else base+=("$tok"); fi
  done
  local r=("${base[@]/OUTSENTINEL/$WORK/$name.ref}")
  ( cd "$WORK" && "$REF" "${r[@]}" >/dev/null 2>&1 ); local rr=$?
  if [ $rr -ne 0 ] || [ ! -f "$WORK/$name.ref" ]; then fail "$name (reference exit=$rr / no output)"; return; fi
  local nt ok=1
  for nt in $THREADS; do
    local c=("${base[@]/OUTSENTINEL/$WORK/$name.t$nt}")
    ( cd "$WORK" && "$CAND" -t "$nt" "${c[@]}" >/dev/null 2>&1 ); local rc=$?
    if [ $rc -ne 0 ]; then fail "$name (candidate -t$nt exit=$rc)"; ok=0; continue; fi
    if ! cmp -s "$WORK/$name.ref" "$WORK/$name.t$nt"; then
      fail "$name (-t$nt differs from reference -t1)"; diff "$WORK/$name.ref" "$WORK/$name.t$nt" | head -4 | sed 's/^/      /'; ok=0
    fi
  done
  [ $ok -eq 1 ] && pass "$name  ($(grep -vc -e '^Match' -e '^Anchor' -e '^Hit' -e '^Bin' -e '^Gene' -e '^$' "$WORK/$name.ref" 2>/dev/null) rows, -t {$THREADS} all == ref)"
}

echo
echo "=== single unstranded file (G/g/b) ==="
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
echo "=== stranded single file: bed6 (e) and cppmatch-stranded (cs) ==="
check e_samestrand   -h e -s s -b c -l s -a s -d g -- @hits.e  @genes_stranded.txt OUT -300 30 21
check e_oppstrand    -h e -s o -b c -l e -a e -d g -- @hits.e  @genes_stranded.txt OUT -300 30 21
check e_both         -h e -s b -b c -l p -a p -d p -- @hits.e  @genes_stranded.txt OUT 0 100 20
check cs_samestrand  -h c -s s -b c -l s -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21
check cs_oppstrand   -h c -s o -b c -l e -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21
check c_unstranded   -h c -s b -b c -l p -a p -d p -- @hits.c  @genes_stranded.txt OUT 0 100 20
check cs_coverage    -v c -h c -s s -b c -a s -d g -- @hits.cs @genes_stranded.txt OUT -300 30 21

echo
echo "=== -p/-m split (plus/minus, G format) ==="
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
check var_bins_frac  -h G -b v -s b -l p -a p -d p    @hits.G @genes_stranded.txt OUT 5000
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
check empty_result   -h G -b c -s b -l p -a p -d p -- @hits.G @genes_stranded.txt OUT 9000000 100 5

echo
echo "=== multi-file (comma-separated) hit inputs, single output -> multi-index column blocks ==="
# These exercise the max_index>1 && outfile.size()==1 branch of the parallel row formatter.
check multi_total    -h G    -b c -s b -l p -a p -d p -- "@hits.G,@hits.G" @genes_stranded.txt OUT 0 100 20
check multi_total_gn -h G    -b c -s b -l p -a s -d g -- "@hits.G,@hits.G" @genes_stranded.txt OUT -500 50 25
check multi_avg      -v a -h G -b c -s b -l p -a s -d g -- "@hits.G,@hits.G" @genes_stranded.txt OUT -500 50 25
check multi_density  -v d -h G -b c -s b -l p -a s -d g -- "@hits.G,@hits.G" @genes_stranded.txt OUT -1000 50 40

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed  (candidate -t {$THREADS} each byte-identical to reference -t 1) ==="
[ $FAIL -eq 0 ]
