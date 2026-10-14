#!/usr/bin/env bash
# Byte-diff harness for bedgraphs2stdBedGraph_streaming.
#
# Proves the streaming script produces byte-identical .bedGraph output vs the
# original across a matrix of synthetic inputs, checks the misaligned / unsorted
# / atomic error paths, and reports a peak-RSS + wall-time comparison.
# Self-contained: needs only perl + coreutils, no cluster or real data.
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
ORIG="$HERE/../bedgraphs2stdBedGraph"
PROTO="$HERE/../bedgraphs2stdBedGraph_streaming"
GEN="$HERE/gen_bedgraphs.pl"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

PASS=0; FAIL=0
pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }
mask_log(){ sed 's/^Initiated .*/Initiated <TS>/' "$1"; }

run_byteident(){
  local name="$1"; shift
  local cdir="$WORK/$name"; local idir="$cdir/input"
  mkdir -p "$idir"
  perl "$GEN" --outdir "$idir" "$@" >/dev/null || { fail "$name (generator)"; return; }
  mkdir -p "$cdir/orig" "$cdir/proto"
  cp "$idir"/* "$cdir/orig"/ 2>/dev/null; cp "$idir"/* "$cdir/proto"/ 2>/dev/null
  ( cd "$cdir/orig"  && perl "$ORIG"  merged ); local ro=$?
  ( cd "$cdir/proto" && perl "$PROTO" merged ); local rp=$?
  if [ $ro -ne 0 ] || [ $rp -ne 0 ]; then fail "$name (nonzero exit orig=$ro proto=$rp)"; return; fi
  if cmp -s "$cdir/orig/merged.bedGraph" "$cdir/proto/merged.bedGraph"; then
    pass "$name  (.bedGraph byte-identical; $(grep -vc '^track' "$cdir/orig/merged.bedGraph") records)"
  else
    fail "$name  (.bedGraph differs)"; diff "$cdir/orig/merged.bedGraph" "$cdir/proto/merged.bedGraph" | head -6 | sed 's/^/      /'
  fi
  if diff -q <(mask_log "$cdir/orig/merged.log") <(mask_log "$cdir/proto/merged.log") >/dev/null; then
    pass "$name  (.log identical modulo timestamp)"
  else
    fail "$name  (.log differs)"
  fi
}

echo "=== byte-identity scenarios ==="
run_byteident int_mixed_base   --files 6  --mode int   --base mixed --seed 1
run_byteident int_all0based    --files 5  --mode int   --base 0     --seed 2
run_byteident onebased_only    --files 3  --mode int   --base 1     --seed 3
run_byteident single_1based    --files 1  --mode int   --base 1     --seed 4
run_byteident float_normalized --files 6  --mode float --base mixed --seed 5
run_byteident many_files_50    --files 50 --mode int   --base mixed --seed 6 --density 0.5
run_byteident sparse_union     --files 8  --mode int   --base mixed --seed 7 --density 0.25
run_byteident with_track_blank --files 4  --mode float --base mixed --seed 8 --track --blank
run_byteident wide_intervals   --files 4  --mode int   --base mixed --seed 9 --width 5
run_byteident chr_cmp_order    --files 4  --mode int   --base mixed --seed 10 --chroms chr1,chr2,chr10,chr11,chrX,chrY,chrM

echo
echo "=== misaligned input: BOTH must die (same end, different start) ==="
mdir="$WORK/misaligned"; mkdir -p "$mdir/orig" "$mdir/proto"
printf 'chr1\t10\t20\t3\n' > "$mdir/orig/a.bedGraph"; printf 'chr1\t12\t20\t4\n' > "$mdir/orig/b.bedGraph"
cp "$mdir/orig"/*.bedGraph "$mdir/proto"/
( cd "$mdir/orig"  && perl "$ORIG"  merged 2>/dev/null ); ro=$?
( cd "$mdir/proto" && perl "$PROTO" merged 2>/dev/null ); rp=$?
if [ $ro -ne 0 ] && [ $rp -ne 0 ]; then pass "misaligned: both exit nonzero (orig=$ro proto=$rp)"
else fail "misaligned: expected both nonzero, got orig=$ro proto=$rp"; fi
if [ ! -f "$mdir/proto/merged.bedGraph" ] && [ -z "$(ls "$mdir/proto"/merged.bedGraph.tmp.* 2>/dev/null)" ]; then
  pass "misaligned: atomic -> no partial .bedGraph and no leftover temp"
else fail "misaligned: partial output or leftover temp present"; fi

echo
echo "=== atomic: a good prior output is NOT clobbered on error ==="
adir="$WORK/atomic"; mkdir -p "$adir"
printf 'chr1\t10\t20\t3\n' > "$adir/a.bedGraph"; printf 'chr1\t12\t20\t4\n' > "$adir/b.bedGraph"
printf 'PRIOR-GOOD\n' > "$adir/merged.bedGraph"
( cd "$adir" && perl "$PROTO" merged 2>/dev/null )
if [ "$(cat "$adir/merged.bedGraph")" = "PRIOR-GOOD" ]; then pass "atomic: prior .bedGraph left intact on error"
else fail "atomic: prior .bedGraph was modified"; fi

echo
echo "=== unsorted-chr input: documented divergence (streaming guards) ==="
udir="$WORK/unsorted"; mkdir -p "$udir/input"
perl "$GEN" --outdir "$udir/input" --files 3 --mode int --base 0 --seed 11 --chroms chr1,chr2,chr10 --unsorted >/dev/null
mkdir -p "$udir/orig" "$udir/proto"; cp "$udir/input"/* "$udir/orig"/; cp "$udir/input"/* "$udir/proto"/
( cd "$udir/orig"  && perl "$ORIG"  merged 2>/dev/null ); ro=$?
( cd "$udir/proto" && perl "$PROTO" merged 2>/dev/null ); rp=$?
if [ $ro -eq 0 ] && [ $rp -ne 0 ]; then pass "unsorted-chr: original succeeds, streaming dies with guard (as designed)"
else fail "unsorted-chr: expected orig=0/proto!=0, got orig=$ro proto=$rp"; fi

echo
echo "=== memory + wall-time comparison on a large input ==="
sample_rss(){
  local t0 t1 pid max hwm rc
  t0=$(date +%s.%N); "$@" & pid=$!; max=0
  while kill -0 "$pid" 2>/dev/null; do
    hwm=$(awk '/VmHWM/{print $2}' "/proc/$pid/status" 2>/dev/null)
    [ -n "$hwm" ] && [ "$hwm" -gt "$max" ] && max=$hwm; sleep 0.02
  done
  wait "$pid"; rc=$?; t1=$(date +%s.%N)
  awk -v m="$max" -v a="$t0" -v b="$t1" 'BEGIN{printf "%d %.2f %d", m, b-a, '"$rc"'}'
}
big="$WORK/big"; mkdir -p "$big/input"
perl "$GEN" --outdir "$big/input" --files 8 --mode int --base mixed --seed 42 \
     --chroms chrBIG --grid 1200000 --spacing-max 1 --density 0.6 >/dev/null
mkdir -p "$big/orig" "$big/proto"; cp "$big/input"/* "$big/orig"/; cp "$big/input"/* "$big/proto"/
echo "  input: $(du -sh "$big/input" | cut -f1) across $(ls "$big/input" | wc -l) files"
read RO_RSS RO_WALL RO_RC < <(cd "$big/orig"  && sample_rss perl "$ORIG"  merged)
read RP_RSS RP_WALL RP_RC < <(cd "$big/proto" && sample_rss perl "$PROTO" merged)
if [ "$RO_RC" -eq 0 ] && [ "$RP_RC" -eq 0 ] && cmp -s "$big/orig/merged.bedGraph" "$big/proto/merged.bedGraph"; then
  pass "big input: byte-identical ($(grep -vc '^track' "$big/orig/merged.bedGraph") records)"
else fail "big input: differ or nonzero (orig_rc=$RO_RC proto_rc=$RP_RC)"; fi
printf "  %-10s %12s %10s\n" "impl" "maxRSS(MB)" "wall(s)"
awk -v r="$RO_RSS" -v w="$RO_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","original", r/1024, w}'
awk -v r="$RP_RSS" -v w="$RP_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","streaming", r/1024, w}'
awk -v o="$RO_RSS" -v p="$RP_RSS" 'BEGIN{ if(p>0) printf "  peak-RSS reduction: %.1fx\n", o/p }'

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
