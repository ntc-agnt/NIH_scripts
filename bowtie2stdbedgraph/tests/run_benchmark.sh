#!/usr/bin/env bash
# Byte-identity + mem/time benchmark for bowtie2stdBedGraph_optimized.pl vs the
# original bowtie2stdBedGraph.pl. Self-contained: needs only perl + coreutils,
# no cluster or real data. (DBI is required by both scripts only for the -r UCSC
# lookup; to stay runnable where DBI is not installed, this harness runs copies
# with `use DBI` swapped for a lazy require. That does not affect any output.)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ORIG_SRC="$HERE/../bowtie2stdBedGraph.pl"
OPT_SRC="$HERE/../bowtie2stdBedGraph_optimized.pl"
GEN="$HERE/gen_bowtie.pl"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
LAZY='BEGIN { eval { require DBI; DBI->import; }; }'
perl -pe "s/^use DBI;.*/$LAZY/" "$ORIG_SRC" > "$W/orig.pl"
perl -pe "s/^use DBI;.*/$LAZY/" "$OPT_SRC"  > "$W/opt.pl"
ORIG="$W/orig.pl"; OPT="$W/opt.pl"
PASS=0; FAIL=0; pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }; fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

echo "=== byte-identity: optimized vs original, across option modes ==="
perl "$GEN" --outfile "$W/small.bowtie" --reads 400000 --sites 80000 --seed 3 >/dev/null
printf 'chr1\t200000000\nchr2\t200000000\nchr3\t200000000\nchr4\t200000000\nchr5\t200000000\nchrX\t200000000\n' > "$W/sizes.txt"
mode(){
  local tag="$1"; shift
  mkdir -p "$W/$tag/o" "$W/$tag/p"
  ( cd "$W/$tag/o" && perl "$ORIG" "$@" "$W/small.bowtie" out ) 2>/dev/null
  ( cd "$W/$tag/p" && perl "$OPT"  "$@" "$W/small.bowtie" out ) 2>/dev/null
  local d=0 n; for f in "$W/$tag/o"/*.bedGraph; do cmp -s "$f" "$W/$tag/p/$(basename "$f")" || d=1; done
  n=$(ls "$W/$tag/o"/*.bedGraph 2>/dev/null | wc -l)
  if [ "$d" -eq 0 ]; then pass "$tag ($n files byte-identical)"; else fail "$tag (differs)"; fi
}
mode default
mode norm     -o n
mode binned   -o b -b 50
mode merged   -o m -s 75
mode all      -o nbmv -b 50 -s 75
mode swap     -x -o n
mode withlist -o mv -s 75 -b 50 -l "$W/sizes.txt"

echo
echo "=== peak memory + wall-time on a larger input (-o nbmv) ==="
sample(){ local t0 t1 pid max tot p rss rc
  t0=$(date +%s.%N); "$@" & pid=$!; max=0
  while kill -0 "$pid" 2>/dev/null; do
    tot=0; for p in $(pgrep -P "$pid" 2>/dev/null) "$pid"; do rss=$(awk '/VmRSS/{print $2}' "/proc/$p/status" 2>/dev/null); tot=$((tot+${rss:-0})); done
    [ "$tot" -gt "$max" ] && max=$tot; sleep 0.03
  done
  wait "$pid"; rc=$?; t1=$(date +%s.%N)
  awk -v m="$max" -v a="$t0" -v b="$t1" 'BEGIN{printf "%d %.2f %d", m, b-a, '"$rc"'}'; }
perl "$GEN" --outfile "$W/big.bowtie" --reads 5000000 --sites 1500000 --seed 11 >/dev/null
echo "  input: $(du -sh "$W/big.bowtie" | cut -f1), 5M reads / 1.5M sites"
mkdir -p "$W/bo" "$W/bp"
read OR OW OC < <(cd "$W/bo" && sample perl "$ORIG" -o nbmv -b 50 -s 75 "$W/big.bowtie" out)
read PR PW PC < <(cd "$W/bp" && sample perl "$OPT"  -o nbmv -b 50 -s 75 "$W/big.bowtie" out)
bi=yes; for f in "$W/bo"/*.bedGraph; do cmp -s "$f" "$W/bp/$(basename "$f")" || bi=NO; done
echo "  byte-identical: $bi  (exit orig=$OC opt=$PC)"
printf "  %-12s %12s %10s\n" impl "peakRSS(MB)" "wall(s)"
awk -v r="$OR" -v w="$OW" 'BEGIN{printf "  %-12s %12.1f %10s\n","original", r/1024, w}'
awk -v r="$PR" -v w="$PW" 'BEGIN{printf "  %-12s %12.1f %10s\n","optimized", r/1024, w}'
awk -v o="$OR" -v p="$PR" 'BEGIN{ if(p>0) printf "  peak-RSS: optimized is %.2fx the original\n", p/o }'
awk -v o="$OW" -v p="$PW" 'BEGIN{ if(p>0) printf "  wall:     optimized is %.2fx the original\n", p/o }'

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
