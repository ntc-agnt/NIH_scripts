#!/usr/bin/env bash
# Validation + benchmark harness for extract_fragments_optimized.pl.
#
# Contract (see README_optimized.md):
#   * bedGraph output  -> BYTE-IDENTICAL to extract_fragments.pl.
#   * BED output       -> SAME SET of lines as the original, and DETERMINISTIC
#                         (byte-stable across repeated optimized runs).
# The harness exercises a mode matrix (-o b|g|a, varied -min/-max/-b, -u 0/1/N,
# -l length list, -D/-M chr handling, bowtie vs SAM input, a same-start/
# different-end tie fixture, and empty/all-filtered inputs), then reports a
# peak-RSS + wall-time comparison on a large input. Self-contained: perl +
# coreutils, no cluster or real data (DBI/-r is not exercised -- needs network).
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
ORIG="$HERE/../extract_fragments.pl"
OPT="$HERE/../extract_fragments_optimized.pl"
GEN="$HERE/gen_fragments.pl"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

# The original extract_fragments.pl does `use DBI` unconditionally (DBI is only
# actually used on the -r UCSC-fetch path, which this harness never exercises).
# Provide a no-op stub DBI so the untouched original compiles on hosts without
# the real DBI module. If a real DBI is installed it takes precedence (PERL5LIB
# is searched after built-in paths for `use`-time resolution is unaffected here,
# and the stub is harmless since -r is never used).
mkdir -p "$WORK/stublib"
printf 'package DBI;\n1;\n' > "$WORK/stublib/DBI.pm"
export PERL5LIB="$WORK/stublib${PERL5LIB:+:$PERL5LIB}"

PASS=0; FAIL=0
pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

# Run both scripts with the same flags on the same input and check the contract
# for whichever outputs the mode produces.
#   $1 name   $2 format(bowtie|sam)   $3.. extra flags + gen options split by '::'
check(){
  local name="$1" fmt="$2" tool_flags="$3" gen_opts="$4"
  local cdir="$WORK/$name"; mkdir -p "$cdir"
  local inp="$cdir/input.txt"
  # shellcheck disable=SC2086
  perl "$GEN" --format "$fmt" --out "$inp" $gen_opts >/dev/null 2>&1 \
      || { fail "$name (generator)"; return; }

  # Each arm runs in its own dir with the SAME output prefix ("out"), so the
  # track-header line (which embeds the prefix) is identical too -> a clean
  # byte comparison of the whole file, header included.
  mkdir -p "$cdir/orig" "$cdir/opt1" "$cdir/opt2"
  # shellcheck disable=SC2086
  ( cd "$cdir/orig" && perl "$ORIG" $tool_flags "$inp" out >stdout 2>&1 ); local ro=$?
  # shellcheck disable=SC2086
  ( cd "$cdir/opt1" && perl "$OPT"  $tool_flags "$inp" out >stdout 2>&1 ); local r1=$?
  # shellcheck disable=SC2086
  ( cd "$cdir/opt2" && perl "$OPT"  $tool_flags "$inp" out >stdout 2>&1 ); local r2=$?
  if [ $ro -ne 0 ] || [ $r1 -ne 0 ] || [ $r2 -ne 0 ]; then
    fail "$name (nonzero exit orig=$ro opt=$r1/$r2)"; return; fi

  # bedGraph: byte-identical (produced for -o g and -o a).
  if [ -f "$cdir/orig/out.bedGraph" ]; then
    if cmp -s "$cdir/orig/out.bedGraph" "$cdir/opt1/out.bedGraph"; then
      pass "$name  bedGraph byte-identical ($(grep -vc '^track' "$cdir/orig/out.bedGraph") bins)"
    else
      fail "$name  bedGraph differs"; diff "$cdir/orig/out.bedGraph" "$cdir/opt1/out.bedGraph" | head -4 | sed 's/^/      /'
    fi
  fi

  # BED: same set of lines + optimized byte-stable (produced for -o b and -o a).
  if [ -f "$cdir/orig/out.bed" ]; then
    if cmp -s <(LC_ALL=C sort "$cdir/orig/out.bed") <(LC_ALL=C sort "$cdir/opt1/out.bed"); then
      pass "$name  BED line-set == original ($(grep -vc '^track' "$cdir/orig/out.bed") lines)"
    else
      fail "$name  BED line-set differs from original"
      diff <(LC_ALL=C sort "$cdir/orig/out.bed") <(LC_ALL=C sort "$cdir/opt1/out.bed") | head -4 | sed 's/^/      /'
    fi
    if cmp -s "$cdir/opt1/out.bed" "$cdir/opt2/out.bed"; then
      pass "$name  BED deterministic (two optimized runs byte-identical)"
    else
      fail "$name  BED NOT deterministic across optimized runs"
    fi
  fi
}

echo "=== mode matrix (bowtie input) ==="
check bowtie_both_default  bowtie "-o a"                       "--pairs 3000 --seed 1"
check bowtie_bed_only      bowtie "-o b"                       "--pairs 3000 --seed 2"
check bowtie_graph_only    bowtie "-o g"                       "--pairs 3000 --seed 3"
check bowtie_minmax        bowtie "-o a -min 80 -max 150"      "--pairs 3000 --seed 4 --minlen 40 --maxlen 260"
check bowtie_binsize       bowtie "-o a -b 50"                 "--pairs 3000 --seed 5"
check bowtie_wideopen      bowtie "-o a -min 1 -max 100000"    "--pairs 3000 --seed 6"

echo
echo "=== dedup (-u) matrix: heavy duplicate + tie content ==="
check dedup_u0             bowtie "-o a -u 0"                  "--pairs 4000 --seed 7 --dup-frac 0.45 --tie-frac 0.25"
check dedup_u1             bowtie "-o a -u 1"                  "--pairs 4000 --seed 7 --dup-frac 0.45 --tie-frac 0.25"
check dedup_u3             bowtie "-o a -u 3"                  "--pairs 4000 --seed 7 --dup-frac 0.45 --tie-frac 0.25"
check dedup_u1_bedonly     bowtie "-o b -u 1"                  "--pairs 4000 --seed 8 --dup-frac 0.45 --tie-frac 0.25"

echo
echo "=== tie-heavy fixture (same start, different end -> determinism) ==="
check tie_heavy            bowtie "-o a"                       "--pairs 5000 --seed 9 --tie-frac 0.6 --dup-frac 0.1 --grid 500"

echo
echo "=== SAM input (-a), incl. line-order swap ==="
check sam_both             sam    "-a -o a"                    "--pairs 3000 --seed 11"
check sam_bed_only         sam    "-a -o b"                    "--pairs 3000 --seed 12"
check sam_dedup_u2         sam    "-a -o a -u 2"               "--pairs 4000 --seed 13 --dup-frac 0.4 --tie-frac 0.3"

echo
echo "=== chr handling: bare names (+chr prepend), -D disable, -M mito, -l length list ==="
check chr_prepend          bowtie "-o a"                       "--pairs 3000 --seed 14 --chroms 1,2,10,MT"
check chr_Dflag            bowtie "-o a -D"                    "--pairs 3000 --seed 14 --chroms 1,2,10,MT"
check chr_mito_M           bowtie "-o a -M MT"                 "--pairs 3000 --seed 15 --chroms 1,2,MT"
# length list: cap the final bin at chrom end (affects bedGraph). Build a list file.
ldir="$WORK/lenlist"; mkdir -p "$ldir"; printf 'chr1\t400\nchr2\t400\n' > "$ldir/sizes.txt"
check length_list          bowtie "-o g -b 50 -l $ldir/sizes.txt" "--pairs 3000 --seed 16 --chroms chr1,chr2 --grid 380"

echo
echo "=== edge cases: all-filtered (empty output) and single-chrom ==="
check all_filtered         bowtie "-o a -min 100000 -max 100001" "--pairs 1500 --seed 17"
check single_chrom         bowtie "-o a"                       "--pairs 2000 --seed 18 --chroms chr1 --grid 300"

echo
echo "=== peak-RSS + wall-time comparison on a large input (-o a, all distinct & in range) ==="
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
big="$WORK/big"; mkdir -p "$big/orig" "$big/opt"
# large grid + no dups/ties + fragment range fully inside the tool's [min,max]
# => every fragment distinct & retained => the %entries / BED path is stressed.
perl "$GEN" --format bowtie --out "$big/in.txt" --pairs 400000 --seed 42 \
     --chroms chr1,chr2,chr3,chr10,chrX --grid 50000000 \
     --minlen 110 --maxlen 190 --dup-frac 0 --tie-frac 0 >/dev/null
echo "  input: $(du -sh "$big/in.txt" | cut -f1), 400000 pairs"
read RO_RSS RO_WALL RO_RC < <(cd "$big/orig" && sample_rss perl "$ORIG" -o a -min 100 -max 200 ../in.txt out)
read RP_RSS RP_WALL RP_RC < <(cd "$big/opt"  && sample_rss perl "$OPT"  -o a -min 100 -max 200 ../in.txt out)
if [ "$RO_RC" -eq 0 ] && [ "$RP_RC" -eq 0 ]; then
  if cmp -s "$big/orig/out.bedGraph" "$big/opt/out.bedGraph" \
     && cmp -s <(LC_ALL=C sort "$big/orig/out.bed") <(LC_ALL=C sort "$big/opt/out.bed"); then
    pass "big input: bedGraph byte-identical + BED set-equal ($(grep -vc '^track' "$big/orig/out.bed") BED lines)"
  else
    fail "big input: output mismatch"
  fi
else
  fail "big input: nonzero exit (orig=$RO_RC opt=$RP_RC)"
fi
printf "  %-10s %12s %10s\n" "impl" "maxRSS(MB)" "wall(s)"
awk -v r="$RO_RSS" -v w="$RO_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","original", r/1024, w}'
awk -v r="$RP_RSS" -v w="$RP_WALL" 'BEGIN{printf "  %-10s %12.1f %10s\n","optimized", r/1024, w}'
awk -v o="$RO_RSS" -v p="$RP_RSS" 'BEGIN{ if(p>0) printf "  peak-RSS reduction: %.2fx\n", o/p }'
awk -v o="$RO_WALL" -v p="$RP_WALL" 'BEGIN{ if(p>0) printf "  wall-time speedup:  %.2fx\n", o/p }'

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
