#!/usr/bin/env bash
# Byte-identity + throughput benchmark for the trim_and_filter optimized variants.
# Verifies trim_and_filter_PE_optimized.pl and trim_and_filter_SE_optimized.pl
# are byte-identical to the originals across all modes, then compares throughput
# (original vs optimized, and vs cutadapt if installed).
# Self-contained: needs only perl + coreutils (+ cutadapt, optional).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
PE_O="$HERE/../trim_and_filter_PE.pl";           PE_P="$HERE/../trim_and_filter_PE_optimized.pl"
SE_O="$HERE/../../trim_and_filter_SE/trim_and_filter_SE.pl"
SE_P="$HERE/../../trim_and_filter_SE/trim_and_filter_SE_optimized.pl"
GEN="$HERE/gen_fastq.pl"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
PASS=0; FAIL=0; pass(){ echo "  PASS: $1"; PASS=$((PASS+1)); }; fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }
wall(){ local t0 t1; t0=$(date +%s.%N); "$@" >/dev/null 2>&1; t1=$(date +%s.%N); awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f",b-a}'; }

echo "=== byte-identity: optimized vs original ==="
perl "$GEN" --out1 "$W/s33_1.fastq" --out2 "$W/s33_2.fastq" --reads 40000 --readlen 50 --qbase 33 --seed 2 >/dev/null
perl "$GEN" --out1 "$W/s64_1.fastq" --out2 "$W/s64_2.fastq" --reads 40000 --readlen 50 --qbase 64 --seed 3 >/dev/null
perl "$GEN" --out  "$W/e33.fastq"   --reads 40000 --readlen 50 --qbase 33 --seed 4 >/dev/null
perl "$GEN" --out  "$W/e64.fastq"   --reads 40000 --readlen 50 --qbase 64 --seed 5 >/dev/null
gzip -kf "$W/s33_1.fastq" "$W/s33_2.fastq" "$W/e33.fastq"
pe_case(){ local tag=$1 in1=$2 in2=$3 qs=$4 a=$5 b=$6 c=$7 d=$8 m=$9; mkdir -p "$W/$tag/o" "$W/$tag/p"
  ( cd "$W/$tag/o" && perl "$PE_O" -1 "$in1" -2 "$in2" -a $a -b $b -c $c -d $d -m $m -q $qs -o s ) 2>/dev/null
  ( cd "$W/$tag/p" && perl "$PE_P" -1 "$in1" -2 "$in2" -a $a -b $b -c $c -d $d -m $m -q $qs -o s ) 2>/dev/null
  local d1=0; for f in "$W/$tag/o"/s.*; do cmp -s "$f" "$W/$tag/p/$(basename "$f")" || d1=1; done
  [ $d1 -eq 0 ] && pass "PE $tag ($(ls "$W/$tag/o"/s.*|wc -l) files identical)" || fail "PE $tag"; }
se_case(){ local tag=$1 in=$2 qs=$3 a=$4 b=$5 m=$6; mkdir -p "$W/$tag/o" "$W/$tag/p"
  ( cd "$W/$tag/o" && perl "$SE_O" -i "$in" -a $a -b $b -m $m -q $qs -o s ) 2>/dev/null
  ( cd "$W/$tag/p" && perl "$SE_P" -i "$in" -a $a -b $b -m $m -q $qs -o s ) 2>/dev/null
  local d1=0; for f in "$W/$tag/o"/s.*; do cmp -s "$f" "$W/$tag/p/$(basename "$f")" || d1=1; done
  [ $d1 -eq 0 ] && pass "SE $tag ($(ls "$W/$tag/o"/s.*|wc -l) files identical)" || fail "SE $tag"; }
pe_case pe_full   "$W/s33_1.fastq"    "$W/s33_2.fastq"    sanger   1 50 1 50 25
pe_case pe_window "$W/s33_1.fastq"    "$W/s33_2.fastq"    sanger   5 45 3 48 30
pe_case pe_ilmn   "$W/s64_1.fastq"    "$W/s64_2.fastq"    illumina 5 45 5 45 25
pe_case pe_gz     "$W/s33_1.fastq.gz" "$W/s33_2.fastq.gz" sanger   5 45 5 45 20
se_case se_full   "$W/e33.fastq"    sanger   1 50 25
se_case se_window "$W/e33.fastq"    sanger   5 45 30
se_case se_ilmn   "$W/e64.fastq"    illumina 5 45 25
se_case se_gz     "$W/e33.fastq.gz" sanger   5 45 20

echo
echo "=== throughput: PE, plain fastq, 2M read pairs (single-thread) ==="
perl "$GEN" --out1 "$W/b1.fastq" --out2 "$W/b2.fastq" --reads 2000000 --readlen 50 --qbase 33 --seed 9 >/dev/null
A=5; B=45; M=20; L=50
mkdir -p "$W/to" "$W/tp"
TO=$(cd "$W/to" && wall perl "$PE_O" -1 "$W/b1.fastq" -2 "$W/b2.fastq" -a $A -b $B -c $A -d $B -m $M -q sanger -o s)
TP=$(cd "$W/tp" && wall perl "$PE_P" -1 "$W/b1.fastq" -2 "$W/b2.fastq" -a $A -b $B -c $A -d $B -m $M -q sanger -o s)
bi=$(cmp -s "$W/to/s.1.trim_${A}_${B}.minQS_${M}.fastq" "$W/tp/s.1.trim_${A}_${B}.minQS_${M}.fastq" && echo YES || echo NO)
printf "  %-28s %8s\n" tool "wall(s)"
printf "  %-28s %8s\n" "trim_and_filter_PE (orig)" "$TO"
printf "  %-28s %8s\n" "trim_and_filter_PE_optimized" "$TP"
awk -v o="$TO" -v p="$TP" 'BEGIN{printf "  optimized speedup: %.2fx (byte-identical output: '"$bi"')\n", o/p}'
if command -v cutadapt >/dev/null 2>&1; then
  FRONT=$((A-1)); TAIL=$((L-B)); AER=$(awk -v m=$M 'BEGIN{printf "%.6f",10^(-m/10)}')
  mkdir -p "$W/tc"
  TC=$(cd "$W/tc" && wall cutadapt -j 1 -u $FRONT -u -$TAIL -U $FRONT -U -$TAIL --max-aer $AER --quality-base 33 -o c1.fastq -p c2.fastq "$W/b1.fastq" "$W/b2.fastq")
  printf "  %-28s %8s\n" "cutadapt -j1 (best-match*)" "$TC"
  echo "  * cutadapt --max-average-error-rate is NOT the same metric as this tool's"
  echo "    mean-Phred filter, so kept-read counts differ; throughput reference only."
else
  echo "  (cutadapt not installed; skipping the cutadapt reference)"
fi

echo
echo "=== SUMMARY: $PASS passed, $FAIL failed ==="
[ $FAIL -eq 0 ]
