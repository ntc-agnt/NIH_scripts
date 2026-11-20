#!/usr/bin/perl
# Synthetic data generator for the make_heatmap golden-master harness.
# Writes a self-consistent fixture set into an output directory: a stranded and
# an unstranded gene list (with overlapping spans, both strands, and a valid
# in-range integer in the description column so it doubles as a -a u anchor),
# and the same canonical hit set serialized into every supported hit format
# (G/g/b/e/c + stranded cppmatch), plus plus/minus split files for -p/-m, a bin
# definition file, and a "messy" file exercising bad/track/CRLF lines. Values
# are a mix of integers and non-round floats so float accumulation is exercised.
# Deterministic for a given --seed. perl + core only.
#
#   --outdir DIR   (required)   --seed S (default 1)
#   --chroms a,b,c (default chr1,chr2,chr3)
#   --ngenes N (default 60)     --nhits N (default 20000)
#   --span-max L (default 4000) --coord-max C (default 200000)
use strict; use warnings;

my %o=(outdir=>undef, seed=>1, chroms=>'chr1,chr2,chr3',
       ngenes=>60, nhits=>20000, 'span-max'=>4000, 'coord-max'=>200000);
while(@ARGV){ my $k=shift @ARGV; $k=~/^--(.+)/ or die "bad arg $k\n";
              die "missing value for --$1\n" unless @ARGV; $o{$1}=shift @ARGV; }
defined $o{outdir} or die "need --outdir DIR\n";
mkdir $o{outdir} unless -d $o{outdir};
srand($o{seed}+0);
my @chroms=split(/,/,$o{chroms});
my @strands=('+','-');

# ---- gene list -------------------------------------------------------------
# Overlapping spans on purpose: multiple genes per region, both strands, so the
# per-hit gene scan and strand splitting are exercised. desc = a valid integer
# anchor inside [start,end] (for -a u) AND human-readable enough for other modes.
open(my $gs,'>',"$o{outdir}/genes_stranded.txt")  or die $!;
open(my $gn,'>',"$o{outdir}/genes_nostrand.txt")  or die $!;
for my $i (1..$o{ngenes}) {
  my $chr=$chroms[int(rand(@chroms))];
  my $start=1+int(rand($o{'coord-max'}));
  my $span=50+int(rand($o{'span-max'}));
  my $end=$start+$span;
  my $anchor=$start+int(rand($end-$start+1));   # valid in-range integer anchor
  my $strand=$strands[int(rand(2))];
  print $gs join("\t","gene$i",$anchor,$chr,$start,$end,$strand),"\n";
  print $gn join("\t","gene$i",$anchor,$chr,$start,$end),"\n";
}
close $gs; close $gn;

# ---- canonical hit set -----------------------------------------------------
# Each record: chr, start0 (0-based), end (1-based, >start0), value, strand.
# Clustered near gene coordinates so plenty actually intersect, plus some spread.
my @hits;
for my $i (1..$o{nhits}) {
  my $chr=$chroms[int(rand(@chroms))];
  my $start0=int(rand($o{'coord-max'}+$o{'span-max'}));
  my $w=1+int(rand(200));                      # feature width
  my $end=$start0+$w;
  # value: ~half integers, ~half non-round floats (exercise float summation)
  my $value = (rand()<0.5) ? (1+int(rand(50))) : sprintf("%.4f", 0.0001+rand(25));
  my $strand=$strands[int(rand(2))];
  push @hits, [$chr,$start0,$end,$value,$strand];
}

# ---- serialize to every format --------------------------------------------
# G: standard bedGraph, 0-based start (tool does start++)
open(my $fG,'>',"$o{outdir}/hits.G") or die $!;
# g: bedGraph, 1-based start (tool uses as-is) -> write start0+1 for same locus
open(my $fg,'>',"$o{outdir}/hits.g") or die $!;
# b: bed3 (tool start++); value forced to 1 by tool
open(my $fb,'>',"$o{outdir}/hits.b") or die $!;
# e: bed6 with strand (tool start++)
open(my $fe,'>',"$o{outdir}/hits.e") or die $!;
# c: cppmatch unstranded: id value chr pstart pend  (physical, no ++)
open(my $fc,'>',"$o{outdir}/hits.c") or die $!;
# cs: cppmatch stranded: + trailing strand word
open(my $fcs,'>',"$o{outdir}/hits.cs") or die $!;
# plus/minus split in G format for -p/-m
open(my $fp,'>',"$o{outdir}/hits_plus.G")  or die $!;
open(my $fm,'>',"$o{outdir}/hits_minus.G") or die $!;
my $id=0;
for my $h (@hits) {
  my ($chr,$s0,$e,$v,$st)=@$h;
  $id++;
  my $sw = ($st eq '+') ? 'plus' : 'minus';
  print $fG  join("\t",$chr,$s0,$e,$v),"\n";
  print $fg  join("\t",$chr,$s0+1,$e,$v),"\n";
  print $fb  join("\t",$chr,$s0,$e),"\n";
  print $fe  join("\t",$chr,$s0,$e,"feat$id",$v,$st),"\n";
  print $fc  join("\t","h$id",$v,$chr,$s0+1,$e),"\n";
  print $fcs join("\t","h$id",$v,$chr,$s0+1,$e,$sw),"\n";
  if($st eq '+'){ print $fp join("\t",$chr,$s0,$e,$v),"\n"; }
  else          { print $fm join("\t",$chr,$s0,$e,$v),"\n"; }
}
close $_ for ($fG,$fg,$fb,$fe,$fc,$fcs,$fp,$fm);

# ---- bin definition file (negative, gaps, variable widths; non-overlapping) -
# NOTE: the original bin_parser's do/while loop treats a trailing newline as an
# empty final record and errors out ("Bin file contains formatting errors"), so
# the bin file must NOT end with a newline. This mirrors the real tool's contract
# (documented in README_optimized.md); we preserve it exactly rather than "fix" it.
open(my $bf,'>',"$o{outdir}/bins.txt") or die $!;
print $bf "-200\t-101\n-100\t-1\n0\t99\n100\t299\n400\t500\n501\t900";
close $bf;

# ---- "messy" G file: bad lines, a track line, CRLF, blank lines ------------
# Exercises the bad-line recursion, "track" suppression, and \r handling. The
# skip messages go to stdout, not the matrix file, so byte-diff of output is
# unaffected; both ref and candidate must handle them identically.
open(my $mz,'>',"$o{outdir}/hits_messy.G") or die $!;
print $mz "track type=bedGraph name=\"x\"\n";
print $mz "$chroms[0]\t100\t200\t5\n";
print $mz "\n";
print $mz "garbage line with no numbers\n";
print $mz "$chroms[0]\t300\t400\t2.5\r\n";     # CRLF
print $mz "$chroms[0]\t50\t60\t7\n";
print $mz "$chroms[0]\tNOTANUM\t60\t7\n";       # bad coord
close $mz;

print "wrote fixtures to $o{outdir} ($o{ngenes} genes, $o{nhits} hits, ",
      scalar(@chroms)," chroms)\n";
