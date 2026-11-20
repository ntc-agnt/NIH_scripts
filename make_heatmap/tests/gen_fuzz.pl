#!/usr/bin/env perl
# Adversarial line generator for the parser-equivalence fuzz. Emits a mix of
# plausible-but-varied hit lines and deliberately malformed ones, to stress the
# manual parser's operator>>-equivalence: wrong field counts, non-numeric values
# in numeric slots, floats where ints are expected, signs, scientific notation,
# space-vs-tab and multi-whitespace separators, leading/trailing whitespace,
# whitespace-only lines, "track" lines, and CRLF. Magnitudes are capped well
# below overflow so the (documented) ERANGE-on-double corner is never hit.
# Deterministic for a given --seed. The harness feeds the SAME file to the
# reference and candidate under several parsers and byte-diffs the output.
use strict; use warnings;
my %o=(out=>undef, lines=>20000, seed=>1, chroms=>'chr1,chr2,chr3');
while(@ARGV){ my $k=shift @ARGV; $k=~/^--(.+)/ or die "bad arg $k\n"; $o{$1}=shift @ARGV; }
defined $o{out} or die "need --out\n";
srand($o{seed}+0);
my @chroms=split(/,/,$o{chroms});

sub ws { my @w=("\t","  ","\t\t"," ","\t "); $w[int(rand(@w))] }   # separator variety
sub an_int { my $v=int(rand(200000)); (rand()<0.15?'-':'').(rand()<0.15?'+':'').$v }
sub a_float { my @f=(sprintf("%.3f",rand(50)), sprintf("%d",1+int(rand(40))),
                     sprintf("%.2e",rand(10)), sprintf("-%.3f",rand(5)),
                     sprintf("%.4f",rand(1))); $f[int(rand(@f))] }
sub junk { my @j=("abc","12x","x12","1.2.3","","NaNN","1e","..","+-3","0x1F","--",
                  "1,000","3.","e5","inf_","chrom"); $j[int(rand(@j))] }
sub field { my $r=rand();
  return an_int()  if $r<0.40;
  return a_float() if $r<0.65;
  return $chroms[int(rand(@chroms))] if $r<0.80;
  return junk();
}

open(my $fh,'>',$o{out}) or die $!;
for my $i (1..$o{lines}) {
  my $r=rand();
  if($r<0.04){ print $fh "track name=\"x$i\" type=bedGraph\n"; next; }
  if($r<0.07){ print $fh "\n"; next; }                       # blank (reader skips)
  if($r<0.10){ print $fh "   \t \n"; next; }                 # whitespace-only
  # choose a plausible-ish layout half the time, random-arity the other half
  my @cols;
  if(rand()<0.5){
    # bedgraph/bed-ish: chr, int, int, [value], [strand]
    @cols=($chroms[int(rand(@chroms))], an_int(), an_int());
    push @cols, a_float() if rand()<0.8;
    push @cols, (("+","-","plus","minus","foo")[int(rand(5))]) if rand()<0.5;
  } else {
    my $n=int(rand(7));                                       # 0..6 arbitrary fields
    @cols=map { field() } 1..$n;
  }
  my $line=join(ws(), @cols);
  $line=(" " x int(rand(3))).$line;                          # leading whitespace sometimes
  my $eol=(rand()<0.1)?"\r\n":"\n";                          # CRLF sometimes
  print $fh $line.$eol;
}
close $fh;
print "wrote $o{lines} fuzz lines to $o{out}\n";
