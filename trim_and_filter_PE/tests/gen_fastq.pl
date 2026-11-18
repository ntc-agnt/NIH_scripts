#!/usr/bin/perl
# Synthetic FASTQ generator for the trim_and_filter benchmark.
# Emits reads with a spread of per-read mean quality so the avg-quality filter
# passes some and fails others. Sanger (Phred+33) by default; --qbase 64 for
# Illumina-scale encoding.
#   PE: --out1 F1 --out2 F2 ;  SE: --out F
use strict; use warnings;
my %o=(out=>undef,out1=>undef,out2=>undef,reads=>1_000_000,readlen=>50,qbase=>33,seed=>1);
while(@ARGV){ my $k=shift @ARGV; $k=~/^--(.+)/ or die "bad arg $k\n"; $o{$1}=shift @ARGV; }
srand($o{seed}+0);
my @B=qw(A C G T);
sub gauss { my($m,$s)=@_; my($u1,$u2)=(rand()||1e-9,rand()); $m + $s*sqrt(-2*log($u1))*cos(2*3.14159265358979*$u2); }
sub one_read {
  my ($fh,$id,$mult)=@_;
  my $L=$o{readlen};
  my $seq=""; $seq.=$B[int(rand(4))] for 1..$L;
  my $mean=8+rand(34);                      # per-read mean Phred in ~[8,42)
  my $q="";
  for(1..$L){ my $ph=int(gauss($mean,4)+0.5); $ph=2 if $ph<2; $ph=41 if $ph>41; $q.=chr($o{qbase}+$ph); }
  print $fh "\@read${id}${mult}\n$seq\n+\n$q\n";
}
if(defined $o{out1} && defined $o{out2}){
  open(my $f1,'>',$o{out1}) or die $!; open(my $f2,'>',$o{out2}) or die $!;
  for my $i (1..$o{reads}){ one_read($f1,$i,"/1"); one_read($f2,$i,"/2"); }
  close $f1; close $f2;
  print "wrote $o{reads} read pairs\n";
} elsif(defined $o{out}){
  open(my $f,'>',$o{out}) or die $!;
  one_read($f,$_,"") for 1..$o{reads};
  close $f;
  print "wrote $o{reads} reads\n";
} else { die "need --out (SE) or --out1/--out2 (PE)\n"; }
