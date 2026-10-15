#!/usr/bin/perl
# Synthetic native-bowtie-format alignment generator for the bowtie2stdBedGraph
# refactor benchmark. Emits lines: name<TAB>strand<TAB>chr<TAB>pos0<TAB>seq<TAB>qual
# (the trailing field provides the tab the parser's regex requires).
#
# Reads sample from a fixed pool of genomic "sites" so that many reads pile onto
# the same positions (realistic deep coverage); the number of distinct sites
# drives the size of the in-memory count hashes (the thing the refactor targets).
#
# Usage: gen_bowtie.pl --outfile F --reads N --sites P --chroms c1,c2,.. \
#                      --readlen L --seed S
use strict; use warnings;
my %o=(outfile=>undef,reads=>1_000_000,sites=>200_000,
       chroms=>'chr1,chr2,chr3,chr4,chr5,chrX',readlen=>36,seed=>1,maxpos=>200_000_000);
while(@ARGV){my $k=shift @ARGV; $k=~/^--(.+)/ or die "bad arg $k\n"; $o{$1}=shift @ARGV;}
defined $o{outfile} or die "need --outfile\n";
srand($o{seed}+0);
my @chroms=split/,/,$o{chroms};
# build the site pool
my @site_chr; my @site_pos;
for(my $i=0;$i<$o{sites};$i++){
  push @site_chr, $chroms[int(rand(@chroms))];
  push @site_pos, 1+int(rand($o{maxpos}));
}
my $seq ="A" x $o{readlen};
my $qual="I" x $o{readlen};
open(my $fh,'>',$o{outfile}) or die "cannot write $o{outfile}: $!\n";
for(my $r=0;$r<$o{reads};$r++){
  my $s=int(rand($o{sites}));
  my $strand = (rand()<0.5) ? '+' : '-';
  print $fh "r$r\t$strand\t$site_chr[$s]\t$site_pos[$s]\t$seq\t$qual\n";
}
close $fh;
print "wrote $o{reads} reads over $o{sites} sites -> $o{outfile}\n";
