#!/usr/bin/perl
# Synthetic bedGraph generator for the bedgraphs2stdBedGraph byte-diff harness.
#
# Produces N input files in --outdir, all describing a shared, sorted,
# non-overlapping coordinate grid per chromosome so that same-interval records
# across files ALIGN (same 0-based start for a given end) and therefore sum.
#
# Each file carries a random SUBSET of the grid (to exercise the union/merge),
# with random values. Files are emitted sorted by (chr in cmp order, start asc)
# -- the precondition the streaming merge relies on and that the ecosystem
# producer (bowtie2stdBedGraph) already satisfies.
#
# A 0-based interval [s,e] is written:
#   .bedGraph  (0-based, offset 0):  "chr \t s   \t e \t val"
#   .bedgraph  (1-based, offset 1):  "chr \t s+1 \t e \t val"   (tool subtracts 1)
# so both map to the same canonical 0-based (s,e).
#
# Usage:
#   gen_bedgraphs.pl --outdir DIR --files N --chroms c1,c2,... --grid G
#       --spacing-max M --density D --mode int|float --base 0|1|mixed
#       --seed S [--track] [--blank]
use strict; use warnings;

my %o = (outdir=>undef, files=>4, chroms=>'chr1,chr2,chr10,chrX,chrM',
         grid=>200, 'spacing-max'=>6, density=>0.7, mode=>'int',
         base=>'mixed', seed=>1, width=>1, track=>0, blank=>0, unsorted=>0);
while(@ARGV){ my $k=shift @ARGV;
  if($k=~/^--(.+)/){ my $n=$1;
    if($n eq 'track'||$n eq 'blank'||$n eq 'unsorted'){ $o{$n}=1; }
    else { $o{$n}=shift @ARGV; }
  } else { die "bad arg $k\n"; }
}
defined $o{outdir} or die "need --outdir\n";
srand($o{seed}+0);
mkdir $o{outdir} unless -d $o{outdir};
my @chroms = split /,/, $o{chroms};

# Build a shared grid of non-overlapping 0-based intervals per chromosome.
my %grid;   # chr -> [ [s,e], ... ] ascending
for my $c (@chroms){
  my @g; my $pos = 1 + int(rand(5));
  for(my $i=0;$i<$o{grid};$i++){
    my $s = $pos;
    my $e = $s + $o{width};
    push @g, [$s,$e];
    $pos = $e + 1 + int(rand($o{'spacing-max'}));  # strictly increasing, non-overlapping
  }
  $grid{$c} = \@g;
}

sub rand_val {
  if($o{mode} eq 'int'){ return 1 + int(rand(50)); }
  # float: a few decimal places, the normalized-bedGraph case
  return sprintf("%.6f", rand()*3.0 + 0.0001) + 0;
}

for(my $f=0; $f<$o{files}; $f++){
  # Decide this file's coordinate base / suffix.
  my $is1 = ($o{base} eq '1') ? 1 : ($o{base} eq '0') ? 0 : ($f % 2);
  my $suffix = $is1 ? 'bedgraph' : 'bedGraph';
  my $name = sprintf("sample_%02d.%s", $f, $suffix);
  my $path = "$o{outdir}/$name";
  open(my $fh,'>',$path) or die "cannot write $path: $!\n";
  print $fh "track type=bedGraph name=\"$name\"\n\n" if $o{track};
  # Real ecosystem inputs are sorted by (chr cmp, start). Emit that way unless
  # --unsorted is requested (to exercise the prototype's sortedness guard).
  my @write_chroms = $o{unsorted} ? @chroms : sort {$a cmp $b} @chroms;
  for my $c (@write_chroms){
    my $blanked=0;
    for my $iv (@{$grid{$c}}){
      next if rand() > $o{density};           # random subset -> union across files
      my ($s,$e) = @$iv;
      my $start_field = $is1 ? $s+1 : $s;      # 1-based start stores s+1
      my $val = rand_val();
      print $fh "$c\t$start_field\t$e\t$val\n";
      if($o{blank} && !$blanked && rand()<0.1){ print $fh "\n"; $blanked=1; }
    }
  }
  close $fh;
}
print "generated $o{files} files in $o{outdir}\n";
