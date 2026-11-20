#!/usr/bin/perl
# Synthetic paired-end alignment generator for the extract_fragments harness.
# Emits bowtie-style or SAM-style PE records (query-name adjacent, mate1 leftmost
# -- the layout extract_fragments assumes) with a controllable spread of fragment
# lengths (some inside [min,max], some outside), deliberate same-start/different-
# end ties (to exercise BED sort determinism), and exact-duplicate fragments (to
# exercise -u dedup). Deterministic for a given --seed. perl + core only.
#
#   --format bowtie|sam   output style (default bowtie)
#   --out FILE            output path (required)
#   --pairs N            number of read pairs (default 2000)
#   --readlen RL         read length of the seq-bearing mate (default 36)
#   --chroms a,b,c       chromosome names (default chr1,chr2,chr10)
#   --grid G             max leftmost start coordinate, 1-based (default 100000)
#   --minlen / --maxlen  fragment-length sampling range (default 50..300; the
#                        tool's own [min,max] filter is separate)
#   --tie-frac F         fraction reusing a prior (chr,start) w/ new end (default 0.15)
#   --dup-frac F         fraction repeating a prior exact (chr,start,end) (default 0.15)
#   --seed S             RNG seed (default 1)
use strict; use warnings;

my %o=(format=>'bowtie', out=>undef, pairs=>2000, readlen=>36,
       chroms=>'chr1,chr2,chr10', grid=>100000, minlen=>50, maxlen=>300,
       'tie-frac'=>0.15, 'dup-frac'=>0.15, seed=>1);
while(@ARGV){ my $k=shift @ARGV; $k=~/^--(.+)/ or die "bad arg $k\n";
              die "missing value for --$1\n" unless @ARGV; $o{$1}=shift @ARGV; }
defined $o{out} or die "need --out FILE\n";
srand($o{seed}+0);
my @chroms=split(/,/,$o{chroms});
my $rl=$o{readlen};
my @B=qw(A C G T);
sub rseq { my $n=shift; my $s=''; $s.=$B[int(rand(4))] for 1..$n; $s }
sub rqual { my $n=shift; 'I' x $n }   # constant qual; unused by the tool

# Build a list of fragments first (so ties/dups can reference earlier picks),
# then emit the reads. Each fragment: [chr, start(1-based leftmost), fraglen].
my @frags;
for my $i (1..$o{pairs}) {
  my $r=rand();
  if(@frags && $r < $o{'dup-frac'}) {                      # exact duplicate
    push @frags, [ @{ $frags[int(rand(@frags))] } ];
  } elsif(@frags && $r < $o{'dup-frac'}+$o{'tie-frac'}) {  # same (chr,start), new len
    my $p=$frags[int(rand(@frags))];
    my $L=$o{minlen}+int(rand($o{maxlen}-$o{minlen}+1));
    $L=$rl if $L<$rl;                                       # fragment must span >= read
    push @frags, [ $p->[0], $p->[1], $L ];
  } else {                                                  # fresh fragment
    my $c=$chroms[int(rand(@chroms))];
    my $start=1+int(rand($o{grid}));
    my $L=$o{minlen}+int(rand($o{maxlen}-$o{minlen}+1));
    $L=$rl if $L<$rl;
    push @frags, [ $c, $start, $L ];
  }
}

open(my $fh,'>',$o{out}) or die "cannot write $o{out}: $!\n";
if($o{format} eq 'sam') {
  print $fh "\@HD\tVN:1.0\tSO:queryname\n";
  print $fh "\@SQ\tSN:$_\tLN:$o{grid}\n" for @chroms;
}
my $id=0;
for my $f (@frags) {
  my ($c,$start,$L)=@$f;
  # fragment total length L = second_location - first_location + readlen
  #   => second_location(1-based) = first_location + L - readlen
  my $first1 = $start;
  my $second1 = $start + $L - $rl;
  my $qn="rp".(++$id);
  my $seq=rseq($rl);
  if($o{format} eq 'bowtie') {
    # bowtie cols: id, strand, chr, 0-based pos, seq, qual
    # mate1 (leftmost) supplies chr+pos; mate2 supplies pos+seq. Tool ++'s both.
    my $f0 = $first1-1; my $s0 = $second1-1;
    print $fh join("\t", $qn, '+', $c, $f0, rseq($rl), rqual($rl)),"\n";
    print $fh join("\t", $qn, '-', $c, $s0, $seq,       rqual($rl)),"\n";
  } else {
    # SAM: forward mate (flag&16==0) is leftmost (chr+pos); reverse mate (flag&16)
    # supplies second pos + seq. POS is 1-based. Randomly swap line order; the
    # tool decides leftmost by the flag, not by line order.
    my $fwd = join("\t", $qn, 99,  $c, $first1,  '40', $rl.'M', '=', $second1, 0, rseq($rl), rqual($rl));
    my $rev = join("\t", $qn, 147, $c, $second1, '40', $rl.'M', '=', $first1,  0, $seq,       rqual($rl));
    if(rand() < 0.5){ print $fh "$fwd\n$rev\n"; } else { print $fh "$rev\n$fwd\n"; }
  }
}
close $fh;
print "wrote $o{pairs} pairs ($o{format}) to $o{out}\n";
