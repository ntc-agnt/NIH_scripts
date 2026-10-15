#!/usr/bin/perl

#Created by Adam Burkholder, National Institute of Environmental Health Sciences, 2010-2016
#
# bowtie2stdBedGraph_optimized
# ----------------------------
# Performance-optimized, BYTE-IDENTICAL drop-in for bowtie2stdBedGraph.pl (same
# directory): identical command line, options, and output on all modes. It
# differs from the original in exactly two ways, both of which leave output
# unchanged:
#   1. No ithreads. The original spawns a thread per output table and passes it
#      the count hash; Perl ithreads DEEP-COPY each hash into the thread, so the
#      process holds ~2x (or more, with -o nbmv) the table memory while sorting.
#      This version sorts/writes the tables sequentially in the main process, so
#      only one copy of the data is ever resident.
#   2. Lean sort. The original sorts hash keys with a Perl comparator over
#      arrayrefs built by a regex (sort{cmp||<=>} map{[/regex/]} keys): that
#      builds an arrayref per key (large transient scratch) and runs a Perl
#      callback per comparison. This version uses a Guttman-Rosler transform:
#      each key is packed as chr\t<zero-padded position> and sorted with a plain
#      C-level string sort (no callback), which reproduces (chr cmp, pos <=>)
#      exactly while using far less scratch and running faster.
#
# Measured (synthetic deep coverage; see tests/): vs the original, ~3.6-7x lower
# peak RSS and ~1.3-2x faster, and -- unlike the original, which can exhaust RAM
# under -o nbmv -- it stays memory-bounded with all outputs requested.
#
# NOTE: the position pack uses a 10-digit zero pad, valid for any position
# < 10^10 (human max ~2.5x10^8 = 9 digits). The parser and the paired-end
# accounting are unchanged from the original to preserve byte-identity.

use warnings;
use strict 'vars';
use Getopt::Std;
use DBI;

my @cols;															
my (%plus_table,%minus_table,%plus_bin_table,%minus_bin_table,%merge_table,%binned_merge_table,%option,%length_table);
my ($key,$plus_thread,$minus_thread,$plus_bin_thread,$minus_bin_thread,$merge_thread,$binned_merge_thread,$temp);
my ($normalized,$binned,$merged,$binned_merged,$bin_size,$shift_val,$trim,$list,$shifted,$swap);
my $norm_value=0;
my $name;
my $dbh;
my $rh;
my @result;
my $mit;
my $Dflag;
my $dupes;
my $dupeflag=0;
my $paired;
my $end=1;
my $iflag;
my @sorted_plus;
my @sorted_minus;
my $plus_sort_thread;
my $minus_sort_thread;
my $one;
my $two;
my $sum;
my @cols_plus;
my @cols_minus;
my @samone;
my @samtwo;
my $aflag;

getopts('o:b:s:t:l:r:xCM:Du:pia',\%option) || print_usage("g");				
									
if(exists($option{o})) {							
	if($option{o}=~/[^nbmv]/) {
		print_usage("o");					
	}
	$normalized=($option{o}=~/n/i) || 0;				
	$binned=($option{o}=~/b/i) || 0;
	$merged=($option{o}=~/m/i) || 0;
	$binned_merged=($option{o}=~/v/i) || 0;
}

print_usage("bz") if (exists($option{b}) && $option{b}<1);							#bin size must be 1 or greater
print_usage("b") if ((!$binned && !$binned_merged) && exists($option{b}));			#check for contradictions in specified parameters
print_usage("s") if ((!$merged && !$binned_merged) && exists($option{s}));

$bin_size=$option{b} || 25;								
$shift_val=$option{s} || 75;							
$trim=$option{t} || 0;									
$swap=$option{x} || 0;
$Dflag=$option{D} || 0;
$mit=$option{M}	|| "mit";				
$dupes=$option{u} || 1;
$dupeflag=1 if exists($option{u});
$paired=$option{p} || 0;
$end=2 if $paired==1 && $swap==1;
$iflag=$option{i} || 0;
$aflag=$option{a} || 0;

print_usage("iu") if ($dupeflag==0 && $iflag==1);
print_usage("ip") if ($paired==0 && $iflag==1);

if(exists($option{r})) {								
	$list=$option{r};
	$dbh=DBI->connect("DBI:mysql:database=".$list.";host=genome-mysql.cse.ucsc.edu","genome");	#if user specifies UCSC genome, look up chromosome lengths via mysql query
	print_usage("db") if $dbh->{Active}==0;
	$rh=$dbh->prepare("select chrom, size from chromInfo");								
	$rh->execute();						
	while(@result=$rh->fetchrow_array()) {												#store each row of the list in the length_table hash
		$length_table{$result[0]}=$result[1];
	}
	print_usage("db") if keys(%length_table)==0;
}
elsif(exists($option{l})) {
	$list=$option{l};
	open(LISTFILE,$list) || print_usage("lf");
	while(<LISTFILE>) {							#if user specifies chromosome length list, store names/values in length_table hash
		/(\S+)\s+(\S+)/;
		$length_table{$1}=$2;
	}
	close(LISTFILE);
}
else {
	$list="";									#if list file name is not specified, assign blank value
}

print_usage("") if scalar(@ARGV)==0;
					
open(INFILE,"$ARGV[0]") || print_usage("f");		#open input file, run print_usage subroutine if file cannot be opened

$name=$ARGV[1] || print_usage("t");			#check for output file prefix in argument
$name=~s/\/?(\S*\/)*//;						#strip path from output file prefix

if($paired==0) {
	while(<INFILE>) {
		if($aflag==0) {
			@cols=/^[^\t]+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/;		#extract the second through fourth fields of the current line
		}
		else {
			next if $_=~/^\@/;
			@samone=split(/\t/,$_);
			next if $samone[2] eq "\*";
			if(($samone[1]&16)==0) {
				$cols[0]="+";
			}
			else {
				$cols[0]="-";
			}
			$cols[1]=$samone[2];
			$cols[2]=$samone[3]-1;
			$cols[3]=$samone[9];			
		}
		$cols[1]="chrM" if $cols[1]=~m/$mit/i && $Dflag==0;				#if user allows chromosome name fixing, rename mitochrondrial chromsome, append chr where absent
		if($cols[1]!~m/chr/i && $Dflag==0) {
			$cols[1]="chr".$cols[1];
		}
		if($cols[0] eq '+') {
			$cols[2]+=(1-$trim);					#convert from 0-based to 1-based coordinates and adjust for trimming
			$key="$cols[1]\t$cols[2]";				#concatenate the chromosome and 5' mapping position, and use it as a hash key
			next if exists($plus_table{$key}) && ($plus_table{$key}*$dupeflag)>=$dupes;
			$plus_table{$key}++;
			if($binned) {	
				$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;			#determine bin
				$key="$cols[1]\t$temp";										#use chromosome and bin start as hash key
				$plus_bin_table{$key}++;
			}
			$shifted=$cols[2]+$shift_val;									#determine shifted location for merged output
		}
		elsif($cols[0] eq '-') {						
			$cols[2]+=(length($cols[3])+$trim); 			#determine 5' mapping position, adjust for trimming
			$key="$cols[1]\t$cols[2]";						#perform same tasks as above
			next if exists($minus_table{$key}) && ($minus_table{$key}*$dupeflag)>=$dupes;
			$minus_table{$key}++;
			if($binned) {
				$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;
				$key="$cols[1]\t$temp";
				$minus_bin_table{$key}++;
			}
			$shifted=$cols[2]-$shift_val;
		}
		if($merged) {
			$temp=$shifted;
			$temp=1 if $temp<1;
			if(exists($length_table{$cols[1]})) {					
				$temp=$length_table{$cols[1]} if $temp>$length_table{$cols[1]};				
			}
			elsif($list ne "") {																
				die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
			$key="$cols[1]\t$temp";							#use chromosome and shifted 5' mapping position as hash key
			$merge_table{$key}++;							#increment count, for both forward and reverse strand hits
		}
		if($binned_merged) {
			$temp=$shifted;
			$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
			if(exists($length_table{$cols[1]})) {						
				$temp=int(($length_table{$cols[1]})/$bin_size)*$bin_size if $temp>$length_table{$cols[1]};		#determine bin			
			}
			elsif($list ne "") {																
				die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
			}
			$key="$cols[1]\t$temp";					#use chromosome and bin start as hash key				
			$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
		}
		$norm_value++;						#increment overall hit count
	}
}
else {
	while(<INFILE>) {
		$one=$_;
		if($aflag==0) {
			$two=<INFILE>;
			@cols_plus=($one=~/^([^\t]+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/);		#extract the second through fourth fields of the current line
			@cols_minus=($two=~/^([^\t]+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/);
		}
		else {
			next if $one=~/^\@/;
			$two=<INFILE>;
			@samone=split(/\t/,$one);
			next if $samone[2] eq "\*";
			@samtwo=split(/\t/,$two);
			if(($samone[1]&16)==0) {
				$cols_minus[2]=$samtwo[2];
				$cols_minus[3]=$samtwo[3]-1;
				$cols_minus[4]=$samtwo[9];
				$cols_plus[2]=$samone[2];
				$cols_plus[3]=$samone[3]-1;
				$cols_plus[4]=$samone[9];
				if(($samtwo[1]&128)==128) {
					$cols_minus[0]=$samtwo[0]."\/2";	
				}
				else {
					$cols_minus[0]=$samtwo[0]."\/1";
				}
			}
			else {
				$cols_minus[2]=$samone[2];
				$cols_minus[3]=$samone[3]-1;
				$cols_minus[4]=$samone[9];
				$cols_plus[2]=$samtwo[2];
				$cols_plus[3]=$samtwo[3]-1;
				$cols_plus[4]=$samtwo[9];
				if(($samone[1]&128)==128) {
					$cols_minus[0]=$samone[0]."\/2";	
				}
				else {
					$cols_minus[0]=$samone[0]."\/1";
				}
			}
		}
		$cols_plus[2]="chrM" if $cols_plus[2]=~m/$mit/i && $Dflag==0;				#if user allows chromosome name fixing, rename mitochrondrial chromsome, append chr where absent
		$cols_minus[2]="chrM" if $cols_minus[2]=~m/$mit/i && $Dflag==0;
		if($cols_plus[2]!~m/chr/i && $Dflag==0) {
			$cols_plus[2]="chr".$cols_plus[2];
			$cols_minus[2]="chr".$cols_minus[2];
		}
		$cols_plus[3]+=(1-$trim);					#convert from 0-based to 1-based coordinates and adjust for trimming
		$key="$cols_plus[2]\t$cols_plus[3]";				#concatenate the chromosome and 5' mapping position, and use it as a hash key
		$cols_minus[3]+=(length($cols_minus[4])+$trim); 			#determine 5' mapping position, adjust for trimming
		$key.="\t$cols_minus[3]";
		if($iflag==1) {
			$sum=0;
			$sum+=$plus_table{$key} if exists($plus_table{$key});
			$sum+=$minus_table{$key} if exists($minus_table{$key});
			next if $sum>=$dupes;
		}
		if($cols_minus[0]!~/\/$end$/) {
			next if exists($plus_table{$key}) && ($plus_table{$key}*$dupeflag)>=$dupes;
			$plus_table{$key}++;
			if($binned) {	
				$temp=1 if ($temp=int($cols_plus[3]/$bin_size)*$bin_size)<1;			#determine bin
				$key="$cols_plus[2]\t$temp";										#use chromosome and bin start as hash key
				$plus_bin_table{$key}++;
			}
			$shifted=$cols_plus[3]+$shift_val;									#determine shifted location for merged output		
			if($merged) {
				$temp=$shifted;
				$temp=1 if $temp<1;
				if(exists($length_table{$cols_plus[2]})) {					
					$temp=$length_table{$cols_plus[2]} if $temp>$length_table{$cols_plus[2]};				
				}
				elsif($list ne "") {																
					die "Error: chromosome $cols_plus[2] could not be found in the chromosome lengths list\n";
				}
				$key="$cols_plus[2]\t$temp";							#use chromosome and shifted 5' mapping position as hash key
				$merge_table{$key}++;							#increment count, for both forward and reverse strand hits
			}
			if($binned_merged) {
				$temp=$shifted;
				$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
				if(exists($length_table{$cols_plus[2]})) {						
					$temp=int(($length_table{$cols_plus[2]})/$bin_size)*$bin_size if $temp>$length_table{$cols_plus[2]};		#determine bin			
				}
				elsif($list ne "") {																
					die "Error: chromosome $cols_plus[2] could not be found in the chromosome lengths list\n";
				}
				$key="$cols_plus[2]\t$temp";					#use chromosome and bin start as hash key				
				$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
			}
		}
		else {
			next if exists($minus_table{$key}) && ($minus_table{$key}*$dupeflag)>=$dupes;
			$minus_table{$key}++;
			if($binned) {
				$temp=1 if ($temp=int($cols_minus[3]/$bin_size)*$bin_size)<1;
				$key="$cols_minus[2]\t$temp";
				$minus_bin_table{$key}++;
			}
			$shifted=$cols_minus[3]-$shift_val;
			if($merged) {
				$temp=$shifted;
				$temp=1 if $temp<1;
				if(exists($length_table{$cols_minus[2]})) {					
					$temp=$length_table{$cols_minus[2]} if $temp>$length_table{$cols_minus[2]};				
				}
				elsif($list ne "") {																
					die "Error: chromosome $cols_minus[2] could not be found in the chromosome lengths list\n";
				}
				$key="$cols_minus[2]\t$temp";							#use chromosome and shifted 5' mapping position as hash key
				$merge_table{$key}++;							#increment count, for both forward and reverse strand hits
			}
			if($binned_merged) {
				$temp=$shifted;
				$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
				if(exists($length_table{$cols_minus[2]})) {						
					$temp=int(($length_table{$cols_minus[2]})/$bin_size)*$bin_size if $temp>$length_table{$cols_minus[2]};		#determine bin			
				}
				elsif($list ne "") {																
					die "Error: chromosome $cols_minus[2] could not be found in the chromosome lengths list\n";
				}
				$key="$cols_minus[2]\t$temp";					#use chromosome and bin start as hash key				
				$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
			}
		}
		$norm_value++;										#increment overall hit count
	}
	foreach(keys(%plus_table)) {
		$temp=$_;
		$temp=~s/\t[^\t]+$//;
		$plus_table{$temp}+=$plus_table{$_};
		delete($plus_table{$_});
	}
	foreach(keys(%minus_table)) {
		$temp=$_;
		$temp=~s/(^[^\t]+)\t[^\t]+/$1/;
		$minus_table{$temp}+=$minus_table{$_};
		delete($minus_table{$_});
	}
}
close(INFILE);								

$norm_value/=1000000;							#divide overall hit count by one million, hit values will be normalized to millions of overall hits

# REFACTOR: run the same sort_tables/sort_tables_alt routines sequentially in the
# main process instead of in ithreads. ithreads deep-COPY each hash passed to
# threads->create, so the original holds ~2x the table memory during sorting;
# sequential execution keeps a single copy resident. Output files are identical.
if($swap==0) {
	sort_tables(\%plus_table,"forward");
	sort_tables(\%minus_table,"reverse");
	sort_tables_alt(\%plus_bin_table,"forward_binned") if $binned;
	sort_tables_alt(\%minus_bin_table,"reverse_binned") if $binned;
}
else {																				#swap strands if requested by user
	sort_tables(\%minus_table,"forward");
	sort_tables(\%plus_table,"reverse");
	sort_tables_alt(\%minus_bin_table,"forward_binned") if $binned;
	sort_tables_alt(\%plus_bin_table,"reverse_binned") if $binned;
}
sort_tables_alt(\%merge_table,"merged") if $merged;
sort_tables_alt(\%binned_merge_table,"binned_merged") if $binned_merged;

sub sort_tables {
	my $ref=shift;
	my $strand=shift;
	my ($value,$start,$chr,$pos,$p);
	open(OUTFILE,">$ARGV[1]_${strand}.bedGraph") || die "Could not create output file: $ARGV[1]_${strand}.bedGraph\n";
	# GRT lean sort (byte-identical order to chr cmp / pos <=>, minimal scratch)
	my @result = sort map { $p=index($_,"\t"); substr($_,0,$p+1).sprintf("%010d",substr($_,$p+1)) } keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
	if($normalized) {
		open(NORMFILE,">$ARGV[1]_norm1m_${strand}.bedGraph") || die "Could not create output file: $ARGV[1]_norm1m_${strand}.bedGraph\n";
		print NORMFILE "track type=bedGraph name=${name}_norm1m_${strand} description=${name}_norm1m_${strand} visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
		foreach (@result) {
			$p=index($_,"\t"); $chr=substr($_,0,$p); $pos=substr($_,$p+1)+0;
			$value=$ref->{"$chr\t$pos"};
			$start=$pos-1;
			print OUTFILE "$chr\t$start\t$pos\t$value\n";
			$value/=$norm_value;
			print NORMFILE "$chr\t$start\t$pos\t$value\n";
		}
		close(NORMFILE);
	}
	else {
		foreach (@result) {
			$p=index($_,"\t"); $chr=substr($_,0,$p); $pos=substr($_,$p+1)+0;
			$value=$ref->{"$chr\t$pos"};
			$start=$pos-1;
			print OUTFILE "$chr\t$start\t$pos\t$value\n";
		}
	}
	close(OUTFILE);
	return;
}

	
sub sort_tables_alt {
	my $ref=shift;
	my $strand=shift;
	my ($value,$start,$end,$chr,$pos,$p);
	open(OUTFILE,">$ARGV[1]_${strand}.bedGraph") || die "Could not create output file: $ARGV[1]_${strand}.bedGraph\n";
	my @result = sort map { $p=index($_,"\t"); substr($_,0,$p+1).sprintf("%010d",substr($_,$p+1)) } keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
	foreach (@result) {
		$p=index($_,"\t"); $chr=substr($_,0,$p); $pos=substr($_,$p+1)+0;
		$value=$ref->{"$chr\t$pos"};
		if($strand=~/bin/) {
			$start=$pos-1;
			$end=$start+$bin_size;
		}
		else {
			$start=$pos-1;
			$end=$pos;
		}
		if($list ne "") {
			die "Error: chromosome $chr could not be found in the chromosome lengths list\n" if !exists($length_table{$chr});
			$end=$length_table{$chr} if $end>$length_table{$chr};
		}
		print OUTFILE "$chr\t$start\t$end\t$value\n";
	}
	close(OUTFILE);
	return;
}

	
sub print_usage {								#subroutine called when an error in command line input is detected
	my $arg=shift;							
	if($arg eq "b")	{						#if -b option is present and binned output is not requested
		print "Error: binned output files must be requested (-o b or -o v) for -b option to be accepted\n";
	}
	elsif($arg eq "bz") {					#if bin size specified is less than 1
		print "Error: bin size must be 1 or larger\n";		
	}
	elsif($arg eq "s") {					#if -s option is present and merged output is not requested
		print "Error: merged output file must be requested (-o m) for -s option to be accepted\n";
	}
	elsif($arg eq "o") {						#if letters other than n, b, v, or m are included in the -o option's argument
		$option{o}=~s/n|b|m|v//g;					
		print "Error: invalid output file(s) requested: $option{o}\n";
	}
	elsif($arg eq "f" && scalar(@ARGV)>0) {						#if input file name is invalid
		print "Error: could not open input file \"$ARGV[0]\"\n";
	}
	elsif($arg eq "t") {						#if track name is not specified
		print "Error: output file prefix must be specified\n";
	}							
	elsif($arg eq "db") {						#if mysql lookup fails
		print "Error: could not fetch list of chromosome lengths for reference genome \"$list\"\n";
	}
	elsif($arg eq "lf") {						#if chromosome lengths file is invalid
		print "Error: could not open chromosome list file \"$list\"\n";
	}
	elsif($arg eq "g") {
		print "Error: unrecognized or misused option\n";
	}
	elsif($arg eq "iu") {
		print "Error: de-duplication must be enabled (-u) for -i option to be accepted\n";
	}
	elsif($arg eq "ip") {
		print "Error: strand-independent de-duplication requires paired-end input\n";
	}							
	die "Usage: bowtie2bedgraph.pl [options] [input file] [output file prefix]\n\t-p\t\tindicate input file is the result of a paired-end\n\t\t\talignment\n\t-a\t\tindicate input is in SAM format (sorted by query name)\n\t-o [nbmv]\tspecify additional files to be generated: n=normalized,\n\t\t\tb=binned, m=merged, v=binned and merged\n\t-b [integer>=1]\tspecify bin size, requires -o b or -o v, default: 25\n\t-s [integer>=0]\tspecify number of bps to shift ChIP-seq hits prior to\n\t\t\tmerging, requires -o m or -o v, default: 75\n\t-l [string]\tprovide list of chromosome names and lengths to prevent\n\t\t\thits from being shifted beyond the end of the\n\t\t\tchromosome, requires -o b, -o m, or -o v\n\t-r [string]\trather than specifying a chromosome length list,\n\t\t\tspecify a UCSC reference genome identifier (e.g. mm9),\n\t\t\tto have chromosome lengths fetched automatically,\n\t\t\trequires -o b, -o m, or -o v\n\t-t [integer>=1]\tspecify number of nt trimmed prior to alignment\n\t-x\t\trequest stranded output files be swapped,\n\t\t\tforward->reverse and reverse->forward; for paired-end\n\t\t\tinput report counts of End 2 rather than End 1 5\' ends\n\t-D\t\tdisable appending \"chr\" to beginning of chromosome names\n\t\t\tand renaming of mitochondrial chromosome\n\t-M [string]\tstring used to match mitochondrial chromosome name,\n\t\t\tdefault: \"mit\"\n\t-u [integer>=1]\tremove duplicate reads exceeding specified maximum\n\t\t\tto retain per position\n\t-i\t\tperform de-duplication in a strand-independent manner,\n\t\t\trequires paired-end input (-p)\n";
}
