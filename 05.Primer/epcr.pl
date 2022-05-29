#!/bin/env perl
# Date: 2022-01-13
# Usage: perl epcr.pl xxx.p3out

use strict;
use warnings;

open IN, "$ARGV[0]";
open OUT, ">$ARGV[1]";

my $id;
my $num;
my $primerL;
my $primerR;
my @epcr;
my $hit;
my $tmL;
my $tmR;
my $gcL;
my $gcR;
my $len;

print OUT "ID	NUM	PrimerL	PrimerR	Hit	TmL	TmR	GCL	GCR	Length\n";

while(<IN>){
	chomp;
	if(/^SEQUENCE_ID=(.+)$/){
		$id = $1;
	}
	if(/PRIMER_LEFT_(\d)_SEQUENCE=(.+)/){
		# 是否需要清空变量?
		$num = $1;
		$primerL = $2;
	}
	if(/PRIMER_RIGHT_\d_SEQUENCE=(.+)/){
		$primerR = $1;
		# 执行e-pcr，并计算返回结果有多少行;
		@epcr = readpipe("re-PCR -s genome.hash -n 1 -g 1 $primerL $primerR 50-1000");
		$hit = @epcr;
		$hit -= 2;
	}	
	if(/PRIMER_LEFT_\d_TM=(.+)/){
		$tmL = $1;
	}
	if(/PRIMER_RIGHT_\d_TM=(.+)/){
		$tmR = $1;
	}
	if(/PRIMER_LEFT_\d_GC_PERCENT=(.+)/){
		$gcL = $1;
	}
	if(/PRIMER_RIGHT_\d_GC_PERCENT=(.+)/){
		$gcR = $1;
	}
	if(/PRIMER_PAIR_\d_PRODUCT_SIZE=(.+)/){
		$len = $1;
		print OUT "$id\t$num\t$primerL\t$primerR\t$hit\t$tmL\t$tmR\t$gcL\t$gcR\t$len\n"
	}
}

close IN;
close OUT;

