#! /usr/bin/perl

use strict;
use warnings;
use List::Util qw/max/;

open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";
my @len;

while(<IN>){
	chomp;
	next if /^SEQ_500bp/;
	my @a=split/\s+/,$_;
	my @b=split(/\,/,$a[5]);
	my @len;
	$len[0]=length($a[4]);
	my $i=1;
	foreach(@b){
		$len[$i]=length($_);
		$i+=1;
	}
	my $max=0;
	foreach(@len){
		if($_>$max){$max=$_}
	}
	if($max>$len[0]){$max=$len[0]}
	print OUT "SEQUENCE_ID=$a[3]\nSEQUENCE_TEMPLATE=$a[0]\nSEQUENCE_TARGET=249,$max\nPRIMER_TASK=generic\nPRIMER_PRODUCT_SIZE_RANGE=80-180\nPRIMER_MIN_SIZE=18\nPRIMER_OPT_SIZE=20\nPRIMER_MAX_SIZE=27\nPRIMER_MIN_GC=35.0\nPRIMER_MAX_GC=65.0\nPRIMER_MIN_TM=50.0\nPRIMER_OPT_TM=57.0\nPRIMER_MAX_TM=65.0\nPRIMER_NUM_RETURN=3\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=/dat1/Project/wpf/tools/primer3-2.4.0/src/primer3_config/\n=\n";
}
close IN;
close OUT;
