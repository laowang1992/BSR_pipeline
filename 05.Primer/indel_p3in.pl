#! /usr/bin/perl

use strict;
use warnings;
use List::Util qw/max/;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --input indelWithSEQ.filter.txt --output p3in --config primer3_config
#
# Required:
#
#	--input <string>			input filename.
#
#	--output <string>			output p3in filename.
#
#	--config <string>			the full path of primer3_config.
#
############################################################


__EOUSAGE__

    ;

my $help_flag;
my $input;
my $output;
my $config;

&GetOptions('help|h' => \$help_flag,
            'input|i=s' => \$input,
            'output|o=s' => \$output,
            'config|c=s' => \$config,
            );

unless ($input && $output && $config) {
	die $usage;
}

open IN,"$input";
open OUT,">$output";
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
	print OUT "SEQUENCE_ID=$a[3]\nSEQUENCE_TEMPLATE=$a[0]\nSEQUENCE_TARGET=249,$max\nPRIMER_TASK=generic\nPRIMER_PRODUCT_SIZE_RANGE=80-180\nPRIMER_MIN_SIZE=18\nPRIMER_OPT_SIZE=20\nPRIMER_MAX_SIZE=27\nPRIMER_MIN_GC=35.0\nPRIMER_MAX_GC=65.0\nPRIMER_MIN_TM=50.0\nPRIMER_OPT_TM=57.0\nPRIMER_MAX_TM=65.0\nPRIMER_NUM_RETURN=3\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=$config\n=\n";
}
close IN;
close OUT;
