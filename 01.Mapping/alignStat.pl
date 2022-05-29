#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Program: deal with STAR result Log.final.out
#
# Usage:  $0 --input <prefix>Log.final.out
#
# Required:
#
#	--input <string>			input filename.
#
############################################################


__EOUSAGE__

    ;

my $help_flag;
my $input;

&GetOptions('help|h' => \$help_flag,
            'input|i=s' => \$input,
            );

unless ($input) {
	die $usage;
}

unless ($input =~ /Log\.final\.out$/) {
	die "Error, input file's name must be like 'PREFIXLog.final.out', dont recognize --input $input ";
}

open IN, "$input";
print STDERR "reading file: $input ...\n";

my $total;
my $uniq;
my $multi;
my $many;
my $um_mis;
my $um_s;
my $um_o;
my $chime;

while(<IN>){
	chomp;
	s/^\s+//;
	if(/^Number of input reads \|\s+(\d+)$/){
		$total = $1;
		#print "$total\n";
	}
	if(/^Uniquely mapped reads number \|\s+(\d+)$/){
		$uniq = $1;
		#print "$uniq\n";
	}
	if(/^Number of reads mapped to multiple loci \|\s+(\d+)$/){
		$multi = $1;
		#print "$multi\n";
	}
	if(/^Number of reads mapped to too many loci \|\s+(\d+)$/){
		$many = $1;
		#print "$many\n";
	}
	if(/^Number of reads unmapped: too many mismatches \|\s+(\d+)$/){
		$um_mis = $1;
		#print "$um_mis\n";
	}
	if(/^Number of reads unmapped: too short \|\s+(\d+)$/){
		$um_s = $1;
		#print "$um_s\n";
	}
	if(/^Number of reads unmapped: other \|\s+(\d+)$/){
		$um_o = $1;
		#print "$um_o\n";
	}
	if(/^Number of chimeric reads \|\s+(\d+)$/){
		$chime = $1;
		#print "$chime\n";
	}
}


my $mu = $multi + $many;
my $um = $um_mis + $um_s + $um_o;
my $uniq_rate = sprintf "%0.2f", $uniq / $total * 100;
my $mu_rate = sprintf "%0.2f", $mu / $total * 100;
my $um_rate = sprintf "%0.2f", $um / $total * 100;
my $chime_rate = sprintf "%0.2f", $chime / $total * 100;

print STDOUT "$total,$uniq,$uniq_rate,$mu,$mu_rate,$um,$um_rate,$chime,$chime_rate\n";

close IN;
