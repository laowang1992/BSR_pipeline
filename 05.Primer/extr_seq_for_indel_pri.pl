#! /usr/bin/perl
#提出indel位点前后各250bp序列
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --input indels.filter.txt --genome genome.fa --output indelWithSEQ.filter.txt
#
# Required:
#
#	--input <string>			input filename.
#
#	--genome <string>			genome fasta file.
#
#	--output <string>			output filename.
#
############################################################


__EOUSAGE__

    ;

my $help_flag;
my $input;
my $genome;
my $output;

&GetOptions('help|h' => \$help_flag,
            'input|i=s' => \$input,
            'genome|g=s' => \$genome,
            'output|o=s' => \$output,
            );

unless ($input && $genome && $output) {
	die $usage;
}

unless ($genome =~ /fa(sta)?$/) {
	die "Error, genome file suffix must be .fa or .fasta, dont recognize --genome $genome ";
}

# 读取genome文件
my %genome;
my $fa = Bio::SeqIO->new (-file =>$genome, -f =>'fasta');
while (my $seq_obj = $fa->next_seq) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	$genome{$id} = ${seq};
	print "$id has been read ...\n";
}


# 输入INDEL信息
open IN,"$input";

# 定义输出文件
open OUT,">$output";

my $m=0;
my $n=0;
my $i=1;

while(<IN>){
	chomp;
	print OUT "SEQ_500bp\t$_\n" if /^CHROM\tPOS\tID\tREF\tALT/;
	next if /^CHROM\tPOS\tID\tREF\tALT/;
	my @a=split/\s+/,$_;
	#print "$a[0]\n";
	my $se=substr($genome{$a[0]},$a[1]-250,500);
	next if $se=~/N+/;
	print OUT "$se\t$_\n";
}

close IN;
close OUT;
