#! /usr/bin/perl
#提出indel位点前后各250bp序列
use strict;
use warnings;

# 输入基因组文件
open IN1,"$ARGV[0]";

# 输入INDEL信息
open IN2,"$ARGV[1]";

# 定义输出文件
open OUT,">$ARGV[2]";

my $m=0;
my $n=0;
my $i=1;
my @chr;
my @seq;
while(<IN1>){
	chomp;
	if($i%2){
		$_=~s/^\>//;
		$chr[$m]=$_;
		$m=$m+1;
	}else{
		$seq[$n]=$_;
		$n=$n+1;
	}
	$i=$i+1;
}
my %genome;
my $b=0;
foreach(@chr){
	$genome{$_}=$seq[$b];
	$b=$b+1;
}
while(<IN2>){
	chomp;
	print OUT "SEQ_500bp\t$_\n" if /^CHROM\tPOS\tID\tREF\tALT/;
	next if /^CHROM\tPOS\tID\tREF\tALT/;
	my @a=split/\s+/,$_;
	#print "$a[0]\n";
	my $se=substr($genome{$a[0]},$a[1]-250,500);
	next if $se=~/N+/;
	print OUT "$se\t$_\n";
}
close IN1;
close IN2;
close OUT;
