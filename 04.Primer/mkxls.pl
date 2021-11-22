open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

print OUT "ID	Left_primer_1	Right_primer_1	Left_TM_1	Right_TM_1	Left_GC_1	Right_GC_1	Length_1	Left_primer_2	Right_primer_2	Left_TM_2	Right_TM_2	Left_GC_2	Right_GC_2	Length_2	Left_primer_3	Right_primer_3	Left_TM_3	Right_TM_3	Left_GC_3	Right_GC_3	Length_3";
while(<IN>){
	chomp;
	if(/SEQUENCE_ID=/){
		$_=~s/SEQUENCE_ID=//;
		print OUT "\n$_\t";
	}
	#
	if(/PRIMER_LEFT_0_SEQUENCE=/){
		$_=~s/PRIMER_LEFT_0_SEQUENCE=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_0_SEQUENCE=/){
		$_=~s/PRIMER_RIGHT_0_SEQUENCE=//;
		print OUT "$_\t";
	}
	if(/PRIMER_LEFT_0_TM=/){
		$_=~s/PRIMER_LEFT_0_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_0_TM=/){
		$_=~s/PRIMER_RIGHT_0_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_LEFT_0_GC_PERCENT=/){
		$_=~s/PRIMER_LEFT_0_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_0_GC_PERCENT=/){
		$_=~s/PRIMER_RIGHT_0_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_PAIR_0_PRODUCT_SIZE=/){
		$_=~s/PRIMER_PAIR_0_PRODUCT_SIZE=//;
		print OUT "$_\t";
	}
	#
	if(/PRIMER_LEFT_1_SEQUENCE=/){
		$_=~s/PRIMER_LEFT_1_SEQUENCE=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_1_SEQUENCE=/){
		$_=~s/PRIMER_RIGHT_1_SEQUENCE=//;
		print OUT "$_\t";
	}	
	if(/PRIMER_LEFT_1_TM=/){
		$_=~s/PRIMER_LEFT_1_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_1_TM=/){
		$_=~s/PRIMER_RIGHT_1_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_LEFT_1_GC_PERCENT=/){
		$_=~s/PRIMER_LEFT_1_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_1_GC_PERCENT=/){
		$_=~s/PRIMER_RIGHT_1_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_PAIR_1_PRODUCT_SIZE=/){
		$_=~s/PRIMER_PAIR_1_PRODUCT_SIZE=//;
		print OUT "$_\t";
	}
	#
	if(/PRIMER_LEFT_2_SEQUENCE=/){
		$_=~s/PRIMER_LEFT_2_SEQUENCE=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_2_SEQUENCE=/){
		$_=~s/PRIMER_RIGHT_2_SEQUENCE=//;
		print OUT "$_\t";
	}
	if(/PRIMER_LEFT_2_TM=/){
		$_=~s/PRIMER_LEFT_2_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_2_TM=/){
		$_=~s/PRIMER_RIGHT_2_TM=//;
		print OUT "$_\t";
	}
	if(/PRIMER_LEFT_2_GC_PERCENT=/){
		$_=~s/PRIMER_LEFT_2_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_RIGHT_2_GC_PERCENT=/){
		$_=~s/PRIMER_RIGHT_2_GC_PERCENT=//;
		print OUT "$_\t";
	}
	if(/PRIMER_PAIR_2_PRODUCT_SIZE=/){
		$_=~s/PRIMER_PAIR_2_PRODUCT_SIZE=//;
		print OUT "$_";
	}
}
close IN;
close OUT;
