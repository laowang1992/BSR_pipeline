#!/bin/bash
. ../script/.conf
bcftools annotate -I '%CHROM\_%POS' -x INFO,^FORMAT/GT ../02.SNP_indel/${filename}.filter.INDELs.vcf | \
	grep -v "##" | \
	sed 's/^#CHROM/CHROM/' > ${filename}.indels.txt

## 这里需进一步对${filename}.indels.txt过滤得到${filename}.indels.filter.txt,
## 比如过滤亲本间相同、亲本杂合、有缺失的位点
## 过滤时不要漏掉表头，CHROM那一行
#head -n 1 ${filename}.indels.txt > ${filename}.indels.filter.txt
#awk '{if(($14=="0/0"&&$17=="1/1")||($14=="1/1"&&$17=="0/0")||($14=="0/0"&&$18=="1/1")||($14=="1/1"&&$18=="0/0"))print $0}' XR_BSA.indels.txt >> ${filename}.indels.filter.txt

perl extr_seq_for_indel_pri.pl --genome ${genome} --input ${filename}.indels.filter.txt --output ${filename}.indelWithSEQ.filter.txt

perl indel_p3in.pl --input ${filename}.indelWithSEQ.filter.txt --output ${filename}.p3in --config /public/home/wangpf/tools/primer3/src/primer3_config/

primer3_core -strict_tags ${filename}.p3in > ${filename}.p3out

perl mkxls.pl --input ${filename}.p3out --output ${filename}.primer.txt

Rscript merge.R ${filename}.indels.filter.txt ${filename}.primer.txt
