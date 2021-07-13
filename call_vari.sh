#!/bin/bash

work_dir=/public/home/wangpf/workspace/LiuZ/202105_bsr
sample=${work_dir}/00.data/samples.txt
genome=${work_dir}/refseq/zs11.genome.fa
gtf=${work_dir}/refseq/zs11.v0.gtf
picard=/public/home/wangpf/tools/picard-tools-2.23.3/picard.jar
gatk=/public/home/wangpf/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
thread=28
filename=LZ_BSR

PATH=/public/home/wangpf/tools/STAR/bin/Linux_x86_64:$PATH
PATH=/public/home/wangpf/tools/fastp/:$PATH
PATH=/public/home/wangpf/tools/samtools-1.9:$PATH
PATH=/public/home/wangpf/tools/FastQC:$PATH
PATH=/public/home/wangpf/tools/R-4.0.2/bin:$PATH
export PATH

# 建立基因组索引
cd ${work_dir}/refseq
#java -jar ${picard} CreateSequenceDictionary R=${genome} O=${genome/fa/dict}
#samtools faidx ${genome}
#STAR --runThreadN ${thread} --runMode genomeGenerate \
#	--genomeDir ./ \
#	--genomeFastaFiles ${genome} \
#	--sjdbGTFfile ${gtf}


cd ${work_dir}
IFS_OLD=$IFS
IFS=$'\n'

for i in $(cat ${sample})
do
	IFS=$'\t'
	i=($i)
	IFS=$IFS_OLD
	
	# 过滤，质控
	cd ${work_dir}/00.data/01.clean_data
	fastp -i ${i[2]} -o ./${i[0]}_1.clean.fastq.gz \
		-I ${i[3]} -O ./${i[0]}_2.clean.fastq.gz \
		--json=./${i[0]}.json --html=${i[0]}.html --report_title="${i[0]} fastp report" \
		--thread=16 --length_required 100
	
	cd ${work_dir}
	
	# 第一次比对
	mkdir -p ${i[0]}
	STAR --runThreadN ${thread} --genomeDir ${work_dir}/refseq/ \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${i[0]}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${i[0]}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${i[0]}/${i[0]}_1
	
	# 根据第一次比对结果，重新建立索引
	STAR --runThreadN ${thread} --runMode genomeGenerate \
		--genomeDir ./${i[0]}/ \
		--genomeFastaFiles ${genome} \
		--sjdbFileChrStartEnd ./${i[0]}/${i[0]}_1SJ.out.tab
	
	# 第二次比对
	STAR --runThreadN ${thread} --genomeDir ./${i[0]} \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${i[0]}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${i[0]}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${i[0]}/${i[0]}_2
done

# 添加标签、排序
cd ${work_dir}
awk '{print $2}' ${sample} | \
parallel -j 6 -I% --max-args 1 \
java -jar ${picard} AddOrReplaceReadGroups \
	I=./%/"%"_2Aligned.out.sam \
	O=./%/"%"_rg_added_sorted.bam \
	SO=coordinate \
	RGID=% \
	RGLB=rna \
	RGPL=illumina \
	RGPU=snpcall \
	RGSM=%


# 去重复
awk '{print $2}' ${sample} | \
parallel -j 6 -I% --max-args 1 \
java -jar ${picard} MarkDuplicates \
	I=./%/"%"_rg_added_sorted.bam \
	O=./%/"%"_dedup.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=./%/"%"_dedup.metrics

# MAPQ同步和reads剪切
awk '{print $2}' ${sample} | \
parallel -j 6 -I% --max-args 1 \
java -jar ${gatk} -T SplitNCigarReads \
	-R ${genome} \
	-I ./%/"%"_dedup.bam \
	-o ./%/"%"_dedup_split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS

# 变异检测
java -jar ${gatk} -T HaplotypeCaller \
	-R ${genome} \
	-I ./LL/LL_dedup_split.bam \
	-I ./RL/RL_dedup_split.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o ./${filename}.gatk.vcf


# 过滤低质量
java -jar ${gatk} \
	-T VariantFiltration \
	-R ${genome} \
	-V ${filename}.gatk.vcf  \
	-window 35 \
	-cluster 3 \
	-filterName FilterFS -filter "FS > 30.0" \
	-filterName FilterQD -filter "QD < 2.0" \
	-o ./${filename}.flt.vcf


grep -vP "\tFilter" ${filename}.flt.vcf > ${filename}.filter.vcf

java -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType SNP \
     -V ${filename}.filter.vcf -o ${filename}.filter.SNPs.vcf

java -jar ${gatk} \
     -R ${genome} -T SelectVariants \
     -selectType INDEL \
     -V ${filename}.filter.vcf -o ${filename}.filter.INDELs.vcf


grep -v "##" ${filename}.filter.SNPs.vcf | grep -v "\./\." | grep -v "1/1.*1/1" | sed 's/^#CHROM/CHROM/' > ./${filename}.filter.SNPs.txt

## 
cd ${work_dir}/00.data/00.raw_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript dataStat.R


