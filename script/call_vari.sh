#!/bin/bash

work_dir=/public/home/wangpf/workspace/LiuZ/202105_bsr

sample=${work_dir}/00.data/samples.txt
genome=${work_dir}/refseq/zs11.genome.fa
gff=${work_dir}/refseq/zs11.v0.gff3
gtf=${work_dir}/refseq/zs11.v0.gtf
sjdbOverhang=149		# sjdbOverhang = readLenght - 1
picard=/public/home/wangpf/tools/picard-tools-2.23.3/picard.jar
gatk=/public/home/wangpf/tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
DISCVRSeq=/public/home/wangpf/tools/DISCVRSeq/DISCVRSeq-1.18.jar
thread=28
jobn=4		# jobs number of picard or gatk
filename=LZ_BSR

PATH=/public/home/wangpf/tools/STAR/bin/Linux_x86_64:$PATH
PATH=/public/home/wangpf/tools/fastp/:$PATH
PATH=/public/home/wangpf/tools/samtools-1.9:$PATH
PATH=/public/home/wangpf/tools/FastQC:$PATH
PATH=/public/home/wangpf/tools/R-4.0.2/bin:$PATH
PATH=/public/home/wangpf/tools/ucsc_utilities:$PATH
PATH=/public/home/wangpf/tools/annovar:$PATH
export PATH

# 建立基因组索引
cd ${work_dir}/refseq
java -jar ${picard} CreateSequenceDictionary R=${genome} O=${genome/fa/dict}
samtools faidx ${genome}
gffread ${gff} -T -o ${gtf}
# STAR索引
STAR --runThreadN ${thread} \
	--runMode genomeGenerate \
	--genomeDir ./ \
	--genomeFastaFiles ${genome} \
	--sjdbGTFfile ${gtf} \
	--sjdbOverhang ${sjdbOverhang}
# annovar索引
gtfToGenePred -genePredExt ${gtf} genome_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile ${genome} genome_refGene.txt --out genome_refGeneMrna.fa

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
		--thread=16 --length_required 50
	
	cd ${work_dir}/01.Mapping
	
	# 第一次比对
	mkdir -p ${i[0]}
	STAR --runMode alignReads \
		--runThreadN ${thread} --genomeDir ${work_dir}/refseq/ \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${i[0]}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${i[0]}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${i[0]}/${i[0]}_1 \
		--outSAMtype BAM
	
	# 根据第一次比对结果，重新建立索引
	STAR --runThreadN ${thread} --runMode genomeGenerate \
		--genomeDir ./${i[0]}/ \
		--genomeFastaFiles ${genome} \
		--sjdbFileChrStartEnd ./${i[0]}/${i[0]}_1SJ.out.tab \
		--sjdbGTFfile ${gtf} \
		--sjdbOverhang ${sjdbOverhang}
	
	# 第二次比对
	STAR --runMode alignReads \
		--runThreadN ${thread} --genomeDir ./${i[0]} \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${i[0]}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${i[0]}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${i[0]}/${i[0]}_2 \
		--outSAMtype BAM
done

# 添加标签、排序
cd ${work_dir}/01.Mapping
awk '{print $2}' ${sample} | \
parallel -j ${jobn} -I% --max-args 1 \
java -jar ${picard} AddOrReplaceReadGroups \
	I=./%/"%"_2Aligned.out.bam \
	O=./%/"%"_rg_added_sorted.bam \
	SO=coordinate \
	RGID=% \
	RGLB=rna \
	RGPL=illumina \
	RGPU=snpcall \
	RGSM=%


# 去重复
awk '{print $2}' ${sample} | \
parallel -j ${jobn} -I% --max-args 1 \
java -jar ${picard} MarkDuplicates \
	I=./%/"%"_rg_added_sorted.bam \
	O=./%/"%"_dedup.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=./%/"%"_dedup.metrics

# MAPQ同步和reads剪切
awk '{print $2}' ${sample} | \
parallel -j ${jobn} -I% --max-args 1 \
java -jar ${gatk} -T SplitNCigarReads \
	-R ${genome} \
	-I ./%/"%"_dedup.bam \
	-o ./%/"%"_dedup_split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS


# 变异检测
cd ${work_dir}/02.SNP_indel

awk '{print $2}' ${sample} | \
parallel -j ${jobn} -I% --max-args 1 \
	java -jar ${gatk} -T HaplotypeCaller \
		-R ${genome} \
		-ERC GVCF \
		-I ../01.Mapping/%/"%"_dedup_split.bam \
		-dontUseSoftClippedBases \
		#-stand_call_conf 20.0 \
		-o ./%.gatk.g.vcf \
		"&>" %.HaplotypeCaller.log

ls *g.vcf > GVCFs.list

# 合并gvcf
java -jar ${gatk} \
     -R ${genome} -T CombineGVCFs \
     -V GVCFs.list \
     -o ${filename}.gatk.g.vcf

# gvcf转vcf
java -jar ${gatk} \
     -R ${genome} -T GenotypeGVCFs \
     -V ${filename}.gatk.g.vcf -o ${filename}.gatk.vcf

# 过滤低质量
java -jar ${gatk} \
	-T VariantFiltration \
	-R ${genome} \
	-V ${filename}.gatk.vcf  \
	-filterName FilterQual -filter "QUAL<30.0" \
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

java -jar ${gatk} \
     -R ${genome} -T VariantsToTable \
     -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
     -V ${filename}.filter.SNPs.vcf -o ../04.Analysis/${filename}.filter.SNPs.table

grep -v "##" ${filename}.filter.SNPs.vcf | sed 's/^#CHROM/CHROM/' > ../04.Analysis/${filename}.filter.SNPs.txt

# vcf QC
java -jar ${DISCVRSeq} VariantQC -O ${filename}.flt.report.html -R ${genome} -V ${filename}.flt.vcf

# annotation
cd ${work_dir}/03.Annotation
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.SNPs.vcf --outfile ./${filename}.filter.SNPs.annovar.input
convert2annovar.pl --format vcf4old ../02.SNP_indel/${filename}.filter.INDELs.vcf --outfile ./${filename}.filter.INDELs.annovar.input
annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${filename}.filter.SNPs.anno --exonsort ./${filename}.filter.SNPs.annovar.input ../refseq
annotate_variation.pl --geneanno --neargene 2000 -buildver genome --dbtype refGene --outfile ./${filename}.filter.INDELs.anno --exonsort ./${filename}.filter.INDELs.annovar.input ../refseq

Rscript AnnoStat.R --snpvar ${filename}.filter.SNPs.anno.variant_function \
	--indelvar ${filename}.filter.INDELs.anno.variant_function \
	--snpex ${filename}.filter.SNPs.anno.exonic_variant_function \
	--indelex ${filename}.filter.INDELs.anno.exonic_variant_function

## 
cd ${work_dir}/00.data/00.raw_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript dataStat.R


