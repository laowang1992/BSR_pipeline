#!/bin/bash

work_dir=/public/home/wangpf/workspace/LiuZ/202105_bsr

sampleInfo=${work_dir}/00.data/samples.txt
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

## 转录组分析
trinity=/public/home/wangpf/tools/trinityrnaseq-v2.9.1		# trinity路径
## 进行定量的结构，必须是存在于gtf文件中的结构
featureType=exon
## 对gene或是transcript定量，gene_id或transcript_id
attrType=gene_id
## 是否是链特异性文库，0 (unstranded), 1 (stranded) and 2 (reversely stranded)
strandSpecific=0
## 差异表达计算方法，无生物学重复：edgeR，有生物学重复：DESeq2
de_method=edgeR
dispersion=0.05
## 富集分析参数
orgdb=org.My.eg.db		# 富集分析OrgDB
species=bna				# either the kegg code, scientific name or the common name of the target species 
de_log2FoldChange=1
de_fdr=0.05
enrich_pvalue=0.05
enrich_qvalue=0.05
pdf=FALSE
#########################################################################


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

cat ${sampleInfo} | while read sample group fq1 fq2
do
	# 过滤，质控
	cd ${work_dir}/00.data/01.clean_data
	fastp -i ${fq1} -o ./${sample}_1.clean.fastq.gz \
		-I ${fq2} -O ./${sample}_2.clean.fastq.gz \
		--json=./${sample}.json --html=${sample}.html --report_title="${sample} fastp report" \
		--thread=${thread} --length_required 50
	
	cd ${work_dir}/01.Mapping
	
	# 第一次比对
	mkdir -p ${sample}
	STAR --runMode alignReads \
		--runThreadN ${thread} --genomeDir ${work_dir}/refseq/ \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${sample}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${sample}/${sample}_1 \
		--outSAMtype BAM Unsorted
	
	# 根据第一次比对结果，重新建立索引
	STAR --runThreadN ${thread} --runMode genomeGenerate \
		--genomeDir ./${sample}/ \
		--genomeFastaFiles ${genome} \
		--sjdbFileChrStartEnd ./${sample}/${sample}_1SJ.out.tab \
		--sjdbGTFfile ${gtf} \
		--sjdbOverhang ${sjdbOverhang}
	
	# 第二次比对
	STAR --runMode alignReads \
		--runThreadN ${thread} --genomeDir ./${sample} \
		--readFilesIn ${work_dir}/00.data/01.clean_data/${sample}_1.clean.fastq.gz ${work_dir}/00.data/01.clean_data/${sample}_2.clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${sample}/${sample}_2 \
		--outSAMtype BAM Unsorted
done

# 添加标签、排序
cd ${work_dir}/01.Mapping
awk '{print $1}' ${sampleInfo} | \
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
awk '{print $1}' ${sampleInfo} | \
parallel -j ${jobn} -I% --max-args 1 \
java -jar ${picard} MarkDuplicates \
	I=./%/"%"_rg_added_sorted.bam \
	O=./%/"%"_dedup.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=./%/"%"_dedup.metrics

# MAPQ同步和reads剪切
awk '{print $1}' ${sampleInfo} | \
parallel -j ${jobn} -I% --max-args 1 \
java -jar ${gatk} -T SplitNCigarReads \
	-R ${genome} \
	-I ./%/"%"_dedup.bam \
	-o ./%/"%"_dedup_split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS

# 统计比对情况
awk '{print $1}' ${sampleInfo} | \
parallel -j ${jobn} -I% --max-args 1 \
samtools flagstat %/%_2Aligned.out.bam ">" %/%_2Aligned.flagstat


# 变异检测
cd ${work_dir}/02.SNP_indel

awk '{print $1}' ${sampleInfo} | \
parallel -j ${jobn} -I% --max-args 1 \
	java -jar ${gatk} -T HaplotypeCaller \
		-R ${genome} \
		-ERC GVCF \
		-I ../01.Mapping/%/"%"_dedup_split.bam \
		-dontUseSoftClippedBases \
		-stand_call_conf 20.0 \
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

## 数据统计
cd ${work_dir}/00.data/00.raw_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *[fastq\|fq].gz
cd ${work_dir}/00.data/01.clean_data
mkdir -p QC
fastqc -o ./QC --nogroup --threads ${thread} *clean.fastq.gz
Rscript dataStat.R

## 比对统计
cd ${work_dir}/01.Mapping
echo -ne "Sample,Total reads,Unique reads,Unique rate(%),Multi-mapping reads,Multi-mapping rate(%),Unmapped reads,Unmapped rate(%),Chimeric reads,Chimeric rate(%)\n" > align_stat.csv
for i in $(cut -f1 ${sampleInfo})
do
cd $i
echo -ne "$i," >> ../align_stat.csv
perl ../alignStat.pl --input ${i}_2Log.final.out >> ../align_stat.csv
cd ..
done

## 转录组分析
# 定量
cd ${work_dir}/06.Quantification
for i in $(cut -f2 ${sampleInfo})
do
Rscript run-featurecounts.R \
	-b ../01.Mapping/${i}/${i}_rg_added_sorted.bam \
	-g ${gtf} -o ${i} \
	--nthread ${thread} \
	--featureType ${featureType} \
	--attrType ${attrType} \
	--strandSpecific ${strandSpecific}
done

# 合并
cd ${work_dir}/07.Merge_result
cut -f2 ${sampleInfo} | sed 's/^/..\/06.Quantification\//' | sed 's/$/.count/' > genes.quant_files.txt
perl script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes

# 差异表达。多为无重复，因此用edgeR
cd ${work_dir}/08.DE_analysis
perl ${trinity}/Analysis/DifferentialExpression/run_DE_analysis.pl \
	--matrix ../07.Merge_result/genes.counts.matrix \
	--method ${de_method} \
	--samples_file ../00.data/samples.txt \
	--contrasts contrasts.txt
	--dispersion ${dispersion}

# upset plot
for i in edgeR.*.dir
do
cd ./${i}
Rscript ../upset.R --de_log2FoldChange ${de_log2FoldChange} --de_padj ${de_fdr}
cd ../
done
# volcano plot
for i in edgeR.*.dir/genes.counts.matrix.*DE_results
do
Rscript Volcano_plot.R \
	--de_result ${i} \
	--padj_cutoff ${de_fdr} \
	--log2FC_cutoff ${de_log2FoldChange}
done

