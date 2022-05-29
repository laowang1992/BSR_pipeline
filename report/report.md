# BSR分析报告
## BSR-seq原理
BSR-seq (bulked segregant RNA-seq)<sup>[1](#ref)</sup>是一种将Bulked‐segregant analysis (BSA)<sup>[2](#ref),[3](#ref)</sup>和转录组测序相结合，快速定位QTL的方法。目标性状有差异的双亲构建的分离群体中分别选取极端表型个体进行等量混合构建两个极端表型bulked RNA pool并进行测序。随后进行变异分析筛选出双亲间SNP位点并分别计算两个bulked RNA pool中每个SNP位点上某一亲本基因型read覆盖深度占该位点总read深度的比值，即SNP index，通过两个bulked RNA pool的SNP index相减即得到ΔSNP index。在基因组所有区域中，目标基因及其连锁的区域由于根据表型受到相反的选择在两个bulked RNA pool中表现出不同的趋势，因此ΔSNP index会显著偏离0附近；另一方面，于目标形状无关的区域则两个bulked RNA pool则表现为相似的变化趋势，因此ΔSNP index会在0附近波动。
## BSR-seq流程图
<img src="./image/QTL-seq分析流程图.png" width="800px" alt="BSR-seq流程图">
## 实验材料
* W2：6170，隐性性状亲本（莲状花朵）；
* W6：V02，显性性状亲本（冠状花朵）；
* W9：F2混池，显性性状（冠状花朵）；
* W11：F2混池，隐性性状（莲状花朵）；

注：每个样本测10G数据。
## 分析流程
### 取样、提取RNA、建库和测序
选取分离群体中极端表型个体各30-50株以及双亲取植物组织，分别提取RNA并将分离群体极端表型个体按照等量原则混合构建bulked RNA pools，然后进行建库和高通量测序。（具体流程应参考实验设计以及公司测序报告）。
### 数据过滤
使用fastp<sup>[4](#ref)</sup>（version: 0.20.0）对raw data进行过滤去除接头序列和低质量reads得到clean data。统计过滤前后total bases、total reads、Q30、Q20、GC content以及有效数据比率（data_stat.csv/txt），同时使用FastQC（version: 0.11.9）对过滤前后的数据进行质量评估（QC/sample_fastqc.html）。
### 比对到参考基因组
使用STAR<sup>[5](#ref)</sup>（version: 2.7.3a）软件将clean data比对到参考基因组上(2-pass mapping mode)，然后使用Picard tools<sup>[6](#ref)</sup>（version: 2.23.2）去除建库过程中产生的PCR重复，每个样本得到一各BAM文件用于后续variations calling。
### SNPs and InDels calling、过滤和注释
变异分析使用Genome Analysis Toolkit，GATK<sup>[7](#ref)</sup>（version: 3.8-0-ge9d806836）完成，首先使用GATK的HaplotypeCaller功能对样本单独分析再使用CombineGVCFs功能合并，随后使用GenotypeGVCFs功能得到SNP和INDEL信息，最后使用VariantFiltration功能按照参数–filterExpression QUAL<50.0 || FS > 30.0 || QD < 2.0过滤原始的变异位点的到可靠的变异信息，使用DISCVRSeq（version：1.18）的VariantQC功能<sup>[8](#ref)</sup>对变异信息进行可视化。

使用ANNOVAR<sup>[9](#ref)</sup>软件对过滤后的可靠变异位点进行注释并统计变异位点的分布情况和突变类型。
### BSR-seq分析
进一步筛选SNP位点，保留QUAL > 100、亲本内纯和且双亲间不同的、亲本及两个混池中覆盖深度 > 4的SNP用于最后的分析。首先计算每个混池的SNP index值，随后计算Delta SNP index并绘图，其中点图按照1Mb区间、200kb步长进行滑窗统计，折线图是由R软件（version: 4.0.2）扩展包QTLseqr<sup>[10](#ref)</sup>（version: 0.7.5.2）按照3Mb窗口大小统计得到，置信区间按照Takagi 等（2013）<sup>[11](#ref)</sup>描述方法计算得到，超出95或99%置信区间的区域视为目标性状候选QTL。
### 引物设计
为进行QTL验证和进一步精细定位，根据变异分析中得到的可靠InDel位点进行全基因组范围的引物设计。

我们首先需要确定引物设计的位点。在vcf文件的第10及以后的列为各样本信息，从中可以得到各样本基因型，对应第9列GT位置，两个数字中间用“ / ”分开，这两个数字表示二倍体sample的基因型。0表示样本中有ref的allele；1表示样本中有variant的allele；2表示有第二个variant的allele。
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	W11	W2	W6	W9
TeORChr1	436921	.	AG	A	84.96	PASS	AC=4;AF=1.00;AN=4;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=60.00;QD=29.43;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,2:2:6:79,6,0	1/1:0,1:1:3:39,3,0	./.:0,0:0:.:0,0,0	./.:0,0:0:.:0,0,0
TeORChr1	616022	.	C	CA	213.73	PASS	AC=4;AF=1.00;AN=4;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=60.00;QD=30.53;SOR=1.609	GT:AD:DP:GQ:PL	1/1:0,3:3:9:107,9,0	1/1:0,4:4:12:141,12,0	./.:0,0:0:.:0,0,0	./.:0,0:0:.:0,0,0
```
在鉴定的INDELs中筛选两个亲本内纯和且在亲本间不同的位点，即两亲本基因型分别为0/0、1/1或者1/1、0/0的位点，再根据INDEL位置提取上下游各250 bp基因组序列，随后使用primer3<sup>[12](#ref)</sup>（version: 2.5.0）进行引物设计（需根据实验室电泳仪器能区分的InDel大小挑选合适的位点，引物根据参考基因组设计，其真实位置和特异性需根据序列比对和实验结果进行验证）。
## 结果
### 数据过滤
过滤前后数据统计结果如下表：

|   sample	|	total_reads.raw	|	total_bases.raw	|	q20_bases.raw	|	q30_bases.raw	|	q20_rate.raw	|	q30_rate.raw	|	read1_mean_length.raw	|	read2_mean_length.raw	|	gc_content.raw	|	total_reads.clean	|	total_bases.clean	|	q20_bases.clean	|	q30_bases.clean	|	q20_rate.clean	|	q30_rate.clean	|	read1_mean_length.clean	|	read2_mean_length.clean	|	gc_content.clean	|	effective.rate  |
|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	----	|	---	|
|   W11	|	90387568	|	13025802859	|	12726623094	|	12085316955	|	0.977032	|	0.927798	|	143	|	144	|	0.437091	|	89235954	|	12923921480	|	12634457106	|	12003918241	|	0.977602	|	0.928814	|	144	|	144	|	0.436853	|	0.992178495 |
|   W2	|	106786040	|	15142890283	|	14811997130	|	14091458542	|	0.978149	|	0.930566	|	141	|	141	|	0.438337	|	105204774	|	15018663559	|	14697395747	|	13988383213	|	0.978609	|	0.9314	|	142	|	142	|	0.438119	|	0.991796366 |
|   W6	|	84814398	|	12234923142	|	11953317370	|	11346170204	|	0.976983	|	0.927359	|	144	|	144	|	0.435221	|	83759784	|	12140080837	|	11866997450	|	11269643712	|	0.977506	|	0.928301	|	144	|	145	|	0.434988	|	0.99224823  |
|   W9	|	81520108	|	11790275827	|	11514960002	|	10925595007	|	0.976649	|	0.926662	|	144	|	144	|	0.435977	|	80530406	|	11700903139	|	11434212885	|	10854501805	|	0.977208	|	0.927664	|	145	|	145	|	0.435754	|	0.992419797 |
### 比对
比对结果统计如下表：

|	Sample	|	Total reads	|	Unique reads	|	Unique rate(%)	|	Multi-mapping reads	|	Multi-mapping rate(%)	|	Unmapped reads	|	Unmapped rate(%)	|	Chimeric reads	|	Chimeric rate(%)	|
|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|	---	|
|	W11	|	44617977	|	25799665	|	57.82	|	908990	|	2.04	|	17909322	|	40.14	|	0	|	0	|
|	W2	|	52602387	|	18429507	|	35.04	|	634169	|	1.21	|	33538711	|	63.76	|	0	|	0	|
|	W6	|	41879892	|	26591178	|	63.49	|	759809	|	1.81	|	14528905	|	34.69	|	0	|	0	|
|	W9	|	40265203	|	24741123	|	61.45	|	804758	|	2	|	14719322	|	36.56	|	0	|	0	|

Unique mapping rate偏低，Unmapped rate偏高。有多种原因可能造成这种情况，包括定位材料与基因组材料亲缘关系较远，测序数据污染等多用因素，需要进一步分析鉴定。
### SNP/InDel calling
原始SNP + InDel共456,193个，过滤后剩341,414个，其中SNPs 310,383个，InDels 31,031个。
### 注释
位于exonic上的SNP数目最多，其次是intronic。其中分布在exonic上的synonymous SNP比nonsynonymous SNP数目更多，且stopgain/stoploss类型则最少。

|	Distribution	|	varType	|	number	|
|	---	|	---	|	---	|
|	downstream	|	SNP	|	11304	|
|	exonic	|	SNP	|	146890	|
|	exonic;splicing	|	SNP	|	5	|
|	intergenic	|	SNP	|	42939	|
|	intronic	|	SNP	|	61036	|
|	ncRNA_exonic	|	SNP	|	4	|
|	splicing	|	SNP	|	439	|
|	upstream	|	SNP	|	12193	|
|	upstream;downstream	|	SNP	|	3890	|
|	UTR3	|	SNP	|	17112	|
|	UTR5	|	SNP	|	14283	|
|	UTR5;UTR3	|	SNP	|	58	|

<img src="./results/03.Annotation\Statistics/Tae_BSR.filter.SNPs.distribution.png" width="600px" alt="SNP分布">

|Type|Number|
|---|---|
|nonsynonymous SNV|54402|
|stopgain|546|
|stoploss|130|
|synonymous SNV|91817|


<img src="./results/03.Annotation\Statistics/Tae_BSR.filter.SNPs.anno.exonic_variant_function.png" width="600px" alt="SNP类型">

InDel则是在intronic上分布最多。且nonframeshift insertion/deletion有更高的比例，只有很少InDel会造成stopgain/stoploss。

|	Distribution	|	varType	|	number	|
|	---	|	---	|	---	|
|	downstream	|	INDEL	|	1764	|
|	exonic	|	INDEL	|	4524	|
|	exonic;splicing	|	INDEL	|	1	|
|	intergenic	|	INDEL	|	3661	|
|	intronic	|	INDEL	|	10255	|
|	ncRNA_exonic	|	INDEL	|	2	|
|	splicing	|	INDEL	|	87	|
|	upstream	|	INDEL	|	2300	|
|	upstream;downstream	|	INDEL	|	670	|
|	UTR3	|	INDEL	|	3229	|
|	UTR5	|	INDEL	|	4475	|
|	UTR5;UTR3	|	INDEL	|	16	|

<img src="./results/03.Annotation\Statistics/Tae_BSR.filter.INDELs.distribution.png" width="600px" alt="InDel分布">

|Type|Number|
|---|---|
|frameshift deletion|727|
|frameshift insertion|911|
|nonframeshift deletion|1558|
|nonframeshift insertion|1265|
|stopgain|52|
|stoploss|12|

<img src="./results/03.Annotation\Statistics/Tae_BSR.filter.INDELs.anno.exonic_variant_function.png" width="600px" alt="INDEL类型">

### BSR分析
对SNP进一步筛选过后，首先统计SNP在染色体上的分布，由下图可以看出，SNP分布并不十分均匀，且与基因在染色体上的分布较为一致。随后分别计算两个混池的SNP index并相减得到Delta SNP index值，在chr2上存在超出99%置信区间的peak，由于部分区域的SNP个数较少造成delta SNP index值波动较大，建议将chr2上多个显著区间合并为一个在后续通过分子标记进行精细定位。同时在chr4上也存在部分区域超出99%置信区间，但超出程度微弱，若在F2分离群体中通过分离比确定为单基因控制，则可以不予考虑后续验证。

<img src="./results/04.Analysis/Tae_BSR.SNP_distribution_histogram.png" width="800px" alt="SNP分布">

![delta SNP index dotplot](./results/04.Analysis/Tae_BSR.delta_SNP_index.png)

![delta SNP index lineplot 99CI](./results/04.Analysis/Tae_BSR.deltaSNPindex.99CI.png)

### Primer设计
对于鉴定到的InDels进行primer设计（PAGE垂直电泳）用于后续精细定位，根据实验室现有PAGE垂直电泳分辨率选择合适差异长度的引物，同时应注意需要根据序列比对和群体验证以确定其特异性。

<div id="ref"></div>

## 参考文献
- 1. Sanzhen Liu, Cheng-Ting Yeh, Ho Man Tang, Dan Nettleton, Patrick S. Schnable. Gene Mapping via Bulked Segregant RNA-Seq (BSR-Seq)[J]. PLoS ONE, 2012, 7(5).
- 2. Giovannoni James J., Wing Rod A., Ganal Martin W., Tanksley Steven D..  Isolation of molecular markers from specific chromosomal intervals using DNA pools from existing mapping populations[J]. Narnia, 1991, 19(23).
- 3. R. W. Michelmore, I. Paran,R. V. Kesseli. Identification of Markers Linked to Disease-Resistance Genes by Bulked Segregant Analysis: A Rapid Method to Detect Markers in Specific Genomic Regions by Using Segregating Populations[J]. Proceedings of the National Academy of Sciences of the United States of America, 1991, 88(21).
- 4. Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ preprocessor[J]. Bioinformatics, 2018, 34(17).
- 5. Dobin Alexander,Davis Carrie A,Schlesinger Felix,Drenkow Jorg,Zaleski Chris,Jha Sonali,Batut Philippe,Chaisson Mark,Gingeras Thomas R. STAR: ultrafast universal RNA-seq aligner.[J]. Bioinformatics (Oxford, England), 2013, 29(1).
- 6. “Picard Toolkit.” 2019. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/; Broad Institute
- 7. Aaron McKenna, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian Cibulskis, Andrew Kernytsky, Kiran Garimella, David Altshuler, Stacey Gabriel, Mark Daly, Mark A. DePristo. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data[J]. Cold Spring Harbor Laboratory Press, 2010, 20(9).
- 8. Yan Melissa Y,Ferguson Betsy,Bimber Benjamin N. VariantQC: a visual quality control report for variant evaluation.[J]. Bioinformatics (Oxford, England), 2019, 35(24).
- 9. Kai Wang, Mingyao Li and Hakon Hakonarson. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data.[J]. Nucleic Acids Research, 2010, 38(16).
- 10. Ben N. Mansfeld, Rebecca Grumet. QTLseqr: An R Package for Bulk Segregant Analysis with Next‐Generation Sequencing[J]. The Plant Genome, 2018, 11(2).
- 11. Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi. QTL ‐seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations[J]. The Plant Journal, 2013, 74(1).
- 12. Untergasser Andreas, Cutcutache Ioana, Koressaar Triinu, Ye Jian, Faircloth Brant C., Remm Maido, Rozen Steven G.. Primer3—new capabilities and interfaces[J]. Narnia, 2012, 40(15).
