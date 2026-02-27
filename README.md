# PacBio HiFi + 单染色体整合

# 全染色体级别 Haplotype 构建与甲基化分型流程（六步法）

------------------------------------------------------------------------

# 原始数据说明与完整性验证

## 数据来源

https://www.ebi.ac.uk/ena/browser/view/ERR14732851

下载文件（Submitted files）：

ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14732851/m54329U_210326_192251.hifi_reads.bam

------------------------------------------------------------------------

## 文件完整性校验

官方 MD5：

c5da0cca759e36d620c5986a845c0e07

本地校验：

``` bash
md5sum m54329U_210326_192251.hifi_reads.bam
```

输出：

c5da0cca759e36d620c5986a845c0e07 m54329U_210326_192251.hifi_reads.bam

结论：

MD5 值一致，说明文件在下载和传输过程中没有任何损坏。

后续将基于：
-   HiFi BAM 中的 MM / ML 标签
-   或 PacBio 官方 methylation 工具

实现 haplotype-resolved methylation 分析。

------------------------------------------------------------------------

# Step1：HiFi 数据 Mapping 与排序

## 0. 安装 pbmm2

## 1. 建立索引（只需运行一次）

``` bash
pbmm2 index hg38.fa hg38.mmi --preset HIFI
```

## 2. 进行比对

``` bash
pbmm2 align /home/student/data-assemble/human_ref/hg38.mmi \
m54329U_210326_192251.hifi_reads.bam \
m54329U_hg38.unsorted.bam \
--preset HIFI \
-j 48
```

## 3. 排序

``` bash
samtools sort -@ 12 -m 8G \
m54329U_hg38.unsorted.bam \
-o m54329U_hg38.aligned.sorted.bam
```

## 4. 仅提取 chr1 测试

``` bash
samtools index m54329U_hg38.aligned.sorted.bam

samtools view -@ 8 -b \
m54329U_hg38.aligned.sorted.bam chr1 \
> m54329U_hg38.aligned.sorted_chr1.bam

samtools index -@ 8 m54329U_hg38.aligned.sorted_chr1.bam
```

## 5. 提取 mapping 前9列

``` bash
samtools view m54329U_hg38.aligned.sorted_chr1.bam | cut -f1-9 > first9_fields_chr1.txt
```
用途：
-   分析 FLAG / MAPQ / CIGAR
-   了解 reads 在 reference 上的对齐消耗情况



------------------------------------------------------------------------

# Step2：DeepVariant 变异检测

## 安装 DeepVariant（必须科学上网）

``` bash
sudo docker pull google/deepvariant:1.9.0
```

## 运行 DeepVariant

``` bash
BAM="/media/jason9/disk2/pacBioStudy/data/chr1_only/m54329U_hg38.aligned.sorted_chr1.bam"
REF="/media/jason9/disk2/pacBioStudy/data/chr1_only/chr1.fa"
OUTDIR="/media/jason9/disk2/pacBioStudy/data/chr1_only/deepvariant_out"
SAMPLE="NA12878"
THREADS=30

sudo docker run --rm \
  -v /:/host \
  google/deepvariant:1.9.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref="/host${REF}" \
  --reads="/host${BAM}" \
  --output_vcf="/host${OUTDIR}/${SAMPLE}.chr1.dv.vcf.gz" \
  --output_gvcf="/host${OUTDIR}/${SAMPLE}.chr1.dv.g.vcf.gz" \
  --num_shards="${THREADS}"
```

## 数据情况分析（非必须）

``` bash
#统计 SNP / INDEL 数量
bcftools stats NA12878.chr1.dv.vcf.gz > chr1.stats.txt
grep -E "number of SNPs|number of indels" chr1.stats.txt

#平均深度
bcftools query -f '[%DP\n]' NA12878.chr1.dv.vcf.gz | awk '{sum+=$1} END {print sum/NR}'

#平均质量
bcftools query -f '[%GQ\n]' NA12878.chr1.dv.vcf.gz | awk '{sum+=$1} END {print sum/NR}'
```

## 过滤高质量杂合位点 （建议使用高质量的hetero，做phasing）

``` bash
bcftools view \
-i 'FILTER="PASS" && GT="het" && FMT/GQ>=30 && FMT/DP>=6' \
-Oz -o NA12878.chr1.dv.het.GQ30.DP6.PASS.vcf.gz \
NA12878.chr1.dv.vcf.gz
```
## NA12878金标数据核对一致性（非必须）
``` bash
#提取数据
IN=/mnt/DataBaseShare/KGP/RawData/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
bcftools view -s NA12878 -Oz -o KGP.NA12878.chr1.vcf.gz "$IN"
bcftools index -t KGP.NA12878.chr1.vcf.gz
#python代码实现比对
cd goldCheck 
bash goldConsistentCheck.sh #自己写代码比
```
------------------------------------------------------------------------

# Step3：WhatsHap Read-backed Phasing

## 安装 WhatsHap

``` bash
pip install whatshap
```
## 运行 WhatsHap(示例代码)
``` bash
whatshap phase \
  --reference hg38.fa \
  --ignore-read-groups \
  -o chr1.phased.vcf.gz \
  chr1.vcf.gz \
  chr1.bam
```

输出：

-   0\|1 / 1\|0 phased genotype\
-   PS（phase set）标签

------------------------------------------------------------------------

# Step4：Reads 定边

## 运行 WhatsHap haplotag(示例代码)
``` bash
whatshap haplotag \
--reference hg38.fa \
--ignore-read-groups \
chr1.phased.vcf.gz \
chr1.bam \
-o chr1.haplotagged.bam
```

为每条 read 添加：

-   HP:i:1 属于 haplotype 1\
-   HP:i:2 属于 haplotype 2\
-   PS:i:xxxxx 所属 phase block

------------------------------------------------------------------------

# Step5：单染色体数据叠加，实现全染色体串接

核心思想：

1.  对每一个 haplotype片段(一对) 找到与单染色体 VCF 重叠 SNP\
2.  计算 hap1 / hap2 与单染色体 allele 的一致比例\
3.  一致则串接

结果：

-   串接所有haplotype片段(一对)
-   跨越中心粒与重复区
-   生成全染色体 hapA / hapB

------------------------------------------------------------------------

# Step6：按 QNAME 分组原始 Reads，实现甲基化分型

根据 haplotag 结果，将原始 HiFi reads 分为：

-   hap1.reads.bam\
-   hap2.reads.bam

随后：

-   提取 MM / ML 标签\
-   构建 haplotype-specific methylome

------------------------------------------------------------------------

# 技术意义

本六步流程实现：

1.  HiFi mapping
2.  高精度变异检测
3.  局部 read-backed phasing
4.  read-level haplotag
5.  单染色体锚定跨 block 串接
6.  haplotype-resolved methylation 分析

该方法可实现 chromosome-scale haplotype methylome 构建。
