## Call variants based on short-read sequencing
### Quality control of short reads using fastp
```
fastp \
-i ${sample_name}_1.fastq.gz \
-I ${sample_name}_2.fastq.gz \
-o ${sample_name}_1.clean.fq.gz \
-O ${sample_name}_2.clean.fq.gz \
-h ${sample_name}.html \
-w 10
```
### Align short reads on the modified RosCfam using bwa-mem2
```
bwa-mem2 mem \
-R "@RG\tID:id_${sample_name}\tSM:${sample_name}\tLB:lib1" -t 8 \
${ref_fa} \
${wk_dir}/${sample_name}.R1.clean.fq.gz \
${wk_dir}/${sample_name}.R2.clean.fq.gz | samtools sort -o ${sample_name}.sorted.bam -
samtools index ${sample_name}.sorted.bam
```
### Calculation of sequencing depth by chromosome-based
```
mosdepth -n --fast-mode --by 100000 ${sample_name} ${wk_dir}/1.bwa2.bam/${sample_name}.sorted.bam
```
### Estimation of sex
Ratio of Y chromosome (single-copy region, chrY1 in this analysis) and autosome is calculated

1.2> ratio >0.8   -> male

0.05 > ratio      -> female
```
chrY1_cov=$(awk '$1=="chrY1" {print $4}' ${sample_name}.mosdepth.summary.txt)
chrX_cov=$(awk '$1=="chrX" {print $4}' ${sample_name}.mosdepth.summary.txt)
total_cov=$(awk '$1=="total" {print $4}' ${sample_name}.mosdepth.summary.txt)
ratioYX=`echo "scale=2;$chrY1_cov / $chrX_cov"|bc`
Yratio=`echo "scale=2;$chrY1_cov / $total_cov"|bc`
```

### Generation of gvcf for each sample
```
configureStrelkaGermlineWorkflow.py \
--bam ${sample_name}.sorted.bam  \
--referenceFasta $ref_fa \
--ploidy sk2_sex.$sample_name.vcf.gz \
--runDir sk2_$sample_name
sk2_$sample_name/runWorkflow.py -m local -j 8
```

### Joint-call gvcf into a cohort of variants (chrY as an example)
```
ls *.genome.S1.chrY.vcf.gz > chrY.gvcf.list
gvcfgenotyper/bin/gvcfgenotyper \
 -f ROS_Cfam1.fa \
 -l chrY.gvcf.list -Ob -o sk2.chrY.bcf
```

### Filtering out bad variants
```
#### Filter variants with 90% call rate and at least 3 alternative genotyping individuals
bcftools view  -i  ' COUNT(FMT/FT=="PASS" | FMT/GT=="0") >= 1048 && COUNT(FMT/FT=="PASS") >3 ' sk2.chrY.bcf |bcftools annotate -x INFO,^FORMAT/GT,FORMAT/FT,FORMAT/DP,FORMAT/GQ > sk2.chrY.09_CR.3_PASS.vcf
#### Filter bad sites (Unpass/ Pass ratio <0.8)
bcftools annotate -x INFO,^FORMAT/FT sk2.chrY.09_CR.3_PASS.vcf |grep -v "^#" |awk '{$1="";$2="";$3="";$4="";$5="";$6="";$7="";$8="";$9="";print $0}' |grep -o -n 'PASS' | cut -d : -f 1 | uniq -c > sk2.chrY.09_CR.3_PASS.PASS.count
bcftools annotate -x INFO,^FORMAT/FT sk2.chrY.09_CR.3_PASS.vcf |grep -v "^#" |awk '{$1="";$2="";$3="";$4="";$5="";$6="";$7="";$8="";$9="";print $0}' |grep -o -n '\.' | cut -d : -f 1 | uniq -c > sk2.chrY.09_CR.3_PASS.dot.count
Rscript  filter_snps.r 1165 sk2.chrY.09_CR.3_PASS.PASS.count sk2.chrY.09_CR.3_PASS.dot.count
grep "^#" sk2.chrY.09_CR.3_PASS.vcf > sk2.chrY.09_CR.3_PASS.vcf.header
grep -v "^#" sk2.chrY.09_CR.3_PASS.vcf | awk 'NR==FNR{l[$0];next;} !(FNR in l)' outlines.txt - |cat sk2.chrY.09_CR.3_PASS.vcf.header - > sk2.chrY.09_CR.3_PASS.filter.vcf
```
