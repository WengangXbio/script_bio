## Short variant discovery (SNVs + Indels)
### pre-processing bam file
```
gatk --java-options "-Xmx4G" AddOrReplaceReadGroups \
-I ${sample_name}.bam \
-O ${sample_name}.sorted.RG.bam \
-PL ILLUMINA \
-LB lib1 \
-PU unit1 \
-SM ${sample_name}

gatk --java-options "-Xmx4G" BaseRecalibrator \
--tmp-dir ./tmp \
-I ${sample_name}.sorted.RG.bam \
-R ROS_Cfam1.fa \
--known-sites RosCfam1.known-sites.bed.gz \
-O ${sample_name}.recal_data.table

gatk --java-options "-Xmx4G" ApplyBQSR \
--tmp-dir ./tmp \
-I ${sample_name}.sorted.RG.bam \
-R ROS_Cfam1.fa \
-O ${sample_name}.sorted.RG.bqsr.bam \
--bqsr-recal-file ${sample_name}.recal_data.table \
--preserve-qscores-less-than 6 \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30
```
### Generate gvcf (diploid model for autosomes)
```
gatk --java-options "-Xmx4g" HaplotypeCaller \
--tmp-dir ./tmp \
-R ROS_Cfam1.fa \
-I ${sample_name}.sorted.RG.bqsr.bam \
-O ${sample_name}.GATK.genome.vcf.gz \
--native-pair-hmm-threads 8 \
-ERC GVCF
```
### Generate gvcf (haploid model for Y chromosome)
```
gatk --java-options "-Xmx4g" HaplotypeCaller \
--tmp-dir ./tmp \
-R ROS_Cfam1.fa \
-I ${sample_name}.sorted.RG.bqsr.bam \
-O ${sample_name}.chrY.GATK.genome.vcf.gz \
--native-pair-hmm-threads 8 \
-ploidy 1 \
-ERC GVCF
```
### joint-call variants
```
gatk --java-options "-Xmx32G" GenomicsDBImport \
--genomicsdb-workspace-path ./genomicsdb_autosome \
--batch-size 50 \
--sample-name-map sample_list \
--consolidate \
--tmp-dir ./tmp

gatk --java-options "-Xmx32G" GenotypeGVCFs \
-R ROS_Cfam1.fa \
-V gendb://genomicsdb_autosome \
-O autosome.vcf.gz \
-all-sites \
--tmp-dir ./tmp
```
### Variants filtering
```
gatk SelectVariants \
    -V autosome.vcf.gz \
    -select-type SNP \
    -O autosome.SNP.vcf.gz
gatk VariantFiltration \
    -V autosome.SNP.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O autosome.SNP.snp_filter.vcf.gz
	
gatk SelectVariants \
    -V autosome.vcf.gz \
    -select-type INDEL \
    -O autosome.INDEL.vcf.gz
	
gatk VariantFiltration \
    -V autosome.INDEL.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O autosome.INDEL.snp_filter.vcf.gz
```
Note: QualByDepth (QD) is the variant confidence; FisherStrand (FS) is the Phred-scaled probability; StrandOddsRatio (SOR) is the strand bias using symmetric odds ratio test; RMSMappingQuality (MQ) is the root mean square mapping quality over all the reads at the site; MappingQualityRankSumTest (MQRankSum) is the u-based z-approximation from the Rank Sum Test for mapping qualities; ReadPosRankSumTest (ReadPosRankSum) is the u-based z-approximation from the Rank Sum Test for site position within reads. 

Details see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
