### Call variants for each sample (gvcf)
```
configureStrelkaGermlineWorkflow.py \
--bam ${wk_dir}/1.bwa2.bam/${sample_name}.sorted.bam  \
--referenceFasta $ref_fa \
--ploidy sk2_sex.$sample_name.vcf.gz \
--runDir sk2_$sample_name
sk2_$sample_name/runWorkflow.py -m local -j 8
```
### Joint-call variants on MSY regions
```
l |awk  '{print $9}' |awk -F'.' '{print $1}' > male.list
for sample_name in `cat male.list`  ; do
bcftools index ${sample_name}.genome.S1.vcf.gz -f
bcftools view ${sample_name}.genome.S1.vcf.gz --regions chrY1,chrY2,chrY3 -Oz -o ${sample_name}.genome.S1.chrY.vcf.gz
done
ls *.genome.S1.chrY.vcf.gz > chrY.gvcf.list
/home/s1874451/schoenebeck_group/WENGANG/WZ_software/gvcfgenotyper/bin/gvcfgenotyper \
 -f ~/schoenebeck_group/WENGANG/Liuyang_cnv/ROS_Cfam/ROS_Cfam1.fa \
 -l chrY.gvcf.list -Ob -o sk2.chrY.bcf
bcftools view  -i  ' COUNT(FMT/FT=="PASS" | FMT/GT=="0") >= 1048 && COUNT(FMT/FT=="PASS") >3 ' sk2.chrY.bcf |bcftools annotate -x INFO,^FORMAT/GT,FORMAT/FT,FORMAT/DP,FORMAT/GQ > sk2.chrY.09_CR.3_PASS.vcf
bcftools view -m2 -M2 -v snps sk2.chrY.09_CR.3_PASS.vcf -Ov -o sk2.chrY.09_CR.3_PASS.bSNPs.vcf
```
### Convert vcf to phylip format
```
subject=sk2.chrY1.09_CR.3_PASS.rm_cnv.rm_pseudo.maf0005.recode.bsnps
plink  --ped ../${subject}.ABR.ped --map ../${subject}.ABR.map --keep keep_ind.txt --recode --out ${subject}.ABR.subset --allow-extra-chr
plink  --ped ${subject}.ABR.subset.ped --map ${subject}.ABR.subset.map --maf 0.005 --recode --out ${subject}.ABR.subset.maf --allow-extra-chr
Rscript ../ped2fasta_haploid_c2.r ${subject}.ABR.subset.maf
module load python/2.7.10
python ../fasta2phylip.py sk2.chrY1.09_CR.3_PASS.rm_cnv.rm_pseudo.maf0005.recode.bsnps.ABR.subset.maf
```
### Run RAxML
```
raxmlHPC -s sk2.chrY1.09_CR.3_PASS.rm_cnv.rm_pseudo.maf0005.recode.bsnps.ABR.subset.maf.phy -n raxml.unroot.ml.n638.out -m ASC_GTRGAMMA -f a -x 123 -p 456 --asc-corr=lewis -N autoMRE
```
