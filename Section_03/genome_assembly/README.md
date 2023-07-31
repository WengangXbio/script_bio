## Genome assembly
### Flye assembly
```
flye --pacbio-raw pacbio_subreads.fa --genome-size 2.5g --out-dir ./flye_assembly --keep-haplotypes
```
### Falcon assembly, unzip and polish
```
fc_run fc_run.cfg
fc_unzip.py fc_unzip.cfg
```
### wtdbg2 assembly
```
wtdbg2.pl -t 16 -x rs -g 2500m -o ./wtdbg2_assembly pacbio_subreads.fa
```
### quickmerge
```
merge_wrapper.py falcon_assembly.fasta flye_assembly.fasta
```
## Genome assessment

### K-mer (merqury)
```
conda activate /home/s1874451/schoenebeck_group/WENGANG/anaconda_env/merqury

meryl k=21 count output LAB339608_DDPL03146-W_HTG5VCCXY_L1_1.meryl  LAB339608_DDPL03146-W_HTG5VCCXY_L1_1.chrY.fq
meryl k=21 count output LAB339608_DDPL03146-W_HTG5VCCXY_L1_2.meryl  LAB339608_DDPL03146-W_HTG5VCCXY_L1_2.chrY.fq
meryl union-sum output LAB339608_DDPL03146-W_HTG5VCCXY_L1.meryl LAB339608_DDPL03146-W_HTG5VCCXY_L1_*.meryl
merqury.sh LAB339608_DDPL03146-W_HTG5VCCXY_L1.meryl ROS_Cfam1_Y.fa chrY_PExrosy > chrY_PExrosy.log &
Rscript $MERQURY/plot/plot_spectra_cn.R -f chrY_PExrosy.ROS_Cfam1_Y.spectra-cn.hist -o chrY_PExrosy.ROS_Cfam1_Y.spectra-cn -m 100

meryl k=19 count output LAB339608_DDPL03146-W_HTG5VCCXY_L1_1.wg.meryl  LAB339608_DDPL03146-W_HTG5VCCXY_L1_1.clean.fq.gz
meryl k=19 count output LAB339608_DDPL03146-W_HTG5VCCXY_L1_2.wg.meryl  LAB339608_DDPL03146-W_HTG5VCCXY_L1_2.clean.fq.gz
meryl union-sum output LAB339608_DDPL03146-W_HTG5VCCXY_L1.wg.meryl *.wg.meryl
merqury.sh LAB339608_DDPL03146-W_HTG5VCCXY_L1.wg.meryl ROS_Cfam1.fa wg_PExRoscfam 
Rscript $MERQURY/plot/plot_spectra_cn.R -f wg_PExRoscfam.ROS_Cfam1.spectra-cn.hist -o wg_PExRoscfam.ROS_Cfam1.spectra-cn -m 100
```
### ABySS short reads assembly
```
abyss-pe name=${sample_name} k=55 in='${sample_name}.R1.clean.fq.gz ${sample_name}.R2.clean.fq.gz' j=10
```
### Mapping short-reads assembly against long-reads assembly and filtering single-copy regions' variants
```
minimap2 -x asm5 ROS_Cfam1.fa denovos/${aa}-8.fa > ${aa}/${aa}_abyss.ROS_Cfam1.paf
awk '($4-$3)/$2 >= 0.95 && $2 >= 1000 {print $0}' ${aa}/${aa}_abyss.ROS_Cfam1.paf | grep "chrY" |cut -f1| sort |uniq -c |awk -F' ' '$1==1 {print $2}' > ${aa}/${aa}_abyss.unique_Y.contigs
~/schoenebeck_group/WENGANG/WZ_software/seqtk/seqtk subseq  denovos/${aa}-8.fa ${aa}/${aa}_abyss.unique_Y.contigs > ${aa}/${aa}_abyss.unique_Y.contigs.fa 
minimap2 -ax asm5  ROS_Cfam1_Y.fa  ${aa}/${aa}_abyss.unique_Y.contigs.fa  |samtools sort -o ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam - 
samtools view ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam |grep chrY |cut -f1 |sort |uniq -c |awk -F ' ' '$1==1 {print $2}' > ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam.mapped.countigs
minimap2 -x asm5  ROS_Cfam1_Y.fa  ${aa}/${aa}_abyss.unique_Y.contigs.fa  > ${aa}/${aa}_abyss.unique_Y.contigs.xRosCfamY.paf
awk -f vlookup.awk ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam.mapped.countigs ${aa}/${aa}_abyss.unique_Y.contigs.xRosCfamY.paf |grep yes |cut -f6,8,9 |sort |uniq > ${aa}/${aa}_abyss.assembled.region
samtools index ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam
bcftools mpileup -O b -o ${aa}/${aa}_abyss.bcf  -f ROS_Cfam1_Y.fa ${aa}/${aa}_abyss.unique_Y.contigs.sorted.bam
bcftools view ${aa}/${aa}_abyss.bcf |awk '$5!="<*>" {print}' > ${aa}/${aa}_abyss.variants.vcf 
bedtools intersect -v -a ${aa}/${aa}_abyss.variants.vcf  -b duplication_region.cnvnator.cnv > ${aa}/${aa}_abyss.variants.rm_cnv.vcf
```
