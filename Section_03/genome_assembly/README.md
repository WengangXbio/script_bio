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
### ABySS short reads assembly
```
abyss-pe name=${sample_name} k=55 in='${sample_name}.R1.clean.fq.gz ${sample_name}.R2.clean.fq.gz' j=10
```
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
