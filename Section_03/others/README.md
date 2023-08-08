### Visulizing optical mapping
```
OMtools --viewrefin ref.fa --viewresin alignment
```
### cnvnator for estimating genome depth using short reads
```
mkdir -p ./3.cnvnator.bed/1.cnvnator_root/
mkdir -p ./3.cnvnator.bed/2.cnvnator_rawcnv/
mkdir -p ./3.cnvnator.bed/0.refChr/

chrs="$(seq -f 'chr%g' 1 38) chrMT chrX chrY1 chrY2 chrY3"
for chr in $chrs; do
        samtools faidx $ref_fa $chr > ./3.cnvnator.bed/0.refChr/$chr.fa
done

cnvnator -root ./3.cnvnator.bed/1.cnvnator_root/${sample_name}.root \
-tree ${sample_name}.sorted.bam \
-chrom $(seq -f 'chr%g' 1 38) chrMT chrY1 chrY2 chrY3
cnvnator -root ./3.cnvnator.bed/1.cnvnator_root/${sample_name}.root \
-his 100 \
-d ./3.cnvnator.bed/0.refChr
cnvnator -root ./3.cnvnator.bed/1.cnvnator_root/${sample_name}.root \
-stat 100
cnvnator -root ./3.cnvnator.bed/1.cnvnator_root/${sample_name}.root \
-partition 100
cnvnator -root ./3.cnvnator.bed/1.cnvnator_root/${sample_name}.root \
-call 100 > ./3.cnvnator.bed/2.cnvnator_rawcnv/${sample_name}.raw.cnv
```

### GC content estimation
```
perl gc_content.pl --fasta genome.fasta --window 1000 --step 100
```
###

### EDTA for TE annotation
```
perl EDTA.pl --genome ref.fa --species others --step all
```
### box plot by "ggpubr" package
```
library(ggpubr)
df=read.table("xxx",head=T,sep="\t")
ggboxplot(df, x = "Type", y = "values",color = "Type", palette = "jco", facet.by = "Comp",add = "jitter",size = 1)+  
rotate_x_text(angle = 45)+
stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(label.y =1.1)
```
