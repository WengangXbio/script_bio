## Iso-Seq annotation
### Processing Iso-Seq raw data
```
ccs movieX.subreads.bam movieX.ccs.bam --min-rq 0.9 --min-passes 1 --report-file ./  
lima --isoseq --dump-clips --peek-guess -j 24 movieX.ccs.bam primers.fasta movieX.fl.bam 
isoseq3 refine movieX.fl.primer_5p--primer_3p.bam primers.fasta movieX.flnc.bam
bamtools convert -format fasta -in movieX.flnc.bam > movieX.flnc.fa
python tama_flnc_polya_cleanup.py -f movieX.flnc.fa -p movieX.tama.flnc
```
### Mapping long reads and converting into annotation (GFF)
```
minimap2 -ax splice -t 10 -uf --secondary=no -C5 -O6,24 -B4 ref.fa SAMPLE.clustered.tama_cleanup.fa > SAMPLE.sam
sort -k 3,3 -k 4,4n SAMPLE.sam > SAMPLE.sorted.sam
conda activate cdna_cupcake
collapse_isoforms_by_sam.py --input SAMPLE.clustered.tama_cleanup.fa -s SAMPLE.sorted.sam --dun-merge-5-shorter -o SAMPLE
```

## RNA-Seq annotation
### 1.Alignment-based
```
hisat2-build RosCfam_1.0.fasta RosCfam_1.0
for sample_name in `cat samplelist`  ; do
fastp \
-i ${sample_name}_1.fastq.gz \
-I ${sample_name}_2.fastq.gz \
-o ${sample_name}_1.clean.fq.gz \
-O ${sample_name}_2.clean.fq.gz \
-h ${sample_name}.html \
-w 10
hisat2 -q -x RosCfam_1.0 -p 10\
 -1 ${sample_name}_1.clean.fq.gz -2 ${sample_name}_2.clean.fq.gz \
 -S ${sample_name}.sam 2> ${sample_name}.hisat2.log
samtools sort -o ${sample_name}.sort.bam ${sample_name}.sam
rm ${sample_name}.sam
done
stringtie -o ${sample_name}.gtf  ${sample_name}.sorted.bam 
stringtie --merge gtf_list --o merge.gtf
```
### 2.denovo-based
```
Trinity --seqType fq --max_memory 20G --left ${sample_name}.1.fastq --right ${sample_name}.2.fastq --CPU 6 --output ${sample_name}
minimap2 -ax splice -C5 ${chry} ${sample_name}_trinity.fasta | samtools sort -o ${sample_name}_trinity.RosCfam1.sort.bam - 
samtools view -bq 55 ${sample_name}_trinity.RosCfam1.sort.bam   > ${sample_name}_trinity.q55.sort.bam
samtools index ${sample_name}_trinity.q55.sort.bam
samtools bam2fq ${sample_name}_trinity.q55.sort.bam  |~/schoenebeck_group/WENGANG/WZ_software/seqkit fq2fa > ${sample_name}_trinity.q55.fa
diamond blastx --db ${prot_db} --query ${sample_name}_trinity.q55.fa  --sensitive -f 6 --out ${sample_name}_trinity.uniprot.result --max-hsps 0 -c1 --strand plus
Rscript extract_uniprot_coding.r ${sample_name}      
java -jar /home/s1874451/schoenebeck_group/WENGANG/WZ_software/picard.jar FilterSamReads  I=${sample_name}_trinity.q55.sort.bam  O=${sample_name}_trinity.q55.c50.uniprot_coding.sort.bam READ_LIST_FILE=${sample_name}_trinity.uniprot_coding.c50.extract FILTER=includeReadList
bedtools bamtobed -bed12 -i ${sample_name}_trinity.q55.c50.uniprot_coding.sort.bam  > ${sample_name}_trinity.q55.c50.uniprot_coding.sort.bed
bedToGenePred ${sample_name}_trinity.q55.c50.uniprot_coding.sort.bed ${sample_name}_trinity.q55.c50.uniprot_coding.sort.genepred
genePredToGtf "file" ${sample_name}_trinity.q55.c50.uniprot_coding.sort.genepred ${sample_name}_trinity.q55.c50.uniprot_coding.gtf
awk '$3=="exon"{print}' ${sample_name}_trinity.q55.c50.uniprot_coding.gtf > ${sample_name}_trinity.q55.c50.uniprot_coding.exon.gtf
agat_convert_sp_gxf2gxf.pl -g ${sample_name}_trinity.q55.c50.uniprot_coding.exon.gtf -o ${sample_name}_trinity.q55.c50.uniprot_coding.agat.gtf
```
## Liftover annotation
```
spaln -W -KP RosY.fa
spaln -Q4 -O0 -M7 -d RosY genes.pep > liftover_RosY.gtf
spaln -W -KD RosY.fa
spaln -Q4 -O0 -M7 -d RosY genes_mrna.fa > liftover_RosY.gtf
```
