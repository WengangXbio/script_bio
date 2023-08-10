### Short-reads-based validation
```
for sample_name in `cat list`  ; do
#sample_name=SRR7120167
samtools  view -L extract_PAB.bed /home/s1874451/Datastore/WZ/Liuyang/1.bwa2.bam/${sample_name}.sorted.bam  -o ${sample_name}.sorted.PAB.bam
bedtools bamtofastq  -i ${sample_name}.sorted.PAB.bam -fq ${sample_name}.sorted.PAB.fastq
awk '{print (NR%4 == 1) ? "@1_" ++i : $0}' ${sample_name}.sorted.PAB.fastq > ${sample_name}.sorted.PAB.rename.fastq
minimap2 -a --splice -O6,24 -B2 -G500  chrY1.PAB.fa ${sample_name}.sorted.PAB.rename.fastq | samtools sort -o ${sample_name}.PAB.sort.bam -
samtools index ${sample_name}.PAB.sort.bam
cp ${sample_name}.PAB.sort.bam*  ~/Datastore/WZ/dropbox/
done
```

### Long-reads-based validation
```
minimap2 -cx map-ont PAB.fa long.fastq.gz >  long.sam
```
