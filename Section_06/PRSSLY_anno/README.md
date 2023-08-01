## Annotation of PRSSLY 
### Based on Iso-Seq data (Mice as an example)
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/093/SRR23170593/SRR23170593_subreads.fastq.gz
zcat SRR23170593_subreads.fastq.gz |sed -n '1~4s/^@/>/p;2~4p' - > SRR23170593_subreads.fa
diamond blastx --db ${prot_db} --query SRR23170593_subreads.fa  --sensitive -f 6 --out SRR23170593_subreads.blast --max-hsps 0 -c1 --strand plus
grep 'PRSSLY' SRR23170593_subreads.blast
### Extract identified subreads named "prssly_potential.fasta"
minimap2 -x map-pb -d mouse_ref.mmi mouse_ref.fa
minimap2 -t 8 -ax splice -uf --secondary=no -C5 -a  mouse_ref.mmi prssly_potential.fasta >  prssly_potential.sam   
samtools sort prssly_potential.sam   > prssly_potential.sort.bam
bedtools bamtobed -bed12 -i prssly_potential.sort.bam  > prssly_potential.sort.bed
### Convert bam to gtf using ucsc toolkits
bedToGenePred prssly_potential.sort.bed prssly_potential.sort.genepred
genePredToGtf "file" prssly_potential.sort.genepred prssly_potential.gtf
```

### Based on RNA-Seq data
```
### Annotating PRSSLY on Y sequences based on protein sequences
spaln -W -KP Ychromosome.fa
spaln -Q4 -O0 -M7 -d Ychromosome PRSSLY.pep > PRSSLY.Ychromosome.spaln.gtf
### Defining exon boundaries by RNA-Seq alignment
hisat2 -q -x ref_genome -p 10\
 -1 ${sample_name}_1.clean.fq.gz -2 ${sample_name}_2.clean.fq.gz \
 -S ${sample_name}.sam 2> ${sample_name}.hisat2.log
samtools sort -o ${sample_name}.sort.bam ${sample_name}.sam
stringtie -o ${sample_name}.gtf  ${sample_name}.sort.bam
### Both spaln-based and RNA-seq-based annotation is displayed on IGV, and PRSSLY annotation is generated manually.
```
