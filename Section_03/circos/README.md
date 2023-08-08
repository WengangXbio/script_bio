## Circos plot
```
circos --conf circos.conf
```
### The following codes are for generating information for circos plot
```
#GC content
#cd /home/s1874451/schoenebeck_group/WENGANG/YYYY/sequence_feature/GC_content
cd ~/Datastore/WZ/Y_paper/GC_content/
conda activate  /home/s1874451/schoenebeck_group/WENGANG/anaconda_env/AGAT/
perl /home/s1874451/schoenebeck_group/WENGANG/WZ_software/GC_content_in_sliding_window/GC_content.pl -fasta chrY1.fa --window 1000 --step 500
perl /home/s1874451/schoenebeck_group/WENGANG/WZ_software/GC_content_in_sliding_window/GC_content.pl -fasta chrY2_1032020-1246049.fa --window 1000 --step 500
perl /home/s1874451/schoenebeck_group/WENGANG/WZ_software/GC_content_in_sliding_window/GC_content.pl -fasta chrY2_1-772310.fa --window 1000 --step 500
perl /home/s1874451/schoenebeck_group/WENGANG/WZ_software/GC_content_in_sliding_window/GC_content.pl -fasta ROS_Cfam1.PAR.fasta --window 1000 --step 500
perl /home/s1874451/schoenebeck_group/WENGANG/WZ_software/GC_content_in_sliding_window/GC_content.pl -fasta chrY3.fa --window 1000 --step 500
awk 'BEGIN{OFS="\t"} {print "chrY1",$2,$3,$4}' chrY1.fa_chrY1.GC_content > chrY1.fa_chrY1.GC_content.txt
awk 'BEGIN{OFS="\t"} {print "chrY2",$2,$3,$4}' chrY2_1-772310.fa_chrY2:1-772310.GC_content > ROS_Cfam1_Y.fa_chrY2_1.GC_content.txt
awk 'BEGIN{OFS="\t"} {print "chrY2",$2+1032020,$3+1032020,$4}' chrY2_1032020-1246049.fa_chrY2:1032020-1246049.GC_content > ROS_Cfam1_Y.fa_chrY2_2.GC_content.txt
awk 'BEGIN{OFS="\t"} {print "chrY3",$2,$3,$4}' chrY3.fa_chrY3.GC_content > chrY3.fa_chrY3.GC_content.txt
awk 'BEGIN{OFS="\t"} {print "chrX",$2,$3,$4}' ROS_Cfam1.PAR.fasta_chrX_PAR.GC_content > ROS_Cfam1.PAR.fasta_chrX_PAR.GC_content.txt
cat ROS_Cfam1.PAR.fasta_chrX_PAR.GC_content.txt chrY1.fa_chrY1.GC_content.txt ROS_Cfam1_Y.fa_chrY2_1.GC_content.txt \
ROS_Cfam1_Y.fa_chrY2_2.GC_content.txt chrY3.fa_chrY3.GC_content.txt > GC_content.circos.txt
awk 'BEGIN{OFS="\t"} {print "chrY1",$2,$3,$4}' chrY1.fa_chrY1.GC_deviation > chrY1.fa_chrY1.GC_deviation.txt
awk 'BEGIN{OFS="\t"} {print "chrY2",$2,$3,$4}' chrY2_1-772310.fa_chrY2:1-772310.GC_deviation > ROS_Cfam1_Y.fa_chrY2_1.GC_deviation.txt
awk 'BEGIN{OFS="\t"} {print "chrY2",$2+1032020,$3+1032020,$4}' chrY2_1032020-1246049.fa_chrY2:1032020-1246049.GC_deviation > ROS_Cfam1_Y.fa_chrY2_2.GC_deviation.txt
awk 'BEGIN{OFS="\t"} {print "chrY3",$2,$3,$4}' chrY3.fa_chrY3.GC_deviation > chrY3.fa_chrY3.GC_deviation.txt
awk 'BEGIN{OFS="\t"} {print "chrX",$2,$3,$4}' ROS_Cfam1.PAR.fasta_chrX_PAR.GC_deviation > ROS_Cfam1.PAR.fasta_chrX_PAR.GC_deviation.txt
cat ROS_Cfam1.PAR.fasta_chrX_PAR.GC_deviation.txt chrY1.fa_chrY1.GC_deviation.txt ROS_Cfam1_Y.fa_chrY2_1.GC_deviation.txt \
ROS_Cfam1_Y.fa_chrY2_2.GC_deviation.txt chrY3.fa_chrY3.GC_deviation.txt > GC_deviation.circos.txt



#TE annotation
cat ROS_Cfam1.fa.mod.EDTA.TEanno.PAR.gff3 ROS_Cfam1.fa.mod.EDTA.TEanno.chrY.gff3 |awk '{print $1,$3,$4,$5,$9}' - |awk -F '[; =]' '{print $1,$2,$3,$4,$8,$14}' > ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf
module load igmm/apps/R/3.6.0
gt=as.matrix(read.table("ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf",head=F))
colnames(gt)=c("chr","type","start","end","ontology","identity")
df=data.frame(chr=0,type=0,start=0,end=0,ontology=0,identity=0)
for(i in 1:nrow(gt)){
if(!length(intersect(c(df[nrow(df),3]:df[nrow(df),4]),c(gt[i,3]:gt[i,4])))){
df=rbind(df,gt[i,])
}else{
if(c(as.numeric(df[nrow(df),4])-as.numeric(df[nrow(df),3]))<c(as.numeric(gt[i,4])-as.numeric(gt[i,3]))){
df[nrow(df),]=gt[i,]
}
}
}
df=df[-1,]
write.table(df,"ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal", sep="\t", row.names=FALSE, quote=FALSE)

module load roslin/blast+/2.9.0
module load igmm/apps/BEDTools/2.27.1
makeblastdb -in NonLTR_new.fa -dbtype nucl
awk 'BEGIN{OFS="\t"} {print $1,$3,$4,$5"_"$2"_"$1"_"$3}' ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal |tail -n +2 > ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.bed
grep -v "Simple_repeat\|Low_complexity" ./RM2_ROS_Cfam1.upY.fa_1609783087/RM2_ROS_Cfam1.upY.fa_1609783087.out |tail -n +4 |awk 'BEGIN{OFS="\t"}{print $5,$6,$7}' > RM2_out_chrY.bed
grep -v "Simple_repeat\|Low_complexity" ./RM2_ROS_Cfam1.PAR.fasta_1652293634/RM2_ROS_Cfam1.PAR.fasta_1652293634.out |tail -n +4 |awk 'BEGIN{OFS="\t"}{print $5,$6,$7}' > RM2_out_PAR.bed
cat RM2_out_PAR.bed  RM2_out_chrY.bed > RM2_out.chrY_PAR.bed
bedtools subtract -a ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.bed -b RM2_out.chrY_PAR.bed -A > ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.bed
awk '$3-$2>100' ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.bed > ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.bed
bedtools getfasta -fi ~/schoenebeck_group/WENGANG/Liuyang_cnv/ROS_Cfam/ROS_Cfam1.fa -bed ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.bed -name > ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.fa
blastn -task blastn-short  -db NonLTR_new.fa -evalue 1e-5 -query ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.fa  -max_hsps 1 -max_target_seqs 1 \
 -out ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.blast_results -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
grep -v "SINEC" ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.blast_results |cut -f1 > LINE.EDTA.list
grep    "SINEC" ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.blast_results |cut -f1 > SINE.EDTA.list
grep -f LINE.EDTA.list ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.bed > LINE.EDTA.bed
grep -f SINE.EDTA.list ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.bed > SINE.EDTA.bed
cat LINE.EDTA.list SINE.EDTA.list |grep -vf - ROS_Cfam1.fa.mod.EDTA.TEanno.chrY_PAR.gff3.mdf.removal.unannotated.length.bed > others.EDTA.bed
cat ./RM2_ROS_Cfam1.PAR.fasta_1652293634/RM2_ROS_Cfam1.PAR.fasta_1652293634.out ./RM2_ROS_Cfam1.upY.fa_1609783087/RM2_ROS_Cfam1.upY.fa_1609783087.out | grep -v "Simple_repeat\|Low_complexity" | tail -n +4 |grep "SINE" |awk 'BEGIN{OFS="\t"} {print $5,$6,$7}' > SINE.RM2.bed
cat ./RM2_ROS_Cfam1.PAR.fasta_1652293634/RM2_ROS_Cfam1.PAR.fasta_1652293634.out ./RM2_ROS_Cfam1.upY.fa_1609783087/RM2_ROS_Cfam1.upY.fa_1609783087.out | grep -v "Simple_repeat\|Low_complexity" | tail -n +4 |grep "LINE" |awk 'BEGIN{OFS="\t"} {print $5,$6,$7}' > LINE.RM2.bed
cat ./RM2_ROS_Cfam1.PAR.fasta_1652293634/RM2_ROS_Cfam1.PAR.fasta_1652293634.out ./RM2_ROS_Cfam1.upY.fa_1609783087/RM2_ROS_Cfam1.upY.fa_1609783087.out | grep -v "Simple_repeat\|Low_complexity" | tail -n +4 |grep -v "SINE\|LINE" |awk 'BEGIN{OFS="\t"} {print $5,$6,$7}' > others.RM2.bed
cat LINE.EDTA.bed LINE.RM2.bed |cut -f1-3 > LINE.bed
cat SINE.EDTA.bed SINE.RM2.bed |cut -f1-3 > SINE.bed
cat others.EDTA.bed others.RM2.bed |cut -f1-3 > others.bed
bedtools intersect -a region_10000.bed -b LINE.bed > region_10000.LINE.bed
bedtools intersect -a region_10000.bed -b SINE.bed > region_10000.SINE.bed
bedtools intersect -a region_10000.bed -b others.bed > region_10000.others.bed

library(tidyr)
line=read.table("region_10000.LINE.bed",head=F)
sine=read.table("region_10000.SINE.bed",head=F)
othe=read.table("region_10000.others.bed",head=F)
region=read.table("region_10000.bed",head=F)
names=unique(region$V4)
length1=c();length2=c();length3=c()
for(i in 1:length(names)){
yw1=line[which(line$V4==as.character(names[i])),]
length1=c(length1,sum(yw1$V3-yw1$V2))
yw2=sine[which(sine$V4==as.character(names[i])),]
length2=c(length2,sum(yw2$V3-yw2$V2))
yw3=othe[which(othe$V4==as.character(names[i])),]
length3=c(length3,sum(yw3$V3-yw3$V2))
}
region$"id=Line"=length1/(region$V3-region$V2)
region$"id=Sine"=length2/(region$V3-region$V2)
region$"id=Other"=length3/(region$V3-region$V2)
out=region[,c(1,2,3,5,6,7)]
write.table(out,"TE_circos.hist.txt", sep="\t", row.names=FALSE,col.names=FALSE, quote=FALSE)

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4","$5","$6}' TE_circos.hist.txt > TE_circos.hist.txt1





#CNVnator
cd /home/s1874451/Datastore/WZ/rosy_1.0_figure1/bn_op_gap
grep "NC_051843.1" LAB339608_DDPL03146-W_HTG5VCCXY_L1.raw.cnv |cut -f2,4 |awk -F '[:-]' 'BEGIN{OFS="\t"} $2<6590800 {print "chrX",$2,$3,$4}' > LAB339608.PAR.raw.cnv.depth
cd /home/s1874451/schoenebeck_group/WENGANG/YYYY/sequence_feature/cnvnator/LAB_for_circus
cat LAB339608.PAR.raw.cnv.depth LAB339608.chrY.raw.cnv.depth > LAB339608.raw.cnv.depth

module load igmm/apps/R/3.6.0
cnv=read.table("LAB339608.raw.cnv.depth",head=F)
cnvX=cnv[which(cnv$V1=="chrX"),]
cnvX1=data.frame(V1="chrX",V2=head(c(1,cnvX$V3+1),-1),V3=cnvX$V2-1,V4=1)
cnvY1=cnv[which(cnv$V1=="chrY1"),]
cnvY11=data.frame(V1="chrY1",V2=head(c(1,cnvY1$V3+1),-1),V3=cnvY1$V2-1,V4=1)
cnvY2=cnv[which(cnv$V1=="chrY2"),]
cnvY21=data.frame(V1="chrY2",V2=head(c(1,cnvY2$V3+1),-1),V3=cnvY2$V2-1,V4=1)
cnvY3=cnv[which(cnv$V1=="chrY3"),]
cnvY31=data.frame(V1="chrY3",V2=head(c(1,cnvY3$V3+1),-1),V3=cnvY3$V2-1,V4=1)
cnva=rbind(cnv,cnvX1,cnvY11,cnvY21,cnvY31)
cnva=cnva[-which(cnva$V3-cnva$V2<500),]
cnva=cnva[with(cnva, order(V1, V2)), ]
write.table(cnva,"LAB339608.raw.cnv.depth.circus", sep="\t", row.names=FALSE, quote=FALSE)


#Similiarity
module load igmm/apps/BEDTools/2.27.1
bedtools getfasta -name -fi ~/schoenebeck_group/WENGANG/Liuyang_cnv/ROS_Cfam/ROS_Cfam1.fa -bed 2000mer.500slide.bed > 2000mer.500slide.fasta
/home/s1874451/schoenebeck_group/WENGANG/WZ_software/minimap2-2.24_x64-linux/minimap2 -x map-hifi -t 8 ROS_Cfam1.exclude_y.fa 2000mer.500slide.fasta > 2000mer.500slide.xRosCfam_excludeY.paf
awk '$6=="chrX" {print $1,($4-$3)/$2}' 2000mer.500slide.xRosCfam_excludeY.paf > 2000mer.500slide.chrX.similarity.txt
awk '$6!="chrX" {print $1,($4-$3)/$2}' 2000mer.500slide.xRosCfam_excludeY.paf > 2000mer.500slide.auto.similarity.txt
head  similarity.stat
tail -n +2 similarity.stat |awk 'BEGIN{OFS="\t"} {print $2,$3-749,$3+750,$4}' > similarity_chrX.circos.txt
tail -n +2 similarity.stat |awk 'BEGIN{OFS="\t"} {print $2,$3-749,$3+750,$5}' > similarity_autosome.circos.txt

module load igmm/apps/BEDTools/2.27.1
bedtools getfasta -name -fi ~/schoenebeck_group/WENGANG/Liuyang_cnv/ROS_Cfam/ROS_Cfam1.fa -bed 5000mer.2000slide.bed > 5000mer.2000slide.fasta
/home/s1874451/schoenebeck_group/WENGANG/WZ_software/minimap2-2.24_x64-linux/minimap2 -x map-hifi -t 8 ROS_Cfam1.exclude_y.fa 5000mer.2000slide.fasta > 5000mer.2000slide.xRosCfam_excludeY.paf
awk '$6=="chrX" {print $1,($4-$3)/$2}' 5000mer.2000slide.xRosCfam_excludeY.paf > 5000mer.2000slide.chrX.similarity.txt
awk '$6!="chrX" {print $1,($4-$3)/$2}' 5000mer.2000slide.xRosCfam_excludeY.paf > 5000mer.2000slide.auto.similarity.txt
tail -n +2 5000mer.2000slide.similarity.stat |awk 'BEGIN{OFS="\t"} {print $2,$3-999,$3+1000,$4}' > 5000mer.2000slide.similarity_chrX.circos.txt
tail -n +2 5000mer.2000slide.similarity.stat |awk 'BEGIN{OFS="\t"} {print $2,$3-999,$3+1000,$5}' > 5000mer.2000slide.similarity_autosome.circos.txt

module load roslin/blast+/2.9.0
makeblastdb -in 2000mer.500slide.fasta -dbtype nucl
blastn  -db 2000mer.500slide.fasta -evalue 1e-8 -query 2000mer.500slide.fasta -out 2000mer_500slide.self.blast_results -outfmt 6
awk '$4>1000&&$1!=$2' 2000mer_500slide.self.blast_results |awk -F'[_"\t"]' 'BEGIN{OFS="\t"} {print $0,$1,$2,$3,$4,$4-$2}' \
|awk '!($13==$15 && $17<2000 && $17>-2000)' |awk '!seen[$1]++' |awk '{print $1,$3*$4/200000}' > 2000mer_500slide.self.blast_results.greater_than_1000.txt
awk '$4>1500&&$1!=$2' 2000mer_500slide.self.blast_results |awk -F'[_"\t"]' 'BEGIN{OFS="\t"} {print $0,$1,$2,$3,$4,$4-$2}' \
|awk '!($13==$15 && $17<2000 && $17>-2000)' |awk '!seen[$1]++' |awk '{print $1,$3*$4/200000}' > 2000mer_500slide.self.blast_results.greater_than_1500.txt
awk '$4>1500&&$1!=$2' 2000mer_500slide.self.blast_results |awk -F'[_"\t"]' 'BEGIN{OFS="\t"} {print $0,$1,$2,$3,$4,$4-$2}' |awk '!($13==$15 && $17<2000 && $17>-2000)' \
|awk 'BEGIN{OFS="\t"} {print $13,$14-499,$14+500,$15,$16-499,$16+499}' > 2000mer_500slide.self.blast_results.greater_than_1500.link.circos.txt
blastn  -db NonLTR_new.fa -evalue 1e-8 -query 2000mer.500slide.fasta -out 2000mer.500slide.NonLTR.blast_results -outfmt 6
awk '$4>1500' 2000mer.500slide.NonLTR.blast_results |awk -F'[_"\t"]' 'BEGIN{OFS="\t"} {print $1,$2-499,$2+500}' |uniq > 2000mer.500slide.NonLTR.reads
grep -vwf 2000mer.500slide.NonLTR.reads 2000mer_500slide.self.blast_results.greater_than_1500.link.circos.txt > 2000mer_500slide.self.blast_results.greater_than_1500.link_non-NonLTR.circos.txt
grep -wf 2000mer.500slide.NonLTR.reads 2000mer_500slide.self.blast_results.greater_than_1500.link.circos.txt > 2000mer_500slide.self.blast_results.greater_than_1500.link_NonLTR.circos.txt


awk '$4>1500&&$1!=$2' 2000mer_500slide.self.blast_results |awk -F'[_"\t"]' 'BEGIN{OFS="\t"} {print $0,$1,$2,$3,$4,$4-$2}' |awk '!($13==$15 && $17<2000 && $17>-2000)' \
|awk '$3>=90' |awk 'BEGIN{OFS="\t"} {print $13,$14-499,$14+500,$15,$16-499,$16+500}' > 2000mer_500slide.self.blast_results.greater_than_1500.link.circos.txt


#Palindromes
cd /home/s1874451/schoenebeck_group/WENGANG/YYYY/sequence_feature/EDTA_TE
cat ./RM2_ROS_Cfam1.PAR.fasta_1652293634/RM2_ROS_Cfam1.PAR.fasta_1652293634.out ./RM2_ROS_Cfam1.upY.fa_1609783087/RM2_ROS_Cfam1.upY.fa_1609783087.out | grep "Simple_repeat\|Low_complexity"  |awk 'BEGIN{OFS="\t"} {print $5,$6,$7}' > ../Palindromes/RM2.Simple_repeat.bed
cd /home/s1874451/schoenebeck_group/WENGANG/YYYY/sequence_feature/Palindromes
cat chrX_PAB_50mer.bed chrY_50mer.bed > 50mer.bed
bedtools subtract -a 50mer.bed -b RM2.Simple_repeat.bed > 50mer.Simple_repeat_removal.bed
```
