### Aligning RNA-Seq reads on the modified RosCfam_1.0
```
for sample_name in `cat samplelist`  ; do
~/schoenebeck_group/WENGANG/WZ_software/fastp \
-i ${sample_name}_1.fastq.gz \
-I ${sample_name}_2.fastq.gz \
-o ${sample_name}_1.clean.fq.gz \
-O ${sample_name}_2.clean.fq.gz \
-h ${sample_name}.html \
-w 10

hisat2 -q -x /home/s1874451/schoenebeck_group/WENGANG/YYYY/rnaseq/ref/GCF_014441545.1.modified -p 10\
 -1 ${sample_name}_1.clean.fq.gz -2 ${sample_name}_2.clean.fq.gz \
 -S ${sample_name}.sam 2> ${sample_name}.hisat2.log
samtools sort -o ${sample_name}.sort.bam ${sample_name}.sam
rm ${sample_name}.sam
done
```
### Modify from Refseq to GFF for featureCounts
```
awk '$3=="gene" {print $9 }' Ros_Cfam.anno.GFF |awk -F '[;,]' '{print $2}' |awk -F '=' '{print $2}' > refseq_anno/Ros_Cfam.anno.gene_id
#gene="GeneID:481097"
for gene in `cat Ros_Cfam.anno.gene_id`  ; do
grep ${gene} ../Ros_Cfam.anno.GFF |awk -v var="$gene" 'BEGIN{OFS="\t"} $3=="exon" {print $1,$4,$5,var,0,$7}' > $gene.bed
sort -k1,1 -k2,2n $gene.bed > $gene.sort.bed
bedtools merge -i $gene.sort.bed  -c 6 -o distinct |awk -v var="$gene" 'BEGIN{OFS="\t"} {print $1,"WZ_mod","exon",$2,$3,0,$4,".","gene_name \""var"\";"}' > $gene.GFF
done
```
### nc_gene (if required)
```
awk '$3=="C_gene_segment" ||$3=="V_gene_segment" || $3=="lnc_RNA" || $3=="pseudogene" {print $9 }' Ros_Cfam.anno.GFF |awk -F '[;,]' '{print $2}' |awk -F '=' '{print $2}' > refseq_anno_ncgene/Ros_Cfam.anno.gene_id
for gene in `cat Ros_Cfam.anno.gene_id`  ; do
grep ${gene} ../Ros_Cfam.anno.GFF |awk -v var="$gene" 'BEGIN{OFS="\t"}  {print $1,$4,$5,var,0,$7}' > $gene.bed
sort -k1,1 -k2,2n $gene.bed > $gene.sort.bed
bedtools merge -i $gene.sort.bed  -c 6 -o distinct |awk -v var="$gene" 'BEGIN{OFS="\t"} {print $1,"WZ_mod","exon",$2,$3,0,$4,".","gene_name \""var"\";"}' > $gene.GFF
done
```
### Count RNA-Seq reads according to coding genes' exon
```
featureCounts -s 0 -t exon -g gene_name -p -B -a Ros_Cfam.anno.exon_collapse.add_Y_v3.GFF -o count_new.txt *.bam
Rscript tpm_calculator.r count_new.txt
cat tpm_calculator.r
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
```
### RNASEQ_chrY_TPM_by_class.tiff (headmap)
```
test_tpm=read.table("K:\\schoenebeck_group\\WZ\\dropbox\\counts.s2_tpm.txt",head=T)
condition=read.table("K:\\schoenebeck_group\\WZ\\dropbox\\sample.condition1",head=F)
test_tpm=test_tpm[c(36610:36645),]
test_tpm=as.matrix(test_tpm[,-c(1:6)])
samples=data.frame(V1=colnames(test_tpm))
yx =left_join(samples,condition,by="V1")
colnames(yx)=c("id","tissue","stage","class")
tpmm=data.frame(tpm=0,class=0,class0=0)
for(i in 1:ncol(test_tpm)){
add=data.frame(tpm=test_tpm[,i],class=paste(yx[i,4],"_",i,sep=""),class0=yx[i,4])
tpmm=rbind(tpmm,add)
}
tpmm=tpmm[-1,]
ggplot(data = tpmm, mapping = aes(x = class, y = tpm,color=class0)) +  geom_boxplot(alpha = 0,lwd=1) + 
scale_color_brewer(palette="Paired")+theme_bw()+ ylim(0,50)+
theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
                        axis.text.y = element_text(colour = "grey20", size = 12),
                        strip.text = element_text(face = "italic"),
                        text = element_text(size = 16))
```
### Other plot 
```
.libPaths('~/schoenebeck_group/WENGANG/R_lib/')
library(dplyr)
library(DESeq2)
library(tidyr)
test=read.table("count_new_tpm.txt",head=T)
rownames(test)=test[,1]
test=as.matrix(test[,-c(1:6)])
test=as.data.frame(t(test))
test$BCORY=test$BCORY1+test$BCORY2
test$CUL4BY=test$CUL4BY_1+test$CUL4BY_4
test$TSPY=test$TSPY_4+test$TSPY_1+test$TSPY_2+test$TSPY_3+test$TSPY_5+test$TSPY_6
test$UBE1Y=test$UBE1Y_1+test$UBE1Y_2
drop <- c("BCORY1","BCORY2","CUL4BY_1","CUL4BY_4","TSPY_4","TSPY_1","TSPY_2","TSPY_3","TSPY_5","TSPY_6","UBE1Y_1","UBE1Y_2")
test = test[,!(names(test) %in% drop)]
test=t(test)
chrY_coor=c(36610:36632)
vsd_chrY=as.data.frame(t(test[chrY_coor, ]))
vsd_chrY$id=rownames(vsd_chrY)
condition=read.table("sample.condition1",head=F)
colnames(condition)=c("id","tissue","development","system")
condition=unite(condition, tissue, -c("id","system")) 
vsd_chrY =left_join(vsd_chrY,condition,by="id")
Ym=matrix(0,nrow=length(colnames(vsd_chrY)[c(1:23)]),ncol=length(unique(vsd_chrY$tissue)))
for(i in 1:length(colnames(vsd_chrY[c(1:23)]))){
for(j in 1:length(unique(vsd_chrY$tissue))){
di=which(vsd_chrY$tissue==unique(vsd_chrY$tissue)[j])
Ym[i,j]=median(vsd_chrY[di,i])
}
}
colnames(Ym)=unique(vsd_chrY$tissue)
rownames(Ym)=colnames(vsd_chrY[c(1:23)])
Ym=log10(Ym+1)
library(pheatmap)
library(seriation)
library(dendextend)
library("RColorBrewer")
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)
#flip tree order
od=c("TMSB4Y","EIF2S3Y","EIF1AY","DDX3Y","USP9Y","UTY","ZFY","KDM5D","AP1S2Y","WWC3Y","TRAPPC2Y","OFD1","CYorf15","TSPY","RBMYL","PRSSY","TETY2","BCORY","UBE1Y","HSFY","CUL4BY","SRY","AMELY")
phtmap <- pheatmap(Ym,col=colors)
col_dend=phtmap[[1]]
col_dend <- rotate(col_dend, order = od)
pheatmap(Ym,col=colors, cluster_rows=as.hclust(col_dend),display_numbers = TRUE)

#expression and expression specification
library(dplyr)
library(DESeq2)
library(mixOmics)
test=read.table("count_new_tpm.txt",head=T)
rownames(test)=test[,1]
condition=read.table("sample.condition1",head=F)
gene_name=test[,1]
test=as.matrix(test[,-c(1:6)])
samples=data.frame(V1=colnames(test))
yx =left_join(samples,condition,by="V1")
colnames(yx)=c("id","tissue","stage","class")
test=as.data.frame(t(test))
test$BCORY=test$BCORY1+test$BCORY2
test$CUL4BY=test$CUL4BY_1+test$CUL4BY_4
test$TSPY=test$TSPY_4+test$TSPY_1+test$TSPY_2+test$TSPY_3+test$TSPY_5+test$TSPY_6
test$UBE1Y=test$UBE1Y_1+test$UBE1Y_2
drop <- c("BCORY1","BCORY2","CUL4BY_1","CUL4BY_4","TSPY_4","TSPY_1","TSPY_2","TSPY_3","TSPY_5","TSPY_6","UBE1Y_1","UBE1Y_2")
test = test[,!(names(test) %in% drop)]
test=t(test)
chrY_coor=c(36610:36632)
YYY=rownames(test)[chrY_coor]
TAU=c();var=c()
for(i in 1:length(YYY)){
geneexp=data.frame(idd=colnames(test),tpm=log10(test[which(rownames(test)==YYY[i]),]+1))
tau=c()
for (j in 1:20){
yx=yx[sample(1:94),]
tissue=c();idd=c()
for(k in 1:nrow(yx)){
if(length(intersect(as.character(yx[k,2]), tissue))==0){
tissue=c(tissue,as.character(yx[k,2]));idd=c(idd,as.character(yx[k,1]))
}
}
idm=data.frame(idd=idd)
gtm=left_join(idm,geneexp,by="idd")
tau=c(tau,(22-(sum(gtm[,2])/max(gtm[,2])))/21)
}
tau=tau[!is.na(tau)]
TAU=c(TAU,mean(tau))
var=c(var,sd(tau))
}
TAU=(TAU-min(TAU))/(max(TAU)-min(TAU))
df=data.frame(Gene=YYY,tau=TAU,SD=var,test="tau")
average=c();sd=c()
for(i in 1:length(YYY)){
geneexp=data.frame(idd=colnames(test),tpm=log10(test[which(rownames(test)==YYY[i]),]+1))
average=c(average,mean(as.numeric(log(geneexp[,2]+1))))
sd=c(sd,sd(as.numeric(log(geneexp[,2]+1)/4.6)))
}
df2=data.frame(Gene=YYY,tau=average/max(average),SD=sd/max(average),test="Expression")
df=rbind(df,df2)
df$Gene=factor(df$Gene, levels=rev(c("TMSB4Y","EIF2S3Y","EIF1AY","DDX3Y","USP9Y","UTY","ZFY","KDM5D","AP1S2Y","WWC3Y","TRAPPC2Y","OFD1","CYorf15","TSPY","RBMYL","PRSSY","TETY2","BCORY","UBE1Y","HSFY","CUL4BY","SRY","AMELY")))
ggplot(df, aes(x=Gene, y=tau,color=test)) +geom_pointrange(aes(ymin=tau-SD, ymax=tau+SD))+ theme_bw()+
theme(axis.text.y = element_text( color="black", size=10))+ coord_flip()+ scale_color_brewer(palette = "Dark2")

ygene=c("AMELY","AP1S2Y","CYorf15","DDX3Y","EIF1AY","EIF2S3Y","HSFY","KDM5D","OFD1","RBMYL","SRY","TMSB4Y","TRAPPC2Y","USP9Y","UTY","WWC3Y","ZFY","BCORY","UBE1Y","CUL4BY","TSPY")
xgene=c("GeneID:100683823","GeneID:611468","GeneID:480853","GeneID:480886","GeneID:119868387","GeneID:119868410","GeneID:111094819","GeneID:491894","GeneID:480841","GeneID:609457","GeneID:492178","GeneID:100684977","GeneID:119868362","GeneID:480885","GeneID:491845","GeneID:491735","GeneID:119868408","GeneID:480880","GeneID:480896","GeneID:492102","GeneID:491893")
df=data.frame(Y_gene=0,X_gene=0,Gene=0)
correlation=c()
for(i in 1:21){
Y_gene=log10(test[which(rownames(test)==ygene[i]),]+1)
X_gene=log10(test[which(rownames(test)==xgene[i]),]+1)
correlation=c(correlation,cor(Y_gene,X_gene))
df1=data.frame(Y_gene,X_gene,Gene=ygene[i])
df=rbind(df,df1)
}
df=df[-1,]
ggplot(df, aes(x=Y_gene, y=X_gene)) +
  geom_point() + 
  geom_smooth(method=lm, fullrange=TRUE)+ facet_wrap(~Gene, scales = "free",ncol=7)
ggsave("Gene_expression_X_Y_corr.png",dpi = 300, width = 20, height =10,limitsize = FALSE)


library(ggplot2)
condition1=condition[which(condition$V3=="Adult"),]
cc=data.frame(V1=colnames(test))
cc=left_join(cc,condition1,by="V1")
test1=test[,-which(cc$V2=="Testis")]
df=data.frame(Y_gene=0,X_gene=0,Gene=0)
correlation=c()
for(i in 1:21){
Y_gene=log10(test1[which(rownames(test1)==ygene[i]),]+1)
X_gene=log10(test1[which(rownames(test1)==xgene[i]),]+1)
correlation=c(correlation,cor(Y_gene,X_gene))
df1=data.frame(Y_gene,X_gene,Gene=ygene[i])
df=rbind(df,df1)
}
df=df[-1,]
ggplot(df, aes(x=Y_gene, y=X_gene)) +
  geom_point() + 
  geom_smooth(method=lm, fullrange=TRUE)+ facet_wrap(~Gene, scales = "free",ncol=7)
 ggsave("Gene_expression_X_Y_corr_ex_testis.png",dpi = 300, width = 20, height = 10,limitsize = FALSE)
udf=data.frame(ygene,correlation)
add_udf=data.frame(ygene=c("PRSSY","TETY2"),correlation=NA)
udf=rbind(udf,add_udf)
udf$ygene=factor(udf$ygene, levels=rev(c("TMSB4Y","EIF2S3Y","EIF1AY","DDX3Y","USP9Y","UTY","ZFY","KDM5D","AP1S2Y","WWC3Y","TRAPPC2Y","OFD1","CYorf15","TSPY","RBMYL","PRSSY","TETY2","BCORY","UBE1Y","HSFY","CUL4BY","SRY","AMELY")))
ggplot(udf, aes(x=ygene, y=correlation)) +   geom_bar(stat="identity", fill="#FFB266",width=0.6)+  theme_minimal()+ coord_flip()

ddff=data.frame(Gene=0,Tissue=0,Chr=0,Expr=0,Expr_sd=0)
uu=data.frame(V1=colnames(test))
condition1=condition[which(condition$V3=="Adult"),]
uu =left_join(uu,condition1,by="V1")
tis=as.character(unique(uu$V2))
tis=tis[!is.na(tis)]
for(i in 1:21){
Y_gene=test[which(rownames(test)==ygene[i]),]+0.01
X_gene=test[which(rownames(test)==xgene[i]),]+0.01
for(j in 1:length(tis)){
if(mean(Y_gene[which(uu$V2==tis[j])])>=1 |mean((X_gene[which(uu$V2==tis[j])]))>=1){
ddff1=data.frame(Gene=ygene[i],Tissue=tis[j],Chr="Y",Expr=mean(Y_gene[which(uu$V2==tis[j])]/X_gene[which(uu$V2==tis[j])]),Expr_sd=sd(Y_gene[which(uu$V2==tis[j])]/mean(X_gene[which(uu$V2==tis[j])])))
ddff2=data.frame(Gene=ygene[i],Tissue=tis[j],Chr="X",Expr=mean(X_gene[which(uu$V2==tis[j])]/X_gene[which(uu$V2==tis[j])]),Expr_sd=sd(X_gene[which(uu$V2==tis[j])]/mean(X_gene[which(uu$V2==tis[j])])))
ddff=rbind(ddff,ddff1,ddff2)
}else{
ddff1=data.frame(Gene=ygene[i],Tissue=tis[j],Chr="Y",Expr=0,Expr_sd=0)
ddff2=data.frame(Gene=ygene[i],Tissue=tis[j],Chr="X",Expr=0,Expr_sd=0)
ddff=rbind(ddff,ddff1,ddff2)
}
}
}
ddff=ddff[-1,]
ddff$Gene <- factor(ddff$Gene, levels =c("TMSB4Y","EIF2S3Y","EIF1AY","DDX3Y","USP9Y","UTY","ZFY","KDM5D","AP1S2Y","WWC3Y","TRAPPC2Y","OFD1","CYorf15","TSPY","RBMYL","BCORY","UBE1Y","HSFY","CUL4BY","SRY","AMELY"))
ggplot(ddff, aes(x=Tissue, y=Expr, fill=Chr)) +  geom_bar(stat="identity",color = "black", position=position_dodge()) + 
 geom_errorbar(aes(ymin=Expr, ymax=Expr+Expr_sd), width=.2,position=position_dodge(.9)) + facet_wrap(~Gene,  ncol=1,scales = 'free') +  scale_y_sqrt()+ theme_bw()+
theme(axis.text.x = element_text( color="black", size=10,angle=90))
ggsave("Gene_expression_X_Y.png",dpi = 300, width = 10, height = 60,limitsize = FALSE)

library(scales)
ddff1=ddff[which(ddff$Chr=="Y"),]
ddff1=ddff1[-which(ddff1$Expr==0),]
ddff1$Expr[which(ddff1$Expr>32)]=32-runif(length(which(ddff1$Expr>32)),min=0,max=20)
ddff1$Expr[which(ddff1$Expr<0.03125)]=0.03125+runif(length(which(ddff1$Expr<0.03125)),min=0,max=0.001)
ddff1$fills="#404040"
ddff1$fills[which(ddff1$Expr<0.6666)]="#66B2FF"
ddff1$fills[which(ddff1$Expr>1.5)]="#FF9999"
ddff1$shapes=0
ddff1$shapes[which(ddff1$Tissue=="Testis")]=15
add_ddff1=data.frame(Gene=c("PRSSY","TETY2","AMELY"),Tissue=NA,Chr="Y",Expr=NA, Expr_sd=NA,fills="#FF9999", shapes=15)####add
ddff1=rbind(ddff1,add_ddff1)
ddff1$Gene=factor(ddff1$Gene, levels=rev(c("TMSB4Y","EIF2S3Y","EIF1AY","DDX3Y","USP9Y","UTY","ZFY","KDM5D","AP1S2Y","WWC3Y","TRAPPC2Y","OFD1","CYorf15","TSPY","RBMYL","PRSSY","TETY2","BCORY","UBE1Y","HSFY","CUL4BY","SRY","AMELY")))
ggplot(ddff1, aes(x=Gene, y=Expr))+ theme_bw()  + coord_flip()+geom_jitter(size=5,width = 0.1,shape = ddff1$shapes,color = ddff1$fills)+
scale_y_continuous(trans='log2',limits = c(0.03125,32),breaks=c(0.03125,0.0625,0.125,0.25,0.5,1,2,4,8,16,32),labels = percent)+
theme(panel.grid.minor = element_blank(),panel.grid.major.y = element_line(color = "#C0C0C0",size = 1,linetype = "dotted"),panel.grid.major.x=element_blank(),axis.text.y = element_text( color="black", size=12))
ggsave("Gene_expression_X_Y_dot.png",dpi = 300, width = 10, height = 10,limitsize = FALSE)
```
