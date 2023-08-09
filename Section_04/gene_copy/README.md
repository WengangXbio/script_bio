## Estimation of MSY gene copy based on cnvnator results
### 1. Convert cnvnator results into .BED format
```
sra=read.table("sra",head=F)
for(j in 1:nrow(sra)){
df=read.table(paste(sra[j,1],".chrY.raw.cnv.depth",sep=""),head=F)
df1=df[df[,1]=="chrY1",]
df2=df[df[,1]=="chrY2",]
df3=df[df[,1]=="chrY3",]
c1s=c(df1[,2][-1],1599364)
c1e=df1[,3]+1
c2s=c(df2[,2][-1],1246049)
c2e=df2[,3]+1
c3s=c(df3[,2][-1],3937623)
c3e=df3[,3]+1
df1a=data.frame(V1=rep("chrY1",length(c1s)),V2=c1e,V3=c1s,V4=1)
df2a=data.frame(V1=rep("chrY2",length(c2s)),V2=c2e,V3=c2s,V4=1)
df3a=data.frame(V1=rep("chrY3",length(c3s)),V2=c3e,V3=c3s,V4=1)
dfm=rbind(df,df1a,df2a,df3a)
dfm=dfm[order(dfm[,1], dfm[,2]), ]
dfm=dfm[-which(dfm[,3]-dfm[,2]<=0),]
write.table(dfm,paste("../gene_copy/",sra[j,1],".depth",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=F)
}
```
### 2. Intersect exons coordinates with depth files
```
for sra in `cat SRA`  ; do
bedtools intersect -a ${sra}.depth -b gene_anno.bed  -wa -wb > ${sra}.exon.depth
done
```
### 3. Estimating gene copy
```
.libPaths('~/schoenebeck_group/WENGANG/R_lib/')
library(ggplot2)
library(ggpubr)
library(tidyr)
sra=read.table("SRA",head=F)
df=read.table(paste(sra[1,1],".exon.depth",sep=""),head=F)
gene=unique(df$V8)
wds=data.frame(gene)
for(j in 1:nrow(sra)){
df=read.table(paste(sra[j,1],".exon.depth",sep=""),head=F)
copy=c()
for(i in 1:length(gene)){
tdf=df[which(df$V8==gene[i]),]
tdf[1,2]=tdf[1,6]
tdf[nrow(tdf),3]=tdf[1,7]
copy=c(copy,sum((tdf$V3-tdf$V2)/sum(tdf$V3-tdf$V2)*tdf$V4))
}
wds=cbind(wds,copy)
}
rownames(wds)=wds$gene
wds=wds[,-1]
wds[14,]=wds[14,]+wds[19,] # BCORY
wds[15,]=wds[15,]+wds[18,] # UBE1Y
wds[23,]=wds[23,]+wds[30,]+wds[31,]+wds[33,] #CUL4BY
wds[24,]=wds[24,]+wds[25,]+wds[26,]+wds[27,]+wds[28,]+wds[32,] #TSPY
wds=wds[-c(18,19,30,25,26,31,33,27,28,32),]
wds=as.data.frame(t(wds))
single_wds=wds[,-c(14,15,21,22,23)]
index=0.5/rowMeans(single_wds)
adj_wds=wds*index
df1=gather(adj_wds, key="Gene", value="Copy", 1:23)
df1$Gene <- factor(df1$Gene, levels=c("TETY2","UTY","DDX3Y","USP9Y","HSFY","EIF1AY","KDM5D","ZFY","EIF2S3Y","WWC3Y","AMELY","AP1S2Y","TMSB4Y","OFD1","TRAPPC2Y","CYorf15","RBMYL","PRSSY","UBE1Y_1","BCORY2","SRY","CUL4BY_1","TSPY_1"))
df1$Copy=log2(df1$Copy*2)+1
df1=df1[-which(df1$Copy<0.5),]
ggboxplot(df1, x = "Gene", y = "Copy", size = 1,add = "jitter")
ggsave("gene_copy_adj.png",dpi = 200, width = 15, height = 5)
write.table(adj_wds,"gene_copy_adj.txt",col.names=T,row.names=T,quote=F,sep="\t")
```
