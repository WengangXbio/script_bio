### 1. Detect SINE sequence on PAR, MSY, and FSX
```
blastn  -evalue 1e-8 -query PAR_MSY_FSX.fa -db SINE.fa -out PAR_MSY_FSX.SINE.blast -outfmt 6 -word_size 11 -reward 1 -penalty -2 -gapopen 5 -gapextend 2
awk '$8-$7>125' PAR_MSY_FSX.SINE.blast |awk 'BEGIN{OFS="\t"} {print $1,$7,$8}'  > PAR_MSY_FSX.SINE.blast.bed
bedtools sort -i PAR_MSY_FSX.SINE.blast.bed > PAR_MSY_FSX.SINE.blast.sorted.bed
bedtools merge -i PAR_MSY_FSX.SINE.blast.sorted.bed > PAR_MSY_FSX.SINE.blast.sorted.merged.bed
```
### 2. Plot the density of SINE around PAB
```
.libPaths('~/schoenebeck_group/WENGANG/R_lib/')     
library(ggplot2)
library(ggpubr)
library(tidyr)
df=read.table("PAR_MSY_FSX.SINE.blast.sorted.merged.bed",head=F)
sine_pab=df[which(df$V1=="chrX"),]
sine_fsx=df[which(df$V1=="FSX_chrX_6590801-8090800"),]
sine_msy=df[which(df$V1=="chrY1"),]
sine_pab$V2=6590800-sine_pab$V2
sine_pab$V3=6590800-sine_pab$V3
sine_fsx_10k=c()
sine_msy_10k=c()
sine_pab_10k=c()
for (i in 1:160){
sta=i*10000-9999;end=i*10000
sine_fsxds=sine_fsx[which(sine_fsx$V2>sta&sine_fsx$V2<end),]
sine_fsx_10k=c(sine_fsx_10k,sum(sine_fsxds$V3 - sine_fsxds$V2))
sine_msyds=sine_msy[which(sine_msy$V2>sta&sine_msy$V2<end),]
sine_msy_10k=c(sine_msy_10k,sum(sine_msyds$V3 - sine_msyds$V2))
sine_pabds=sine_pab[which(sine_pab$V2>sta&sine_pab$V2<end),]
sine_pab_10k=c(sine_pab_10k,sum(sine_pabds$V2 - sine_pabds$V3))
}
group_fsx_msy=c()
for (k in 1:ceiling(length(sine_fsx_10k)/10)){
group_fsx_msy=c(group_fsx_msy,rep(k,10))
}
group_pab=c()
for (k in 1:ceiling(length(sine_pab_10k)/10)){
group_pab=c(group_pab,rep(-k,10))
}
df2=data.frame(FSX=sine_fsx_10k,MSY=sine_msy_10k,coor=group_fsx_msy)
df2_long <- gather(df2, key = "chr", value = "Value", FSX, MSY)
df3_long=data.frame(coor=group_pab,chr="PAB",Value=sine_pab_10k)
dfall=rbind(df2_long,df3_long)
dfall$Value=dfall$Value/10000
ggboxplot(dfall, x = "coor", y = "Value", size = 2,color="chr")+ scale_color_manual(values=c(  "#009999","#D6604D","#0F0F0F"))
ggsave("SINE_all_distribution.png",dpi = 300, width = 15, height = 10)
```
