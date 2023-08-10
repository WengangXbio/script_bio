### 1. Call variants on males and females
```
ls *.PAB.vcf.gz >  gvcf.list
gvcfgenotyper \
 -f ROS_Cfam1.fa \
 -l gvcf.list -Ob -o All_PAB.sk2.chrY.bcf
```
### 2. Filter PAB variants
```
 bcftools view  -i  ' COUNT(FMT/FT=="PASS" | FMT/GT=="0/0") >= 716&& COUNT(FMT/FT=="PASS") >1 ' All_PAB.sk2.chrY.bcf |bcftools annotate -x INFO,^FORMAT/GT,FORMAT/FT,FORMAT/DP,FORMAT/GQ > All_PAB.08_CR.1_PASS.vcf; 
  bcftools annotate -x INFO,^FORMAT/GT,FORMAT/FT,FORMAT/DP,FORMAT/GQ All_PAB.08_CR.1_PASS.vcf |grep -v "^#" > All_PAB.08_CR.1_PASS.body
grep "^#" All_PAB.08_CR.1_PASS.vcf > All_PAB.08_CR.1_PASS.header
 bcftools annotate -x INFO,^FORMAT/FT All_PAB.08_CR.1_PASS.vcf |grep -v "^#" \
 |awk '{$1="";$2="";$3="";$4="";$5="";$6="";$7="";$8="";$9="";print $0}' \
 |grep -o -n 'PASS' | cut -d : -f 1 | uniq -c > All_PAB.08_CR.1_PASS.PASS.count
bcftools annotate -x INFO,^FORMAT/FT All_PAB.08_CR.1_PASS.vcf |grep -v "^#" \
 |awk '{$1="";$2="";$3="";$4="";$5="";$6="";$7="";$8="";$9="";print $0}' \
 |grep -o -n '\.' | cut -d : -f 1 | uniq -c > All_PAB.08_CR.1_PASS.dot.count
Rscript  filter_snps.r 895 All_PAB.08_CR.1_PASS.PASS.count All_PAB.08_CR.1_PASS.dot.count
awk 'NR==FNR{l[$0];next;} !(FNR in l)' outlines.txt All_PAB.08_CR.1_PASS.body > All_PAB.08_CR.1_PASS.filter.body
cat All_PAB.08_CR.1_PASS.header All_PAB.08_CR.1_PASS.filter.body > All_PAB.08_CR.1_PASS.filter.vcf
```
### 3. Filter bSNPs and convert vcf to plink
```
bcftools view -m2 -M2 All_PAB.08_CR.1_PASS.filter.vcf -Ov -o All_PAB.08_CR.1_PASS.filter.biallele.vcf
plink --vcf All_PAB.08_CR.1_PASS.filter.biallele.vcf --recode12 --out All_PAB.08_CR.1_PASS.filter.biallele
plink --ped All_PAB.08_CR.1_PASS.filter.biallele.ped --map All_PAB.08_CR.1_PASS.filter.biallele.map --pca --out All_PAB.08_CR.1_PASS.filter.biallele --mind --cow
awk 'BEGIN{OFS="\t"}{print $1,$4,$3,$4}' All_PAB.08_CR.1_PASS.filter.biallele.map > All_PAB.08_CR.1_PASS.filter.biallele.map1
```
### 4. Call male- and female-specific variants
```
awk -F'\t' '$4=="Dog"' All_PAB.ind_inf.txt |grep Male |awk -F'\t' '$5=="keep"' |awk 'BEGIN{OFS="\t"} {print $1,$1}' > Dog_male.txt
awk -F'\t' '$4=="Dog"' All_PAB.ind_inf.txt |grep Female |awk -F'\t' '$5=="keep"' |awk 'BEGIN{OFS="\t"} {print $1,$1}' > Dog_female.txt
plink --ped All_PAB.08_CR.1_PASS.filter.biallele.ped --map All_PAB.08_CR.1_PASS.filter.biallele.map1 --keep Dog_male.txt --mind --out Dog_male --freq
plink --ped All_PAB.08_CR.1_PASS.filter.biallele.ped --map All_PAB.08_CR.1_PASS.filter.biallele.map1 --keep Dog_female.txt --mind --out Dog_female --freq
paste Dog_male.frq Dog_female.frq |awk '{print $2,$5,$11}' |awk '$2>=0.01||$3>=0.01' |awk '$2==0||$3==0' |column -t >  Dog.results
```
### 5. Plot in R
```
library(ggplot2)
library(tidyr)
library(dplyr)
df=read.table("Dog.results",head=F)
df$V3=-(3+log10(df$V3*2))
df$V2=3+log10(df$V2*2)
df=gather(df, key="sex", value="MAF", 2:3)
df=df[-which(df$MAF=="Inf"),]
df=df[-which(df$MAF=="-Inf"),]
ggplot(df, aes(x=V1, y=MAF,color=sex)) +
  geom_segment( aes(x=V1, xend=V1, y=0, yend=MAF),size=1.5) +
  geom_point( size=3.5, shape=21, stroke=2) +xlim(6500000,6600000)+theme_bw()+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +ylim(-3.5,3.5)
ggsave("PAB_0.1Mb.png",dpi = 300, width = 15, height = 10)
```
