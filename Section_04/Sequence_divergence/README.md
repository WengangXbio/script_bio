### Sequence divergence using kaks_calculator 2.0
#### 1. Aligning mRNA sequences using MEGA7 and output as fasta format;
#### 2. Convert .fasta into .clustal format using ALTER (https://www.sing-group.org/ALTER/);
#### 3. Convert .clustal into .axt format using kakscalculator2 script:
```
kakscalculator2/bin/AXTConvertor sry_Canidae_cds.fas.alter.aln sry_Canidae_cds.fas.alter.axt
```
#### 4. Calculate KaKs ratio using kaks_calculator (window version)

### Calculate the divergence time of two tips (species) based on a time-scaled phylogeny
```
.libPaths('~/schoenebeck_group/WENGANG/R_lib/')
library(ape)
tree <- read.tree("Canidae.nwk")
pairwise_distances <- cophenetic(tree)
divergence_time=c()
spe1=c()
spe2=c()
for (i in 1:(nrow(pairwise_distances) - 1)) {
  for (j in (i + 1):nrow(pairwise_distances)) {
    divergence_time=c(divergence_time,pairwise_distances[i,j]/2)
	spe1=c(spe1,colnames(pairwise_distances)[i])
	spe2=c(spe2,rownames(pairwise_distances)[j])	
  }
}
data.frame(spe1,spe2,divergence_time)
```
