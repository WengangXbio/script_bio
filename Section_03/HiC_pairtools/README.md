### 1.Parse .sam into pairs
```
pairtools parse -o 190521_E00365_0391.pairs.gz -c ROS_Cfam1.chrom.sizes \
  --drop-sam --drop-seq --output-stats 190521_E00365_0391.stats \
  --assembly roscfam --no-flip \
  --add-columns mapq \
  --walks-policy mask \
190521_E00365_0391.sam
```

### 2.Sorting and deduplication
```
pairtools sort --nproc 5 -o 190521_E00365_0391.sorted.pairs.gz 190521_E00365_0391.pairs.gz --tmpdir=./
pairtools dedup \
    --max-mismatch 3 \
    --mark-dups \
    --output \
        >( pairtools split \
            --output-pairs 190521_E00365_0391.nodups.pairs.gz \
            --output-sam 190521_E00365_0391.nodups.bam \
         ) \
    --output-unmapped \
        >( pairtools split \
            --output-pairs 190521_E00365_0391.unmapped.pairs.gz \
            --output-sam 190521_E00365_0391.unmapped.bam \
         ) \
    --output-dups \
        >( pairtools split \
            --output-pairs 190521_E00365_0391.dups.pairs.gz \
            --output-sam 190521_E00365_0391.dups.bam \
         ) \
    --output-stats 190521_E00365_0391.dedup.stats \
190521_E00365_0391.sorted.pairs.gz
```
### 3.Select qualified pairs
```
gzip -dc 190521_E00365_0391.nodups.pairs.gz | grep -v "#" | head -n 10
gzip -dc 190521_E00365_0391.dups.pairs.gz | grep -v "#" | head -n 10
pairtools select "mapq1>0 and mapq2>0" 190521_E00365_0391.nodups.pairs.gz -o 190521_E00365_0391.nodups.UU.pairs.gz
```
### 4. convert .pairs into .mcool
```
conda activate /home/s1874451/schoenebeck_group/WENGANG/anaconda_env/cooler
cooler cload pairs \
-c1 2 -p1 3 -c2 4 -p2 5 \
--assembly roscfam \
ROS_Cfam1.chrom.sizes:10000 \
190521_E00365_0391.nodups.UU.pairs.gz \
190521_E00365_0391.10000.cool

cooler zoomify \
--nproc 5 \
--out 190521_E00365_0391.1000000.mcool \
--resolutions 1000000,2000000 \
--balance \
190521_E00365_0391.1000000.cool
````
