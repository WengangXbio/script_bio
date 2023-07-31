### Flye assembly
```
flye --pacbio-raw pacbio_subreads.fa --genome-size 2.5g --out-dir ./flye_assembly --keep-haplotypes
```
### Falcon assembly, unzip and polish
```
fc_run fc_run.cfg
fc_unzip.py fc_unzip.cfg
```
### wtdbg2 assembly
```
wtdbg2.pl -t 16 -x rs -g 2500m -o ./wtdbg2_assembly pacbio_subreads.fa
```
