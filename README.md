# Intro






# Protein Folding

## Pograms used
ColabFold v1.5.2 (AlphaFold v2.3.1)
FoldSeek (8-ef4e960)
ESMfold v1.0.3
Python 3.11
SignalP - 4.1

## OS
CentOS Linux 7

## Prediction pipeline
### SignalP
```
signalp -fasta /path/to/genome_or_effector.fasta -mature
```


### ColabFold
1: MMSeqs2 alignment to ColabFold database

```
colabfold_search genome_or_effector.fasata /path/to/ColabFold/DataBase /path/to/out/dir/msas -s 6 --threads 64 --db-load-mode 1
```


2: Protein fold prediction

```
colabfold_batch --pair-mode 'unpaired_paired' \
--num-recycle 3 \
--num-models 3 \
--stop-at-score 100 \
--zip \
path/to/msas/1 \
path/to/out
```

### pLDDT score and pTM analyses

$${\color{!!! add this script later !!!}Red}$$




### ESMFold
Use: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb with default settings

### Foldseek
1: Make database
```
foldseek createdb /path/to/pdb/files /path/to/db/out/folder
```

2: all vs all search 
```
db=/path/to/db/out/folder
cd /path/to/desired/result/folder
foldseek easy-search $db $db aln tmp --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,prob,pident,alntmscore,qtmscore,ttmscore,u
```










