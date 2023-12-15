# Intro






# Protein Folding

## Pograms used
ColabFold v1.5.2 (AlphaFold v2.3.1)
FoldSeek (8-ef4e960)
ESMfold v1.0.3
Python 3.11

## OS
CentOS Linux 7

## Prediction pipeline

### ColabFold
1: MMSeqs2 alignment to ColabFold database

```
colabfold_search in.fasata /path/to/ColabFold/DataBase /path/to/out/dir/msas -s 6 --threads 64 --db-load-mode 1
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

```
python ProteinFolding/script.py

```












