# Intro

# Protein Folding

## Programs used
ColabFold v1.5.2 (AlphaFold v2.3.1)

FoldSeek (8-ef4e960)

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

Adjust settings like work directory and score strictness in ProteinFolding/AnalyseMTEandpLDDT.py. Then run:

```
python3 ProteinFolding/AnalyseMTEandpLDDT.py
```


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

# 1-Class network creation

## Programmes Used
R 4.2.1

Gephi 0.10.1

## OS
MacOS Ventura

## Workflow
### Network creation
1_layer_make_array.R requires 1 input file, containing the expression data for all timepoints of each gene, and the associated attributes to be visualised on the network.

1_layer_make_array_into_gexf.R accepts the output files of make_array.R and creates a Gephi compatible gexf file for visualising the network.

To create the visualisations seen in the paper, first apply Force Atlas, followed by Fruchterman Reingold (area of 10,000).

# 2-Class network creation and analysis

## Programs used
R 4.2.2

Gephi 0.10

## OS
Windows 11

## Workflow
### Network creation

2_layer_make_array.R requires 2 input files, containing the expression data for all timepoints of each gene, and the associated attributes to be visualised on the network, for genes in two distinct classes (i.e. nematode TFs and effectors).

2_layer_make_array_into_gexf.R accepts the output files of make_array.R and creates a Gephi compatible gexf file for visualising the network.

To create the visualisations seen in the paper, apply Fruchterman Reingold (area of 100,000), followed by Network Splitter 3d (https://gephi.org/plugins/#/plugin/network-splitter-3d). 

For the effector-TF network, [z] was set equal to the degree of the node for TFs, and to 0 for effectors.

The paramaters applied were then:

Z-Maximum Level: equal to the largest [z] in the network

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

For all effector-plant network [z] was set to 80 for all effectors, and to 0 for all plant genes.

The parameters applied were then:

Z-Maximum Level: 80

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

### Adding GO Slim attributes
The list of GO terms and GO slim terms assosciated with plant genes in the network was downloaded from TAIR bulk GO annotation search, functional categorisation and download tool (https://v2.arabidopsis.org/tools/bulk/go/). This was inputted into GO_Slim_attribute_applier_1.R, which outputs a list of genes associated with a particular term. Presence or absence of these terms in the list of genes of interest (i.e. plant genes in the network) was then added to the data as a binary attribute using GO_Slim_attribute_applier_2.R. 

### Analysis of stress or defence genes
Bootstrapping to identify enrichment or depletion of gene sets (i.e. immune/defence-related genes) in each supercluster (as defined by Siddique et al. 2022) was done using Attribute_analysis.R. 
