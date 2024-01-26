# Intro






# Protein Folding

## Programs used
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

Adjust settings like work directory and score strictness in ProteinFolding/AnalyseMTEandpLDDT.py. Then run:

```
python3 ProteinFolding/AnalyseMTEandpLDDT.py
```


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

# 1-Class network creation

## Programmes Used
R 4.2.1

Gephi 0.10.1

## OS
MacOS Ventura

## Workflow
### Network creation
Requires 1 file, containing the expression data for all time points of the genes to be considered and the associated attributes to be visualised on the network.

Load this file into final_with_negpos.R, then run make_array_into_gexf.R with the output. To improve computation speeds for larger datasets, final_with_negpos_partition.R can be applied instead of final_with_negpos.R, followed by compiler.R, before running make_array_into_gexf.R. 

This will create a gexf file containing a visualisation of the network that can be opened in Gephi. To create the visualisations seen in the paper, the layout applied was first Force Atlas, followed by Fruchterman Reingold with an area of 10,000.

# 2-Class network creation and analysis

## Programs used
R 4.2.2

Gephi 0.10

## OS
Windows 11

## Workflow
### Network creation
Requires 2 files, containing the expression data across all time points of the genes to be considered in both gene classes

Load these files into 2_class_network_creator.R, then run 2_class_gexf_creator on the output

This will create a gexf file containing a visualisation of the network that can be opened in Gephi. To create the visualisations seen in the paper, the layout applied was first Fruchterman Reingold with an area of 100,000, followed by Network Splitter 3d (https://gephi.org/plugins/#/plugin/network-splitter-3d). 

For the effector-TF network, [z] was set equal to the degree of the node for TFs, and to 0 for effectors. 
The paramaters applied were then:

Z-Maximum Level: 50 (equal to the largest [z] in the network)

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

For all effector-plant network visualisations, [z] was set to 80 for all effectors, and to 0 for all plant genes. 
The parameters applied were then:

Z-Maximum Level: 80

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

### Adding GO Slim attributes
The list of GO terms and GO slim terms assosciated with plant genes in the network was downloaded from TAIR, and the subset of genes with assosciated terms of interest identified by GO_Slim_attribute_applier_1.R. Presence or absence in this subset was then applied to the data as a binary attribute by GO_Slim_attribute_applier_2.R. 

### Analysis of stress or defence genes
Bootstrapping to analyse an over or underabundance of edges between effector superclusters and genes involved in immunity or defense was done using Attribute_analysis.R. 




