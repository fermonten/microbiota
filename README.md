# Pipelines for microbiota analysis

This repository contains code for the analysis of 16S rRNA sequences as part of a Master's thesis. Three scripts were used in the final analysis: one script for QIIME2; one for PICRUSt2 and another one for all the downstream analysis in R.

## QIIME2
The script qiime.sh was used for the processing of raw sequences, merging, quality control, taxonomic classification and phylogeny building. The outputs produced are the ASV table, taxonomy table, phylogenetic tree. 

## PICRUSt2
The script picrust.sh was used to generate enzyme and pathway data from the ASV table and the representative sequences. The outputs produced are EC, KEGG and MetaCyc pathway tables.

## R downstream analysis
The microeco.R file conatins all the downstream analysis performed in R, with the microeco package and others.
