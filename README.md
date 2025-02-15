# mipsa-analysis

This pipeline is used for analyzing immunoprecipitation experiments from [MIPSA](https://www.nature.com/articles/s41551-022-00925-y) screening. 

There are 2 files needed for the analysis from fastq.gz files to fold change calculation. 
1. MIPSA_preprocessing_alignment.sh: Process the raw fastq.gz files from NGS to count matrix for fold change analysis.
2. MIPSA_fc_analysis.Rmd: Read in a count file and study design matrix to perform fold change analysis between samples and negative control (beads only).

Created by Puwanat Sangkhapreecha

Contact: psangkh1@jhu.edu
