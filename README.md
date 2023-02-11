# mipsa-analysis

This pipeline contains 3 files needed for the analysis from fastq files to fold change files at both DNA barcodes and their aggregates. 
1. cutadapt.bowtie.sh: Process the raw fastq files from NGS and perform alignment.
2. alignment.R: Read in sam output files and create count files at both DNA barcode and protein level.
3. fc.analysis.R: Read in the count files and study design matrix to perform fold change analysis between samples and negative control (beads only).

Created by Puwanat Sangkhapreecha

Contact: psangkh1@jhu.edu
