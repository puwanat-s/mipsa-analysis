## mipsa-analysis

This pipeline is used for analyzing immunoprecipitation experiments from [MIPSA](https://www.nature.com/articles/s41551-022-00925-y) screening. 

There are 2 files needed for the analysis. 
1. MIPSA_preprocessing_alignment.sh: Process the raw fastq.gz files from NGS, perform Bowtie alignment, and create a count matrix.
2. MIPSA_fc_analysis.Rmd: Perform fold change analysis between samples and negative control (beads only).

## Usage <a name="usage"></a>

### Installation

```
git clone https://github.com/puwanat-s/mipsa-analysis
```

### Prerequisite

For the Larman lab, you need a Rockfish account to run the bash script. 

For the bash script to run properly,

1. fastq.gz files need to be inside a folder name "data"
2. Name of fastq.gz files need to follow this pattern: MIPSA_run-type_of_library-sample_name-replicate number. For example, MIPSA22-UltimateHumanORF-beads_only-1


Created by Puwanat Sangkhapreecha
