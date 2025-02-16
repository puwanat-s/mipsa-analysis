## MIPSA analysis 

This pipeline is used for analyzing immunoprecipitation experiments from [MIPSA](https://www.nature.com/articles/s41551-022-00925-y) screening. 

There are 2 files needed for the analysis. 
1. `MIPSA_preprocessing_alignment.sh`: Process the raw fastq.gz files from NGS, perform Bowtie alignment, and create a count matrix.
2. `MIPSA_fc_analysis.Rmd`: Perform fold change analysis between samples and negative control (beads only).

## Usage <a name="usage"></a>

### Installation

```
git clone https://github.com/puwanat-s/mipsa-analysis
```

### Prerequisite

For the bash script to run properly,

1. fastq.gz files need to be inside a folder name `data`
2. Name of fastq.gz files need to follow this pattern: `MIPSA22-horf-beads-1_S10_R1_001.fastq.gz`. If not, you can change line 86 in the bash script to print the desired sample name and replicate number
3. Change `line 8` in the bash script to your email address
4. The number of cpus can be adjusted proportionally to the number of fastq.gz files. Change `line 5 and 16` to do so.
5. The preprocessing and alignment run very fast, so the requested time was set for 2 hrs. The time limit can be changed as needed in `line 6`.
6. The current dictionary being used is Ultimate Human ORFeome #3, constructed from Oxford Nanopore Sequencing with Minimap2 alignment. If you need to use a different dictionary, please build Bowtie indices from Bowtie v.1 and create .fai file using samtools. Store them in a folder name `dictionary`

### Quick start

In the directory on Rockfish containing MIPSA_preprocessing_alignment.sh, data, and dictionary, run the following commands:

```
chmod +x MIPSA_preprocessing_alignment.sh
sbatch MIPSA_preprocessing_alignment.sh
```

There should be 2 output folders:
1. `bowtie_output` which contains bam files and `bowtie_results.tsv` which collects the alignment stats from Bowtie
2. `count_matrix` which contains `merged_barcode_counts.tsv` and `study_design.tsv`

Create a new project in Rstudio and copy the Rmd file into the folder. This can be done locally or on Rockfish as suggested in the [user guide](https://www.arch.jhu.edu/guide/). 

The Rmd file contains necessary libraries (e.g., `tidyverse` and `edgeR`) for differential analysis to identify statistically significant enriched targets. It can also generate plots, such as sample-vs-sample scatter plot, MDS plot, and volcano plot. Details can be found in the file. Please feel free to add visualizations for your own investigation! 

You can take a look at `MIPSA_fc_analysis.html` for outputs from the codes. Here are some visualizations:


