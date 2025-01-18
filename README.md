# RNAseq-analysis
Workflow for the analysis of RNA seq data (differential gene expression, alternative splicing)

# General information
This repository provides an example of a data analysis pipeline for **differential gene expression (DGE)** and **alternative splicing (AS)** analysis using short-read bulk RNAseq data.
The dataset used in this example is derived from the SRA study [SRP242571](https://www.ncbi.nlm.nih.gov/sra/?term=SRP242571). 
The R version and package versions used für this analysis are summarized in the [package_versions.txt](https://github.com/stefan-meinke/RNAseq-analysis/blob/main/results/package_versions.txt) file.

## Workflow overview
1. Pre-processing: raw RNA-seq reads were processed using [fastp](https://github.com/OpenGene/fastp) for adapter trimming and quality filtering
2. Alignment: reads were aligned to the reference genome GRCh38.111 (Ensembl) using [STAR](https://github.com/alexdobin/STAR) with default settings.
3. DGE analysis:
   - STAR-derived count files (*ReadsPerGene.out.tab.gz) were used as input.
   - Differential gene expression analysis was performed using the [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) package
4. Alternative splicing analysis:
   - rMATS-turbo ([Wang et al. (2024) Nat Protoc](https://www.nature.com/articles/s41596-023-00944-2#citeas), [rmats-turbo](https://github.com/Xinglab/rmats-turbo)) derived *MATS.JC.txt result files were used as input
   
... To be continued ...
