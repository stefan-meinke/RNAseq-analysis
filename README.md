# RNAseq-analysis
Workflow for the analysis of RNA seq data (differential gene expression, alternative splicing)

# General information
This repository represents an example of a data analysis pipeline for differential gene expression (DGE) and alternative splicing (AS) analysis using short-read bulk RNAseq data.
In this example, the SRA study SRP264954 was analyzed. The files were pre-processed using fastp for adapter trimming and quality filtering, and reads were aligned to the reference genome GRCh38.111 (Ensembl) using STAR (default settings). For DGE analysis STAR-derived count files (*ReadsPerGene.out.tab.gz) files were used. 
For DGE the package DESeq2 is used. For further information see: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

