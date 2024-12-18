# ----------------------------------------------------------------------------------------------------------- #
# Script for loading the counts files und performing differential gene expression (DGE) analysis using DESeq2 #
# ----------------------------------------------------------------------------------------------------------- #


# install and load necessary libraries
required_libraries <- c("devtools",
                        "BiocManager",
                        "readr", 
                        "writexl",
                        "ReportingTools",
                        "tidyverse", 
                        "magrittr", 
                        "ggplot2", 
                        "ggsignif",
                        "gridExtra",
                        "VennDiagram",
                        "ggvenn",
                        "ggVennDiagram",
                        "dplyr", 
                        "ggtranscript", 
                        "patchwork", 
                        "readxl",
                        "ggbeeswarm",
                        "stringr",
                        "knitr",
                        "kableExtra",
                        "scales",
                        "pheatmap",
                        "purrr",
                        "VennDiagram",
                        "eulerr",
                        "UpSetR",
                        "clusterProfiler",
                        "DESeq2",
                        "edgeR",
                        "rtracklayer",
                        "biomaRt",
                        "org.Hs.eg.db",
                        "GenomicFeatures",
                        "openxlsx",
                        "apeglm")

# Check and install/load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


for (package in required_libraries) {
  if (package == "ggtranscript") {
    # Install ggtranscript from GitHub
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("dzhang32/ggtranscript")
  } else {
    # Check if the package is already installed
    if (!requireNamespace(package, quietly = TRUE)) {
      # Install Bioconductor packages
      if (package %in% BiocManager::available()) {
        BiocManager::install(package)
      } else {
        # Install CRAN packages
        install.packages(package, dependencies = TRUE)
      }
    }
  }
  # Load the package
  library(package, character.only = TRUE)
}



# source the helper functions
source("scripts/helper_functions.R")

# Define the root directory of the project
project_dir <- getwd()

# Define the directory for raw read count files
raw_data_dir <- file.path(project_dir, "data", "raw")






# -------------------------------------------------- #
# Load and prepare the counts files for DGE analysis #
# -------------------------------------------------- #

# Load the read count files
load_reads_files(raw_data_dir)


# join all counts dataframes, select only the "gene" and "counts_secondStrand" column
dfs <- lapply(dfs, function(x) {
  first_col_name <- colnames(x)[1]
  last_col_name <- colnames(x)[ncol(x)]
  
  # Select the first and last columns using the column names
  selected_cols <- dplyr::select(x, geneID = !!first_col_name, Sample = !!last_col_name)
})

# save the dataframe names in a variable
df_names <- names(dfs)

# rename all "counts_secondStrand" columns according to the dataframe name
dfs <- map2(dfs, df_names, ~setNames(.x, c("geneID", .y)))

# Perform left joins using reduce()
counts_raw <- purrr::reduce(dfs, left_join, by = "geneID")





# ---------------------------------------------------------------------------------------------- #
# Use biomaRt package to get the gene names and description of genes based on the Ensembl geneID #
# ---------------------------------------------------------------------------------------------- #

# Specify the Ensembl dataset and attributes you want to retrieve
ensembl_dataset <- "hsapiens_gene_ensembl"

# define the ensembl attributes to retrieve
ensembl_attributes <- c("ensembl_gene_id", "external_gene_name", "description") 

# Create a biomart object to access the Ensembl database
ensembl <- useMart("ensembl", dataset = ensembl_dataset)#, optionally use: host = "https://useast.ensembl.org/")

# Define the gene IDs you want to convert
ensembl_gene_ids <- counts_raw$geneID

# Retrieve gene names using the Ensembl gene IDs
gene_ID <- getBM(attributes = ensembl_attributes, filters = "ensembl_gene_id", values = ensembl_gene_ids, mart = ensembl) %>%
  dplyr::rename(geneID = ensembl_gene_id) %>%
  mutate(description = str_extract(description, ".*(?=\\[)")) %>%
  dplyr::rename(Gene = external_gene_name)


# add "Gene" column to the counts_raw dataframe by joining the counts_raw dataframe with the gene_ID df
counts_raw %<>%
  left_join(gene_ID, by = "geneID") %>% 
  relocate(Gene)


# save the raw counts dataframe
# write_xlsx(counts_raw, "data/processed/counts_raw.xlsx")


# ---------------------------------------------------------- #
# Perform differential gene expression analysis using DESeq2 #
# ---------------------------------------------------------- #


# create a read counts matrix
counts_matrix <- counts_raw %>% 
  dplyr::select(geneID,
                matches("SRR")) %>% 
  column_to_rownames("geneID")


# load the meta data (SraRunTable.csv)
SraRunTable <- read_csv("data/SraRunTable.csv")

# create the sampleTable dataframe
sampleTable <- SraRunTable %>% 
  dplyr::select(Sample = Run, group = treatment) %>% 
  column_to_rownames("Sample")

groups <- unique(sampleTable$group)

sampleTable$group <- factor(sampleTable$group, levels=groups)

# ensure that the order of columns in read_counts matches the colData rownames
counts_matrix <- counts_matrix[, rownames(sampleTable)]

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sampleTable,
                              design = ~ group)

# Normalize data
# read for count normalization: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/02_DGE_count_normalization.html
dds <- DESeq(dds)

# filter only for genes with read counts >= 10 
dds <- dds[rowSums(counts(dds)) >= 10, ]

res_names <- resultsNames(dds)

group_contrasts <- res_names[grep("^group_", res_names)]

extract_groups <- function(contrast_string) {
  # Remove the "group_" prefix
  clean_string <- gsub("^group_", "", contrast_string)
  # Split the string by "_vs_" to separate the two groups
  groups <- unlist(strsplit(clean_string, "_vs_"))
  return(groups)
}

group_list <- lapply(group_contrasts, extract_groups)

DESeq_results <- list()

# Loop through the list and build contrasts dynamically
for (i in seq_along(group_list)) {
  
  # Extract the groups for the current contrast
  baseline <- group_list[[i]][2]  # baseline is typically the second element
  treatment <- group_list[[i]][1]  # treatment (e.g., 3.weeks.RMPI or 3.weeks.MM)
  
  # Build the contrast vector
  contrast <- c("group", treatment, baseline)
  
  # Get the results
  res <- results(dds, contrast = contrast)
  
  # Create a coef name based on the res_names entry
  coef_name <- group_contrasts[i]
  
  # Shrink the log fold changes using lfcShrink
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # Convert the result to a data frame and add a column with the contrast name
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("geneID") %>%
    mutate(group = paste(treatment, "vs", baseline))
  
  # Store the results in the list
  DESeq_results[[coef_name]] <- res_df
}


# Combine all DESeq_results dataframes into one dataframe
DESeq_results <- bind_rows(DESeq_results)

# remove NA values, join with gene_ID
DESeq_results %<>% 
  drop_na %>% 
  left_join(gene_ID, by = "geneID") %>%
  relocate(Gene)


# Filter DESeq results for significance
DESeq_results_sig_padj <- DESeq_results %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))


# save the result files
# write_xlsx(DESeq_results, "data/processed/DESeq_results.xlsx")
# write_xlsx(DESeq_results_sig_padj, "data/processed/DESeq_results_sig_padj.xlsx")

