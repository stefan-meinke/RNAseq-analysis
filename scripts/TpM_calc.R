# ------------------------------------------------------------------------------ #
# Script for calculating gene expression values in Transcripts per Million (TpM) #
# ------------------------------------------------------------------------------ #

# To get gene expression values in Transcripts per Million (TpM) use edgeR's built-in functions to calculate 
# rpkm, FPKM (Fragments Per Kilobase Million) values and subsequently convert them to TpM.


# install and load all necessary libraries
required_libraries <- c("devtools",
                        "BiocManager",
                        "readxl",
                        "edgeR",
                        "GenomicFeatures")

# Check and install/load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


for (package in required_libraries) {
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
  # Load the package
  library(package, character.only = TRUE)
}



# load the counts_raw.xlsx and make the counts_matrix
counts_raw <- read_xlsx("data/processed/counts_raw.xlsx")

counts_matrix <- counts_raw %>% 
  dplyr::select(geneID,
                matches("SRR")) %>% 
  column_to_rownames("geneID")


# load the gene_ID dataframe
gene_ID <- read_xlsx("data/gene_ID.xlsx")

# load the sampleTable
sampleTable <- read_xlsx("data/processed/sampleTable.xlsx")



# load the gene annotation file (gtf) file 
gtf_url <- "ftp://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"

# create TxDb object
txdb <- makeTxDbFromGFF(gtf_url, format = "gtf")

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum them
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exonic.gene.sizes <- as.data.frame(exonic.gene.sizes)

# organizes the data, renaming columns and creating a data frame for gene lengths.
length_df <- exonic.gene.sizes %>% 
  dplyr::rename(length = exonic.gene.sizes)
length_df <- merge(counts_matrix, length_df, by = 0, all = FALSE)
rownames(length_df) <- length_df$Row.names

#set up the data for differential gene expression analysis using the DGEList function, filtering by expression levels, and calculating normalization factors.
dge <- DGEList(counts=counts_matrix, genes=data.frame(Length=length_df), group = sampleTable$group)
# dge <- dge[filterByExpr(dge), , keep.lib.sizes=FALSE] # not necessarily needed for TpM calculation
dge <- calcNormFactors(dge)

# use edgeR's built-in functions to calculate rpkm, FPKM (Fragments Per Kilobase Million) values and subsequently convert them to TPM.
fpkm <- rpkm(dge, gene.length = dge$genes$Length.length, normalized.lib.sizes = TRUE)
tpm <- fpkm / colSums(fpkm) * 1e6

# if wanted, convert to log2 tpm
log_tpm <- log2(tpm)



# adapt the tpm dataframe 
tpm <- tpm %>% 
  as.data.frame() %>% 
  filter(if_all(where(is.numeric), ~ . >= 1)) %>% 
  rownames_to_column("geneID") %>% 
  left_join(gene_ID, by = "geneID") %>% 
  relocate(Gene)

tpm_long <- tpm %>% 
  pivot_longer(cols=-c(Gene,geneID,description),names_to = "Sample", values_to = "TpM") %>% 
  left_join(sampleTable %>% rownames_to_column("Sample"), by = c("Sample")) %>% 
  dplyr::select(Gene, geneID, Sample, TpM, group)
