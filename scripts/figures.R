# ----------------------------------------------------------- #
# Script for visualizing differential gene expression results #
# ----------------------------------------------------------- #

required_libraries <- c("devtools",
                        "BiocManager",
                        "readr", 
                        "tidyverse", 
                        "magrittr", 
                        "ggplot2", 
                        "ggsignif",
                        "gridExtra",
                        "VennDiagram",
                        "ggvenn",
                        "ggVennDiagram",
                        "dplyr", 
                        "patchwork", 
                        "readxl",
                        "ggbeeswarm",
                        "stringr",
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
                        "PoiClaClu") # Poisson distance plot

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



# Define a set of colors and a ggplot theme
my_colors <- c(upregulated = "#D6604D", downregulated = "#4393C3")

my_theme <-   theme(line = element_line(color = "black"),
                    text = element_text(size = 10, color = "black"),
                    panel.background = element_blank(),
                    panel.border = element_rect(color = "black", fill = NA),
                    strip.background = element_blank(),
                    strip.text.x = element_text(face = "bold"),
                    axis.ticks = element_line(color = "black"),
                    axis.text = element_text(color = "black", size = 10),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.position = "right",
                    legend.text = element_text(size = 10),
                    legend.title = element_blank())



# load DGE results
DESeq_results <- read_xlsx("results/DESeq_results.xlsx")

DESeq_results_sig <- DESeq_results %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))


# load TpM table
tpm_long <- read_xlsx("results/tpm.xlsx")

# set the groups as factor and specify the levels
tpm_long$group <- factor(tpm_long$group, levels = c("baseline", "3 weeks MM", "3 weeks RMPI"))

# load the read counts table
counts <- read_xlsx("data/processed/counts_raw.xlsx")

sampleTable <- read_xlsx("data/processed/sampleTable.xlsx")




# -------------- #
# Distance plots #
# -------------- #

# generate distance plots (and also PCA plot) using variance-stabilized read counts

counts %<>% 
  dplyr::select(geneID, matches("SRR")) %>% 
  column_to_rownames("geneID")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleTable %>% column_to_rownames("Sample"),
                              design = ~ group)

# Normalize data
dds <- DESeq(dds)

# filter only for genes with read counts >=10
dds <- dds[rowSums(counts(dds)) >= 10, ]

vsd = varianceStabilizingTransformation(object = dds, 
                                        blind = TRUE,
                                        fitType = "parametric")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL

# Open a PDF device
# pdf("results/figures/sample_distance_heatmap.pdf", width = 4.6, height = 3.3)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(20))

# Close the PDF device
# dev.off()

### based on Poisson Distance using raw counts (not normalized)
poisd <- PoissonDistance((t(counts(dds))))

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- dds$Sample
colnames(samplePoisDistMatrix) <- NULL



# pdf("results/figures/sample_poisson_distance_heatmap.pdf", width = 4.6, height = 3.3)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(20))

# dev.off()





# ------------------------------------------------------------------------------------------ #
# generate boxplots visualizing gene expression levels in TpM of different genes of interest #
# ------------------------------------------------------------------------------------------ #

genes_of_interest <- c("NIBAN1", "PRDM16")
Gene_df <- tpm_long %>% 
  filter(Gene %in% genes_of_interest)

# Generate pairwise comparisons
comparisons <- combn(unique(as.character(Gene_df$group)), 2, simplify = FALSE)

# make a list object to store the boxplots
tpm_plots <- list()

for(gene in Gene_df$Gene){
  
  tmp_df <- Gene_df %>% 
    filter(Gene == gene)
  
  p <- ggplot(tmp_df, aes(group, TpM)) +
    geom_boxplot(aes(fill = group), alpha = .3) +
    geom_beeswarm(shape = 21, size = 4, position = position_dodge(), aes(fill = group), color="black") +
    geom_signif(comparisons = comparisons,
                map_signif_level=TRUE,
                step_increase = 0.1,
                textsize = 4,
                test = "t.test",
                margin_top = 0.1,
                vjust = 0) +
    ggtitle(gene) +
    labs(x="") +
    my_theme +
    theme(legend.position = "None")
  
  tpm_plots[[gene]] <- p
  
  print(p)
}


# Display all plots in a grid
grid.arrange(grobs = tpm_plots)

# save the plots
for(i in seq_along(tpm_plots)){
  plot_name <- names(tpm_plots)[i]
  plot <- tpm_plots[[i]]
  file_name <- paste0("results/figures/",plot_name,".pdf")
  
  ggsave(file_name, plot, height = 3.3, width = 3.1)
}




# ------------ #
# Venn Diagram #
# ------------ #

# make a list containing each DESeq results group as a single dataframe
single_dfs <- list()

unique_groups <- unique(DESeq_results_sig$group)

for (i in unique_groups){
  filtered_data <- DESeq_results_sig %>% 
    filter(group == i)
  
  single_dfs[[i]] <- filtered_data 
}


# regulated_genes <- split(DESeq_results_sig$Gene, DESeq_results_sig$group)


# Dynamically generate the VennDiag list
VennDiag <- setNames(
  lapply(single_dfs, function(df) unique(df$Gene)), 
  unique_groups
)


# Create the eulerr object
eul <- euler(VennDiag, shape = "ellipse", proportional = TRUE)

# Generate the Venn diagram
venn.plot <- plot(eul, #fill = c("grey90", "grey50", "grey70"),
                  quantities = list(type = "counts"),
                  cex = 1)

# Open a PDF device
pdf("results/figures/venn_diagram.pdf", height = 3.1, width = 6)

# Print the plot to the PDF
print(venn.plot)

# Close the PDF device
dev.off()






# ---------- #
# Upset Plot #
# ---------- #

# Pre-allocate the list
upset_list <- list()

# Loop through the unique combinations and extract Ensembl IDs
for (i in seq_along(unique_combinations)) {
  # Extract the current combination name
  comparison_name <- unique_combinations[i]
  
  # Extract the corresponding Ensembl IDs (this assumes the order is correct)
  DESeq_results_sig_padj_filtered <- DESeq_results_sig_padj %>% filter(group == comparison_name)
  gene_ids <- unique(DESeq_results_sig_padj_filtered$geneID)
  
  # Assign to the list
  upset_list[[comparison_name]] <- gene_ids
}

upset_plot <- upset(fromList(upset_list),
                    order.by = "freq",
                    sets.x.label = "number of genes",
                    text.scale = 1.5,
                    point.size = 2)

print(upset_plot)

# Create a named vector for each set to indicate membership
all_genes <- unique(unlist(upset_list))
all_genes <- all_genes[!is.na(all_genes)]
upset_data <- as.data.frame(
  sapply(upset_list, function(x) all_genes %in% x)
)
row.names(upset_data) <- all_genes

# Extract intersections using a custom function
extract_intersections <- function(upset_data) {
  comb_matrix <- as.matrix(upset_data)
  set_names <- colnames(comb_matrix)
  intersect_list <- list()
  
  # Iterate through each combination
  for (i in 1:nrow(comb_matrix)) {
    included_sets <- set_names[comb_matrix[i, ] == TRUE]
    intersection_name <- paste(included_sets, collapse = " & ")
    intersect_list[[intersection_name]] <- c(intersect_list[[intersection_name]], rownames(comb_matrix)[i])
  }
  
  return(intersect_list)
}

# Get intersections
intersect_genes <- extract_intersections(upset_data)


# Create a named vector for mapping geneIDs with gene names
gene_id_to_name <- setNames(gene_ID$Gene, gene_ID$geneID)

# Function to map gene IDs to gene names
map_gene_names <- function(gene_ids) {
  gene_names <- gene_id_to_name[gene_ids]
  return(gene_names)
}

# Apply the function to the intersect_genes list
intersect_genes_with_names <- lapply(intersect_genes, map_gene_names)




# -------- #
# PCA plot #
# -------- #

# perform PCA analysis using variance-stablized counts
pcaData <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

pcaData %<>%
  mutate(group = factor(group, levels = c("baseline", "3 weeks MM", "3 weeks RMPI")))

PCAplot <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size = 4, shape = 21, aes(fill = group)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_label_repel(aes(label = Sample), size = 3, box.padding = 0.5) +
  # stat_ellipse() +
  # coord_fixed(ratio = 1) +
  my_theme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor =element_blank())


# save the PCA plot
# ggsave("results/figures/PCA.pdf", PCAplot, height = 2.8, width = 4.7)




# -------------------------------------- #
# number of up- and down-regulated genes #
# -------------------------------------- #

DGE_numbers <- DESeq_results_sig %>% 
  mutate(direction = case_when(
    log2FoldChange > 0 ~ "upregulated",
    log2FoldChange < 0 ~ "downregulated",
    log2FoldChange == 0 ~ "unchanged")) %>% 
  group_by(group, direction) %>% 
  dplyr::summarize(n = n())

DGE_numbers$direction <- factor(DGE_numbers$direction, level = c("downregulated", "upregulated"))


DGE_numbers_plot <- ggplot(DGE_numbers, aes(x = group, y = n)) +
  geom_bar(stat = "identity", aes(fill = direction), color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#4393C3","#D6604D")) +
  labs(x = "",
       fill = "") +
  ggtitle("Number of regulated genes") +
  my_theme +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom")


# ggsave("results/figures/DGE_numbers.pdf", DGE_numbers_plot, height = 2.2, width = 4.7)









# dumbell plot showing top regulated genes 
# heatmap
# heatmap of shared genes
gene_groups <- DESeq_results_sig %>%
  select(Gene, group) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = group, values_from = present, values_fill = 0)

pheatmap(as.matrix(gene_groups[,-1]), cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = FALSE,
         main = "Regulated Genes Shared Across Groups")
# volcano plot
# MA plot
# new script for enrichment analysis



# -------------------------------------------- #
# R and library versions used in this analysis #
# -------------------------------------------- #
# Get R version
r_version <- R.version.string

# Get package versions
package_versions <- sapply(required_libraries, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    paste(pkg, packageVersion(pkg), sep = ": ")
  } else {
    paste(pkg, "NOT INSTALLED", sep = ": ")
  }
})

# Save version information to a file
writeLines(c(
  paste("R Version:", r_version),
  "Package Versions:",
  paste(package_versions, collapse = "\n")
), "results/package_versions.txt")
