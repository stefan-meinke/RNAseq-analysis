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




