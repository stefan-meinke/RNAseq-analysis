# ------------------------------------------------------------------ #
# This R file contains the functions used in the DGE analysis script #
# ------------------------------------------------------------------ #



########################################
#### Load gene expression (GE) data ####
########################################
load_reads_files <- function(dir, col_names = c("gene", "counts_unstranded", "counts_firstStrand", "counts_secondStrand"), verbose = TRUE){
  
  # Check if the directory exists
  if (!dir.exists(dir)) {
    stop("The directory does not exist: ", dir)
  }
  
  # create the files list
  files_list <- list.files(path = dir, pattern = "*.ReadsPerGene.out.tab", full.names = TRUE, recursive = TRUE)
  
  # load *.tab files containing the read counts
  dfs <- lapply(files_list, function(x) {
    if (verbose) message("Loading file: ", basename(x))
    
    # Error handling for reading files
    tryCatch({
      read_tsv(x, skip = 4, col_names = col_names)
    }, error = function(e) {
      warning("Failed to load file: ", basename(x), " - ", e$message)
      NULL  # Return NULL if file fails to load
    })
  })
  
  # adapt the sample names
  dfs <- setNames(dfs, gsub("\\_ReadsPerGene.out.tab$", "", basename(files_list)))
  
  # export the list object to the global environment
  assign("dfs", dfs, envir = .GlobalEnv)
  
  if (verbose) message("Successfully loaded ", length(dfs), " files.")
}



