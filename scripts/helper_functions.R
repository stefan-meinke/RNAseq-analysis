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



##################################
#### Load rMATS-turbo results ####
##################################
load_rmats_data <- function(dir, filter_data = TRUE){
  # List subdirectories within the rMATS-turbo directory
  subdirs <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
  
  # create lists to store results
  results_list <- list()        # To store data for each subdirectory
  results_combined_list <- list() # To store combined data for each subdirectory
  
  # Process each subdirectory
  for (subdir in subdirs) {
    # List all files in the subdirectory
    files_list <- list.files(
      path = subdir, 
      pattern = "*.MATS.JC.txt", 
      full.names = TRUE, 
      recursive = TRUE
    )
    
    # Read all files into a list of data frames
    dfs <- lapply(files_list, function(x) read_tsv(x, col_names = TRUE))
    
    # Extract the sample names and assign them as names to the list
    names(dfs) <- gsub("\\.MATS\\.JC\\.txt$", "", basename(files_list))
    
    # Optionally filter the data frames
    if (filter_data) {
      dfs <- lapply(dfs, function(df) {
        df %>%
          filter(
            sapply(IJC_SAMPLE_1, sum_comma_separated) >= 10 | 
              sapply(IJC_SAMPLE_2, sum_comma_separated) >= 10
          ) %>%
          filter(FDR < 0.01 & abs(IncLevelDifference) > 0.15)
      })
    } else {
      dfs <- lapply(dfs, function(df) {
        df %>%
          filter(
            sapply(IJC_SAMPLE_1, sum_comma_separated) >= 10 | 
              sapply(IJC_SAMPLE_2, sum_comma_separated) >= 10
          )
      })
    }
    
    # Add the type of alternative splicing event as a column (e.g., SE, RI)
    dfs <- imap(dfs, ~ mutate(.x, Type = .y))
    
    # combine all dataframes for this subdirectory
    AS_combined <- bind_rows(dfs)
    
    return(AS_combined)
    # # Perform final processing
    # AS_combined <- AS_combined %>%
    #   select(-matches("^ID\\.\\.\\.")) %>% # Remove any columns matching ID... pattern
    #   rename(
    #     ID = ID...1, 
    #     geneID = GeneID, 
    #     Gene = geneSymbol
    #   )
    
    # Save results in the lists
    subdir_name <- basename(subdir)
    results_list[[subdir_name]] <- dfs  # Individual files
    results_combined_list[[subdir_name]] <- AS_combined  # Combined data
  }
}

#############################################################
#### Function to sum rMATS IncLevel values for filtering ####
#############################################################
sum_comma_separated <- function(x) {
  # Split the string by commas
  values <- strsplit(as.character(x), ",")[[1]]
  
  # Convert the resulting character vector to numeric
  numeric_values <- as.numeric(values)
  
  # Sum the numeric values
  sum_value <- sum(numeric_values, na.rm = TRUE)
  
  return(sum_value)
}
