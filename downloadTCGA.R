library(TCGAbiolinks)


# 1. Run the entire script to declare all functions.
# 2. Step through 'main()' and change the parameters as required.

main <- function() {
  # change to absolute path of where you cloned the repository
  out_dir <- "/home/john/Desktop/cancer-pathway-interactive-visualizer"
  out_dir <- "/home/john/Desktop/cancer_vis"
  
  setwd(out_dir)
  
  # download tumor data
  # change project and sample_types as required
  data_all_primary_tumor <- get_gdc_data("TCGA-CHOL", "Primary Tumor", download = TRUE) # (project, sample_types, download = TRUE, customdir = NULL) 
  # IF YOU ALREADY DOWNLOADED set download = FALSE
  # data_all <-  get_gdc_data("TCGA-CHOL", "Solid Tissue Normal", download = FALSE)
  
  # subset for our personal use...2 samples
  first_five_sample_ids_tumor <- colnames(data_all_primary_tumor)[1:2]
  data_all_primary_tumor_subset <- data_all_primary_tumor[, first_five_sample_ids_tumor]
  
  # process (remove gene version numbers, ) save as CSV
  assay_tumor_data_subset = process_and_write_csv(data_all_primary_tumor_subset, "tumor_rawcounts.csv", TRUE)
  
  
  # do the same for normal tissue
  data_all_primary_normal <- get_gdc_data("TCGA-CHOL", "Solid Tissue Normal")
  
  first_five_sample_ids_normal <- colnames(data_all_primary_normal)[1:2]
  
  # After preparing the data with GDCprepare, you correctly subset the data to include only the first five samples. Ensure that data_all is structured as expected (typically a SummarizedExperiment object):
  data_all_primary_normal_subset <- data_all_primary_normal[, first_five_sample_ids_normal]
  
  # write the raw counts and return matrix
  assay_normal_data_subset = process_and_write_csv(data_all_primary_normal_subset, "normal_rawcounts.csv", TRUE)
 
  
  # for other functions used in the shiny app you can use main.R
}


# Task: get gene tumor OR/AND normal data - if not already downoaded 
get_gdc_data <- function(project, sample_types, download = TRUE, customdir = NULL) {
  # Construct the query
  query <- GDCquery(project = project, 
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "STAR - Counts",
                    sample.type = sample_types)
  
  # Download data if requested
  if (download) {
    GDCdownload(query, method = "api", files.per.chunk = 5)
  }
  
  # Prepare the data
  if (is.null(customdir)) {
    data_all <- GDCprepare(query)
  } else {
    data_all <- GDCprepare(query, directory = customdir)  
  }
  
  return(data_all)
}


# The result is written to a CSV file with the specified filename.
process_and_write_csv <- function(data_subset, filename, writerawcounts=FALSE) {
  # browser()
  # Convert to a regular data frame or matrix
  assay_data_subset <- assay(data_subset)
  #---^ what the user should be giving us -- the rawcounts
  # then we need to do the below preproc and then can get the log folds ourselves
  
  # Write to CSV raw counts
  if (writerawcounts) {
    write.csv(assay_data_subset, file = filename, row.names = TRUE)
  }
  head(assay_data_subset)
  
  # Remove version numbers from Ensembl IDs
  rownames(assay_data_subset) <- gsub("\\..*", "", rownames(assay_data_subset))
  
  # Averaging duplicates: when there are multiple versions of the same gene
  assay_data_subset <- aggregate(assay_data_subset, by = list(rownames(assay_data_subset)), FUN = mean)
  
  # Rename the rownames
  rownames(assay_data_subset) <- assay_data_subset$Group.1
  assay_data_subset <- assay_data_subset[-1]
  
  
  # Assuming assay_data_subset is your data frame or matrix
  
  # Check for duplicate row names
  duplicate_row_names <- rownames(assay_data_subset)[duplicated(rownames(assay_data_subset))]
  
  # See if there are any duplicates
  if (length(duplicate_row_names) > 0) {
    print("There are duplicate row names:")
    print(duplicate_row_names)
  } else {
    print("There are no duplicate row names.")
  }
  
  return(assay_data_subset)
}








