library(shiny)
library(KEGGREST)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(pathview)
library(TCGAbiolinks)
library(Biostrings) # get seq from ensembl
library(biomaRt)
library(SummarizedExperiment)
library(DESeq2)
# Define server logic required to run the application
# server <- function(input, output) {
#   # The server part is empty as we are just testing the UI dropdown
# }

# To Do 
# allow pathway change
# the input gene fileS should have all possible genes
#   tell user what genes from their file are in the pathway
# make control example data have all genes - so that it can fit with any user and pathway data
# 
# what data can we expect exactly
# 
# Issues: 
# bug: slecting pathway  40, 53 
#Warning: Error in DESeqDataSet: NA values are not allowed in the count matrix

#


server <- function(input, output, session) {
  # browser()
  # Load the CSV file
  pathways_data <- read.csv("data/hsa_pathways_with_genes.csv", stringsAsFactors = FALSE, header = FALSE, col.names = c('pathwayids','pathwaynames','pathwayidsnames'))
  
  # Check the structure of pathways_data
  print(dim(pathways_data))
  print(head(pathways_data))
  
  # Apply the sub function
  # pathway_ids <- sub(":.*", "", pathways_data$Pathway)
  pathway_ids <- pathways_data$pathwayids
  
  
  # Check the result of the sub function
  print(length(pathway_ids))
  print(head(pathway_ids))
  
  # Call this reactive expression in the UI
  # Assuming pathway_ids is correctly extracted
  observe({
    # updateSelectInput(session, "pathwayId", choices = pathways_data$pathwayidsnames)
    updateSelectInput(session, "pathwayId", choices = pathways_data$pathwayidsnames, selected = NULL)
    })
  # then after they select it split at the : to get the pathway id
  
  # Reactive expression to fetch gene names for pathway
  pathway_data <- reactive({
    req(input$pathwayId)  # Ensure that a pathway ID is selected
    message('tying pathway_data')
    
    # will need to get to the left of : later when using full pathway name/title
    # also replace other instances of input$pathwayId   
    # browser()
    
    # while we have the 3 column file form  c('pathwayids','pathwaynames','pathwayidsnames')
    # Retrieve the selected pathway name
    selected_name <- input$pathwayId
    
    # Find the corresponding pathway ID
    selected_id <- pathways_data$pathwayids[which(pathways_data$pathwayidsnames == selected_name)]
    
    # browser()
    gene_names <- get_gene_names_for_pathway(selected_id)  
    # Assuming gene_names is your data frame
    colnames(gene_names)[2] <- selected_name
    # Extracting the second column (now named 'selected_name')
    # selected_name_column <- gene_names$selected_name
    selected_name_column <- gene_names[, 2]
    
    # gene_names <- gene_names[, 0:1]
    # Extract the first column
    first_column <- gene_names[, 1, drop = FALSE]
    
    # Get the row names
    row_names <- rownames(gene_names)
    
    # Combine row names and the first column into a new data frame
    gene_names <- data.frame(RowNames = row_names, FirstColumn = first_column)
    gene_names <- first_column
    
    # Assuming selected_name_column is your character vector
    selected_name_df <- setNames(data.frame(selected_name_column), selected_name)
    
    
    list(gene_names=gene_names, selected_name_column=selected_name_df, selected_id=selected_id)
    # get_gene_names_for_pathway(input$pathwayId)  
    
    # return(selected_id)
  })  
  
  
  # Output for gene names table
  # If pathway_data could be empty or NULL at times (like before a file is uploaded or processed), ensure your code can handle such cases:
  output$geneNames <- renderTable({
    # browser()
    # pathway_data()
    if (is.null(pathway_data()$gene_names) || nrow(pathway_data()$gene_names) == 0) {
      return()  # Return early if data is not available
    }
    # head(pathway_data()$gene_names, n = 7)  # Adjust 'n' as needed
    # head(pathway_data()$selected_name_column, n = 7)  # Adjust 'n' as needed
    pathway_data()$selected_name_column  # Adjust 'n' as needed
    
  })
  
  # Helper Function to get gene names for a pathway
  get_gene_names_for_pathway <- function(pathwayid) {
    # browser()
    message('in getgenenames for pathway?')
    # Attempt to fetch data from KEGG
    result <- tryCatch({
      pathway_info <- keggGet(pathwayid)
      
      # browser()
      
      # Extract gene information
      gene_info <- pathway_info[[1]][["GENE"]]
      
      # Process to extract gene names and descriptions
      genes <- gene_info[seq(2, length(gene_info), by = 2)]
      
      # Extract just the gene names (up to the semicolon)
      gene_names <- sapply(strsplit(genes, ";"), function(x) x[1])
      
      # (might do automatically) find the type of gene name used in the users file so we can translate it to ensemble
      
      # Map gene symbols to Ensembl IDs
      mapped_ids <- mapIds(org.Hs.eg.db, 
                           keys = gene_names, 
                           column = "ENSEMBL", 
                           keytype = "SYMBOL",
                           multiVals = "first")
      # browser()
      # Return the gene names as a data frame for table output
      data.frame(GeneNames = mapped_ids, GeneDetails = genes)
      
    }, error = function(e) {
      # Handle error (e.g., KEGG site down)
      message("Error in fetching data from KEGG: ", e$message)
      
      # Return an empty data frame or a data frame with an error message
      return(data.frame(GeneNames = "Data unavailable due to KEGG server issue"))
    })
    # return(result[["GeneNames"]])
    return(result)
  }

  # done getting pathway details
  
  ### $$$ get missing genes from user/example data
  missingGenesInControl <- reactiveVal()
  missingGenesInExperiment <- reactiveVal()
  
  updateMissingGenes <- function() {
    req(deaResults())
    
    # browser()
    # Get gene names from pathway_data
    pathway_genes <- pathway_data()$gene_names[["GeneNames"]]
    
    # Find missing genes in control data
    missing_control_genes <- pathway_genes[!pathway_genes %in% rownames(geneData())]
    
    # Find missing genes in experimental data
    missing_experiment_genes <- pathway_genes[!pathway_genes %in% rownames(geneData2())]
    
    # Update the reactive values
    missingGenesInControl(missing_control_genes)
    missingGenesInExperiment(missing_experiment_genes)
  }
  
  
  observe({
    # Assuming pathway_data, geneData, and geneData2 are reactive sources
    updateMissingGenes()
  })
  
  # output$missingControlGeneList <- renderText({
  #   paste("Missing Genes in Control:", paste(missingGenesInControl(), collapse = ", "))
  # })
  # Render missing genes in control data
  
  # output$missingControlGeneList <- renderText({
  #   missing_control_genes <- missingGenesInControl()
  #   if (length(missing_control_genes) > 0) {
  #     paste("Missing Genes in Control:", paste(missing_control_genes, collapse = ", "))
  #   } else {
  #     "No missing genes in Control data."
  #   }
  # })
  
  # output$missingExperimentGeneList <- renderText({
  #   paste("Missing Genes in Experiment:", paste(missingGenesInExperiment(), collapse = ", "))
  # })
  # Render missing genes in experimental data
  # output$missingExperimentGeneList <- renderText({
  #   missing_experiment_genes <- missingGenesInExperiment()
  #   if (length(missing_experiment_genes) > 0) {
  #     paste("Missing Genes in Experiment:", paste(missing_experiment_genes, collapse = ", "))
  #   } else {
  #     "No missing genes in Experimental data."
  #   }
  # })
  
  
  
  ### $$$ 
  
  
  # Reactive expression for DEA
  deaResults <- reactive({
    req(geneData(), geneData2())
    print('in dea')
    # browser()
    
    #   ctl         exper
    
    # geneDataX return the ensembleIds and expression 
    
    # fix ensemble -- check if their data has decimal and remove if they do -- ours does not so skip for now
    # rownames(geneData()) <- gsub("\\..*", "", rownames(geneData()))
    # rownames(geneData2()) <- gsub("\\..*", "", rownames(geneData2()))
    
    # rmove the irrelavent genes from the expression set for the user Genedata
    expression_selected_ctl <- geneData()[pathway_data()$gene_names[["GeneNames"]], ]
    expression_selected_exp <- geneData2()[pathway_data()$gene_names[["GeneNames"]], ]
    
    expression_selected_ctl <- na.omit(expression_selected_ctl)
    expression_selected_exp <- na.omit(expression_selected_exp)
    
    
    # having NA in expression_selected_ctl  means that genes from pathway are not in the user/example geneData ?
    
    
    
    
    # gene_expression_sample1 <- expression_selected[, 1:2]  # limit the samples if needed
    
    # if  SummarizedExperiment, ExpressionSet,  then we need to convert using assay   this is fine for now for us sicne we have matrix alreadyy
    # class(expression_selected_exp) # If it's a data.frame, you cannot use assay() directly on it.
    
    # expression_data_matrix_ctl <- assay(expression_selected_ctl)
    # expression_data_matrix_exp <- assay(expression_selected_exp)
    
    # join the matrices and other steps to DEA
    # Get the column names from both matrices
    melanoma_names <- colnames(expression_selected_exp)
    # normal_names <- colnames(normal_expression_matrix)
    kirc_names <- colnames(expression_selected_ctl)
    
    # Combine the column names to create sample_names
    sample_names <- c(melanoma_names, kirc_names)
    
    
    sample_info <- data.frame(
      condition = c(rep("melanoma", 2), rep("normal", 2)),
      row.names = sample_names
    )
    ##$$$$    THIS    THESE MUST MATCH FOR 
    
    combined_expression_data <- cbind(expression_selected_exp, expression_selected_ctl)
    combined_expression_data <- combined_expression_data[, sample_names]   ## Order the columns in the combined_expression_data matrix
    
    dds <- DESeqDataSetFromMatrix(
      countData = combined_expression_data,
      colData = sample_info,
      design = ~ condition
    )
    
    # Run the DEA
    dds <- DESeq(dds)
    
    # Get results
    results <- results(dds)
    
    # 5 Extracting Log Fold Changes:
    # You'll want to extract the log fold changes from the results for visualization in pathview.
    log_fold_changes <- results$log2FoldChange
    names(log_fold_changes) <- rownames(results)
    
    # also return the genes in the pathway that are not in either one of the control or exp data
    
    return(list(log_fold_changes=log_fold_changes, names_log_fold_changes=names(log_fold_changes)))
    
  })
  
  pathway_diagram <- reactive({
    # Ensure file1 is uploaded
    # req(input$file1)
    req(deaResults())
    # browser()
    
    # Define the output directory
    outputDir <- file.path(getwd()) #, "www"
    
    gc()
    # most memory intensive point: 1330 MB 1000 ms
    # Generate pathway diagram
    pathview(gene.data = deaResults()$log_fold_changes,
             # pathway.id = input$pathwayId,
             pathway.id = pathway_data()$selected_id,
             species = "hsa",
             out.dir = outputDir,  # tempdir() in production pros and cons to this
             gene.idtype = "ensembl",
             png = TRUE)
    
    # Construct the dynamic file name
    imageName <- paste0(pathway_data()$selected_id, ".pathview.png")
    
    # Construct the file path for the image
    imagePath <- file.path(outputDir, imageName)
    
    # Return path to the generated image
    list(image = list(src = imagePath))
    
    # Return path to the generated image
    # list(image = list(src = file.path(tempdir(), "pathway_diagram.png")))
  })
  
  # Define a reactive expression for image readiness
  isImageReady <- reactive({
    # Your logic to determine if the image is ready
    # Example: Check if a necessary input is provided
    !is.null(pathway_diagram())
  })
  
  output$pathwayUI <- renderUI({
    if (isImageReady()) {
      imageOutput("pathwayOutput")
    }
  })
  
  # Render the pathway diagram  KEEP THIS
  output$pathwayOutput <- renderImage({
    # browser()
    pathway_diagram()$image  # Assuming pathway_diagram() returns the image information
  }, deleteFile = FALSE)
  
  # output both path diagram and gene list in  the same to keep them block level w each other 
  # Render the pathway diagram
  # output$dynamicContent <- renderUI({
  #   # Assuming pathway_diagram() returns the path to the image
  #   image_path <- pathway_diagram()$image_path
  #   
  #   # Create an image element
  #   image_output <- tags$img(src = image_path, style = "display: block; max-width: 100%; height: auto;")
  #   
  #   # Text outputs for missing genes
  #   missing_control_output <- missingGenesInControl()
  #   missing_experiment_output <- missingGenesInExperiment()
  #   
  #   # Combine image and text outputs in a list
  #   tagList(
  #     image_output,
  #     HTML(paste("Missing Genes in Control:", missing_control_output)),
  #     HTML(paste("Missing Genes in Experiment:", missing_experiment_output))
  #   )
  # }) # , deleteFile = FALSE)
  

    #   # # Reactive value to store example control data status
  exampleControlDataLoaded <- reactiveVal(FALSE)
  # # # Reactive value to store the example data
  # reactiveExampleData <- reactiveVal(NULL)
  
  # Reactive value to store gene data
  geneData <- reactiveVal()
  
  # disease data This reactive expression will process the file when it's uploaded
  # geneData <- reactive({
  observeEvent(input$file1, {
    # browser()
    req(input$file1)
    
    # when we load the example data may want to reverse this to remove user data and load example?
    if (exampleControlDataLoaded()) {
      geneData(NULL)
      geneData2(NULL)
      exampleControlDataLoaded(FALSE)
    }

    # browser()
    # When no file is uploaded, return NULL

    inFile <- input$file1
    data_all <- preprocess_gene_data(inFile)
    # Fallback for an empty or non-existent file
    if (is.null(input$file1)) {
      # return(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
      geneData(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
      
    } else {
      # data_all <- read.csv(input$file1$datapath)
      # colnames(data_all)[1] <- "Gene Names"
      # return(data_all)
      geneData(data_all)
      
    }
  })
  
  geneData2 <- reactiveVal()
  
  # control data  This reactive expression will process the file when it's uploaded
  # geneData2 <- reactive({
  observeEvent(input$file2, {
    # browser()
    req(input$file2)

    if (exampleControlDataLoaded()) {
      geneData(NULL)
      geneData2(NULL)
      exampleControlDataLoaded(FALSE)
    }
    
    inFile <- input$file2
    data_all <- preprocess_gene_data(inFile)
    
    # Fallback for an empty or non-existent file
    if (is.null(input$file2)) {
      # return(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
      geneData2(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
      
    } else {
      # data_all <- read.csv(input$file2$datapath)
      # colnames(data_all)[1] <- "Gene Names"
      # return(data_all)
      geneData2(data_all)
    }
  })
  
  processFile <- function(inFile) {
    # browser()
    tryCatch({
      # If inFile is a file uploaded via Shiny
      if (is.list(inFile)) {
        data <- read.csv(inFile$datapath, row.names = 1)
      } else {
        # If inFile is just a file path
        data <- read.csv(inFile, row.names = 1)
      }
      return(data)
    }, error = function(e) {
      # Error handling
      cat("Error in reading file: ", e$message, "\n")
      return(NULL)
    })
  }
  
  
  preprocess_gene_data <- function(inFile) {
    # browser()
    
    # Read the uploaded file
    # data_all <- read.csv(inFile$datapath, row.names = 1)
    data_all <- processFile(inFile)
    
    
    print(class(data_all))  # "data.frame"
    # data_all <- assay(data_all)  # cant do that
    
    # Create a new column with modified gene IDs (without decimals)
    data_all$GeneID <- gsub("\\..*", "", rownames(data_all))
    
    # Averaging duplicates
    data_all <- aggregate(. ~ GeneID, data = data_all, FUN = mean)
    
    # Set the modified gene IDs as row names
    rownames(data_all) <- data_all$GeneID
    
    # Remove the temporary GeneID column
    data_all$GeneID <- NULL
    
    return(data_all)

  }
  
  observeEvent(input$loadExample, {
    # browser()
    # Load the example data
  
    # just make this load into input$file1 input$file2
      # exampleData <- read.csv("kirc_expression_matrix.csv")
  
    # Paths to the example files
    path_to_file1 <- 'data/primary_tumor_test_rawcounts.csv'
    path_to_file2 <- 'data/primary_normal_test_rawcounts.csv'
    
    # Process the data as if it was uploaded
    file1gd <- processGeneDataFile(path_to_file1) 
    file2gd <- processGeneDataFile(path_to_file2)
    
    # Store the loaded data in a reactive value
    geneData(file1gd)
    geneData2(file2gd)

    # Indicate that the example data is loaded
    exampleControlDataLoaded(TRUE)
  })
  
  processGeneDataFile <- function(inputFile) {
    # browser()
    # Fallback for an empty or non-existent file
    if (is.null(inputFile)) {
      return(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
      
    } else {
      data_all <- preprocess_gene_data(inputFile)
      return(data_all)
    }
  }
  
}
    
    
  
  
    # # Averaging duplicates
    # data_all <- aggregate(data_all, by = list(rownames(data_all)), FUN = mean)
    # 
    # # Rename the rownames
    # rownames(data_all) <- data_all$Group.1
    # data_all <- data_all[-1]
    # 
    # return(data_all)
  
  
  
  
  # experimental data   debug had to fix Warning: Error in <Anonymous>: error in evaluating the argument 'x' in selecting a method for function 'head': 
  # by making geneData return empty dataframe by default 
  # output$controlTable <- renderTable({
  #   # Use browser() to inspect the data
  #   # browser()
  #   data_to_display <- geneData()
  #   if (is.data.frame(data_to_display)) {
  #     head(data_to_display, n = 5)
  #   } else {
  #     return(data.frame(Message = "Data is not available"))
  #   }
  # })
  # # control group  Render the table in the UI
  # output$controlTable <- renderTable({
  #   # geneData()
  #   # browser()
  #   head(geneData(), n = 5)  # Adjust 'n' as needed
  #   
  # })
  
  # # React to the button click
  # load user disease data or example data---v
  # Combined renderTable for both example and uploaded data
  # output$experimentalTable <- renderTable({
  #   if (!is.null(input$file2) && input$file2$size > 0) {
  #     # If a file is uploaded, display user-uploaded data
  #     head(geneData2(), n = 5)  # Adjust 'n' as needed
  #   } else {
  #     # Display example data
  #     head(reactiveExampleData(), n = 5)  # Adjust 'n' as needed
  #   }
  # })
    # } else {
    #   # Optional: Placeholder if no data is loaded
    #   return(data.frame(Message = "No data loaded"))
    # }
  
  # output$combinedCtlExpTable <- renderTable({  
  #   # This will trigger 'deaResults' to run and use its output here
  #   deaResults()
  # })
  
  # dont need to dispay anything after DEA
  # if anything display the pathway genes that are not in one or both of the uploaded data
  # output$combinedCtlExpTable <- renderTable({  
  #   # Check if 'deaResults()' returns a valid dataframe
  #   result <- deaResults()
  #   if (is.data.frame(result)) {
  #     return(result)
  #   } else {
  #     return(data.frame(Message = "No data available"))
  #   }
  # })
  
  

  # 
  # # Output for gene names table
  # output$exampleData <- renderTable({
  #   # Check if example data is loaded
  #   if (exampleControlDataLoaded()) {
  #     # Display head of the example data
  #     head(reactiveExampleData(), n = 5)  # Adjust 'n' as needed
  #   } else {
  #     # Optional: Return a placeholder or message if example data is not loaded
  #     return(data.frame(Message = "Example data not loaded"))
  #   }
  # })
  

