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
# 
# 
# 

server <- function(input, output, session) {
  
  # Reactive expression to fetch gene names for pathway
  pathway_data <- reactive({
    message('tying pathway_data')
    req(input$pathwayId)  # Ensure that a pathway ID is selected
    
    # browser()
    get_gene_names_for_pathway(input$pathwayId)
    
  })  
  # Output for gene names table
  output$geneNames <- renderTable({
    # pathway_data()
    head(pathway_data(), n = 5)  # Adjust 'n' as needed
    
  })
  


    # Reactive expression for DEA
  deaResults <- reactive({
    req(geneData(), geneData2())
    print('in dea')
    
    browser()
    #   ctl         exper
    
    # geneDataX return the ensembleIds and expression 
    
    # fix ensemble -- check if their data has decimal and remove if they do -- ours does not so skip for now
    # rownames(geneData()) <- gsub("\\..*", "", rownames(geneData()))
    # rownames(geneData2()) <- gsub("\\..*", "", rownames(geneData2()))
    
    # rmove the irrelavent genes from the expression set for the user Genedata
    expression_selected_ctl <- geneData()[pathway_data()[["GeneNames"]], ]
    expression_selected_exp <- geneData2()[pathway_data()[["GeneNames"]], ]
    
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
    
    # # 6 Visualization with Pathview:
    # #   Now, use the log_fold_changes for visualization:
    # pathview(gene.data = log_fold_changes, 
    #          pathway.id = "hsa04010", 
    #          species = "hsa", 
    #          gene.idtype = "ensembl",
    #          png = TRUE)
    
  })
  
  # generate pathway diagram
  pathway_diagram <- reactive({
    # Ensure file1 is uploaded
    # req(input$file1)
    req(deaResults())
    # browser()
    
    # Define the output directory
    outputDir <- file.path(getwd()) #, "www"
    
    # Generate pathway diagram
    pathview(gene.data = deaResults()$log_fold_changes,
             pathway.id = input$pathwayId,
             species = "hsa",
             out.dir = outputDir,  # tempdir() in production pros and cons to this
             gene.idtype = "ensembl",
             png = TRUE)
    
    # Construct the dynamic file name
    imageName <- paste0(input$pathwayId, ".pathview.png")
    
    # Construct the file path for the image
    imagePath <- file.path(outputDir, imageName)
    
    # Return path to the generated image
    list(image = list(src = imagePath))
    
    # Return path to the generated image
    # list(image = list(src = file.path(tempdir(), "pathway_diagram.png")))
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
      data.frame(GeneNames = mapped_ids)
      
    }, error = function(e) {
      # Handle error (e.g., KEGG site down)
      message("Error in fetching data from KEGG: ", e$message)
      
      # Return an empty data frame or a data frame with an error message
      return(data.frame(GeneNames = "Data unavailable due to KEGG server issue"))
    })
    
    # return(result[["GeneNames"]])
    return(result)
    
  }
  
  # Render the pathway diagram in UI
  output$pathwayOutput <- renderImage({
    pathway_diagram()$image
  }, deleteFile = FALSE)

  # Output the gene names
  # output$geneNames <- renderTable({
  #   pathway_data()$gene_names
  # })

  # Output the mapped IDs
  # output$mappedIds <- renderTable({
  #   pathway_data()$mapped_ids
  # })
  
  
  
    # # Reactive value to store example control data status
  exampleControlDataLoaded <- reactiveVal(FALSE)
  # # Reactive value to store the example data
  reactiveExampleData <- reactiveVal(NULL)
  

  # Reset example data when a file is uploaded
  observeEvent(input$file2, {
    if (!is.null(input$file2)) {
      # Reset example data and its loaded status
      reactiveExampleData(NULL)
      exampleControlDataLoaded(FALSE)
    }
  })
  
  # disease data This reactive expression will process the file when it's uploaded
  geneData <- reactive({
    req(input$file1)
    # browser()
    # When no file is uploaded, return NULL

    inFile <- input$file1
    data_all <- preprocess_gene_data(inFile)
    # Fallback for an empty or non-existent file
    if (is.null(input$file1)) {
      return(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
    } else {
      # data_all <- read.csv(input$file1$datapath)
      # colnames(data_all)[1] <- "Gene Names"
      return(data_all)
    }
  })
  
  # control data  This reactive expression will process the file when it's uploaded
  geneData2 <- reactive({
    req(input$file2)
    
    inFile <- input$file2
    data_all <- preprocess_gene_data(inFile)
    
    # Fallback for an empty or non-existent file
    if (is.null(input$file2)) {
      return(data.frame(GeneNames = character(0), stringsAsFactors = FALSE))
    } else {
      # data_all <- read.csv(input$file2$datapath)
      # colnames(data_all)[1] <- "Gene Names"
      return(data_all)
    }
  })
  
  preprocess_gene_data <- function(inFile) {
    # browser()
    # Read the uploaded file
    data_all <- read.csv(inFile$datapath, row.names = 1)
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
  output$controlTable <- renderTable({
    # Use browser() to inspect the data
    # browser()
    data_to_display <- geneData()
    if (is.data.frame(data_to_display)) {
      head(data_to_display, n = 5)
    } else {
      return(data.frame(Message = "Data is not available"))
    }
  })
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
  output$experimentalTable <- renderTable({
    if (!is.null(input$file2) && input$file2$size > 0) {
      # If a file is uploaded, display user-uploaded data
      head(geneData2(), n = 5)  # Adjust 'n' as needed
    } else {
      # Display example data
      head(reactiveExampleData(), n = 5)  # Adjust 'n' as needed
    }
    # } else {
    #   # Optional: Placeholder if no data is loaded
    #   return(data.frame(Message = "No data loaded"))
    # }
  })
  
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
  
  
  
  observeEvent(input$loadExample, {
    # browser()
    # Load the example data
    exampleData <- read.csv("kirc_expression_matrix.csv")

    # Store the loaded data in a reactive value
    reactiveExampleData(exampleData)

    # Indicate that the example data is loaded
    exampleControlDataLoaded(TRUE)
  })
}

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
  
  # load user disease data or example data----^
  
 
# # pathway_diagram is a reactive expression that generates a new pathway diagram based on user input.
# # renderImage is used to display the generated pathway diagram in the UI.
# 
# get_gene_names_for_pathway <- function(pathwayid) {
#   # Retrieve pathway information
#   pathway_info <- keggGet(pathwayid)
#   
#   # Extract gene information
#   gene_info <- pathway_info[[1]][["GENE"]]
#   
#   # Process to extract gene names
#   genes <- gene_info[seq(2, length(gene_info), by = 2)]
#   gene_names <- sapply(strsplit(genes, ";"), function(x) x[1])
#   
#   # Map gene symbols to Ensembl IDs
#   mapped_ids <- mapIds(org.Hs.eg.db, 
#                        keys = gene_names, 
#                        column = "ENSEMBL", 
#                        keytype = "SYMBOL",
#                        multiVals = "first")
#   
#   return(list(gene_names = gene_names, mapped_ids = mapped_ids))
# }
# 
# # Example usage
# pathway_genes <- get_gene_names_for_pathway("hsa04010")
# head(pathway_genes$gene_names)
# head(pathway_genes$mapped_ids)

# load example data button: need to load sample expression data 
