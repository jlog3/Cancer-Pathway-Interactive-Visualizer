# ui.R

library(shiny)

# Run Your App:
# use this to run the app -- the other way does not update the cached and will run old versions
# shiny::runApp("/home/john/Desktop/cancer_vis/shiny_app/")

ui <- fluidPage(
  titlePanel("Cancer Pathway Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for selecting a pathway ID
      selectInput("pathwayId", "Select a Pathway ID", choices = c("hsa04010", "hsa04910")),
      
      # File input for gene expression data
      fileInput('file1', 'Choose Tumor Expression Data File',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values,text/plain',
                  '.csv')
      ),
      
      # File input for control group data
      fileInput('file2', 'Choose Control Group Data File',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values,text/plain',
                  '.csv')
      ),
      
      # Button to load example control data
      actionButton("loadExample", "Load Example Control Data")
      
      # Add other input controls here if needed
    ),
    
    mainPanel(
      

      # after Control input
      # Pathway Genes That are not in control data 
      
      # Table for displaying gene names
      tableOutput("geneNames"),
      
      # Disease data  Table for displaying gene names
      tableOutput("controlTable"),
      
      # Control / Example ctl data: Table for displaying mapped IDs
      tableOutput("experimentalTable"),
      
      # combined
      tableOutput("combinedCtlExpTable"),
      
      # Display the pathway diagram
      imageOutput("pathwayOutput"),
      
    )
  )
)



# To use maybe
# pathways <- keggList("pathway", "hsa")
# print(pathways)
# pathway_ids <- names(pathways)
# pathway_names <- pathways


# Old Verions
# ui <- fluidPage(
#   titlePanel("Cancer Pathway Visualization"),
#   sidebarLayout(
#     sidebarPanel(
#       # Add input controls here if needed
#     ),
#     mainPanel(
#       img(src = "hsa04010.pathview.png", height = "600px"),
#       # Add other UI elements here like tables or text
#     )
#     # Input for selecting a pathway ID
#     selectInput("pathwayId", "Select a Pathway ID", choices = c("hsa04010", "hsa04910", ...)),
#
#     # Display results
#     tableOutput("geneNames"),
#     tableOutput("mappedIds")
#   )
# )


# Do:
# Options
# <dropdown> gene data
# <dropdown> 
# after finished selecting option...
# Process gene data and make expression data
# 
# Load example data button: need to load sample expression data 
# 



# Notes
# can have pathway categories  by cancer type to later see/ana connectedness between paths
# 
# No Need for Pre-downloaded Pathways: pathview automatically downloads the required pathway data from KEGG when you specify a pathway ID. Therefore, there's no need to pre-download and store these pathways on your server for production use.
# 
# user might want to gve their own control data but we have some 
# 

# Scraps
# Set the working directory to the parent directory of MyShinyApp in R.
# Use runApp("shiny_app") to start your app.
# out_dir <- "/home/john/Desktop/cancer_vis"
# setwd(out_dir)
# shinyApp(ui = ui, server = server)

