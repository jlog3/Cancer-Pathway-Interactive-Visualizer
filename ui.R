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
      # selectInput("pathwayId", "Select a Pathway ID", choices = c("hsa04010", "hsa04910")),
      selectInput("pathwayId", "Select a Pathway ID", choices = NULL),  # popultae from file in server

      
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
      # Use a div to contain the table with absolute positioning
      # div(style = "position: absolute; left: -49%; top: 410px; ",
      #     tableOutput("geneNames")
      # ),
      # Include custom CSS
      tags$head(
        tags$style(HTML("
            .custom-table-style {
                width: 50%; /* Make the table full-width */
                margin-bottom: 20px; /* Same as form.well */
                position: absolute;
                left: -49%; 
                top: 410px;
                /* Add other styling as needed to match form.well */
            }
        "))
      ),
      
      # ... rest of your UI elements ...
      
      # Apply the custom class to your table
      div(class = "custom-table-style", 
          tableOutput("geneNames")
      ),
      
      # Place the uiOutput outside of the div with the absolute positioning
      uiOutput("pathwayUI")
    )
    
  )
)

      #  Can make this all of the gene names in the pathway a drop down with hyperlinks to gene references
      # Consider wrapping tables in their own container\
      # make position absolute so plot gen doesnt mov it 
    
          # tableOutput("controlTable"),
          # tableOutput("experimentalTable"),
          # tableOutput("combinedCtlExpTable")
      
      
      # # Table for displaying gene names
      # tableOutput("geneNames"),
      # 
      # # cotnrl data  Table for displaying gene names
      # tableOutput("controlTable"),
      # 
      # # Control / Example ctl data: Table for displaying mapped IDs
      # tableOutput("experimentalTable"),
      # 
      # # combined
      # tableOutput("combinedCtlExpTable"),
      
      # Display the pathway diagram
      # imageOutput("pathwayOutput"),
      # Dynamic UI for the pathway diagram
      


# Do:
# get selectze to work
# <dropdown> gene data
# <dropdown> 
# after finished selecting option...
# Process gene data and make expression data
# 
# Load example data button: need to load sample expression data 

      
# selectizeInput("pathwayId", "Select a Pathway", choices = NULL, options = list(maxOptions = 10000)),
# selectizeInput("pathwayId", "Select a Pathway", choices = NULL),
# selectizeInput("pathwayId", "Select a Pathway ID", choices = NULL, options = list(create = TRUE)),
# selectizeInput("pathwayId", "Select a Pathway ID", choices = c("Option 1", "Option 2")),
# selectInput("pathwayId", "Select a Pathway ID", choices = NULL, multiple = FALSE, selected=NULL, selectize = TRUE),  # popultae from file in server
# selectizeInput("pathwayId", "Select a Pathway ID", choices = NULL, selected=NULL),




# Notes
# can have pathway categories  by cancer type to later see/ana connectedness between paths
  # here somewhere  https://chat.openai.com/share/2ef33a3b-ba9c-4106-bc69-b1d7eeac3320
  #  But it's not simple to do for now this is fine to  Use all of the genes that have data for them