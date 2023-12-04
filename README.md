## Getting Started

To get started with this R/Shiny app, follow these steps:

1. **Create a New Directory**: Open your terminal or file explorer and create a new directory where you want to set up the app. You can do this with the following command:

   ```bash
   mkdir cancer-pathway-interactive-visualizer
   ```
   **Clone Repository**
   ```bash
   cd cancer-pathway-interactive-visualizer
   git clone https://github.com/jlog3/Cancer-Pathway-Interactive-Visualizer.git cancer-pathway-interactive-visualizer
   ```
2. **Install Required R Packages**:
```R
# Install required R packages
install.packages("shiny", dependencies = TRUE)     # shiny (version >= 1.8.0)
install.packages("pathview", dependencies = TRUE)  # pathview (version >= 1.42.0)
install.packages("KEGGREST", dependencies = TRUE)  # KEGGREST (version >= 1.42.0)
install.packages("TCGAbiolinks", dependencies = TRUE)  # TCGAbiolinks (version >= 2.30.0)
install.packages("Biostrings", dependencies = TRUE)   # Biostrings (version >= 2.70.1)
install.packages("biomaRt", dependencies = TRUE)      # biomaRt (version >= 2.58.0)
install.packages("SummarizedExperiment", dependencies = TRUE)  # SummarizedExperiment (version >= 1.32.0)
install.packages("DESeq2", dependencies = TRUE)       # DESeq2 (version >= 1.42.0)
install.packages("AnnotationDbi", dependencies = TRUE) # AnnotationDbi (version >= 1.64.1)
install.packages("org.Hs.eg.db", dependencies = TRUE)  # org.Hs.eg.db (version >= 3.18.0)
```

3. **Run the Shiny App**:
After installing the packages, you can run the Shiny app using the following command in R:
```R
shiny::runApp("/path/to/shiny_app/")
```
Click "Load Example Data" and wait for plot to update. 


## User Data Upload


#### Considerations: 


## Background

## Future Development
