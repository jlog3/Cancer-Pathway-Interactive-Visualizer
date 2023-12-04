## Getting Started

To get started with this R/Shiny app, follow these steps:

1. **Create a New Directory**: Open your terminal or file explorer and create a new directory where you want to set up the app. You can do this with the following command:

   ```bash
   mkdir cancer-pathway-interactive-visualizer
   ```
   **Clone Repository**
   ```bash
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
Two CSV files are required for plotting. 

- Experiment/Tumor Data
- Control/Healthy Data

One of these files should look like this:
![Alt Text](Image_URL)

Rows are ensembl gene ids. 
Columns are sample names.

Each file/group must contain at least two samples (columns). 
This allows you to make meaningful comparisons and identify genes that are differentially expressed between the conditions and gain insights into the molecular differences associated with the disease or condition of interest.

#### Data acquisition:
Go to https://portal.gdc.cancer.gov/projects
Under the “Program” dropdown select “TCGA”
Choose your project of interest.  
Follow the R script  downloadTCGA.R
Follow directions for changing the project and sample type as you require. 

#### Data considerations: 
Tissue-Specific Expression:
-For accurate pathway analysis, the data should ideally come from the same tissue or cell type that is under study in the pathway.
-Gene expression can vary significantly between different tissues or cell types. Certain genes may be highly expressed in one tissue but not in others.
Stage of Disease:
-The stage of the disease can also affect gene expression profiles. Users should consider whether their data reflects the disease stage relevant to their research.
Inter-Patient Variability:
-Variability between patients can impact gene expression profiles. If the goal is to draw broader conclusions, it’s beneficial to use data that represents a cross-section of the patient population affected by the specific cancer.
Consistency in Data Sources:
-Consistency in the source of data is important. Mixing data from different tissues, disease stages, or patient profiles without proper controls can lead to misleading conclusions.

You are encouraged to submit data containing the genes in the selected pathway, but the app is flexible and will tell you which genes in the pathway are not found in your data.


## Background

## Future Development
