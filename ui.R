# install and load the packages

# install.packages("shiny")
# install.packages("shinyFiles")
# install.packages("shinythemes")
# install.packages("shinyWidgets")
# install.packages("ggalluvial")
# install.packages("DT")
# install.packages("scales")
# install.packages("ggrepel")
# install.packages("tidyr")
# 
# install.packages("BiocManager")
# BiocManager::install(version="3.15")
# library(BiocManager)
# BiocManager::install("GeomxTools")
# BiocManager::install("SpatialExperiment")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("clusterProfiler")
# BiocManager::install("SpatialDecon")
# 
# install.packages("devtools")
# library(devtools)
# devtools::install_github("DavisLaboratory/standR")
# 
# library(shiny)
# library(shinyFiles)
# library(shinythemes)
# library(GeomxTools)
# library(shinyWidgets)
# library(ggalluvial)
# library(DT)
# library(scales)
# library(SpatialExperiment)
# library(standR)
# library(edgeR)
# library(limma)
# library(ggrepel)
# library(clusterProfiler)
# library(enrichplot)
# library(SpatialDecon)
# library(tidyr)

# BiocManager::install("org.Hs.eg.db")
# library(org.Hs.eg.db)
# BiocManager::install("org.Mm.eg.db")
# library(org.Mm.eg.db)

library(shiny)
library(shinyFiles)
library(shinythemes)
library(shinyWidgets)
library(GeomxTools)
library(ggalluvial)
library(DT)
library(scales)
library(SpatialExperiment)
library(standR)
library(edgeR)
library(limma)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(SpatialDecon)
library(tidyr)
library(tibble)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(Seurat)
library(patchwork)
library(dplyr)


ui <- navbarPage(
    title = "SpaXSeq",
    theme = shinytheme("flatly"),
    tabPanel(
        title = "Nanostring GeoMx DSP platform",
        navlistPanel(
            widths = c(2, 10),
            "GeoMx-NGS RNA",
            tabPanel(
                title = "Upload Data",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Upload Data"),
                        helpText("Upload data files generated from Nanostring’s GeoMx DSP platform"),
                        strong("Select DCC Folder :"), p(''),
                        shinyDirButton('DCC_dir', 'Browse...', 'Please select a folder containing all DCC files', FALSE, class = "btn-primary"),
                        br(),
                        br(),
                        strong("Select PKC Folder :"), p(''),
                        shinyDirButton('PKC_dir', 'Browse...', 'Please select a folder containing all PKC files', FALSE, class = "btn-primary"),
                        br(),
                        br(),
                        strong("Select Annotation File :"), p(''),
                        shinyFilesButton('anno_file', 'Browse...', 'Please select a file in .xlsx or .xls format', FALSE, class = "btn-primary"),
                        br(),
                        br(),
                        actionButton("upload_button", "Upload"),
                        br(),
                        h3("Research Objective"),
                        helpText("Main experimental questions of interest"),
                        pickerInput(inputId = "factors_of_interest_biology", "Select two key biological metadata variables of concern :", 
                                    choices = "", multiple = TRUE),
                        p(''),
                        actionButton("biology_button", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Study Design"),
                        helpText("Describe the different types of samples"),
                        strong("Select categorical metadata variables displayed on the X-axis :"), p(''),
                        pickerInput(inputId = "factors_of_interest",
                                    choices = "", multiple = TRUE),
                        plotOutput("alluvium_plot"),
                        h3("Study Overview"),
                        helpText("Describe the dimensions of the raw count matrix in probe level"),
                        dataTableOutput("summary_table"),
                    )
                )      
            ),
            "Data Preprocessing",
            tabPanel(
                title = "Segment QC",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Segment QC"),
                        helpText("Assess sequencing quality and adequate tissue sampling for every ROI/AOI segment"),
                        strong("Select parameter cutoffs :"), p(''),
                        numericInput("Reads", "Minimum number of raw reads",
                                    min = 500, max = 10000, value = 1000, step = 100),
                        numericInput("Trimmed", "Minimum percentage of reads trimmed (%)",
                                    min = 30, max = 100, value = 80, step = 5),
                        numericInput("Stitched", "Minimum percentage of reads stitched (%)",
                                    min = 30, max = 100, value = 80, step = 5),
                        numericInput("Aligned", "Minimum percentage of reads aligned (%)",
                                    min = 30, max = 100, value = 80, step = 5),
                        numericInput("Saturation", "Minimum sequencing saturation (%)",
                                    min = 0, max = 100, value = 50, step = 5),
                        numericInput("NegativeCount", "Minimum negative probe count geomean",
                                    min = 0, max = 50, value = 10, step = 1),
                        numericInput("NTCCount", "Maximum No Template Control (NTC) count",
                                    min = 0, max = 9000, value = 60, step = 10),
                        numericInput("Nuclei", "Minimum nuclei count",
                                    min = 0, max = 500, value = 200, step = 10),
                        numericInput("Area", "Minimum surface area (μm²)",
                                    min = 500, max = 20000, value = 16000, step = 100),
                        br(),
                        actionButton("segmentQC_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Segment QC Statistic"),
                        helpText("Describe the distributions of the segments for the different QC parameters"),
                        strong("Select a categorical metadata variable displayed on the X-axis :"), p(''),
                        pickerInput(inputId = "factor_of_interest_segqc", choices = ""),
                        tabsetPanel(type = "tabs",
                                    tabPanel("Raw reads", 
                                             plotOutput("raw_reads"),
                                             sliderInput("Reads_2", "Change parameter cutoffs :",
                                                         min = 500, max = 10000, value = 1000, step = 100, width = "800px")),
                                    tabPanel("% Trimmed", 
                                             plotOutput("reads_trimmed"),
                                             sliderInput("Trimmed_2", "Change parameter cutoffs :",
                                                         min = 30, max = 100, post = "%", value = 80, step = 5, width = "800px")),
                                    tabPanel("% Stitched", 
                                             plotOutput("reads_stitched"),
                                             sliderInput("Stitched_2", "Change parameter cutoffs :",
                                                         min = 30, max = 100, post = "%", value = 80, step = 5, width = "800px")),
                                    tabPanel("% Aligned", 
                                             plotOutput("reads_aligned"),
                                             sliderInput("Aligned_2", "Change parameter cutoffs :",
                                                         min = 30, max = 100, post = "%", value = 80, step = 5, width = "800px")),
                                    tabPanel("% Saturation", 
                                             plotOutput("sequencing_saturation"),
                                             sliderInput("Saturation_2", "Change parameter cutoffs :",
                                                         min = 0, max = 100, post = "%", value = 50, step = 5, width = "800px")),
                                    tabPanel("Negative count", 
                                             dropdownButton(
                                                 selectInput("panel", "Switch panel used :",
                                                             choices = ""),
                                                 icon = icon("gear"),
                                                 width = "300px",
                                                 tooltip = tooltipOptions(title = "Click to switch panel used")
                                             ),
                                             plotOutput("negative_control_counts"),
                                             sliderInput("NegativeCount_2", "Change parameter cutoffs :",
                                                         min = 0, max = 50, value = 10, step = 1, width = "800px")),
                                    tabPanel("NTC", 
                                             dataTableOutput("NTC_count")),
                                    tabPanel("Nuclei", 
                                             plotOutput("nuclei"),
                                             sliderInput("Nuclei_2", "Change parameter cutoffs :",
                                                         min = 0, max = 500, value = 200, step = 10, width = "800px")),
                                    tabPanel("Area", 
                                             plotOutput("area"),
                                             sliderInput("Area_2", "Change parameter cutoffs :",
                                                         min = 500, max = 20000, post = "μm²", value = 16000, step = 100, width = "800px"))),
                        h3("Segment QC Result"),
                        helpText("Describe the number of the segments after performing segment QC"),
                        dataTableOutput("summary_seg")
                    )
                )
            ),
            tabPanel(
                title = "Probe QC",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Probe QC"),
                        helpText("An outlier probe removal process"),
                        strong("Remove global outliers :"), p(''),
                        numericInput("ProbeRatio", "(geomean probe in all segments / geomean probes within given gene) less than",
                                     min = 0, max = 1, value = 0.1, step = 0.1),
                        numericInput("FailGrubbs", "Fail Grubb’s test in at least % of the segments",
                                     min = 0, max = 100, value = 20, step = 5),
                        strong("Remove local outliers :"), p(''),
                        switchInput(inputId = "LocalOutliers", value = TRUE, onLabel = "Yes", offLabel = "No"),
                        br(),
                        actionButton("probeQC_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Probe QC Result"),
                        helpText("Describe the number of the probes after performing Probe QC"),
                        dataTableOutput("summary_probe"),
                        h3("Study Overview"),
                        helpText("Describe the dimensions of the collapsed count matrix in gene level"),
                        dataTableOutput("summary_table_gene")
                    )
                )
            ),
            tabPanel(
                title = "Filtering",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Filtering"),
                        helpText("Filter out either segments and/or genes with abnormally low signal based on LOQ"),
                        strong("Select segment filtering cutoffs :"), p(''),
                        numericInput("seg_thr", "Gene detection threshold in each segment (%)",
                                    min = 0, max = 100, value = 10, step = 1),
                        strong("Select gene filtering cutoffs :"), p(''),
                        numericInput("gene_thr", "Global gene detection rate threshold (%)",
                                    min = 0, max = 100, value = 10, step = 1),
                        br(),
                        actionButton("filter_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Filtering Statistic"),
                        tabsetPanel(type = "tabs",
                                    tabPanel("Segment filtering", 
                                             h3("Segment Gene Detection"),
                                             helpText("Describe the distribution of segments with respect to their percentage of genes detected"),
                                             dropdownButton(
                                                 selectInput(inputId = "factor_of_interest_filter", 
                                                             "Select a grouping metadata variable :", 
                                                             choices = ""),
                                                 icon = icon("gear"),
                                                 width = "300px",
                                                 tooltip = tooltipOptions(title = "Click to change grouping metadata variable :")
                                             ),
                                             plotOutput("seg_filter_plot")),
                                    tabPanel("Gene filtering", 
                                             h3("Gene Detection Rate"),
                                             helpText("Describe the total number of genes detected in different percentages of segments"),
                                             plotOutput("gene_filter"))),
                        h3("Filtering Result"),
                        helpText("Describe the dimensions of the collapsed count matrix in gene level after performing filtering"),
                        dataTableOutput("summary_table_filter")
                    )
                )
            ),
            tabPanel(
                title = "Normalization",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Normalization"),
                        helpText("Rectify technical variations"),
                        pickerInput(inputId = "nor_method", "Select normalization method :",
                                    choices = c("TMM", "RPKM", "TPM", "CPM", "upperquartile", "sizefactor"),
                                    selected = "TMM"),
                        helpText("Note: RPKM and TPM require gene length information"),
                        h3("Batch correction"),
                        helpText("Remove batch effect introduced by the different slides using RUV4 method"),
                        numericInput("k_best", "Select the number k of factors to remove :",
                                    min = 0, max = 10, value = 2, step = 1),
                        em("Hint: The optimal k is the smallest value at which 
                            the observed technical variation is no longer present"),
                        p(""),
                        actionButton("nor_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Comparison of Results"),
                        helpText("Describe the effects of normalization"),
                        tabsetPanel(type = "tabs",
                                    tabPanel("RLE",
                                             h3("Raw Data"),
                                             plotOutput("nor_before"),
                                             h3("Normalized Data"),
                                             plotOutput("nor_after"),
                                             h3("Batch Corrected Data"),
                                             plotOutput("nor_batch")
                                    ),
                                    tabPanel("PCA",
                                             h3("Raw Data"),
                                             plotOutput("pca_before"),
                                             h3("Normalized Data"),
                                             plotOutput("pca_after"),
                                             h3("Batch Corrected Data"),
                                             plotOutput("pca_batch")
                                    )
                        )
                    )
                )
            ),
            "Statistical Analysis",
            tabPanel(
                title = "Dimension Reduction",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Dimension Reduction"),
                        helpText("Identify broad patterns in the high-plex expression data"),
                        strong("Select plot parameters :"), p(''),
                        pickerInput("factor_of_interest_color", "Color by",
                                    choices = ""),
                        pickerInput("factor_of_interest_shape", "Shape by",
                                    choices = ""),
                        actionButton("dim_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Dimension Reduction Results"),
                        helpText("Describe the structure of data"),
                        tabsetPanel(type = "tabs",
                                    tabPanel("PCA", plotOutput("pca_plot", height = "800px")),
                                    tabPanel("MDS",plotOutput("mds", height = "800px")),
                                    tabPanel("UMAP", plotOutput("umap", height = "800px")),
                                    tabPanel("t-SNE",plotOutput("tsne", height = "800px")))
                     )
                )
            ),
            tabPanel(
                title = "Differential Expression",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Differential Expression"),
                        helpText("Identify genes that exhibit significant expression differences between distinct sample groups"),
                        strong("Select pairwise comparison :"), p(''),
                        pickerInput("contrast_1", "Group 1",
                                    choices = ""),
                        pickerInput("contrast_2", "Group 2",
                                    choices = ""),
                        strong("Select significance parameters :"), p(''),
                        numericInput("p_value", "Adjusted P-value threshold",
                                     min = 0, max = 1, value = 0.05, step = 0.001),
                        numericInput("logfc", "Log2 fold change threshold",
                                     min = 0, max = 10, value = 1, step = 0.1),
                        actionButton("expre_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Differential Expression Result"),
                        tabsetPanel(type = "tabs",
                                    tabPanel("Summary",
                                             h3("Difference Overview"),
                                             helpText("Describe the number of differential genes"),
                                             dataTableOutput("summary_de_table"),
                                             h3("Differential Expression Genes"),
                                             helpText("Describe statistical significance and magnitude for each differential gene"),
                                             dataTableOutput("de_table")
                                             ),
                                    tabPanel("Visualization",
                                             h3("Visualization Plot"),
                                             helpText("Visualize differential gene expression patterns"),
                                             tabsetPanel(type = "tabs",
                                                         tabPanel("Volcano", plotOutput("volcano", height = "800px")),
                                                         tabPanel("MA", plotOutput("ma", height = "800px")),
                                                         tabPanel("Heatmap",plotOutput("heatmap", height = "800px")))
                                             ),
                                    tabPanel("Genes of interest",
                                             h3(),
                                             helpText(""),
                                             textInput(inputId = "gene_name",
                                                       label = "Enter gene symbol(s) separated with /, i.e. NLGN1/A2M",
                                                       width = "800px"),
                                             plotOutput("violin_gene_of_interest")))
                    )
                )
            ),
            "Downstream Analysis",
            tabPanel(
                title = "Gene Set Enrichment Analysis",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("GO GSEA"),
                        helpText("Interpret DGE data by identifying predefined functional GO gene sets"),
                        strong("Select species :"), p(''),
                        pickerInput("species",
                                    choices = c("Human","Mouse"),
                                    selected = "Human"),
                        strong("Select significance parameter :"), p(''),
                        numericInput("p_value_gsea", "Adjusted P-value threshold",
                                     min = 0, max = 1, value = 0.05, step = 0.001),
                        actionButton("gsea_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("GO GSEA Result"),
                        tabsetPanel(type = "tabs",
                                    tabPanel("Summary",
                                             h3("Enriched GO terms"),
                                             helpText("Describe significantly enriched GO terms"),
                                             dataTableOutput("gsea_table")),
                                    tabPanel("Visualization",
                                             h3("Visualization Plot"),
                                             helpText("Visualize and interpret enrichment results"),
                                             tabsetPanel(type = "tabs",
                                                         tabPanel("Bar Plot", plotOutput("bar_plot"))))
                                                         # tabPanel("MA", plotOutput("test")),
                                                         # tabPanel("Heatmap",plotOutput("test_2"))))
                        )
                    )
                )
            ),
            tabPanel(
                title = "Cellular Deconvolution",
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        h3("Cellular Deconvolution"),
                        helpText("Estimate the proportions of different cell types in samples"),
                        strong("Select a pre-specified cell profile matrix :"), p(''),
                        pickerInput("species_2", "Species",
                                    choices = c("Human","Mouse"),
                                    selected = "Human"),
                        pickerInput("age_group", "Age group",
                                    choices = c("Adult", "Fetal", "Fetal/E14.5", "Fetal/E9.5-13.5", "Neonatal", "COVID-Infected")),
                        textInput("matrix_name", "Name of profile matrix (e.g. Gut_HCA"),
                        em("Details: complete list of matrices can be found on the"),
                        a(href = "https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/NewProfileMatrices",
                          em("CellProfileLibrary GitHub Page")),
                        br(),
                        actionButton("dec_submit", "Submit")
                    ),
                    mainPanel(
                        width = 9,
                        h3("Cellular Deconvolution Result"),
                        helpText(""),
                        plotOutput("dec_heatmap"),
                        plotOutput("dec_plot")
                    )
                )
            )
        )
    ),
    tabPanel(
        title = "10X Genomics Visium platform")
)












