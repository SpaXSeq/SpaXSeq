
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
            strong("Select DCC folder :"), p(''),
            shinyDirButton('DCC_dir', 'Browse...', 'Please select a folder containing all DCC files', FALSE, class = "btn-primary"),
            br(),
            br(),
            strong("Select PKC folder :"), p(''),
            shinyDirButton('PKC_dir', 'Browse...', 'Please select a folder containing all PKC files', FALSE, class = "btn-primary"),
            br(),
            br(),
            strong("Select annotation file :"), p(''),
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
            DT::dataTableOutput("summary_table")
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
                                 DT::dataTableOutput("NTC_count")),
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
            DT::dataTableOutput("summary_seg")
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
            DT::dataTableOutput("summary_probe"),
            h3("Study Overview"),
            helpText("Describe the dimensions of the collapsed count matrix in gene level"),
            DT::dataTableOutput("summary_table_gene")
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
                                   tooltip = tooltipOptions(title = "Click to change grouping metadata variable")
                                 ),
                                 plotOutput("seg_filter_plot")),
                        tabPanel("Gene filtering",
                                 h3("Gene Detection Rate"),
                                 helpText("Describe the total number of genes detected in different percentages of segments"),
                                 plotOutput("gene_filter_plot"))),
            h3("Filtering Result"),
            helpText("Describe the dimensions of the collapsed count matrix in gene level after performing filtering"),
            DT::dataTableOutput("summary_table_filter")
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
            br(),
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
                                 plotOutput("pca_batch")))
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
            br(),
            actionButton("dim_submit", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("Visualization"),
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
            strong("Select significance thresholds :"), p(''),
            numericInput("p_value", "Adjusted P-value",
                         min = 0, max = 1, value = 0.05, step = 0.001),
            numericInput("logfc", "Log2FC of Group1/Group2",
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
                                 DT::dataTableOutput("summary_de_table"),
                                 h3("Differential Expression Genes"),
                                 helpText("Describe statistical significance and magnitude for each differential gene"),
                                 DT::dataTableOutput("de_table")),
                        tabPanel("Visualization",
                                 h3("Visualization Plot"),
                                 helpText("Visualize differential gene expression patterns"),
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Volcano", plotOutput("volcano", height = "800px")),
                                             tabPanel("MA", plotOutput("ma", height = "800px")),
                                             tabPanel("Heatmap", plotOutput("heatmap", height = "800px")),
                                             tabPanel("Gene of interest",
                                                      p(""),
                                                      strong("Enter gene symbol :"), p(''),
                                                      textInput(inputId = "gene_name",
                                                                label = NULL,
                                                                placeholder = "e.g. NLGN1"),
                                                      plotOutput("violin_gene_of_interest", width = "400px")))))
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
            helpText("Interpret DEG data by identifying predefined functional GO gene sets"),
            pickerInput("species", "Select species :",
                        choices = c("Human","Mouse"),
                        selected = "Human"),
            pickerInput("subontology", "Select subontology :",
                        choices = c("MF","BP","CC","ALL"),
                        selected = "BP"),
            numericInput("p_value_gsea", "Select P-value threshold :",
                         min = 0, max = 1, value = 0.05, step = 0.001),
            br(),
            actionButton("gsea_submit", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("GO GSEA Result"),
            tabsetPanel(type = "tabs",
                        tabPanel("Summary",
                                 h3("Enriched GO terms"),
                                 helpText("Describe significantly enriched GO terms"),
                                 DT::dataTableOutput("gsea_table")),
                        tabPanel("Visualization",
                                 h3("Visualization Plot"),
                                 helpText("Visualize functional enrichment results"),
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Dot Plot",
                                                      dropdownButton(
                                                        numericInput("category_dot", "Select the number of terms displayed :",
                                                                     min = 1, max = 100, value = 20, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("dot_plot", height = "800px")),
                                             tabPanel("Network Plot",
                                                      dropdownButton(
                                                        numericInput("category_network", "Select the number of terms displayed :",
                                                                     min = 1, max = 50, value = 5, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("network", height = "700px")),
                                             tabPanel("Ridgeline plot",
                                                      dropdownButton(
                                                        numericInput("category_ridgeline", "Select the number of terms displayed :",
                                                                     min = 1, max = 100, value = 20, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("ridgeline_plot", height = "800px")),
                                             tabPanel("Enrichment score",
                                                      dropdownButton(
                                                        numericInput("category_score", "Select the number of terms displayed :",
                                                                     min = 1, max = 20, value = 3, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("score_plot")))))
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
            strong("Select a pre-defined cell profile matrix :"), p(''),
            pickerInput("species_2", "Species",
                        choices = c("Human","Mouse"),
                        selected = "Human"),
            pickerInput("age_group", "Age group",
                        choices = c("Adult", "Fetal", "Fetal/E14.5", "Fetal/E9.5-13.5", "Neonatal", "COVID-Infected")),
            textInput("matrix_name",
                      "Name of profile matrix",
                      placeholder = "e.g. Gut_HCA"),
            em("Details: a complete list of matrices can be found on the"),
            a(href = "https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/NewProfileMatrices",
              em("CellProfileLibrary GitHub Page")),
            br(),
            br(),
            actionButton("dec_submit", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("Cellular Deconvolution Result"),
            tabsetPanel(type = "tabs",
                        tabPanel("Overview",
                                 h3("Heatmap of Cell Abundance"),
                                 helpText("Describe the relative abundance of cell types for each segment"),
                                 plotOutput("dec_heatmap", height = "800px")),
                        tabPanel("Contrast",
                                 h3("Differential Abundance of Cell Types"),
                                 helpText("Describe the proportion distribution of cell types in each segment of pairwise comparision"),
                                 plotOutput("dec_barplot", height = "800px")))
          )
        )
      )
    )
  ),
  tabPanel(
    title = "10X Genomics Visium platform",
    navlistPanel(
      widths = c(2, 10),
      "10X Visium RNA",
      tabPanel(
        title = "Upload Data",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            h3("Upload Data"),
            helpText("Upload data files generated from 10X Genomics Visium platform"),
            strong("Select folder :"), p(''),
            shinyDirButton('folder', 'Browse...', 'Please select a folder containing the H5 file and the image data', FALSE, class = "btn-primary"),
            br(),
            br(),
            actionButton("upload_button_seurat", "Upload")
          ),
          mainPanel(
            width = 9,
            h3("Study Overview"),
            helpText("Describe the dimensions of the raw count matrix in spot level"),
            DT::dataTableOutput("summary_table_seurat")
          )
        )
      ),
      "Data Preprocessing",
      tabPanel(
        title = "Filtering & Normalization",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            h3("Filtering & Normalization"),
            helpText("Filter out low-quality spots and normalization using sctransform method"),
            strong("Select filtering cutoffs :"), p(''),
            numericInput("nCount_min","minimum UMIs per spot",
                         min = 0, max = 100000, value = 0, step = 500),
            numericInput("nCount_max","maximum UMIs per spot",
                         min = 0, max = 100000, value = 100000, step = 500),
            numericInput("nFeature_min", "minimum genes per spot",
                         min = 0, max = 10000, value = 0, step = 100),
            numericInput("nFeature_max", "maximum genes per spot",
                         min = 0, max = 10000, value = 10000, step = 100),
            br(),
            actionButton("filter_submit_seurat", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("Filtering Statistic"),
            helpText("Describe the distribution of molecular counts and gene counts detected across spots"),
            tabsetPanel(type = "tabs",
                        tabPanel("UMI",
                                 plotOutput("count_filter_seurat", height = "500px"),
                                 sliderInput("nCount_2", "Change parameter cutoffs :",
                                             min = 0, max = 100000, value = c(0,100000), step = 500, width = "800px")),
                        tabPanel("Gene",
                                 plotOutput("gene_filter_seurat", height = "500px"),
                                 sliderInput("nFeature_2", "Change parameter cutoffs :",
                                             min = 0, max = 10000, value = c(0,10000), step = 100, width = "800px"))),
            h3("Filtering Result"),
            helpText("Describe the dimensions of the normalized count matrix in spot level after performing filtering"),
            DT::dataTableOutput("summary_table_filter_seurat")
          )
        )
      ),
      "Statistical Analysis",
      tabPanel(
        title = "Dimension Reduction",
        h3("Clustering Analysis"),
        helpText("Group spots with similar expression profiles into distinct clusters"),
        tabsetPanel(type = "tabs",
                    tabPanel("UMAP", plotOutput("umap_seurat", height = "500px")),
                    tabPanel("Heatmap", plotOutput("heatmap_seurat", height = "800px"))
        )
      ),
      tabPanel(
        title = "Differential expression",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            h3("Differential Expression"),
            helpText("Identify spatially variable genes"),
            strong("Select pairwise comparison :"), p(''),
            numericInput("cluster_1", "Cluster 1",
                         min = 0, max = 30, value = 0, step = 1),
            numericInput("cluster_2", "Cluster 2",
                         min = 0, max = 30, value = 1, step = 1),
            strong("Select significance thresholds :"), p(''),
            numericInput("p_value_seurat", "Adjusted P-value",
                         min = 0, max = 1, value = 0.05, step = 0.001),
            numericInput("logfc_seurat", "Log2FC of Cluster1/Cluster2",
                         min = 0, max = 10, value = 1, step = 0.1),
            actionButton("expre_submit_seurat", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("Differential Expression Result"),
            tabsetPanel(type = "tabs",
                        tabPanel("Summary",
                                 h3("Difference Overview"),
                                 helpText("Describe the number of differential genes"),
                                 DT::dataTableOutput("summary_de_table_seurat"),
                                 h3("Differential Expression Genes"),
                                 helpText("Describe statistical significance and magnitude for each differential gene"),
                                 DT::dataTableOutput("de_table_seurat")
                        ),
                        tabPanel("Visualization",
                                 h3("Gene of interest"),
                                 helpText("Visualize differential gene expression patterns"),
                                 strong("Enter gene symbol :"), p(''),
                                 textInput(inputId = "gene_name_seurat",
                                           label = NULL,
                                           placeholder = "e.g. NLGN1"),
                                 plotOutput("gene_of_interest_seurat",height = "500px")))
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
            helpText("Interpret DEG data by identifying predefined functional GO gene sets"),
            pickerInput("species_seurat", "Select species :",
                        choices = c("Human","Mouse"),
                        selected = "Human"),
            pickerInput("subontology_seurat", "Select subontology :",
                        choices = c("MF","BP","CC","ALL"),
                        selected = "BP"),
            numericInput("p_value_gsea_seurat", "Select P-value threshold :",
                         min = 0, max = 1, value = 0.05, step = 0.001),
            br(),
            actionButton("gsea_submit_seurat", "Submit")
          ),
          mainPanel(
            width = 9,
            h3("GO GSEA Result"),
            tabsetPanel(type = "tabs",
                        tabPanel("Summary",
                                 h3("Enriched GO terms"),
                                 helpText("Describe significantly enriched GO terms"),
                                 DT::dataTableOutput("gsea_table_seurat")),
                        tabPanel("Visualization",
                                 h3("Visualization Plot"),
                                 helpText("Visualize functional enrichment results"),
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Dot Plot",
                                                      dropdownButton(
                                                        numericInput("category_dot_seurat", "Select the number of terms displayed :",
                                                                     min = 1, max = 100, value = 20, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("dot_plot_seurat", height = "800px")),
                                             tabPanel("Network Plot",
                                                      dropdownButton(
                                                        numericInput("category_network_seurat", "Select the number of terms displayed :",
                                                                     min = 1, max = 50, value = 5, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("network_seurat", height = "700px")),
                                             tabPanel("Ridgeline plot",
                                                      dropdownButton(
                                                        numericInput("category_ridgeline_seurat", "Select the number of terms displayed :",
                                                                     min = 1, max = 100, value = 20, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("ridgeline_plot_seurat", height = "800px")),
                                             tabPanel("Enrichment score",
                                                      dropdownButton(
                                                        numericInput("category_score_seurat", "Select the number of terms displayed :",
                                                                     min = 1, max = 20, value = 3, step = 1),
                                                        icon = icon("gear"),
                                                        width = "300px",
                                                        tooltip = tooltipOptions(title = "Click to change the number of terms")
                                                      ),
                                                      plotOutput("score_plot_seurat")))))

          )
        )
      )
    )
  )
)
