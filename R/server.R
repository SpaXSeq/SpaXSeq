
# define function for visualizing study design using alluvial plot
plot_alluvium <- function(object, factors_of_interest)
{
  pData(object) %>% dplyr::select(factors_of_interest) %>% ggalluvial::to_lodes_form() %>%
    ggplot(aes(x = x, stratum = stratum, alluvium = alluvium, fill = stratum, label = stratum)) +
    ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "NA") +
    ggalluvial::geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Metadata") +
    ylab("Frequency")
}


# define function for visualizing segment QC using violin plot
QC_violin <- function(object = NULL, fill_by = NULL, annotation = NULL, thr = NULL, scale_trans = NULL)
{
  plt <- ggplot(object, aes_string(x = paste0("unlist(`", fill_by, "`)"),
                                   y = paste0("unlist(`", annotation, "`)"),
                                   fill = paste0("unlist(`", fill_by, "`)"))) +
    geom_violin(scale = "width") +
    geom_jitter(shape=20, size = 1.8, position = position_jitter(width = 0.2)) +
    geom_hline(yintercept = thr, lty = "dashed", color = "black") +
    theme_classic() + guides(fill = "none") +
    labs(x = fill_by, y = "values of segments", title = annotation)
  if(!is.null(scale_trans)){
    plt <- plt + scale_y_continuous(trans = scale_trans)
  }
  plt
}

server <- function(input, output, session) {

  ####################### GeoMx workflow ##########################

  ### Upload Data Tab ###

  # get file paths
  volumes = getVolumes()
  shinyDirChoose(input, 'DCC_dir', roots = volumes())
  shinyDirChoose(input, 'PKC_dir', roots = volumes())
  shinyFileChoose(input, 'anno_file', roots = volumes())


  # define file paths
  DCCFiles <- reactive({
    dir(parseDirPath(roots = volumes(), input$DCC_dir), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
  })

  PKCFiles <- reactive({
    dir(parseDirPath(roots = volumes(), input$PKC_dir), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
  })

  SampleAnnotationFile <- reactive({
    as.character(parseFilePaths(roots = volumes(), input$anno_file)$datapath)
  })


  # load three file types to creat GeomMxset object
  demoData_raw <- reactive({
    input$upload_button
    isolate(readNanoStringGeoMxSet(dccFiles = DCCFiles(),
                                   pkcFiles = PKCFiles(),
                                   phenoDataFile = SampleAnnotationFile(),
                                   phenoDataSheet = "Sheet1",
                                   phenoDataDccColName = "Sample_ID"))
  })


  # obtain modules used
  modules <- reactive({
    gsub(".pkc", "", annotation(demoData_raw()))
  })


  # update the options according to annotation file
  observeEvent(input$upload_button,{
    updatePickerInput(session, inputId = "factors_of_interest_biology",
                      choices = colnames(pData(demoData_raw())))
    updatePickerInput(session, inputId = "factors_of_interest",
                      choices = colnames(pData(demoData_raw())))
    updatePickerInput(session, inputId = "factor_of_interest_segqc",
                      choices = colnames(pData(demoData_raw())))
    updateSelectInput(session, inputId = "panel",
                      choices = modules())
    updateSelectInput(session, inputId = "factor_of_interest_filter",
                      choices = colnames(pData(demoData_raw())))
    updatePickerInput(session, inputId = "factor_of_interest_color",
                      choices = colnames(pData(demoData_raw())))
    updatePickerInput(session, inputId = "factor_of_interest_shape",
                      choices = colnames(pData(demoData_raw())))
  })


  # study overview table: the dimensions of dataset
  output$summary_table <- DT::renderDataTable({
    req(input$upload_button)
    datatable(data.frame(Features = nrow(fData(demoData_raw())),
                         Samples = nrow(pData(demoData_raw()))))
  })


  # study design alluvial plot: visualize sample types
  output$alluvium_plot <- renderPlot({
    req(input$factors_of_interest)
    plot_alluvium(demoData_raw(), input$factors_of_interest)
  })


  ### Segment QC Tab ###

  # shift expression counts with a value of 0 to 1
  demoData <- reactive({
    shiftCountsOne(demoData_raw(), useDALogic = TRUE)
  })


  # add segment QC flags according to parameter cutoffs
  demoData_seg_qc <- reactive({
    input$segmentQC_submit
    isolate(setSegmentQCFlags(demoData(), qcCutoffs = list(minSegmentReads = input$Reads, percentTrimmed = input$Trimmed, percentStitched = input$Stitched,
                                                           percentAligned = input$Aligned, percentSaturation = input$Saturation, minNegativeCount = input$NegativeCount,
                                                           maxNTCCount = input$NTCCount, minNuclei = input$Nuclei, minArea = input$Area)))
  })


  # segment QC statistic violin plots: visualize the distribution of segments
  output$raw_reads <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()), input$factor_of_interest_segqc, "Raw", input$Reads_2)})

  output$reads_trimmed <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "Trimmed (%)", input$Trimmed_2)})

  output$reads_stitched <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "Stitched (%)", input$Stitched_2)})

  output$reads_aligned <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "Aligned (%)", input$Aligned_2)})

  output$sequencing_saturation <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "Saturated (%)", input$Saturation_2)})

  output$nuclei <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "nuclei", input$Nuclei_2, scale_trans = "log10")})

  output$area <- renderPlot({
    req(input$factor_of_interest_segqc)
    QC_violin(sData(demoData_seg_qc()),input$factor_of_interest_segqc, "area", input$Area_2, scale_trans = "log10")})

  output$NTC_count <- DT::renderDataTable({
    req(input$factor_of_interest_segqc)
    datatable(as.data.frame(table(NTC_Count = sData(demoData_seg_qc())$NTC)), colnames = c("NTC Count", "Number of Segments"), caption = "NTC Count Summary")})

  output$negative_control_counts <- renderPlot({
    req(input$factor_of_interest_segqc)
    negativeGeoMeans <- esBy(negativeControlSubset(demoData()),
                             GROUP = "Module",
                             FUN = function(x) { assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") })
    negCols <- paste0("NegGeoMean_", colnames(negativeGeoMeans))
    p_data_demo <- pData(demoData())
    p_data_demo[, negCols] <- negativeGeoMeans
    QC_violin(p_data_demo,input$factor_of_interest_segqc, paste0("NegGeoMean_", input$panel), input$NegativeCount_2, scale_trans = "log10")})


  # segment QC result table: summary the number of the segments
  output$summary_seg <- DT::renderDataTable({
    req(input$factor_of_interest_segqc)
    QCResults <- protocolData(demoData_seg_qc())[["QCFlags"]]
    QC_Summary <- data.frame(Pass = colSums(!QCResults), Warning = colSums(QCResults))
    QCResults$QCStatus <- apply(QCResults, 1L, function(x) {ifelse(sum(x) == 0L, "PASS", "WARNING")})
    QC_Summary["TOTAL", ] <- c(sum(QCResults[, "QCStatus"] == "PASS"), sum(QCResults[, "QCStatus"] == "WARNING"))
    datatable(QC_Summary, caption = "QC Summary Table for each segment")
  })


  # remove low-performing segments
  demoData_seg_qc_pass <- reactive({
    req(input$segmentQC_submit)
    QCResults <- protocolData(demoData_seg_qc())[["QCFlags"]]
    QCResults$QCStatus <- apply(QCResults, 1L, function(x) {ifelse(sum(x) == 0L, "PASS", "WARNING")})
    demoData_seg_qc()[, QCResults$QCStatus == "PASS"]
  })


  ### Probe QC Tab ###

  # add probe QC flags according to parameter cutoffs
  demoData_pro_qc <- reactive({
    input$probeQC_submit
    isolate(setBioProbeQCFlags(demoData_seg_qc_pass(),
                               qcCutoffs = list(minProbeRatio = input$ProbeRatio, percentFailGrubbs = input$FailGrubbs),
                               removeLocalOutliers = input$LocalOutliers))
  })


  # probe QC result table: summary the number of the probes
  output$summary_probe <- DT::renderDataTable({
    req(input$segmentQC_submit)
    ProbeQCResults <- fData(demoData_pro_qc())[["QCFlags"]]
    qc_df <- data.frame(Pass = sum(!ProbeQCResults$LowProbeRatio & !ProbeQCResults$GlobalGrubbsOutlier),
                        LowProbeRatio = sum(ProbeQCResults$LowProbeRatio),
                        GlobalGrubbsOutlier = sum(ProbeQCResults$GlobalGrubbsOutlier),
                        LocalGrubbsOutlier = sum(rowSums(ProbeQCResults[, -2:-1]) > 0 & !ProbeQCResults$GlobalGrubbsOutlier))
    datatable(qc_df, caption = "QC Summary Table for each probe")
  })


  # remove low-performing probes
  demoData_pro_qc_pass <- reactive({
    req(input$probeQC_submit)
    subset(demoData_pro_qc(), fData(demoData_pro_qc())[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
             fData(demoData_pro_qc())[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
  })


  # aggregate probe counts to gene level
  target_demoData <- reactive({
    aggregateCounts(demoData_pro_qc_pass())
  })


  # study overview table: the dimensions of dataset
  output$summary_table_gene <- DT::renderDataTable({
    req(input$probeQC_submit)
    datatable(data.frame(Features = nrow(fData(target_demoData())),
                         Samples = nrow(pData(target_demoData()))))
  })


  ### Filtering Tab ###

  # calculate LOQ for each module used
  LOQ_Mat <- reactive({
    LOQ <- data.frame(row.names = colnames(target_demoData()))
    for(module in modules()) {
      vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
      LOQ[, module] <- pmax(2, pData(target_demoData())[, vars[1]] * pData(target_demoData())[, vars[2]] ^ 2)
    }
    LOQ_Mat <- c()
    for(module in modules()) {
      ind <- fData(target_demoData())$Module == module
      Mat_i <- t(esApply(target_demoData()[ind, ], MARGIN = 1, FUN = function(x) { x > LOQ[, module]}))
      LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
    }
    LOQ_Mat <- LOQ_Mat[fData(target_demoData())$TargetName, ]
    LOQ_Mat
  })


  # segment filtering: calculate gene detection rate in each segment based on LOQ
  p_data <- reactive({
    p_data <- pData(target_demoData())
    p_data$GenesDetected <- colSums(LOQ_Mat(), na.rm = TRUE)
    p_data$GeneDetectionRate <- p_data$GenesDetected / nrow(target_demoData())
    p_data$DetectionThreshold <- cut(p_data$GeneDetectionRate,
                                     breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
                                     labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
    p_data
  })


  # segment filtering: segment gene detection barplot
  output$seg_filter_plot <- renderPlot({
    req(input$probeQC_submit)
    ggplot(p_data(), aes(x = DetectionThreshold)) +
      geom_bar(aes_string(fill = paste0("unlist(`", input$factor_of_interest_filter, "`)"))) +
      geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
      theme_bw() +
      scale_x_discrete(drop=FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Gene Detection Rate", y = "Number of Segments", fill = input$factor_of_interest_filter)
  })


  # filter out low-signal segments
  target_demoData_seg <- reactive({
    input$filter_submit
    isolate(target_demoData()[, p_data()$GeneDetectionRate >= input$seg_thr/100])
  })


  # Gene filtering: calculate global gene detection rate based on LOQ
  f_data <- reactive({
    LOQ_Mat_gene <- LOQ_Mat()[, colnames(target_demoData_seg())]
    f_data <- fData(target_demoData_seg())
    f_data$DetectedSegments <- rowSums(LOQ_Mat_gene, na.rm = TRUE)
    f_data$DetectionRate <- f_data$DetectedSegments / nrow(pData(target_demoData_seg()))
    f_data
  })


  # Gene filtering: gene detection rate barplot
  plot_detect <- reactive({
    plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
    plot_detect$Number <- unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                        function(x) {sum(f_data()$DetectionRate >= x)}))
    plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData_seg()))
    rownames(plot_detect) <- plot_detect$Freq
    plot_detect
  })

  output$gene_filter_plot <- renderPlot({
    req(input$probeQC_submit)
    ggplot(plot_detect(), aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
                vjust = 1.6, color = "black", size = 4) +
      scale_fill_gradient2(low = "orange2", mid = "lightblue", high = "dodgerblue3",
                           midpoint = 0.65, limits = c(0,1), labels = scales::percent) +
      theme_bw() +
      scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = expansion(mult = c(0, 0))) +
      labs(x = "% of Segments", y = "Genes Detected, % of Panel > LOQ")
  })


  # filter out low-signal genes
  target_demoData_gene <- reactive({
    input$filter_submit
    isolate(target_demoData_seg()[f_data()$DetectionRate >= input$gene_thr/100,])
  })


  # filtering result table: the dimensions of dataset
  output$summary_table_filter <- DT::renderDataTable({
    req(input$filter_submit)
    datatable(data.frame(Features = nrow(fData(target_demoData_gene())),
                         Samples = nrow(pData(target_demoData_gene()))))
  })


  ### Normalization Tab ###

  # coerce to SpatialExperiment object to use functions in standR package
  spe <- reactive({
    req(input$filter_submit)
    spe <- as.SpatialExperiment(target_demoData_gene(), normData = "exprs", forceRaw = TRUE)
    assayNames(spe) <- "counts"
    raw_count <- assay(spe, "counts")
    assay(spe, "logcounts") <- edgeR::cpm(raw_count, log = TRUE)
    col_data <- colData(spe)
    col_data$biology <- paste0(col_data[,input$factors_of_interest_biology[1]], "_", col_data[,input$factors_of_interest_biology[2]])
    colData(spe) <- col_data
    spe
  })


  # update the options according to biological variables of interest
  observeEvent(input$filter_submit,{
    updatePickerInput(session, inputId = "contrast_1",
                      choices = unique(colData(spe())$biology))
    updatePickerInput(session, inputId = "contrast_2",
                      choices = unique(colData(spe())$biology))
  })


  # perform normalization
  spe_norm <- reactive({
    input$nor_submit
    isolate(geomxNorm(spe(), method = input$nor_method))
  })


  # perform batch correction using RUV4 method
  spe_norm_ncg <- reactive({
    findNCGs(spe_norm(), batch_name = "slide name", top_n = 300)
  })

  spe_ruv <- reactive({
    input$nor_submit
    isolate(geomxBatchCorrection(spe_norm_ncg(), factors = "biology", NCGs = metadata(spe_norm_ncg())$NCGs, k = input$k_best))
  })


  # visualization
  # RLE: Raw data
  output$nor_before <- renderPlot({
    req(input$filter_submit)
    plotRLExpr(spe(), ordannots = "slide name", color = `slide name`) + ggtitle("Raw")
  })

  # RLE: Normalized data
  output$nor_after <- renderPlot({
    req(input$filter_submit)
    input$nor_submit
    isolate(plotRLExpr(spe_norm(), assay = 2, ordannots = "slide name", color = `slide name`) + ggtitle(input$nor_method))
  })

  # RLE: batch corrected data
  output$nor_batch <- renderPlot({
    req(input$filter_submit)
    input$nor_submit
    isolate(plotRLExpr(spe_ruv(), assay = 2, ordannots = "slide name", color = `slide name`) + ggtitle(paste0("RUV4 k = ", input$k_best)))
  })


  # visualization
  # PCA: Raw data
  output$pca_before <- renderPlot({
    req(input$filter_submit)
    drawPCA(spe(), assay = 1,
            color = `slide name`,
            shape = colData(spe())[,input$factors_of_interest_biology[2]]) +
      labs(shape = input$factors_of_interest_biology[2]) +
      ggtitle("Raw")
  })

  # PCA: Normalized data
  spe_norm_pca <- reactive({
    set.seed(100)
    scater::runPCA(spe_norm())
  })

  pca_results_norm <- reactive({
    reducedDim(spe_norm_pca(), "PCA")
  })

  output$pca_after <- renderPlot({
    req(input$filter_submit)
    input$nor_submit
    isolate(drawPCA(spe_norm_pca(), precomputed = pca_results_norm(),
                    color = `slide name`,
                    shape = colData(spe_norm_pca())[,input$factors_of_interest_biology[2]]) +
              labs(shape = input$factors_of_interest_biology[2]) +
              ggtitle(input$nor_method))
  })

  # PCA: batch corrected data
  spe_batch_pca <- reactive({
    scater::runPCA(spe_ruv())
  })

  pca_results_batch <- reactive({
    reducedDim(spe_batch_pca(), "PCA")
  })

  output$pca_batch <- renderPlot({
    req(input$filter_submit)
    input$nor_submit
    isolate(drawPCA(spe_batch_pca(), precomputed = pca_results_batch(),
                    color = `slide name`,
                    shape = colData(spe_batch_pca())[,input$factors_of_interest_biology[2]]) +
              labs(shape = input$factors_of_interest_biology[2]) +
              ggtitle(paste0("RUV4 k = ", input$k_best)))
  })


  ### Dimension Reduction Tab ###

  # PCA
  output$pca_plot <- renderPlot({
    req(input$nor_submit)
    drawPCA(spe_batch_pca(), precomputed = pca_results_batch(),
            color = colData(spe_batch_pca())[,input$factor_of_interest_color],
            shape = colData(spe_batch_pca())[,input$factor_of_interest_shape]) +
      labs(color = input$factor_of_interest_color, shape = input$factor_of_interest_shape)
  })


  # MDS
  output$mds <- renderPlot({
    req(input$nor_submit)
    standR::plotMDS(spe_batch_pca(), assay = 2,
                    color = colData(spe_batch_pca())[,input$factor_of_interest_color],
                    shape = colData(spe_batch_pca())[,input$factor_of_interest_shape]) +
      labs(color = input$factor_of_interest_color, shape = input$factor_of_interest_shape)
  })


  # UMAP
  spe_batch_umap <- reactive({
    set.seed(100)
    scater::runUMAP(spe_batch_pca(), dimred = "PCA")
  })

  output$umap <- renderPlot({
    req(input$nor_submit)
    plotDR(spe_batch_umap(), dimred = "UMAP",
           color = colData(spe_batch_umap())[,input$factor_of_interest_color],
           shape = colData(spe_batch_umap())[,input$factor_of_interest_shape]) +
      labs(color = input$factor_of_interest_color, shape = input$factor_of_interest_shape)
  })


  # TSNE
  spe_batch_tsne <- reactive({
    set.seed(100)
    scater::runTSNE(spe_batch_pca(), dimred = "PCA")})

  output$tsne <- renderPlot({
    req(input$nor_submit)
    plotDR(spe_batch_tsne(), dimred = "TSNE",
           color = colData(spe_batch_tsne())[,input$factor_of_interest_color],
           shape = colData(spe_batch_tsne())[,input$factor_of_interest_shape]) +
      labs(color = input$factor_of_interest_color, shape = input$factor_of_interest_shape)
  })


  ### Differential Expression Tab ###

  # coerce to DGEList object to use limma-voom pipeline
  dge_all <- reactive({
    dge <- SE2DGEList(spe_ruv())
    calcNormFactors(dge)
  })


  # establish a design matrix
  design <- reactive({
    weight = ""
    for (k in 1:input$k_best) {
      ruv_k = paste0(" + ruv_W", k)
      weight = paste0(weight, ruv_k)
    }
    design <- model.matrix(as.formula(paste0("~0 + biology", weight)), data = colData(spe_ruv()))
    colnames(design) <- gsub("^biology","",colnames(design))
    colnames(design) <- gsub(" ","_",colnames(design))
    design
  })


  # establish pairwise comparison
  contr_number <- reactive({
    input$expre_submit
    isolate(paste(gsub(" ","_",input$contrast_1), "-", gsub(" ","_",input$contrast_2)))
  })

  contr_matrix <- reactive({
    req(input$expre_submit)
    makeContrasts(contrasts = contr_number(), levels = colnames(design()))
  })


  # perform differential expression analysis
  efit <- reactive({
    req(input$expre_submit)
    v <- voom(dge_all(), design())
    fit <- lmFit(v, design())
    fit_contrast <- contrasts.fit(fit, contrasts = contr_matrix())
    eBayes(fit_contrast)
  })


  # extract differential expression genes table
  de_genes_table <- reactive({
    input$expre_submit
    isolate(topTable(efit(), sort.by = "P", n = Inf, p.value = input$p_value, lfc = input$logfc))
  })

  genes_table <- reactive({
    topTable(efit(), sort.by = "P", n = Inf)
  })


  # difference overview table: summary the number of differential genes
  output$summary_de_table <- DT::renderDataTable({
    req(input$expre_submit)
    datatable(data.frame(Significant_genes = nrow(de_genes_table()),
                         Up_regulated_genes = length(which(de_genes_table()$logFC > 0)),
                         Down_regulated_genes = length(which(de_genes_table()$logFC < 0)),
                         non_significant_genes = nrow(genes_table()) - nrow(de_genes_table())),
              colnames = c("Significant genes", "Up-regulated genes", "Down-regulated genes", "Non-significant genes"))
  })


  # significantly differential expression genes table
  output$de_table <- DT::renderDataTable({
    req(input$expre_submit)
    datatable(de_genes_table()[,c("logFC", "AveExpr", "P.Value", "adj.P.Val")])
  })


  # subset count matrix according to contrasts
  spe_ruv_contrast <- reactive({
    spe_ruv()[de_genes_table()$TargetName,colData(spe_ruv())$biology == input$contrast_1 |
                colData(spe_ruv())$biology == input$contrast_2]
  })


  # visualization
  # Volcano plot
  output$volcano <- renderPlot({
    req(input$expre_submit)
    ggplot(genes_table(), aes(x=logFC, y=-log10(adj.P.Val))) +
      geom_point(aes(colour = "a")) +
      geom_point(data = subset(de_genes_table(), logFC < 0), aes(colour = "b")) +
      geom_point(data = subset(de_genes_table(), logFC > 0), aes(colour = "c")) +
      labs(title = "Volcano", x= "Log2FoldChange", y = "-Log10(p.adj)") +
      theme_bw() +
      geom_vline(xintercept = -input$logfc, linetype = "dashed", color = "grey50", size = 0.5) +
      geom_vline(xintercept = input$logfc, linetype = "dashed", color = "grey50", size = 0.5) +
      geom_hline(yintercept = -log10(input$p_value), linetype = "dashed", color = "grey50", size = 0.5) +
      geom_label_repel(data= subset(de_genes_table(), logFC < 0),
                       aes(label=row.names(subset(de_genes_table(), logFC < 0)), colour="b"),
                       show.legend = FALSE, size=2.5) +
      geom_label_repel(data= subset(de_genes_table(), logFC > 0),
                       aes(label=row.names(subset(de_genes_table(), logFC > 0)), colour="c"),
                       show.legend = FALSE, size=2.5) +
      scale_colour_manual(values=c("grey60", "dodgerblue", "coral2"), labels=c("Not DE","Down","Up"), name="Significance")
  })


  # MA plot
  output$ma <- renderPlot({
    req(input$expre_submit)
    ggplot(genes_table(), aes(x=AveExpr, y=logFC)) +
      geom_point(aes(colour = "a")) +
      geom_point(data = subset(de_genes_table(), logFC < 0), aes(colour = "b")) +
      geom_point(data = subset(de_genes_table(), logFC > 0), aes(colour = "c")) +
      labs(title = "MA", x= "Average Expression (log2CPM)", y = "Log2FoldChange") +
      theme_bw() +
      geom_hline(yintercept = -input$logfc, linetype = "dashed", color = "grey50", size = 0.5) +
      geom_hline(yintercept = input$logfc, linetype = "dashed", color = "grey50", size = 0.5) +
      geom_label_repel(data= subset(de_genes_table(), logFC < 0),
                       aes(label=row.names(subset(de_genes_table(), logFC < 0)), colour="b"),
                       show.legend = FALSE, size=2.5) +
      geom_label_repel(data= subset(de_genes_table(), logFC > 0),
                       aes(label=row.names(subset(de_genes_table(), logFC > 0)), colour="c"),
                       show.legend = FALSE, size=2.5) +
      scale_colour_manual(values=c("grey60", "dodgerblue", "coral2"), labels=c("Not DE","Down","Up"), name="Significance")
  })


  # Heatmap plot
  output$heatmap <- renderPlot({
    req(input$expre_submit)
    pheatmap(assay(spe_ruv_contrast(), "logcounts"),
             color = colorRampPalette(c('dodgerblue','white','coral2'))(100),
             border_color = NA,
             scale = "row",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             cutree_rows = 2,
             cutree_cols = 2,
             legend = TRUE,
             clustering_method = "average",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             show_rownames = FALSE,
             show_colnames = FALSE,
             annotation_col = as.data.frame(colData(spe_ruv_contrast())[,c(input$factors_of_interest_biology[1],input$factors_of_interest_biology[2])]),
             annotation_legend = TRUE,
             annotation_names_col = TRUE)
  })


  # Genes of interest
  output$violin_gene_of_interest <- renderPlot({
    req(input$expre_submit)
    ggplot(as.data.frame(colData(spe_ruv_contrast())), aes(x = biology, fill = biology,
                                                           y = as.numeric(assay(spe_ruv_contrast(), "logcounts")[input$gene_name,]))) +
      geom_violin() +
      geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
      theme_bw() +
      labs(x = "Groups", y = "Expression", title = input$gene_name) +
      scale_y_continuous(trans = "log2") +
      scale_fill_manual(values = c("dodgerblue", "coral2"), name = "Groups")
  })


  ### Gene Set Enrichment Analysis Tab ###

  # prepare a ranked genelist
  geneList <- reactive({
    d = de_genes_table()[,c("TargetName", "logFC")]
    geneList = d[,2]
    names(geneList) = as.character(d[,1])
    sort(geneList, decreasing = TRUE)
  })


  # perform GO GSEA
  go_gsea <- reactive({
    req(input$gsea_submit)
    if (input$species == "Human") {
      gseGO(geneList = geneList(),
            OrgDb = org.Hs.eg.db,
            keyType = "SYMBOL",
            ont = input$subontology,
            minGSSize = 5,
            maxGSSize = 500,
            pvalueCutoff = input$p_value_gsea,
            verbose = FALSE)
    } else {
      gseGO(geneList = geneList(),
            OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL",
            ont = input$subontology,
            minGSSize = 5,
            maxGSSize = 500,
            pvalueCutoff = input$p_value_gsea,
            verbose = FALSE)
    }
  })


  # significantly enriched GO terms table
  output$gsea_table <- DT::renderDataTable({
    req(input$gsea_submit)
    datatable(go_gsea()@result)
  })


  # visualization
  # dot plot
  output$dot_plot <- renderPlot({
    req(input$gsea_submit)
    dotplot(go_gsea(), showCategory = input$category_dot)
  })


  # network plot
  output$network <- renderPlot({
    req(input$gsea_submit)
    cnetplot(go_gsea(), showCategory = input$category_network, foldChange = geneList(), circular = TRUE, colorEdge = TRUE)
  })


  # ridgeline plot
  output$ridgeline_plot <- renderPlot({
    req(input$gsea_submit)
    ridgeplot(go_gsea(), showCategory = input$category_ridgeline)
  })


  # running score and preranked list plot
  output$score_plot <- renderPlot({
    req(input$gsea_submit)
    gseaplot2(go_gsea(), geneSetID = 1:input$category_score)
  })


  ### Cellular Deconvolution Tab ###

  # download a cell profile matrix
  profile_matrix <- reactive({
    req(input$dec_submit)
    input$dec_submit
    isolate(download_profile_matrix(species = input$species_2,
                                    age_group = input$age_group,
                                    matrixname = input$matrix_name))
  })


  # prepare deconvolution function data
  # obtain negative control probe names
  neg_probes <- reactive({
    negativeProbefData <- subset(fData(target_demoData_seg()), CodeClass == "Negative")
    unique(negativeProbefData$TargetName)
  })


  # include negative control probes in the normalized expression data
  target_demoData_neg <- reactive({
    target_demoData_neg <- target_demoData_seg()[f_data()$DetectionRate >= input$gene_thr/100 |
                                                   f_data()$TargetName %in% neg_probes(), ]
    pData(target_demoData_neg)$biology = paste0(pData(target_demoData_neg)[,input$factors_of_interest_biology[1]],
                                                "_", pData(target_demoData_neg)[,input$factors_of_interest_biology[2]])
    spe <- as.SpatialExperiment(target_demoData_neg, normData = "exprs", forceRaw = TRUE)
    assayNames(spe) <- "counts"
    raw_count <- assay(spe, "counts")
    assay(spe, "logcounts") <- edgeR::cpm(raw_count, log = TRUE)
    spe_norm <- geomxNorm(spe, method = input$nor_method)
    spe_norm <- findNCGs(spe_norm, batch_name = "slide name", top_n = 300)
    spe_ruv <- geomxBatchCorrection(spe_norm, factors = "biology", NCGs = metadata(spe_norm)$NCGs, k = input$k_best)
    norm_count <- assay(spe_ruv, "logcounts")
    assayDataElement(target_demoData_neg, elt = "exprs_norm") = norm_count
    target_demoData_neg
  })


  # Perform basic deconvolution
  res <- reactive({
    runspatialdecon(object = target_demoData_neg(),
                    norm_elt = "exprs_norm",
                    raw_elt = "exprs",
                    X = profile_matrix(),
                    align_genes = TRUE)
  })


  # visualization
  # heatmap plot
  output$dec_heatmap <- renderPlot({
    req(input$dec_submit)
    pheatmap(t(pData(res())$beta),
             color = colorRampPalette(c("white", "#FED98E", "#FE9929", "#D95F0E", "#993404"))(100),
             border_color = NA,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             legend = TRUE,
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = pData(target_demoData_neg())[,c(input$factors_of_interest_biology[1],input$factors_of_interest_biology[2])],
             annotation_legend = TRUE,
             annotation_names_col = TRUE)
  })


  # bar plot
  output$dec_barplot <- renderPlot({
    req(input$dec_submit)
    samples_subset <- colnames(target_demoData_neg())[pData(target_demoData_neg())$biology %in% c(input$contrast_1, input$contrast_2)]
    subset_prop <- t(res()$prop_of_all)[,samples_subset]
    group = pData(target_demoData_neg())[colnames(subset_prop),"biology"]
    colnames(subset_prop) = paste0(group, "_", colnames(subset_prop))
    subset_prop %>%
      as.data.frame() %>%
      rownames_to_column("CellTypes") %>%
      gather(samples, prop, -CellTypes) %>%
      ggplot(aes(samples, prop, fill = CellTypes)) +
      geom_bar(stat = "identity", position = "stack", color = "black", width = .7) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "bottom")
  })


  ####################### Visium workflow ##########################

  ### Upload Data Tab ###

  # Get file path
  shinyDirChoose(input, 'folder', roots = volumes())


  # creat seurat object containing both the expression data and the associated image data.
  tissue <- reactive({
    req(input$upload_button_seurat)
    Load10X_Spatial(data.dir = parseDirPath(roots = volumes(), input$folder))
  })


  # study overview table: the dimensions of dataset
  output$summary_table_seurat <- DT::renderDataTable({
    req(input$upload_button_seurat)
    datatable(data.frame(Features = nrow(tissue()),
                         Samples = ncol(tissue())))
  })


  ### Filtering & Normalization Tab ###

  # the distribution of molecular counts
  output$count_filter_seurat <- renderPlot({
    req(input$upload_button_seurat)
    plot1 <- VlnPlot(tissue(), features = "nCount_Spatial", pt.size = 0.1) + NoLegend() +
      geom_hline(yintercept = input$nCount_2[1], lty = "dashed", color = "black") +
      geom_hline(yintercept = input$nCount_2[2], lty = "dashed", color = "black")
    plot2 <- SpatialFeaturePlot(tissue(), features = "nCount_Spatial") + theme(legend.position = "right")
    wrap_plots(plot1, plot2)
  })


  # the distribution of gene counts
  output$gene_filter_seurat <- renderPlot({
    req(input$upload_button_seurat)
    plot3 <- VlnPlot(tissue(), features = "nFeature_Spatial", pt.size = 0.1) + NoLegend() +
      geom_hline(yintercept = input$nFeature_2[1], lty = "dashed", color = "black") +
      geom_hline(yintercept = input$nFeature_2[2], lty = "dashed", color = "black")
    plot4 <- SpatialFeaturePlot(tissue(), features = "nFeature_Spatial") + theme(legend.position = "right")
    wrap_plots(plot3, plot4)
  })


  # filter out low-quality spots
  tissue_filter <- reactive({
    req(input$filter_submit_seurat)
    input$filter_submit_seurat
    isolate(subset(tissue(), subset = nFeature_Spatial > input$nFeature_min & nFeature_Spatial < input$nFeature_max &
             nCount_Spatial > input$nCount_min & nCount_Spatial < input$nCount_max))
  })


  # filtering result table: the dimensions of dataset
  output$summary_table_filter_seurat <- DT::renderDataTable({
    req(input$filter_submit_seurat)
    datatable(data.frame(Features = nrow(tissue_filter()),
                         Samples = ncol(tissue_filter())))
  })


  # normalization
  tissue_nor <- reactive({
    SCTransform(tissue_filter(), assay = "Spatial", verbose = FALSE)
  })


  ### Dimension Reduction Tab ###

  # clustering analysis
  tissue_umap <- reactive({
    tissue_dim <- RunPCA(tissue_nor(), assay = "SCT", verbose = FALSE)
    tissue_dim <- FindNeighbors(tissue_dim, reduction = "pca", dims = 1:30)
    tissue_dim <- FindClusters(tissue_dim, verbose = FALSE)
    RunUMAP(tissue_dim, reduction = "pca", dims = 1:30)
  })


  # visualization
  # UMAP plot
  output$umap_seurat <- renderPlot({
    req(input$filter_submit_seurat)
    p1 <- DimPlot(tissue_umap(), reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(tissue_umap(), label = TRUE, label.size = 3)
    p1 + p2
  })


  # Heatmap
  tissue_markers <- reactive({
    FindAllMarkers(tissue_umap(), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
  })

  topgene <- reactive({
    tissue_markers() %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  })

  output$heatmap_seurat <- renderPlot({
    DoHeatmap(tissue_umap(), features = topgene()$gene, size = 3)
  })


  ### Differential expression ###

  # perform differential expression analysis
  de_markers <- reactive({
    FindMarkers(tissue_umap(), ident.1 = input$cluster_1, ident.2 = input$cluster_2,
                logfc.threshold = input$logfc_seurat, return.thresh = input$p_value_seurat)
  })


  # difference overview table: the number of differential genes
  output$summary_de_table_seurat <- DT::renderDataTable({
    datatable(data.frame(Significant_genes = nrow(de_markers()),
                         Up_regulated_genes = length(which(de_markers()$avg_log2FC > 0)),
                         Down_regulated_genes = length(which(de_markers()$avg_log2FC < 0)),
                         non_significant_genes = nrow(tissue_umap()) - nrow(de_markers())),
              colnames = c("Significant genes", "Up-regulated genes", "Down-regulated genes", "Non-significant genes"))
  })


  # differential expression genes table
  output$de_table_seurat <- DT::renderDataTable({
    datatable(de_markers(), rownames = TRUE)
  })


  # gene of interest
  output$gene_of_interest_seurat <- renderPlot({
    SpatialFeaturePlot(tissue_umap(), features = input$gene_name_seurat, alpha = c(0.1, 1))
  })


  ### Gene Set Enrichment Analysis Tab ###
  # prepare genelist
  geneList_seurat <- reactive({
    d = data.frame(gene = rownames(de_markers()), logfc = de_markers()$avg_log2FC)
    geneList = d[,2]
    names(geneList) = as.character(d[,1])
    sort(geneList, decreasing = TRUE)
  })

  # perform GO GSEA
  go_gsea_seurat <- reactive({
    req(input$gsea_submit_seurat)
    if (input$species_seurat == "Human") {
      gseGO(geneList = geneList_seurat(),
            OrgDb = org.Hs.eg.db,
            keyType = "SYMBOL",
            ont = "ALL",
            minGSSize = 5,
            maxGSSize = 500,
            pvalueCutoff = input$p_value_gsea_seurat,
            verbose = FALSE)
    } else {
      gseGO(geneList = geneList_seurat(),
            OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL",
            ont = "ALL",
            minGSSize = 5,
            maxGSSize = 500,
            pvalueCutoff = input$p_value_gsea_seurat,
            verbose = FALSE)
    }
  })


  # GSEA Overview
  output$gsea_table_seurat <- DT::renderDataTable({
    datatable(go_gsea_seurat()@result)
  })


  # visualization
  output$dot_plot_seurat <- renderPlot({
    dotplot(go_gsea_seurat(), showCategory = 20)
  })


  output$network_seurat <- renderPlot({
    p1 <- cnetplot(go_gsea_seurat(), foldChange = geneList_seurat())
    p2 <- cnetplot(go_gsea_seurat(), categorySize = "pvalue", foldChange = geneList_seurat())
    p3 <- cnetplot(go_gsea_seurat(), foldChange = geneList_seurat(), circular = TRUE, colorEdge = TRUE)
    cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
  })

}
