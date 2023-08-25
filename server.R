
# define function for visualizing study design
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

# define function for visualizing segment QC
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
    
    # Get file path
    volumes = getVolumes()
    shinyDirChoose(input, 'DCC_dir', roots = volumes())
    shinyDirChoose(input, 'PKC_dir', roots = volumes())
    shinyFileChoose(input, 'anno_file', roots = volumes())
    
    # define file path
    DCCFiles <- reactive({
        dir(parseDirPath(roots = volumes(), input$DCC_dir), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
    })
    
    PKCFiles <- reactive({
        dir(parseDirPath(roots = volumes(), input$PKC_dir), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
    })
    
    SampleAnnotationFile <- reactive({
        as.character(parseFilePaths(roots = volumes(), input$anno_file)$datapath)
    })
    
    # creat GeomMxset object
    demoData_raw <- reactive({
        input$upload_button
        isolate(readNanoStringGeoMxSet(dccFiles = DCCFiles(),
                                       pkcFiles = PKCFiles(),
                                       phenoDataFile = SampleAnnotationFile(),
                                       phenoDataSheet = "Sheet1",
                                       phenoDataDccColName = "Sample_ID"))})
    
    # modules used
    modules <- reactive({
        gsub(".pkc", "", annotation(demoData_raw()))})
    
    # update the options
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
    
    # study overview
    output$summary_table <- renderDataTable({
        req(input$upload_button)
        datatable(data.frame(Features = nrow(fData(demoData_raw())),
                             Samples = nrow(pData(demoData_raw()))))
    })
    
    # study design
    output$alluvium_plot <- renderPlot({
        req(input$factors_of_interest)
        plot_alluvium(demoData_raw(), input$factors_of_interest)})
    
    # Shift expression counts with a value of 0 to 1
    demoData <- reactive({
        shiftCountsOne(demoData_raw(), useDALogic = TRUE)})
    
    # Segment QC
    demoData_seg_qc <- reactive({
        input$segmentQC_submit
        isolate(setSegmentQCFlags(demoData(), qcCutoffs = list(minSegmentReads = input$Reads, percentTrimmed = input$Trimmed, percentStitched = input$Stitched, 
                                                               percentAligned = input$Aligned, percentSaturation = input$Saturation, minNegativeCount = input$NegativeCount, 
                                                               maxNTCCount = input$NTCCount, minNuclei = input$Nuclei, minArea = input$Area)))})
    # Segment QC Statistic plot
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
    
    output$NTC_count <- renderDataTable({
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
    
    # Segment QC Result
    output$summary_seg <- renderDataTable({
        req(input$factor_of_interest_segqc)
        QCResults <- protocolData(demoData_seg_qc())[["QCFlags"]]
        QC_Summary <- data.frame(Pass = colSums(!QCResults), Warning = colSums(QCResults))
        QCResults$QCStatus <- apply(QCResults, 1L, function(x) {ifelse(sum(x) == 0L, "PASS", "WARNING")})
        QC_Summary["TOTAL", ] <- c(sum(QCResults[, "QCStatus"] == "PASS"), sum(QCResults[, "QCStatus"] == "WARNING"))
        datatable(QC_Summary, caption = "QC Summary Table for each Segment")})
    
    # remove low-performing segments
    demoData_seg_qc_pass <- reactive({
        req(input$segmentQC_submit)
        QCResults <- protocolData(demoData_seg_qc())[["QCFlags"]]
        QCResults$QCStatus <- apply(QCResults, 1L, function(x) {ifelse(sum(x) == 0L, "PASS", "WARNING")})
        demoData_seg_qc()[, QCResults$QCStatus == "PASS"]
    })
    
    # Probe QC
    demoData_pro_qc <- reactive({
        input$probeQC_submit
        isolate(setBioProbeQCFlags(demoData_seg_qc_pass(), 
                                   qcCutoffs = list(minProbeRatio = input$ProbeRatio, percentFailGrubbs = input$FailGrubbs), 
                                   removeLocalOutliers = input$LocalOutliers))})
    
    # Probe QC Result
    output$summary_probe <- renderDataTable({
        req(input$segmentQC_submit)
        ProbeQCResults <- fData(demoData_pro_qc())[["QCFlags"]]
        qc_df <- data.frame(Pass = sum(!ProbeQCResults$LowProbeRatio & !ProbeQCResults$GlobalGrubbsOutlier),
                            LowProbeRatio = sum(ProbeQCResults$LowProbeRatio),
                            GlobalGrubbsOutlier = sum(ProbeQCResults$GlobalGrubbsOutlier),
                            LocalGrubbsOutlier = sum(rowSums(ProbeQCResults[, -2:-1]) > 0 & !ProbeQCResults$GlobalGrubbsOutlier))
        datatable(qc_df, caption = "QC Summary Table for each probe")})
    
    # remove low-performing probes
    demoData_pro_qc_pass <- reactive({
        req(input$probeQC_submit)
        subset(demoData_pro_qc(), fData(demoData_pro_qc())[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
                                  fData(demoData_pro_qc())[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)})
    
    # collapse to targets
    target_demoData <- reactive({
        aggregateCounts(demoData_pro_qc_pass())})
    
    # study overview
    output$summary_table_gene <- renderDataTable({
        req(input$probeQC_submit)
        datatable(data.frame(Features = nrow(fData(target_demoData())),
                             Samples = nrow(pData(target_demoData()))))
    })
    
    # Calculate LOQ per panel used
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
    
    # Segment filtering
    p_data <- reactive({
        p_data <- pData(target_demoData())
        p_data$GenesDetected <- colSums(LOQ_Mat(), na.rm = TRUE)
        p_data$GeneDetectionRate <- p_data$GenesDetected / nrow(target_demoData())
        p_data$DetectionThreshold <- cut(p_data$GeneDetectionRate,
                                         breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
                                         labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
        p_data
    })
    
    # Segment Gene Detection
    output$seg_filter_plot <- renderPlot({
        req(input$probeQC_submit)
        ggplot(p_data(), aes(x = DetectionThreshold)) +
            geom_bar(aes_string(fill = paste0("unlist(`", input$factor_of_interest_filter, "`)"))) +
            geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
            theme_bw() +
            scale_x_discrete(drop=FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            labs(x = "Gene Detection Rate", y = "Number of Segments", fill = input$factor_of_interest_filter)})
    
    # filter out segments
    target_demoData_seg <- reactive({
        input$filter_submit
        isolate(target_demoData()[, p_data()$GeneDetectionRate >= input$seg_thr/100])})
    
    # Gene filtering
    f_data <- reactive({
        LOQ_Mat_gene <- LOQ_Mat()[, colnames(target_demoData_seg())]
        f_data <- fData(target_demoData_seg())
        f_data$DetectedSegments <- rowSums(LOQ_Mat_gene, na.rm = TRUE)
        f_data$DetectionRate <- f_data$DetectedSegments / nrow(pData(target_demoData_seg()))
        f_data
    })
    
    plot_detect <- reactive({
        plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
        plot_detect$Number <- unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                            function(x) {sum(f_data()$DetectionRate >= x)}))
        plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData_seg()))
        rownames(plot_detect) <- plot_detect$Freq
        plot_detect
    })
    
    # Gene Detection Rate
    output$gene_filter <- renderPlot({
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
    
    # filter out genes
    target_demoData_gene <- reactive({
        input$filter_submit
        isolate(target_demoData_seg()[f_data()$DetectionRate >= input$gene_thr/100,])})
    
    # Filtering Result
    output$summary_table_filter <- renderDataTable({
        req(input$filter_submit)
        datatable(data.frame(Features = nrow(fData(target_demoData_gene())),
                             Samples = nrow(pData(target_demoData_gene()))))
    })
    
    # coerce to SpatialExperiment object
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
    
    # update the options
    observeEvent(input$filter_submit,{
        updatePickerInput(session, inputId = "contrast_1",
                          choices = unique(colData(spe())$biology))
        updatePickerInput(session, inputId = "contrast_2",
                          choices = unique(colData(spe())$biology))
    })
    
    # normalization
    spe_norm <- reactive({
        input$nor_submit
        isolate(geomxNorm(spe(), method = input$nor_method))})
    
    # Batch correction
    spe_norm_ncg <- reactive({
        findNCGs(spe_norm(), batch_name = "slide name", top_n = 300)
    })
    
    spe_ruv <- reactive({
        input$nor_submit
        isolate(geomxBatchCorrection(spe_norm_ncg(), factors = "biology", NCGs = metadata(spe_norm_ncg())$NCGs, k = input$k_best))
    })
    
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
        scater::runUMAP(spe_batch_pca(), dimred = "PCA")})
    
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
    
    # Create DGEList object 
    dge_all <- reactive({
        dge <- SE2DGEList(spe_ruv())
        calcNormFactors(dge)
    })
    
    # Establish a design matrix
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
    
    # pairwise comparison
    contr_number <- reactive({
      input$expre_submit
      isolate(paste(gsub(" ","_",input$contrast_1), "-", gsub(" ","_",input$contrast_2)))
    })
    
    contr_matrix <- reactive({
      req(input$expre_submit)
      makeContrasts(contrasts = contr_number(), levels = colnames(design()))
    })
    
    # Differential Expression
    efit <- reactive({
        req(input$expre_submit)
        v <- voom(dge_all(), design()) 
        fit <- lmFit(v, design())
        fit_contrast <- contrasts.fit(fit, contrasts = contr_matrix())
        eBayes(fit_contrast)
    })

    # Difference Overview
    de_genes_table <- reactive({
      input$expre_submit
      isolate(topTable(efit(), sort.by = "P", n = Inf, p.value = input$p_value, lfc = input$logfc))
    })
    
    genes_table <- reactive({
      topTable(efit(), sort.by = "P", n = Inf)
    })
    
    output$summary_de_table <- renderDataTable({
      req(input$expre_submit)
      input$expre_submit
      isolate(datatable(data.frame(Significant_genes = nrow(de_genes_table()),
                             Up_regulated_genes = length(which(de_genes_table()$logFC > 0)),
                             Down_regulated_genes = length(which(de_genes_table()$logFC < 0)),
                             non_significant_genes = nrow(genes_table()) - nrow(de_genes_table())),
                             colnames = c("Significant genes", "Up-regulated genes", "Down-regulated genes", "Non-significant genes"),
                             caption = paste(input$contrast_1,"vs.",input$contrast_2)))
    })
    
    # Differential Expression Gene
    output$de_table <- renderDataTable({
      req(input$expre_submit)
      input$expre_submit
      isolate(datatable(de_genes_table()[,c("logFC", "AveExpr", "P.Value", "adj.P.Val")], caption = paste(input$contrast_1,"vs.",input$contrast_2)))
    })
    
    # Volcano
    output$volcano <- renderPlot({
        req(input$expre_submit)
        ggplot(genes_table(), aes(x=logFC, y=-log10(adj.P.Val))) + 
            geom_point(aes(colour = "a")) + 
            geom_point(data = subset(de_genes_table(), logFC < 0), aes(colour = "b")) + 
            geom_point(data = subset(de_genes_table(), logFC > 0), aes(colour = "c")) +
            labs(title = paste(input$contrast_1,"vs.",input$contrast_2), x= "Log2FoldChange", y = "-Log10(p.adj)") +
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

    # MA
    output$ma <- renderPlot({
        req(input$expre_submit)
        ggplot(genes_table(), aes(x=AveExpr, y=logFC)) + 
            geom_point(aes(colour = "a")) + 
            geom_point(data = subset(de_genes_table(), logFC < 0), aes(colour = "b")) + 
            geom_point(data = subset(de_genes_table(), logFC > 0), aes(colour = "c")) +
            labs(title = paste(input$contrast_1,"vs.",input$contrast_2), x= "Average Expression (log2CPM)", y = "Log2FoldChange") +
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
    
    # Heatmap
    
    # heatmap_table <- reactive({
    #   
    #   spe_ruv()
    #   assay(spe, "logcounts")
    #   
    #   assay(spe_ruv(), "logcounts")[de_genes_table()$TargetName[1:50], colData(spe_ruv())$biology == contrast_1|
    #     
    #   colData(spe_ruv())[colData(spe_ruv())$biology == contrast_1 | colData(spe_ruv())$biology == contrast_2,c("Tissue","Timepoint")]
    # })
    # test2 = colData(spe_ruv)[colData(spe_ruv)$biology == contrast_1 | colData(spe_ruv)$biology == contrast_2,c("Tissue","Timepoint")]
    # test2 = as.data.frame(test2)
    # test = assay(spe_ruv, "logcounts")[de_genes_table$TargetName[1:50], colData(spe_ruv)$biology == contrast_1|
    #                                     +                                   colData(spe_ruv)$biology == contrast_2]
    
    # output$heatmap <- renderPlot({
    #   req(input$expre_submit)
    #   pheatmap(data(),  
    #            color = colorRampPalette(c('dodgerblue','white','coral2'))(100),
    #            border_color = NA,  
    #            scale = "row", 
    #            cluster_rows = TRUE, 
    #            cluster_cols = TRUE, 
    #            angle_col = '45',
    #            cutree_rows = 2,
    #            cutree_cols = 2,
    #            legend = TRUE, 
    #            show_rownames = TRUE, 
    #            show_colnames = FALSE, 
    #            annotation_col = anno_col(),
    #            annotation_legend = TRUE,
    #            annotation_names_col = TRUE)
    # })
    # 
    
    # Genes of interest
    # 
    # output$violin_gene_of_interest <- renderPlot({
    #   req(input$expre_submit)
    #   
    #   
    # })
    
    # prepare genelist
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
            
            # BiocManager::install("org.Hs.eg.db")
            # library(org.Hs.eg.db)
            
            gseGO(geneList = geneList(),
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "ALL",
                  minGSSize = 5,
                  maxGSSize = 500,
                  pvalueCutoff = input$p_value_gsea,
                  verbose = FALSE)
        } else {
            
            # BiocManager::install("org.Mm.eg.db")
            # library(org.Mm.eg.db)
            
            gseGO(geneList = geneList(),
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "ALL",
                  minGSSize = 5,
                  maxGSSize = 500,
                  pvalueCutoff = input$p_value_gsea,
                  verbose = FALSE)
        }
    })
    
    # GSEA Overview
    output$gsea_table <- renderDataTable({
      req(input$gsea_submit)
      datatable(go_gsea()@result,
                caption = paste(input$contrast_1,"vs.",input$contrast_2))
    })
    
    # visualization
    output$dot_plot <- renderPlot({
      req(input$gsea_submit)
      dotplot(go_gsea(), showCategory = 20)
    })
    
    output$network <- renderPlot({
      req(input$gsea_submit)
      p1 <- cnetplot(go_gsea(), foldChange = geneList())
      p2 <- cnetplot(go_gsea(), categorySize = "pvalue", foldChange = geneList())
      p3 <- cnetplot(go_gsea(), foldChange = geneList(), circular = TRUE, colorEdge = TRUE) 
      cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
    })
    
    # download cell profile matrix
    profile_matrix <- reactive({
        req(input$dec_submit)
        download_profile_matrix(species = input$species_2,
                                age_group = input$age_group,
                                matrixname = input$matrix_name)
    })
    
    # prepare deconvolution
    neg_probes <- reactive({
      negativeProbefData <- subset(fData(target_demoData_seg()), CodeClass == "Negative")
      unique(negativeProbefData$TargetName)
    })
    
    target_demoData_neg <- reactive({
      target_demoData_neg <- target_demoData_seg()[f_data()$DetectionRate >= input$gene_thr/100 |
                                                   f_data()$TargetName %in% neg_probes(), ]
      pData(target_demoData_neg)$biology = paste0(pData(target_demoData_neg)[,input$factors_of_interest_biology[1]],
                                                  "_", pData(target_demoData_neg)[,input$factors_of_interest_biology[2]])
      target_demoData_neg
      
    })
    
    norm_count <- reactive({
        spe <- as.SpatialExperiment(target_demoData_neg(), normData = "exprs", forceRaw = TRUE)
        assayNames(spe) <- "counts"
        raw_count <- assay(spe, "counts")
        assay(spe, "logcounts") <- edgeR::cpm(raw_count, log = TRUE)
        # col_data <- colData(spe)
        # col_data$biology <- paste0(col_data[,input$factors_of_interest_biology[1]], "_", col_data[,input$factors_of_interest_biology[2]])
        # colData(spe) <- col_data
        spe_norm <- geomxNorm(spe, method = input$nor_method)
        spe_norm <- findNCGs(spe_norm, batch_name = "slide name", top_n = 300)
        spe_ruv <- geomxBatchCorrection(spe_norm, factors = "biology", NCGs = metadata(spe_norm)$NCGs, k = input$k_best)
        assay(spe_ruv, "logcounts")
    })
    
    
    target_demoData_neg_norm <- reactive({
      target_demoData_neg_norm <- target_demoData_neg()
      assayDataElement(target_demoData_neg_norm , elt = "exprs_norm" ) = norm_count()
      target_demoData_neg_norm
      
    })
      
    bg <- reactive({
      derive_GeoMx_background(norm = target_demoData_neg_norm()@assayData$exprs_norm,
                              probepool = fData(target_demoData_neg_norm())$Module,
                              negnames = neg_probes())
    })
      
    # Performing deconvolution
    res <- reactive({
      target_demoData_neg_filter = target_demoData_neg_norm()[,colnames(target_demoData_neg_norm())[pData(target_demoData_neg_norm())$biology %in%  c(input$contrast_1, input$contrast_2)]]
      runspatialdecon(object = target_demoData_neg_filter,
                      norm_elt = "exprs_norm",
                      raw_elt = "exprs",
                      X = profile_matrix(),
                      align_genes = TRUE)
      
    })
     
    # visualization
    # heatmap
    output$dec_heatmap <- renderPlot({
      req(input$dec_submit)
      heatmap(t(res()$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))
      
    })
    
    # Get file path
    shinyDirChoose(input, 'folder', roots = volumes())
    
    # creat seurat object
    tissue <- reactive({
      Load10X_Spatial(data.dir = parseDirPath(roots = volumes(), input$folder))
    })
    
    # study overview
    output$summary_table_seurat <- renderDataTable({
      req(input$upload_button_seurat)
      datatable(data.frame(Features = nrow(tissue()),
                           Samples = ncol(tissue())))
      
    })
    
    tissue_filter <- reactive({
      subset(tissue(), subset = nFeature_Spatial > input$nFeature[1] & nFeature_Spatial < input$nFeature[2] & 
               nCount_Spatial > input$nCount[1] & nCount_Spatial < input$nCount[2])
    })
    
    output$count_filter <- renderPlot({
      req(input$upload_button_seurat)
      plot1 <- VlnPlot(tissue_filter(), features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
      plot2 <- SpatialFeaturePlot(tissue_filter(), features = "nCount_Spatial") + theme(legend.position = "right")
      wrap_plots(plot1, plot2)
    })
    
    output$gene_filter <- renderPlot({
      req(input$upload_button_seurat)
      plot3 <- VlnPlot(tissue_filter(), features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
      plot4 <- SpatialFeaturePlot(tissue_filter(), features = "nFeature_Spatial") + theme(legend.position = "right")
      wrap_plots(plot3, plot4)
    })
    
    output$summary_table_filter_seurat <- renderDataTable({
      req(input$filter_submit_seurat)
      datatable(data.frame(Features = nrow(tissue_filter()),
                           Samples = ncol(tissue_filter())))
    })
    
    # normalization
    tissue_nor <- reactive({
      SCTransform(tissue_filter(), assay = "Spatial", verbose = FALSE)
    })
    
    # clustering analysis
    tissue_umap <- reactive({
      tissue_dim <- RunPCA(tissue_nor(), assay = "SCT", verbose = FALSE)
      tissue_dim <- FindNeighbors(tissue_dim, reduction = "pca", dims = 1:30)
      tissue_dim <- FindClusters(tissue_dim, verbose = FALSE)
      RunUMAP(tissue_dim, reduction = "pca", dims = 1:30)
    })
    
    output$umap_seurat <- renderPlot({
      req(input$filter_submit_seurat)
      p1 <- DimPlot(tissue_umap(), reduction = "umap", label = TRUE)
      p2 <- SpatialDimPlot(tissue_umap(), label = TRUE, label.size = 3)
      p1 + p2
    })
    
    tissue_markers <- reactive({
      FindAllMarkers(tissue_umap(), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
    })
    
    topgene <- reactive({
      tissue_markers() %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    })
    
    output$heatmap_seurat <- renderPlot({
      DoHeatmap(tissue_umap(), features = topgene()$gene, size = 3)
    })
    
    de_markers <- reactive({
      FindMarkers(tissue_umap(), ident.1 = input$cluster_1, ident.2 = input$cluster_2,
                  logfc.threshold = input$logfc_seurat, return.thresh = input$p_value_seurat)
    })
    
    output$summary_de_table_seurat <- renderDataTable({
      datatable(data.frame(Significant_genes = nrow(de_markers()),
                           Up_regulated_genes = length(which(de_markers()$avg_log2FC > 0)),
                           Down_regulated_genes = length(which(de_markers()$avg_log2FC < 0)),
                           non_significant_genes = nrow(tissue_umap()) - nrow(de_markers())),
                colnames = c("Significant genes", "Up-regulated genes", "Down-regulated genes", "Non-significant genes"))
    })
    
    output$de_table_seurat <- renderDataTable({
      datatable(de_markers(), rownames = TRUE)
    })
    
    output$gene_of_interest_seurat <- renderPlot({
      SpatialFeaturePlot(tissue_umap(), features = input$gene_name_seurat, alpha = c(0.1, 1)) 
    })
    
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
        
        # BiocManager::install("org.Hs.eg.db")
        # library(org.Hs.eg.db)
        
        gseGO(geneList = geneList_seurat(),
              OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "ALL",
              minGSSize = 5,
              maxGSSize = 500,
              pvalueCutoff = input$p_value_gsea_seurat,
              verbose = FALSE)
      } else {
        
        # BiocManager::install("org.Mm.eg.db")
        # library(org.Mm.eg.db)
        
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
    output$gsea_table_seurat <- renderDataTable({
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
