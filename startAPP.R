startSpaXSeq <- function() {
  if (interactive()) {
    # Load packages
    require(shiny)
    require(shinyFiles)
    require(shinythemes)
    require(shinyWidgets)
    require(GeomxTools)
    require(ggalluvial)
    require(DT)
    require(scales)
    require(SpatialExperiment)
    require(standR)
    require(edgeR)
    require(limma)
    require(ggrepel)
    require(clusterProfiler)
    require(enrichplot)
    require(SpatialDecon)
    require(tidyr)
    require(tibble)
    require(org.Mm.eg.db)
    require(org.Hs.eg.db)
    require(Seurat)
    require(patchwork)
    require(dplyr)
    require(pheatmap)
    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)
  }
}