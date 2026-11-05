source("config.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
Idents(obj) <- "celltype"
genes <- c("TNFRSF13C","TMEM156","HK2","PFKFB3","NLRP3","SLC16A7","EP300","BDNF")
p <- DotPlot(obj, features = genes, group.by = "celltype", split.by = "condition", dot.scale = 7) +
  RotatedAxis() + theme(axis.text.x = element_text(size = 10)) +
  ggtitle("Hub and key lactylation-related genes")
ggsave(file.path(FIG_DIR, "09c_keygene_dotplot.png"), p, width = 8.5, height = 5.5, dpi = 200)
cat("ok\n")
