source("config.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(Seurat); library(dplyr)})
obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
Idents(obj) <- "celltype"
# 仅在高变基因上找 marker：足够做注释质检，速度约为全基因的 1/18
feats <- VariableFeatures(obj)
cat("特征数:", length(feats), "| 细胞类型:", paste(levels(Idents(obj)), collapse = ", "), "\n")
mk <- FindAllMarkers(obj, features = feats, only.pos = TRUE,
                     logfc.threshold = 0.25, min.pct = 0.25, verbose = FALSE)
write.csv(mk, file.path(RES_DIR, "09c_celltype_markers.csv"), row.names = FALSE)
cat("marker 记录:", nrow(mk), "\n")
print(mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% ungroup() %>% select(cluster, gene, avg_log2FC, pct.1, pct.2))
cat("MARKERS_DONE\n")
