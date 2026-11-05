# ============================================================
# 09b_cluster.R — GSE190452 snRNA 聚类与皮层细胞类型注释
# 输入：results/09_sc_filtered.rds（QC 后 64,667 核，含 condition/batch）
# 输出：
#   results/09b_sc_annotated.rds        注释后对象（含 celltype）
#   results/09b_cluster_scores.csv      每簇 × 每类 marker 平均分
#   results/09b_celltype_proportions.csv 细胞类型构成（按样本/条件）
#   figures/09b_elbow.png  09b_umap_cluster.png  09b_umap_celltype.png
#   figures/09b_umap_condition.png  09b_umap_batch.png  09b_umap_sample.png
#   figures/09b_marker_dotplot.png
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(dplyr); library(patchwork)
})
options(width = 120)

obj <- readRDS(file.path(RES_DIR, "09_sc_filtered.rds"))
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
cat("载入:", ncol(obj), "核 x", nrow(obj), "基因\n")

## ---- 1. 标准预处理 ----
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

png(file.path(FIG_DIR, "09b_elbow.png"), width = 1400, height = 900, res = 140)
print(ElbowPlot(obj, ndims = 50))
dev.off()

PCS <- 30
obj <- RunUMAP(obj, dims = 1:PCS, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:PCS, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
ncl <- length(levels(obj))
cat("簇数(0.4):", ncl, "\n")

## ---- 2. marker 打分注释 ----
marker_sets <- list(
  ExNeuron        = c("SLC17A7","SATB2","CAMK2A","NRGN","RORB","CUX2","TBR1"),
  InNeuron        = c("GAD1","GAD2","PVALB","SST","VIP","LHX6"),
  Astrocyte       = c("AQP4","SLC1A3","SLC1A2","ALDH1L1"),
  Oligodendrocyte = c("MBP","PLP1","MOG","MAL","MAG"),
  OPC             = c("PDGFRA","OLIG1","OLIG2","VCAN","CSPG4"),
  Microglia       = c("C1QA","C1QB","C1QC","CX3CR1","P2RY12","AIF1"),
  Endothelial     = c("CLDN5","PECAM1","FLT1","VWF"),
  Pericyte_SMC    = c("PDGFRB","RGS5","NOTCH3","ACTA2")
)
expr <- GetAssayData(obj, layer = "data")
cl_vec <- as.character(Idents(obj))

# 每个细胞的每类平均表达 -> 按簇汇总
score_cell <- sapply(marker_sets, function(gs) {
  gs2 <- intersect(gs, rownames(expr))
  Matrix::colMeans(expr[gs2, , drop = FALSE])
})
score_cell <- as.data.frame(score_cell)
score_cell$cluster <- cl_vec
cluster_scores <- score_cell %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean), n_cells = n(), .groups = "drop")
set_names <- names(marker_sets)
assign_map <- setNames(
  set_names[apply(as.matrix(cluster_scores[, set_names]), 1, which.max)],
  cluster_scores$cluster
)
# 手动覆盖（首轮结果出来后按需修改）
OVERRIDE <- c()
if (length(OVERRIDE) > 0) for (cl in names(OVERRIDE)) assign_map[[cl]] <- OVERRIDE[[cl]]

cat("\n== 簇 × marker 分数 ==\n"); print(as.data.frame(cluster_scores), digits = 2)
cat("\n簇 -> 类型：\n"); print(assign_map)

obj$celltype <- unname(assign_map[cl_vec])
write.csv(as.data.frame(cluster_scores), file.path(RES_DIR, "09b_cluster_scores.csv"), row.names = FALSE)

## ---- 3. 构成表 ----
prop_tbl <- obj@meta.data %>%
  count(sample, condition, batch, celltype, name = "n") %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.csv(prop_tbl, file.path(RES_DIR, "09b_celltype_proportions.csv"), row.names = FALSE)
cat("\n== 条件平均构成 ==\n")
print(prop_tbl %>% group_by(condition, celltype) %>%
        summarise(mean_frac = mean(frac), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = condition, values_from = mean_frac))

## ---- 4. 图 ----
pal <- scales::hue_pal()(ncl)
p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("Cell types")
p2 <- DimPlot(obj, group.by = "condition") + ggtitle("Condition")
p3 <- DimPlot(obj, group.by = "batch")  + ggtitle("Batch")
p4 <- DimPlot(obj, group.by = "sample") + ggtitle("Sample")
ggsave(file.path(FIG_DIR, "09b_umap_celltype.png"),  p1, width = 7.5, height = 6, dpi = 200)
ggsave(file.path(FIG_DIR, "09b_umap_condition.png"), p2, width = 7.5, height = 6, dpi = 200)
ggsave(file.path(FIG_DIR, "09b_umap_batch.png"),     p3, width = 7.5, height = 6, dpi = 200)
ggsave(file.path(FIG_DIR, "09b_umap_sample.png"),    p4, width = 9.5, height = 6, dpi = 200)

png(file.path(FIG_DIR, "09b_umap_cluster.png"), width = 1600, height = 1100, res = 160)
print(DimPlot(obj, group.by = "seurat_clusters", label = TRUE) + ggtitle("Seurat clusters (res 0.4)"))
dev.off()

goi <- unique(unlist(marker_sets)); goi <- intersect(goi, rownames(obj))
p <- DotPlot(obj, features = goi, group.by = "celltype", dot.scale = 6) +
  RotatedAxis() + coord_flip() + ggtitle("Canonical cortical markers")
ggsave(file.path(FIG_DIR, "09b_marker_dotplot.png"), p, width = 9, height = 8.5, dpi = 200)

saveRDS(obj, file.path(RES_DIR, "09b_sc_annotated.rds"))
cat("\n09b 完成：已保存 09b_sc_annotated.rds\n")
