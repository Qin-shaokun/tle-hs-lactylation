# ============================================================
# 09c_localize.R — hub 基因与乳酸化基因的细胞类型定位
# 输入：results/09b_sc_annotated.rds
# 输出：
#   results/09c_hub_localization.csv          hub 基因 × 细胞类型 × 条件
#   results/09c_lactylation_localization.csv  49 乳酸化基因同上
#   results/09c_celltype_markers.csv          FindAllMarkers（每类型 vs 其余）
#   figures/09c_hub_featureplot.png  09c_hub_dotplot.png  09c_lactylation_dotplot.png
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})
options(width = 120)

obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "celltype"
expr <- GetAssayData(obj, layer = "data")
ct   <- obj$celltype
cond <- obj$condition
cat("载入:", ncol(obj), "核 | 类型:", paste(unique(as.character(ct)), collapse = ", "), "\n")

hub <- read.csv(file.path(RES_DIR, "04_hub_genes.csv"), stringsAsFactors = FALSE)$hub_gene
gs  <- read.csv(GENESET_CSV, stringsAsFactors = FALSE)
lac <- unique(gs$gene_symbol)
cat("hub 基因:", paste(hub, collapse = ", "), "| 乳酸化基因:", length(lac), "\n")

## ---- 1. 通用定位统计函数 ----
localize <- function(genes) {
  genes <- intersect(genes, rownames(expr))
  miss  <- setdiff(unique(c(hub, lac)), genes)
  if (length(miss)) cat("注意：表达矩阵中缺失", paste(miss, collapse = ", "), "\n")
  out <- list()
  for (g in genes) {
    e   <- expr[g, ]
    pos <- e > 0
    for (k in unique(as.character(ct))) {
      kpos <- ct == k
      if (sum(kpos) < 50) next
      for (cc in c("TLE", "Control")) {
        idx <- kpos & cond == cc
        n <- sum(idx); p <- sum(pos[idx]); me <- if (n > 0) mean(e[idx]) else NA
        out[[length(out) + 1]] <- data.frame(
          gene = g, celltype = k, condition = cc,
          n_cells = n, n_expr = p, pct_expr = p / n,
          mean_expr_all = me, stringsAsFactors = FALSE
        )
      }
      # 组内 Wilcoxon（两类均有 >=50 细胞才做）
      i1 <- kpos & cond == "TLE"; i0 <- kpos & cond == "Control"
      if (sum(i1) >= 50 && sum(i0) >= 50 && sum(pos[i1]) > 5 && sum(pos[i0]) > 5) {
        w <- suppressWarnings(wilcox.test(e[i1], e[i0], exact = FALSE))
        pv <- w$p.value
        for (cc in c("TLE", "Control")) {
          j <- length(out) - if (cc == "Control") 0 else 1
          out[[j]]$wilcox_p_vs_other <- pv
        }
      }
    }
    # 跨类型分布（不分组）
    for (k in unique(as.character(ct))) {
      kpos <- ct == k
      out[[length(out) + 1]] <- data.frame(
        gene = g, celltype = k, condition = "All",
        n_cells = sum(kpos), n_expr = sum(pos[kpos]),
        pct_expr = mean(pos[kpos]),
        mean_expr_all = mean(e[kpos]), stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

cat("\n== hub 基因定位 ==\n")
hub_loc <- localize(hub)
write.csv(hub_loc, file.path(RES_DIR, "09c_hub_localization.csv"), row.names = FALSE)
print(hub_loc %>% filter(condition == "All") %>%
        select(gene, celltype, pct_expr, mean_expr_all) %>%
        arrange(gene, desc(pct_expr)), digits = 2)

cat("\n== 乳酸化基因定位（49 个）==\n")
lac_loc <- localize(lac)
write.csv(lac_loc, file.path(RES_DIR, "09c_lactylation_localization.csv"), row.names = FALSE)
cat("乳酸化定位表：", nrow(lac_loc), "行\n")

## ---- 2. 图 ----
# hub FeaturePlot（分条件）
ggsave(file.path(FIG_DIR, "09c_hub_featureplot.png"),
       FeaturePlot(obj, features = hub, split.by = "condition", ncol = 2, order = TRUE) &
         theme(title = element_text(size = 10)),
       width = 11, height = 8.5, dpi = 200)

ggsave(file.path(FIG_DIR, "09c_hub_dotplot.png"),
       DotPlot(obj, features = hub, group.by = "celltype", split.by = "condition") +
         RotatedAxis() + ggtitle("Hub genes: TNFRSF13C / TMEM156"),
       width = 8, height = 5.5, dpi = 200)

ggsave(file.path(FIG_DIR, "09c_lactylation_dotplot.png"),
       DotPlot(obj, features = intersect(lac, rownames(obj)), group.by = "celltype",
               dot.scale = 5) + RotatedAxis() + coord_flip() +
         ggtitle("Lactylation-related genes by cell type"),
       width = 9.5, height = 14, dpi = 200)

## ---- 3. FindAllMarkers（探索性，耗时最长放最后）----
cat("\n== FindAllMarkers ==\n")
Idents(obj) <- "celltype"
mk <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.25,
                     min.pct = 0.25, verbose = FALSE)
mk$cluster_name <- as.character(mk$cluster)
write.csv(mk, file.path(RES_DIR, "09c_celltype_markers.csv"), row.names = FALSE)
cat("marker 记录:", nrow(mk), "\n")
print(mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3) %>% ungroup())

cat("\n09C_DONE\n")
