# 08b_kegg.R — KEGG ORA + GSEA 专用（带重试；08 主脚本里 KEGG 超时时的兜底）
# 输入：results/02_deg_main.rds；输出：08_enrich_KEGG.csv / 08_gsea_KEGG.csv
#       figures/08_enrich_KEGG.png|pdf, 08_gsea_KEGG.png|pdf
source("config.R", encoding = "UTF-8")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
})
options(timeout = 600)

d   <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
deg <- d$deg
sig <- subset(deg, abs(logFC) >= THRESH$lfc & adj.P.Val <= THRESH$fdr)
cat("显著 DEG：", nrow(sig), "个\n")

bm_sig <- bitr(sig$gene, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("ENTREZ 映射：", nrow(bm_sig), "/", nrow(sig), "\n")

## ---- 重试封装 ----
retry <- function(expr, n_try = 4, wait = 30, tag = "") {
  for (i in seq_len(n_try)) {
    r <- tryCatch(expr, error = function(e) { e })
    if (!inherits(r, "error")) return(r)
    message(sprintf("[%s] 第 %d/%d 次失败：%s", tag, i, n_try,
                    conditionMessage(r)))
    if (i < n_try) Sys.sleep(wait)
  }
  NULL
}

## ---- KEGG ORA ----
ekegg <- retry(enrichKEGG(bm_sig$ENTREZID, organism = "hsa",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.1),
               tag = "ORA")
if (!is.null(ekegg)) {
  ekegg <- setReadable(ekegg, org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg), file.path(RES_DIR, "08_enrich_KEGG.csv"),
            row.names = FALSE)
  cat("KEGG ORA 显著通路：", sum(ekegg@result$p.adjust < 0.05), "\n")

  p_k <- dotplot(ekegg, showCategory = 15) +
    ggtitle("KEGG pathway enrichment of DEGs") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(FIG_DIR, "08_enrich_KEGG.png"), p_k,
         width = 7.2, height = 6.2, dpi = PLOT_THEME$save_dpi)
  ggsave(file.path(FIG_DIR, "08_enrich_KEGG.pdf"), p_k,
         width = 7.2, height = 6.2)
} else {
  message("KEGG ORA 最终失败，跳过绘图")
}

## ---- KEGG GSEA（全基因 t 值排序）----
rk <- deg[order(deg$t, decreasing = TRUE), ]
bm_all <- bitr(rk$gene, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
bm_all <- merge(bm_all, rk[, c("gene", "t")], by.x = "SYMBOL", by.y = "gene")
bm_all <- bm_all[!duplicated(bm_all$ENTREZID), ]
stat_vec <- bm_all$t; names(stat_vec) <- bm_all$ENTREZID
stat_vec <- sort(stat_vec, decreasing = TRUE)

gse_k <- retry(gseKEGG(stat_vec, organism = "hsa", pAdjustMethod = "BH",
                       verbose = FALSE, seed = TRUE, eps = 0),
               tag = "GSEA")
if (!is.null(gse_k)) {
  write.csv(as.data.frame(gse_k), file.path(RES_DIR, "08_gsea_KEGG.csv"),
            row.names = FALSE)
  cat("KEGG GSEA 显著通路：", sum(gse_k@result$p.adjust < 0.05), "\n")

  ## 焦点通路条形图（代谢/突触/神经炎症相关）
  gk <- as.data.frame(gse_k)
  kw <- c("Oxidative phosphorylation", "Ribosome", "Proteasome",
          "Insulin secretion", "Glycolysis", "Citrate", "Pyruvate",
          "Carbon metabolism", "Synapse", "Long-term potentiation",
          "Long-term depression", "Axon guidance", "Neurotrophin",
          "Cytokine", "NF-kappa", "Complement", "Phagosome", "Lysosome",
          "HIF-1", "AMPK", "mTOR")
  hit <- gk[grepl(paste(kw, collapse = "|"), gk$Description, ignore.case = TRUE) &
              gk$p.adjust < 0.25, ]
  if (nrow(hit) > 0) {
    top <- head(hit[order(-abs(hit$NES)), ], 12)
    top$dir <- ifelse(top$NES > 0, "Up in TLE-HS", "Down in TLE-HS")
    p3 <- ggplot(top, aes(reorder(Description, NES), NES, fill = dir)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("Up in TLE-HS" = "#C94F4F",
                                   "Down in TLE-HS" = "#7FA6C9"), name = "") +
      labs(x = "", y = "Normalized enrichment score (NES)",
           title = "KEGG GSEA: metabolism, translation, synapse and inflammation") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
    ggsave(file.path(FIG_DIR, "08_gsea_KEGG.png"), p3,
           width = 8, height = 5.6, dpi = PLOT_THEME$save_dpi)
    ggsave(file.path(FIG_DIR, "08_gsea_KEGG.pdf"), p3, width = 8, height = 5.6)
    cat("焦点通路：", nrow(top), "条\n")
  }
} else {
  message("KEGG GSEA 最终失败，跳过")
}

message("08b 完成")
