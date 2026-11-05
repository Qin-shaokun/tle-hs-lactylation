# ============================================================
# 08_enrichment.R — GO / KEGG / GSEA 富集（2026-11-04 对接 02_deg_main）
# ORA：显著 DEG；GSEA：全基因 t 值排序（gseGO 离线 + gseKEGG 尽力而为）
# 输出：results/08_enrich_GO.csv, 08_enrich_KEGG.csv, 08_gsea_GO.csv,
#       08_gsea_focus.csv；figures/08_enrich_GO.png, 08_enrich_KEGG.png,
#       08_gsea_focus.png
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
})

d   <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
deg <- d$deg   # data.frame: gene, logFC, t, adj.P.Val, ...
sig <- subset(deg, abs(logFC) >= THRESH$lfc & adj.P.Val <= THRESH$fdr)
cat("显著 DEG：", nrow(sig), "个\n")

bm_sig <- bitr(sig$gene, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("ENTREZ 映射：", nrow(bm_sig), "/", nrow(sig), "\n")

## ---- ORA：GO BP + KEGG ----
ego <- enrichGO(bm_sig$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", readable = TRUE, universe = NULL)
write.csv(as.data.frame(ego), file.path(RES_DIR, "08_enrich_GO.csv"), row.names = FALSE)
cat("GO BP 显著通路：", sum(ego@result$p.adjust < 0.05), "\n")

ekegg <- tryCatch(enrichKEGG(bm_sig$ENTREZID, organism = "hsa",
                             pAdjustMethod = "BH"),
                  error = function(e) { message("KEGG ORA 失败：", conditionMessage(e)); NULL })
if (!is.null(ekegg)) {
  ekegg <- setReadable(ekegg, org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg), file.path(RES_DIR, "08_enrich_KEGG.csv"), row.names = FALSE)
  cat("KEGG 显著通路：", sum(ekegg@result$p.adjust < 0.05), "\n")
}

p_go <- dotplot(ego, showCategory = 15) + ggtitle("GO biological process enrichment of DEGs") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIG_DIR, "08_enrich_GO.png"), p_go,
       width = 8.5, height = 7, dpi = PLOT_THEME$save_dpi)
if (!is.null(ekegg)) {
  p_k <- dotplot(ekegg, showCategory = 15) + ggtitle("KEGG pathway enrichment of DEGs") +
    theme_minimal(base_size = 11)
  ggsave(file.path(FIG_DIR, "08_enrich_KEGG.png"), p_k,
         width = 8.5, height = 7, dpi = PLOT_THEME$save_dpi)
}

## ---- GSEA：全基因 t 排序（修复原版错位 bug）----
rk <- deg[order(deg$t, decreasing = TRUE), ]
bm_all <- bitr(rk$gene, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
bm_all <- bm_all[!duplicated(bm_all$SYMBOL), ]
idx <- match(bm_all$SYMBOL, rk$gene)
stat_vec <- setNames(rk$t[idx], bm_all$ENTREZID)
stat_vec <- sort(stat_vec, decreasing = TRUE)
cat("GSEA 排序基因：", length(stat_vec), "\n")

gse_go <- gseGO(stat_vec, OrgDb = org.Hs.eg.db, ont = "BP",
                minGSSize = 15, maxGSSize = 800,
                pAdjustMethod = "BH", eps = 0, verbose = FALSE)
res_go <- as.data.frame(gse_go)
if (nrow(res_go) > 0) res_go <- res_go[order(res_go[["p.adjust"]]), ]
write.csv(res_go, file.path(RES_DIR, "08_gsea_GO.csv"), row.names = FALSE)

## 重点通路快查（乳酸化/糖酵解/突触/神经炎症/乳酸相关）
kw <- "lactat|lactyl|glycolys|synap|long.term|neuroinflam|hippocamp|axon|oxidative phosphorylation|HIF"
focus_idx <- unique(c(grep(kw, res_go$Description, ignore.case = TRUE),
                      grep(kw, res_go$ID, ignore.case = TRUE)))
focus_tbl <- if (length(focus_idx)) res_go[focus_idx, c("ID", "Description", "NES", "p.adjust", "pvalue")] else res_go[NULL, ]
focus_tbl <- focus_tbl[order(-abs(focus_tbl$NES)), ]
write.csv(focus_tbl, file.path(RES_DIR, "08_gsea_focus.csv"), row.names = FALSE)
cat("\n== GSEA 焦点通路（按 |NES| 排序，前 20）==\n")
print(head(focus_tbl, 20))

## 焦点通路 NES 条形图
if (nrow(focus_tbl) > 0) {
  top <- head(focus_tbl, 20)
  top$dir <- ifelse(top$NES > 0, "Up in TLE-HS", "Down in TLE-HS")
  p3 <- ggplot(top, aes(reorder(Description, NES), NES, fill = dir)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("Up in TLE-HS" = "#C94F4F", "Down in TLE-HS" = "#7FA6C9")) +
    labs(x = "", y = "Normalized enrichment score (NES)",
         title = "GSEA of focused pathways") +
    theme_minimal(base_size = 10)
  ggsave(file.path(FIG_DIR, "08_gsea_focus.png"), p3,
         width = 9, height = 7, dpi = PLOT_THEME$save_dpi)
}

## KEGG GSEA（需联网，尽力而为）
gse_k <- tryCatch(gseKEGG(stat_vec, organism = "hsa", pAdjustMethod = "BH",
                          verbose = FALSE),
                  error = function(e) { message("KEGG GSEA 失败：", conditionMessage(e)); NULL })
if (!is.null(gse_k)) {
  write.csv(as.data.frame(gse_k), file.path(RES_DIR, "08_gsea_KEGG.csv"), row.names = FALSE)
}
message("08 完成")
