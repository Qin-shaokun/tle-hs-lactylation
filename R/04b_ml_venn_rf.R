# 04b_ml_venn_rf.R — Fig 3 缺失面板：三算法韦恩图 + RF 重要性图
# 只用已保存的 04_ml.rds，不重算
source("config.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(ggplot2))

ml <- readRDS(file.path(RES_DIR, "04_ml.rds"))
gs <- read.csv(GENESET_CSV)

lasso_genes <- ml$lasso
svm_genes   <- head(ml$svm_rank, 10)
rf_genes    <- ml$rf_genes
hub         <- ml$hub

n <- function(x) length(x)
A <- n(lasso_genes); B <- n(svm_genes); C <- n(rf_genes)
AB <- length(intersect(lasso_genes, svm_genes))
AC <- length(intersect(lasso_genes, rf_genes))
BC <- length(intersect(svm_genes, rf_genes))
ABC <- length(Reduce(intersect, list(lasso_genes, svm_genes, rf_genes)))
cat("LASSO:", A, " SVM-top10:", B, " RF-top15:", C,
    " | 两两交集:", AB, AC, BC, " 三交集:", ABC, "\n")

## ============ (A) 三算法韦恩图（base R 手绘）============
draw_venn <- function() {
  ## 上：韦恩图本体
  par(fig = c(0, 1, 0.17, 1), mar = c(0, 0, 2.8, 0), xpd = FALSE)
  plot(NA, xlim = c(-1.62, 1.62), ylim = c(-1.55, 1.35), asp = 1,
       axes = FALSE, xlab = "", ylab = "")
  title(main = "Consensus of three feature-selection algorithms", cex.main = 1.05)
  col_a <- rgb(0.23, 0.51, 0.96, 0.45)   # LASSO 蓝
  col_b <- rgb(0.93, 0.47, 0.20, 0.45)   # SVM-RFE 橙
  col_c <- rgb(0.20, 0.72, 0.38, 0.45)   # RF 绿
  symbols(c(-0.58, 0.58, 0), c(0.28, 0.28, -0.52), circles = rep(0.90, 3),
          add = TRUE, inches = FALSE,
          bg = c(col_a, col_b, col_c), fg = c("steelblue", "darkorange", "forestgreen"))
  lab <- function(x, y, s, col = "black", font = 1, cex = 1.05)
    text(x, y, s, col = col, font = font, cex = cex)
  lab(-0.92, 0.82, A, "steelblue", 2, 1.15)
  lab( 0.92, 0.82, B, "darkorange", 2, 1.15)
  lab( 0.00, -1.12, C, "forestgreen", 2, 1.15)
  lab( 0.00, 0.58, AB)          # LASSO∩SVM only
  lab(-0.40, -0.02, AC - ABC)   # LASSO∩RF only
  lab( 0.40, -0.02, BC - ABC)   # SVM∩RF only
  lab( 0.00, 0.06, ABC, "firebrick", 2, 1.25)

  ## 下：图例 + hub 说明（独立区域，避免裁切/压字）
  par(fig = c(0, 1, 0, 0.19), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("top", bty = "n", cex = 0.82, ncol = 3, xpd = NA,
         legend = c(sprintf("LASSO (n = %d)", A),
                    sprintf("SVM-RFE top10 (n = %d)", B),
                    sprintf("Random forest top15 (n = %d)", C)),
         text.col = c("steelblue", "darkorange", "forestgreen"))
  text(0.5, 0.22, sprintf("Two-algorithm consensus (LASSO ∩ RF): %s",
                         paste(hub, collapse = ", ")),
       col = "firebrick", font = 2, cex = 0.85, xpd = NA)
}

png(file.path(FIG_DIR, "04_venn.png"), 2100, 1800, res = PLOT_THEME$save_dpi)
draw_venn(); dev.off()
pdf(file.path(FIG_DIR, "04_venn.pdf"), width = 7, height = 6)
draw_venn(); dev.off()

## ============ (B) RF 重要性 Top15 ============
sc <- read.csv(file.path(RES_DIR, "04_ml_scores.csv"))
rf_top <- sc[sc$gene %in% rf_genes, ]
rf_top <- rf_top[order(rf_top$rf_importance), ]
rf_top$is_hub   <- rf_top$gene %in% hub
rf_top$in_lact  <- rf_top$gene %in% gs$gene_symbol

p <- ggplot(rf_top, aes(rf_importance, reorder(gene, rf_importance))) +
  geom_col(aes(fill = ifelse(is_hub, "Hub gene (LASSO ∩ RF)",
                             ifelse(in_lact, "Lactylation gene set", "Other"))),
           width = 0.68) +
  scale_fill_manual(values = c("Hub gene (LASSO ∩ RF)" = "#C0392B",
                               "Lactylation gene set" = "#8E44AD",
                               "Other" = "grey62"), name = "") +
  geom_text(aes(label = sprintf("%.2f", rf_importance)), hjust = -0.15,
            size = 2.9, color = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Random forest importance (decrease in Gini)",
       y = NULL,
       title = "Top 15 genes by random forest importance") +
  theme_bw(PLOT_THEME$base_size) +
  theme(legend.position = "top",
        axis.text.y = element_text(face = ifelse(rf_top$gene[order(rf_top$rf_importance)] %in% hub,
                                                 "bold", "plain")))
ggsave(file.path(FIG_DIR, "04_rf_importance.png"), p,
       width = 6.2, height = 5.2, dpi = PLOT_THEME$save_dpi)
ggsave(file.path(FIG_DIR, "04_rf_importance.pdf"), p,
       width = 6.2, height = 5.2)

## 输出交集明细，供图注核对
cat("\nLASSO∩RF（两算法共识，即 hub 来源）：",
    paste(intersect(lasso_genes, rf_genes), collapse = ", "), "\n")
cat("SVM-RFE top10：", paste(svm_genes, collapse = ", "), "\n")
cat("RF top15 中属于乳酸化基因集的：",
    paste(intersect(rf_genes, gs$gene_symbol), collapse = ", "), "\n")

message("04b 完成")
