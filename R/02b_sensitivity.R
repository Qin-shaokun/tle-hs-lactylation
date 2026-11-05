# 02b_sensitivity.R — S1 敏感性分析图（复用 02 已保存结果，不重算）
# 面板：A DEG 计数 / B 主模型 vs 校正脑区 logFC / C 主模型 vs 成人子集 / D 三模型 DEG 交集
source("config.R", encoding = "UTF-8")

mn <- read.csv(file.path(RES_DIR, "02_deg_all.csv"))
rg <- read.csv(file.path(RES_DIR, "02_deg_sens_region.csv"))
ad <- read.csv(file.path(RES_DIR, "02_deg_sens_adult.csv"))

sig_of <- function(d) d$gene[abs(d$logFC) >= THRESH$lfc & d$adj.P.Val <= THRESH$fdr]
sm <- sig_of(mn); sr <- sig_of(rg); sa <- sig_of(ad)
dir_of <- function(d, g) sign(d$logFC[match(g, d$gene)])

## 统计量
pair_stats <- function(a, b, ga, gb) {
  m <- merge(a[, c("gene", "logFC")], b[, c("gene", "logFC")], by = "gene")
  sh <- intersect(ga, gb)
  mg <- merge(data.frame(gene = sh, lfc_a = a$logFC[match(sh, a$gene)]),
              data.frame(gene = sh, lfc_b = b$logFC[match(sh, b$gene)]), by = "gene")
  ## 方向一致性：以主模型 DEG 为基准，看其在敏感性模型中是否同向
  base <- intersect(ga, b$gene)
  la <- a$logFC[match(base, a$gene)]; lb <- b$logFC[match(base, b$gene)]
  list(rho_all  = cor(m$logFC.x, m$logFC.y, method = "spearman"),
       rho_deg  = cor(mg$lfc_a, mg$lfc_b, method = "spearman"),
       n_shared = nrow(mg),
       concord  = 100 * mean(sign(la) == sign(lb)))
}
st_r <- pair_stats(mn, rg, sm, sr)
st_a <- pair_stats(mn, ad, sm, sa)

cat(sprintf("region: rho_all=%.3f rho_sharedDEG=%.3f shared=%d concord=%.1f%%\n",
            st_r$rho_all, st_r$rho_deg, st_r$n_shared, st_r$concord))
cat(sprintf("adult : rho_all=%.3f rho_sharedDEG=%.3f shared=%d concord=%.1f%%\n",
            st_a$rho_all, st_a$rho_deg, st_a$n_shared, st_a$concord))

## ============ 绘图 ============
draw_s1 <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.2, 1.2), cex.lab = 1.0,
      cex.axis = 0.9, mgp = c(2.4, 0.7, 0))

  ## ---- A：DEG 计数 ----
  cnt_up <- c(sum(mn$logFC[mn$gene %in% sm] > 0), sum(rg$logFC[rg$gene %in% sr] > 0),
              sum(ad$logFC[ad$gene %in% sa] > 0))
  cnt_dn <- c(sum(mn$logFC[mn$gene %in% sm] < 0), sum(rg$logFC[rg$gene %in% sr] < 0),
              sum(ad$logFC[ad$gene %in% sa] < 0))
  mm <- rbind(Up = cnt_up, Down = cnt_dn)
  bp <- barplot(mm, col = c("#C0392B", "#2471A3"), border = NA,
                ylim = c(0, max(colSums(mm)) * 1.28),
                names.arg = c("Main\n(n = 78)", "Region-adjusted\n(n = 78)",
                              "Adults only\n(n = 59)"),
                ylab = "Number of DEGs", main = "A  DEG counts across models",
                cex.main = 1.0, font.main = 2, las = 1)
  text(bp, colSums(mm) + max(colSums(mm)) * 0.04, labels = colSums(mm),
       font = 2, cex = 0.95, xpd = NA)
  legend("topright", c("Up in TLE-HS", "Down in TLE-HS"), fill = c("#C0392B", "#2471A3"),
         border = NA, bty = "n", cex = 0.82, inset = c(0.01, 0.02))

  ## ---- 散点函数 ----
  scatter <- function(b, st, xlab, ylab, main, note) {
    m <- merge(mn[, c("gene", "logFC")], b[, c("gene", "logFC")], by = "gene")
    shared <- intersect(sm, intersect(sm, b$gene[abs(b$logFC) >= THRESH$lfc &
                                                   b$adj.P.Val <= THRESH$fdr]))
    is_sh <- m$gene %in% shared
    rng <- range(c(m$logFC.x, m$logFC.y), na.rm = TRUE)
    plot(m$logFC.x[!is_sh], m$logFC.y[!is_sh], pch = 16, cex = 0.35,
         col = rgb(0.6, 0.6, 0.6, 0.25), xlim = rng, ylim = rng,
         xlab = xlab, ylab = ylab, main = main, cex.main = 1.0, font.main = 2)
    points(m$logFC.x[is_sh], m$logFC.y[is_sh], pch = 16, cex = 0.45,
           col = rgb(0.75, 0.22, 0.17, 0.55))
    abline(0, 1, lty = 2, col = "grey30", lwd = 1.2)
    abline(h = c(-THRESH$lfc, THRESH$lfc), v = c(-THRESH$lfc, THRESH$lfc),
           lty = 3, col = "grey70")
    legend("bottomright", bty = "n", cex = 0.78, pch = 16,
           col = c(rgb(0.6, 0.6, 0.6, 0.5), rgb(0.75, 0.22, 0.17, 0.8)),
           legend = c(sprintf("Other genes (n = %d)", sum(!is_sh)),
                      sprintf("Shared DEGs (n = %d)", sum(is_sh))))
    legend("topleft", bty = "o", bg = "white", box.col = "grey75",
           cex = 0.7, legend = strsplit(note, "\n")[[1]])
  }

  scatter(rg, st_r,
          "log2FC, main model (age + sex)",
          "log2FC, region-adjusted model",
          "B  Region confounding check",
          sprintf("Spearman rho = %.3f (all genes)\nSpearman rho = %.3f (shared DEGs)\nMain DEGs same direction = %.0f%%",
                  st_r$rho_all, st_r$rho_deg, st_r$concord))

  scatter(ad, st_a,
          "log2FC, main model (age + sex)",
          "log2FC, adults only (>=18 y)",
          "C  Adult-subset sensitivity",
          sprintf("Spearman rho = %.3f (all genes)\nSpearman rho = %.3f (shared DEGs)\nMain DEGs same direction = %.0f%%",
                  st_a$rho_all, st_a$rho_deg, st_a$concord))

  ## ---- D：三模型 DEG 交集（手绘韦恩）----
  A <- length(sm); B <- length(sr); C <- length(sa)
  AB <- length(intersect(sm, sr)); AC <- length(intersect(sm, sa))
  BC <- length(intersect(sr, sa)); ABC <- length(Reduce(intersect, list(sm, sr, sa)))
  plot(NA, xlim = c(-1.55, 1.55), ylim = c(-1.5, 1.25), asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = "D  Overlap of DEG sets", cex.main = 1.0,
       font.main = 2)
  symbols(c(-0.55, 0.55, 0), c(0.22, 0.22, -0.48), circles = rep(0.88, 3),
          add = TRUE, inches = FALSE,
          bg = c(rgb(0.23, 0.51, 0.96, 0.35), rgb(0.20, 0.72, 0.38, 0.35),
                 rgb(0.93, 0.47, 0.20, 0.35)),
          fg = c("steelblue", "forestgreen", "darkorange"))
  txt <- function(x, y, s, cex = 0.95, font = 1) text(x, y, s, cex = cex, font = font)
  txt(-0.88, 0.75, A, 1.05, 2); txt(0.88, 0.75, B, 1.05, 2)
  txt(0.00, -1.02, C, 1.05, 2)
  txt(0.00, 0.52, AB - ABC); txt(-0.38, -0.10, AC - ABC); txt(0.38, -0.10, BC - ABC)
  txt(0.00, 0.02, ABC, 1.1, 2)
  legend("bottom", bty = "n", cex = 0.72, ncol = 3, xpd = NA,
         legend = c(sprintf("Main (%d)", A), sprintf("+Region (%d)", B),
                    sprintf("Adults (%d)", C)),
         text.col = c("steelblue", "forestgreen", "darkorange"))
}

png(file.path(FIG_DIR, "S1_sensitivity.png"), 2600, 2300, res = PLOT_THEME$save_dpi)
draw_s1(); dev.off()
pdf(file.path(FIG_DIR, "S1_sensitivity.pdf"), width = 8.6, height = 7.6)
draw_s1(); dev.off()

## 明细表，供补充材料
out <- data.frame(
  comparison = c("Main vs region-adjusted", "Main vs adults only"),
  rho_all_genes = c(st_r$rho_all, st_a$rho_all),
  rho_shared_DEG = c(st_r$rho_deg, st_a$rho_deg),
  n_shared_DEG = c(st_r$n_shared, st_a$n_shared),
  direction_concordance_pct = c(st_r$concord, st_a$concord))
write.csv(out, file.path(RES_DIR, "02b_sensitivity_stats.csv"), row.names = FALSE)
print(out)
message("02b 完成")
