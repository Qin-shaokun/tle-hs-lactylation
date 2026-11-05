# ============================================================
# 05b_nomogram_firth.R — Firth 校正列线图（2026-11-04）
# 背景：原 05 脚本 lrm 因完全分离导致系数 SE 巨大，列线图不可用
# 方案：logistf（Firth 惩罚似然）拟合 → 按 nomogram 数学定义自绘
#   每个预测变量的点数 = 100 * (beta_i*x_i - min贡献) / 总贡献极差
#   总点数 → 线性预测子 → 概率
# 输出：figures/05_nomogram.pdf / .png（覆盖旧图），
#       results/05_nomogram_model.rds
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(logistf); library(pROC)})

ml  <- readRDS(file.path(RES_DIR, "04_ml.rds"))
hub <- ml$hub
stopifnot(length(hub) >= 1)

d   <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
dat <- data.frame(t(d$sub_expr[hub, , drop = FALSE]),
                  y = as.numeric(d$sub_pd$condition == GSE256068$case_label))
cat("样本：", nrow(dat), " 事件：", sum(dat$y), "\n")

## ---- Firth logistic（完全分离数据需放宽迭代）----
f <- logistf(y ~ ., data = dat, firth = TRUE,
             control   = logistf.control(maxit = 500, maxstep = 3),
             plcontrol = logistpl.control(maxit = 200, maxstep = 5))
se <- sqrt(diag(f$var))
tab <- data.frame(
  coef         = round(coef(f), 4),
  OR           = round(exp(coef(f)), 3),
  ci_low_wald  = round(exp(coef(f) - 1.96 * se), 3),
  ci_high_wald = round(exp(coef(f) + 1.96 * se), 3),
  p_wald       = signif(f$prob, 3))
cat("注：Firth 系数 / Wald CI / Wald p（profile-PL CI 在本数据未收敛，不报）\n")
print(tab)

## 样本内 AUC（Firth 拟合概率）
pf   <- f$predict
r_tr <- roc(dat$y, pf, quiet = TRUE)
cat("训练集 AUC（Firth，样本内，仅作参考）：", round(as.numeric(auc(r_tr)), 3), "\n")

## ---- 点数换算（rms 惯例：效应最大的变量占满 100 点）----
beta      <- coef(f)[-1]
intercept <- coef(f)[1]
contrib   <- lapply(hub, function(g) beta[g] * dat[[g]]); names(contrib) <- hub
cmin      <- sapply(contrib, min); cmax <- sapply(contrib, max)
span      <- cmax - cmin
k         <- 100 / max(span)              # 点数/单位线性预测子
pts_of    <- function(g, x) (beta[g] * x - cmin[[g]]) * k
Tmin      <- sum(cmin)                    # 总贡献下界
TPmax     <- sum(span) * k                # 总点数上限
lp_of_pts <- function(P) intercept + Tmin + P / k
p_of_pts  <- function(P) 1 / (1 + exp(-lp_of_pts(P)))

cat("\n点数标定：k =", round(k, 3), " 总点数上限 =", round(TPmax, 1),
    " 概率全距", signif(p_of_pts(0), 2), "~", signif(p_of_pts(TPmax), 2), "\n")

## ---- 绘图 ----
draw_nomogram <- function() {
  n_lines <- length(hub) + 3              # 点数 + 基因 + 总点数 + 概率
  yy <- seq(n_lines, 1)
  par(mar = c(2.2, 10, 2.4, 2.5), xpd = TRUE)
  plot(NA, xlim = c(-0.30 * TPmax, TPmax * 1.02), ylim = c(0.45, n_lines + 0.55),
       axes = FALSE, xlab = "", ylab = "",
       main = "Nomogram for P(TLE-HS) - Firth logistic")

  axis_line <- function(y, at, labels, title, cex = 0.85, stagger = FALSE) {
    segments(0, y, max(at), y, lwd = 1.3, col = "grey30")
    for (i in seq_along(at)) {
      segments(at[i], y - 0.10, at[i], y + 0.10, lwd = 1, col = "grey30")
      dy <- if (stagger && i %% 2 == 0) 0.58 else 0.30
      text(at[i], y + dy, labels[i], cex = cex, adj = c(0.5, 0))
    }
    text(-0.035 * TPmax, y, title, cex = 0.95, adj = c(1, 0.5), font = 2)
  }

  ## 1) Points 0-100
  axis_line(yy[1], seq(0, 100, 20), as.character(seq(0, 100, 20)), "Points")

  ## 2) 基因轴（各自 0 ~ 该基因最大贡献点数）
  for (i in seq_along(hub)) {
    g   <- hub[i]
    brk <- pretty(range(dat[[g]]), n = 4)
    brk <- brk[brk >= min(dat[[g]]) & brk <= max(dat[[g]])]
    axis_line(yy[1 + i], pts_of(g, brk), formatC(brk, format = "fg", digits = 2),
              paste0(g, "\nlog2 expression"))
  }

  ## 3) Total Points
  tp <- seq(0, floor(TPmax / 20) * 20, by = 20)
  axis_line(yy[1 + length(hub) + 1], tp, as.character(tp), "Total points")

  ## 4) 概率轴
  cand <- c(0.05, 0.2, 0.5, 0.8, 0.95)
  p_lo <- p_of_pts(0); p_hi <- p_of_pts(TPmax)
  cand <- cand[cand > p_lo & cand < p_hi]
  p_pts <- k * (log(cand / (1 - cand)) - intercept - Tmin)
  axis_line(yy[n_lines], p_pts, sprintf("%.2f", cand), "P(TLE-HS)",
            cex = 0.8, stagger = TRUE)

  mtext("Sum points, then read probability on the bottom axis", side = 1,
        line = 1.2, cex = 0.8, col = "grey35")
}

png(file.path(FIG_DIR, "05_nomogram.png"), 2100, 1500, res = PLOT_THEME$save_dpi)
draw_nomogram(); dev.off()
pdf(file.path(FIG_DIR, "05_nomogram.pdf"), width = 9, height = 6.4)
draw_nomogram(); dev.off()

saveRDS(list(model = f, hub = hub, method = "Firth (logistf)",
             coef = coef(f), wald_ci = cbind(coef(f) - 1.96 * se, coef(f) + 1.96 * se),
             p_wald = f$prob, auc_insample = as.numeric(auc(r_tr)),
             calib = list(Tmin = Tmin, TPmax = TPmax, k = k, intercept = intercept)),
        file.path(RES_DIR, "05_nomogram_model.rds"))
message("05b 完成：Firth 列线图已写入 figures/05_nomogram.pdf 与 .png")
