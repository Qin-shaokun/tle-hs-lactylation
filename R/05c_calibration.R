# ============================================================
# 05c_calibration.R — 列线图校准（2026-11-04）
# 背景：Firth 模型 AUC=1.0、概率近乎二值，用常规 glm 再校准求斜率
#       同样遭遇完全分离（斜率 17.7 → 校正 22.1，无意义）。
# 方案：
#   1) 校准截距/斜率改用 Firth 再校准（logistf(y ~ logit(p))）→ 有限估计
#   2) 主图用十分位校准（观测率 vs 平均预测 + 95% CI），不用平滑曲线
#   3) 副图画预测概率分布（病例/对照），直观说明"近乎完全分离"
#   4) 报 Brier、AUC、Hosmer-Lemeshow（10 组）
# 输出：figures/05_calibration.pdf/.png、results/05_calibration.csv
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(logistf); library(pROC)})

ml  <- readRDS(file.path(RES_DIR, "04_ml.rds"))
hub <- ml$hub
d   <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
dat <- data.frame(t(d$sub_expr[hub, , drop = FALSE]),
                  y = as.numeric(d$sub_pd$condition == GSE256068$case_label))
n <- nrow(dat); set.seed(42)

f  <- readRDS(file.path(RES_DIR, "05_nomogram_model.rds"))$model
p  <- f$predict
lp <- log(p / (1 - p))
X  <- as.matrix(dat[hub])

## ---- 表观校准（Firth 再校准，避免分离导致的无穷斜率）----
g0 <- logistf(y ~ lp, data = dat, firth = TRUE, pl = FALSE)
a0 <- unname(coef(g0)[1]); b0 <- unname(coef(g0)[2])

brier <- mean((p - dat$y)^2)
auc   <- as.numeric(auc(roc(dat$y, p, quiet = TRUE)))

## ---- Hosmer-Lemeshow（10 组，按预测概率等分）----
ord  <- order(p)
grp  <- cut(seq_len(n), breaks = 10, labels = FALSE)
obs  <- tapply(dat$y[ord], grp, sum)
exp_ <- tapply(p[ord],     grp, sum)
nn   <- as.numeric(table(grp))
hl   <- sum((obs - exp_)^2 / (exp_ * (1 - exp_ / nn)))
hl_p <- 1 - pchisq(hl, df = length(nn) - 2)

p_g   <- tapply(p[ord], grp, mean)
o_g   <- obs / nn
se_g  <- sqrt(o_g * (1 - o_g) / nn)

cat(sprintf("表观校准（Firth）：截距 %.2f  斜率 %.2f\n", a0, b0))
cat(sprintf("Brier %.4f  AUC %.3f  HL chi2 %.2f (df %d, p %.3f)\n",
            brier, auc, hl, length(nn) - 2, hl_p))

## ---- bootstrap 乐观度校正（模型 + 再校准均用 Firth）----
B <- 1000
oa <- ob <- rep(NA_real_, B); ok <- 0
fit_ctrl <- logistf.control(maxit = 500, maxstep = 3)
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  fb  <- tryCatch(logistf(y ~ ., data = dat[idx, ], firth = TRUE,
                          control = fit_ctrl, pl = FALSE), error = function(e) NULL)
  if (is.null(fb)) next
  lp_te <- as.numeric(cbind(1, X) %*% coef(fb))          # 原始样本
  lp_tr <- as.numeric(fb$linear.predictors)               # boot 内
  gt <- tryCatch(logistf(dat$y[idx] ~ lp_tr, firth = TRUE, pl = FALSE), error = function(e) NULL)
  ge <- tryCatch(logistf(dat$y      ~ lp_te, firth = TRUE, pl = FALSE), error = function(e) NULL)
  if (is.null(gt) || is.null(ge)) next
  oa[ok + 1] <- coef(gt)[1] - coef(ge)[1]
  ob[ok + 1] <- coef(gt)[2] - coef(ge)[2]
  ok <- ok + 1
}
opt_a <- mean(oa[seq_len(ok)]); opt_b <- mean(ob[seq_len(ok)])
a_c <- a0 - opt_a; b_c <- b0 - opt_b
cat(sprintf("bootstrap %d/%d；optimism 截距 %.2f 斜率 %.2f → 校正后 %.2f / %.2f\n",
            ok, B, opt_a, opt_b, a_c, b_c))

## ---- 图：左十分位校准，右预测概率分布 ----
png(file.path(FIG_DIR, "05_calibration.png"), 2400, 1250, res = PLOT_THEME$save_dpi)
pdf(file.path(FIG_DIR, "05_calibration.pdf"), width = 10, height = 5.4)
draw_cal <- function() {
  par(mfrow = c(1, 2), mar = c(4.4, 4.8, 3.0, 1.8), mgp = c(2.6, 0.7, 0))

  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Mean predicted probability", ylab = "Observed proportion",
       main = "Calibration by decile")
  abline(0, 1, lty = 2, col = "grey55", lwd = 1.5)
  segments(p_g, pmax(o_g - 1.96 * se_g, 0), p_g, pmin(o_g + 1.96 * se_g, 1),
           col = "grey40", lwd = 1.4)
  lines(p_g, o_g, col = "steelblue", lwd = 1.6)
  points(p_g, o_g, pch = 19, cex = 1.2)
  legend("topleft", bty = "n", cex = 0.85, lty = c(2, 1), lwd = c(1.5, 1.6),
         col = c("grey55", "steelblue"), legend = c("Ideal", "Observed"))
  text(0.98, 0.32, sprintf("Brier = %.4f\nAUC = %.2f\nH-L p = %.2f\nCalibration slope %.2f\n(Firth; corrected %.2f)",
                           brier, auc, hl_p, b0, b_c),
       adj = c(1, 1), cex = 0.8, col = "grey25")

  ## 右：预测概率分布（说明近乎完全分离）
  brks <- seq(0, 1, 0.05)
  h0 <- hist(p[dat$y == 0], breaks = brks, plot = FALSE)
  h1 <- hist(p[dat$y == 1], breaks = brks, plot = FALSE)
  ymax <- max(c(h0$counts, h1$counts)) * 1.15
  plot(h0, col = rgb(0.2, 0.45, 0.7, 0.55), border = "white", xlim = c(0, 1),
       ylim = c(0, ymax), xlab = "Predicted probability of TLE-HS",
       ylab = "Number of samples", main = "Predicted probability distribution")
  plot(h1, col = rgb(0.8, 0.25, 0.2, 0.55), border = "white", add = TRUE)
  legend("topleft", bty = "n", cex = 0.85, fill = c(rgb(0.2, 0.45, 0.7, 0.55),
         rgb(0.8, 0.25, 0.2, 0.55)), border = "white",
         legend = c("Control (n = 14)", "TLE-HS (n = 64)"))
}
draw_cal(); dev.off()
draw_cal(); dev.off()

out <- data.frame(
  metric = c("calib_intercept_firth", "calib_slope_firth",
             "calib_intercept_corrected", "calib_slope_corrected",
             "optimism_intercept", "optimism_slope",
             "brier", "auc_insample", "hl_chisq", "hl_df", "hl_p", "bootstrap_ok"),
  value  = round(c(a0, b0, a_c, b_c, opt_a, opt_b,
                   brier, auc, hl, length(nn) - 2, hl_p, ok), 4))
write.csv(out, file.path(RES_DIR, "05_calibration.csv"), row.names = FALSE)
print(out)
message("05c 完成：校准图已写入 figures/05_calibration.pdf 与 .png")
