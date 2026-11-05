# ============================================================
# 05_validation_roc.R — 外部验证 + 列线图（2026-11-04 重写）
# 验证集：
#   V1 GSE190451：3 TLE vs 3 非癫痫对照，脑组织 bulk RNA-seq，
#      官方 TPM 矩阵（同平台 GPL24676）
#   V2 GSE71058：12 例 mTLE 齿状回颗粒细胞（LCM），HS 5 vs noHS 7，
#      官方 raw counts（按患者聚合后 log2CPM）→ 硬化程度方向性验证
# 训练模型：hub gene 多基因 logistic（GSE256068 主分析）
# 输出：results/05_validation.csv, 05_nomogram_model.rds；
#       figures/05_roc_valid.png (Fig 4A), S5_roc_severity.png (Fig S5)
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(pROC); library(ggplot2); library(rms)})

ml  <- readRDS(file.path(RES_DIR, "04_ml.rds"))
hub <- ml$hub
stopifnot(length(hub) >= 1)
cat("Hub genes:", paste(hub, collapse = ", "), "\n")

val_tbl <- data.frame()

## ============ V1: GSE190451（病例-对照 AUC）============
v1 <- read.delim(gzfile(file.path(DATA_DIR, "GSE190451_TLE_samples.txt.gz")),
                 check.names = FALSE)
stopifnot(all(hub %in% v1$Gene_name))
emat1 <- as.matrix(v1[, paste0("Normal_", 1:3), drop = FALSE])
emat1 <- cbind(emat1, as.matrix(v1[, paste0("TLE_", 1:3)]))
rownames(emat1) <- v1$Gene_name
grp1 <- factor(c(rep("control", 3), rep("case", 3)), levels = c("control", "case"))
emat1 <- log2(emat1 + 1)   # TPM → log2(TPM+1)，与训练集可比较
cat("\n[V1 GSE190451] log2(TPM+1) 后取值范围:", round(min(emat1), 2), "~", round(max(emat1), 2), "\n")

for (g in hub) {
  r <- roc(grp1, emat1[g, ], quiet = TRUE, levels = c("control", "case"),
           direction = "auto")
  val_tbl <- rbind(val_tbl, data.frame(
    dataset = "GSE190451", gene = g, n_case = 3, n_ctrl = 3,
    auc = as.numeric(auc(r)),
    ci_low = as.numeric(ci.auc(r))[1], ci_high = as.numeric(ci.auc(r))[3],
    logFC_valid = mean(emat1[g, 4:6]) - mean(emat1[g, 1:3]), p_ttest = NA))
}
## 联合 panel AUC（线性组合，训练集系数不可靠时改用等权均值）
panel_v1 <- roc(grp1, as.vector(colMeans(emat1[hub, , drop = FALSE])), quiet = TRUE,
                levels = c("control", "case"), direction = "auto")
val_tbl <- rbind(val_tbl, data.frame(
  dataset = "GSE190451", gene = "Panel(mean)", n_case = 3, n_ctrl = 3,
  auc = as.numeric(auc(panel_v1)),
  ci_low = NA, ci_high = NA, p_ttest = NA,
  logFC_valid = NA))
cat("[V1] 单基因与 panel AUC：\n"); print(val_tbl[val_tbl$dataset == "GSE190451", ])

## ============ V2: GSE71058（HS vs noHS 方向性）============
cnt <- read.delim(gzfile(file.path(DATA_DIR, "GSE71058_rawcounts.txt.gz")),
                  check.names = FALSE)
## 同名基因符号多行：按符号聚合求和
agg_gene <- rowsum(as.matrix(cnt[, -1]), group = cnt[[1]], reorder = FALSE)
cnt <- agg_gene
rownames(cnt) <- rownames(agg_gene)
## 患者分组（来自 series matrix characteristics，与列顺序一一对应）
pid  <- sub("(rep|cap)[0-9]+$", "", colnames(cnt))
diag_map <- c(dukeepi1513 = "HS", dukeepi1930 = "HS", dukeepi3949 = "HS",
              dukeepi1046 = "HS", dukeepi2638 = "HS",
              dukeepi3169 = "noHS", dukeepi4000 = "noHS", dukeepi2213 = "noHS",
              dukeepi2953 = "noHS", dukeepi1895 = "noHS", dukeepi2948 = "noHS",
              dukeepi3956 = "noHS")
stopifnot(all(pid %in% names(diag_map)))
## 按患者聚合（replicate counts 求和）
agg <- t(sapply(split(seq_len(ncol(cnt)), pid), function(ii)
  rowSums(cnt[, ii, drop = FALSE])))
grp2 <- factor(ifelse(diag_map[rownames(agg)] == "HS", "HS", "noHS"),
               levels = c("noHS", "HS"))
cat("\n[V2 GSE71058] 患者：", table(grp2), "（noHS, HS）\n")
## log2 CPM（基因×患者）
aggT <- t(agg)                                    # 基因 × 患者
lcpm2 <- log2(sweep(aggT, 2, colSums(aggT), "/") * 1e6 + 1)
stopifnot(all(hub %in% rownames(lcpm2)))

for (g in hub) {
  tt <- tryCatch(t.test(lcpm2[g, grp2 == "HS"], lcpm2[g, grp2 == "noHS"]),
                 error = function(e) NULL)
  r  <- roc(grp2, lcpm2[g, ], quiet = TRUE,
            levels = c("noHS", "HS"), direction = "auto")
  val_tbl <- rbind(val_tbl, data.frame(
    dataset = "GSE71058", gene = g, n_case = sum(grp2 == "HS"),
    n_ctrl = sum(grp2 == "noHS"),
    auc = as.numeric(auc(r)),
    ci_low = NA, ci_high = NA,
    logFC_valid = mean(lcpm2[g, grp2 == "HS"]) - mean(lcpm2[g, grp2 == "noHS"]),
    p_ttest = if (!is.null(tt)) tt$p.value else NA))
}
cat("[V2] HS vs noHS 方向性验证：\n")
print(val_tbl[val_tbl$dataset == "GSE71058", ])

write.csv(val_tbl, file.path(RES_DIR, "05_validation.csv"), row.names = FALSE)

## ============ Fig 4A：GSE190451 病例-对照验证（单panel）============
## ============ Fig S5：GSE71058 严重度验证（单panel，补充图）==========
## 2026-11-04 按导师意见拆分：原双panel合并图改为两个独立文件
roc_v1 <- lapply(hub, function(g) roc(grp1, emat1[g, ], quiet = TRUE,
                                      levels = c("control", "case"), direction = "auto"))
roc_v2 <- lapply(hub, function(g) roc(grp2, lcpm2[g, ], quiet = TRUE,
                                      levels = c("noHS", "HS"), direction = "auto"))
cols_g <- c("firebrick", "darkorange", "purple")[seq_along(hub)]

draw_v1 <- function() {
  par(mar = c(4.8, 4.8, 3.6, 1.2), mgp = c(2.7, 0.7, 0))
  plot(roc_v1[[1]], col = cols_g[1], legacy.axes = TRUE, lwd = 2,
       main = "GSE190451: TLE (n = 3) vs control (n = 3)")
  for (i in seq_along(hub)[-1]) plot(roc_v1[[i]], add = TRUE, col = cols_g[i], lwd = 2)
  plot(panel_v1, add = TRUE, col = "steelblue", lty = 2, lwd = 2)
  abline(0, 1, col = "grey70", lty = 3)
  legend("bottomright", bty = "n", cex = 0.85, lwd = 2, lty = c(rep(1, length(hub)), 2),
         col = c(cols_g, "steelblue"),
         legend = c(sprintf("%s (AUC = %.3f)", hub, sapply(roc_v1, function(r) as.numeric(auc(r)))),
                    sprintf("Two-gene panel (AUC = %.3f)", as.numeric(auc(panel_v1)))))
}

draw_v2 <- function() {
  par(mar = c(4.8, 4.8, 3.6, 1.2), mgp = c(2.7, 0.7, 0))
  plot(roc_v2[[1]], col = cols_g[1], legacy.axes = TRUE, lwd = 2,
       main = "GSE71058: HS (n = 5) vs noHS (n = 7)")
  for (i in seq_along(hub)[-1]) plot(roc_v2[[i]], add = TRUE, col = cols_g[i], lwd = 2)
  abline(0, 1, col = "grey70", lty = 3)
  legend("bottomright", bty = "n", cex = 0.85, lwd = 2, col = cols_g,
         legend = sprintf("%s (AUC = %.3f)", hub, sapply(roc_v2, function(r) as.numeric(auc(r)))))
}

png(file.path(FIG_DIR, "05_roc_valid.png"), 1500, 1500, res = PLOT_THEME$save_dpi)
draw_v1(); dev.off()
pdf(file.path(FIG_DIR, "05_roc_valid.pdf"), width = 5, height = 5)
draw_v1(); dev.off()

png(file.path(FIG_DIR, "S5_roc_severity.png"), 1500, 1500, res = PLOT_THEME$save_dpi)
draw_v2(); dev.off()
pdf(file.path(FIG_DIR, "S5_roc_severity.pdf"), width = 5, height = 5)
draw_v2(); dev.off()

## ============ 训练集列线图 → 已移至 05b（Firth）+ 05c（校准）============
## 原 lrm 因完全分离不可用，此处不再绘图，避免覆盖 Firth 版
message("05 完成（列线图见 05b_nomogram_firth.R，校准见 05c_calibration.R）")
