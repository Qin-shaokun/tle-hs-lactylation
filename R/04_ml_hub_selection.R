# ============================================================
# 04_ml_hub_selection.R — 机器学习筛选 hub gene（LASSO + SVM-RFE + 随机森林）
# 【2026-11-04 重写版】对接 02_deg_main.rds / 03_wgcna.rds
# 候选池：显著 DEG ∩ TLE 相关 WGCNA 模块（非 grey）；
#         乳酸化基因集成员单独标注（解释性亮点，不作为硬过滤）
# 输出：results/04_hub_genes.csv, 04_ml_scores.csv；
#       figures/04_lasso.png, 04_roc_train.pdf/png
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(glmnet); library(e1071); library(randomForest)
  library(pROC); library(ggplot2)
})

d <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
w <- readRDS(file.path(RES_DIR, "03_wgcna.rds"))
gs <- read.csv(GENESET_CSV)

## ---- 1. 候选池 ----
deg_pool <- d$deg$gene[abs(d$deg$logFC) >= 0.585 & d$deg$adj.P.Val <= 0.05]

## 与 TLE 相关性排序的非 grey 模块
mt <- w$mm[, "TLE", drop = FALSE]
mt <- cbind.data.frame(module = sub("^ME", "", rownames(mt)), r = mt[, 1])
mt <- mt[mt$module != "grey", ]
mt <- mt[order(-abs(mt$r)), ]
top_mods <- head(mt$module, 2)
cat("与 TLE 最相关模块：", paste(top_mods, collapse = ", "),
    "（r =", paste(sprintf("%.2f", head(mt$r, 2)), collapse = ", "), "）\n")

mod_pool <- w$modules$gene[w$modules$module %in% top_mods]
candidates <- intersect(deg_pool, mod_pool)
if (length(candidates) < 10) {   # 兜底：DEG ∩ 全部非 grey 模块
  candidates <- intersect(deg_pool, w$modules$gene[w$modules$module != "grey"])
}
if (length(candidates) < 10) candidates <- deg_pool
cat("候选池（DEG ∩ 模块）：", length(candidates), "个；",
    "其中乳酸化基因集成员：", sum(candidates %in% gs$gene_symbol), "个\n")

X <- t(d$sub_expr[candidates, ])
y <- as.numeric(d$sub_pd$condition == GSE256068$case_label)  # TLE-HS=1, 对照=0

## ---- 2. LASSO ----
cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
lasso_genes <- coef(cv, s = THRESH$lasso_lambda) |>
  (\(m) rownames(m)[as.numeric(m[, 1]) != 0])() |> setdiff("(Intercept)")
cat("LASSO 保留：", length(lasso_genes), "\n")

png(file.path(FIG_DIR, "04_lasso.png"), 2000, 1000, res = PLOT_THEME$save_dpi)
par(mfrow = c(1, 2))
plot(cv$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(cv$lambda.min), lty = 2, col = "red")
plot(cv)
dev.off()

## ---- 3. SVM-RFE（标准实现：线性核，|w|^2 排序，每轮剔除 10%）----
svmrfe <- function(X, y, frac_drop = 0.1) {
  feats <- colnames(X); n <- length(feats); ranks <- c()
  while (length(feats) > 2) {
    m <- svm(X[, feats, drop = FALSE], factor(y), scale = TRUE,
             kernel = "linear", cost = 1)
    w_vec <- setNames(as.vector(t(m$coefs) %*% m$SV), colnames(m$SV))  # 权重向量（保留特征名！）
    score <- w_vec^2
    n_drop <- max(1, floor(length(feats) * frac_drop))
    rk <- names(sort(score, decreasing = FALSE))[seq_len(n_drop)]  # 权重最小者先剔
    ranks <- c(ranks, rk); feats <- setdiff(feats, rk)
  }
  c(feats, ranks)   # 排名越前越重要
}
svm_rank <- tryCatch(svmrfe(X, y), error = function(e) {
  warning("SVM-RFE 失败：", conditionMessage(e)); character(0) })
svm_genes <- head(svm_rank, 10)
cat("SVM-RFE top10：", paste(svm_genes, collapse = ", "), "\n")

## ---- 4. 随机森林 ----
rf <- randomForest(X, factor(y), ntree = THRESH$rf_ntree, importance = TRUE)
imp <- sort(importance(rf)[, 1], decreasing = TRUE)
rf_genes <- names(head(imp, 15))
cat("RF top15：", paste(rf_genes, collapse = ", "), "\n")

## ---- 5. 交集 → hub gene ----
hub3 <- Reduce(intersect, list(lasso_genes, svm_genes, rf_genes))
hub2 <- setdiff(Reduce(intersect, list(lasso_genes, rf_genes)), hub3)
hub  <- c(hub3, hub2)
if (length(hub) < 2) {   # 最终兜底：三算法各自排名前列的并集取 RF∩LASSO
  hub <- intersect(head(svm_rank, 20), rf_genes)
}
hub <- head(unique(hub), 10)

scores <- data.frame(
  gene = candidates,
  in_lasso = candidates %in% lasso_genes,
  svm_rank = match(candidates, svm_rank),   # NA = 未进 top
  rf_importance = imp[match(candidates, names(imp))],
  in_lactylation_geneset = candidates %in% gs$gene_symbol,
  module = w$modules$module[match(candidates, w$modules$gene)],
  stringsAsFactors = FALSE)
scores <- scores[order(-scores$in_lasso, is.na(scores$svm_rank),
                       -replace(scores$rf_importance, is.na(scores$rf_importance), 0)), ]
write.csv(scores, file.path(RES_DIR, "04_ml_scores.csv"), row.names = FALSE)
write.csv(data.frame(hub_gene = hub,
                     in_lasso = hub %in% lasso_genes,
                     in_svm_top10 = hub %in% svm_genes,
                     in_rf_top15 = hub %in% rf_genes,
                     in_lactylation_geneset = hub %in% gs$gene_symbol),
          file.path(RES_DIR, "04_hub_genes.csv"), row.names = FALSE)
cat("\n★ Hub genes：", paste(hub, collapse = ", "), "\n")
cat("其中乳酸化基因集成员：", paste(intersect(hub, gs$gene_symbol), collapse = ", "), "\n")

## ---- 6. 训练集 ROC（单基因，多 panel 并排，避免循环覆盖）----
auc_tbl <- data.frame()
for (g in hub) {
  r <- roc(y, d$sub_expr[g, ], quiet = TRUE)
  auc_tbl <- rbind(auc_tbl, data.frame(gene = g, AUC = as.numeric(auc(r))))
}
png(file.path(FIG_DIR, "04_roc_train.png"), 1000 * length(hub), 1000,
    res = PLOT_THEME$save_dpi)
par(mfrow = c(1, length(hub)))
for (i in seq_along(hub)) {
  r <- roc(y, d$sub_expr[hub[i], ], quiet = TRUE)
  plot(r, main = paste0(hub[i], "  AUC=", round(as.numeric(auc(r)), 3)),
       col = "steelblue", legacy.axes = TRUE)
}
dev.off()
pdf(file.path(FIG_DIR, "04_roc_train.pdf"), width = 5 * length(hub), height = 5)
par(mfrow = c(1, length(hub)))
for (i in seq_along(hub)) {
  r <- roc(y, d$sub_expr[hub[i], ], quiet = TRUE)
  plot(r, main = paste0(hub[i], "  AUC=", round(as.numeric(auc(r)), 3)),
       col = "steelblue", legacy.axes = TRUE)
}
dev.off()
write.csv(auc_tbl, file.path(RES_DIR, "04_hub_auc.csv"), row.names = FALSE)
print(auc_tbl)

saveRDS(list(hub = hub, lasso = lasso_genes, svm_rank = svm_rank,
             rf_genes = rf_genes, candidates = candidates, auc = auc_tbl),
        file.path(RES_DIR, "04_ml.rds"))
message("04 完成")
