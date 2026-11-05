# ============================================================
# 02_deg_main.R — 主分析：TLE-HS 颞叶 vs 对照皮层（GSE256068）
# 【2026-11-03 修正版】
#   原"海马子集"设计不成立（海马样本全是对照），改为：
#     TLE-HS(64, 颞叶) vs ControlCortex(14, 额/颞/顶叶混合)
#   主模型：~ age + sex + condition
#   敏感性 1：~ age + sex + region + condition（校正脑区混杂）
#   敏感性 2：成人子集（≥18 岁）
# 矩阵格式：空格分隔、字段带引号、首行为样本编号（对应 metadata$description）
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(pheatmap)
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})

## ---- 1. 读入表达矩阵 ----
f <- file.path(DATA_DIR, "GSE256068_processed_data.csv.gz")
stopifnot(file.exists(f))
con <- gzfile(f, "r"); hdr_line <- readLines(con, n = 1, warn = FALSE); close(con)
samples <- scan(text = hdr_line, what = "", quiet = TRUE)
cat("表头样本数:", length(samples), "\n")

dt <- fread(f, skip = 1, header = FALSE, sep = " ", quote = "\"", data.table = FALSE)
genes_ens <- dt[[1]]
mat <- as.matrix(dt[, -1]); storage.mode(mat) <- "double"
rownames(mat) <- genes_ens
stopifnot(ncol(mat) == length(samples))
colnames(mat) <- samples
cat("矩阵：", nrow(mat), "行 ×", ncol(mat), "列；取值范围 ",
    round(min(mat), 2), "~", round(max(mat), 2), "\n")

## ---- 2. Ensembl → Symbol ----
ens <- sub("\\..*$", "", rownames(mat))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens),
          keytype = "ENSEMBL", columns = "SYMBOL"))
map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$ENSEMBL), ]
idx <- match(ens, map$ENSEMBL)
keep <- !is.na(idx)
mat <- mat[keep, ]; sym <- map$SYMBOL[idx[keep]]
## 同 symbol 多探针：保留平均表达最高的
ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
mat <- mat[ord, ]; sym <- sym[ord]
dup <- duplicated(sym)
mat <- mat[!dup, ]; rownames(mat) <- sym[!dup]
cat("转 symbol 后：", nrow(mat), "个基因\n")

## ---- 3. 样本表 ----
pd <- read.csv(file.path(DATA_DIR, "GSE256068_metadata_raw.csv"),
               check.names = FALSE, fileEncoding = "UTF-8")
names(pd) <- sub(":ch1$", "", names(pd))
pd$group  <- pd$`disease state`
pd$age    <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pd$age)))
pd$sex    <- factor(pd$Sex)
pd$region <- factor(pd$tissue)
pd$sample_id <- pd$description          # 与矩阵列名对应
stopifnot(all(samples %in% pd$sample_id))

pd <- pd[match(samples, pd$sample_id), ]  # 按矩阵列顺序对齐
sel <- pd$group %in% GSE256068$groups_keep
sub_expr <- mat[, sel, drop = FALSE]
sub_pd <- pd[sel, ]
sub_pd$condition <- factor(sub_pd$group,
                           levels = c(GSE256068$ctrl_label, GSE256068$case_label))
sub_pd$region <- droplevels(sub_pd$region)
cat("\n分析样本：", ncol(sub_expr), "例\n"); print(table(sub_pd$condition, sub_pd$region))
cat("年龄："); print(tapply(sub_pd$age, sub_pd$condition,
                            function(x) round(c(n = length(x), mean = mean(x), sd = sd(x)), 1)))

## ---- 4. limma ----
run_limma <- function(expr, pd, fml, coef_name) {
  design <- model.matrix(as.formula(fml), data = pd)
  fit <- eBayes(lmFit(expr, design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt
}
coef_main <- paste0("condition", GSE256068$case_label)
deg <- run_limma(sub_expr, sub_pd, "~ age + sex + condition", coef_main)
sig <- subset(deg, abs(logFC) >= THRESH$lfc & adj.P.Val <= THRESH$fdr)
cat("\n【主模型】显著 DEG：", nrow(sig), "（上调 ", sum(sig$logFC > 0),
    " / 下调 ", sum(sig$logFC < 0), "）\n")
write.csv(deg, file.path(RES_DIR, "02_deg_all.csv"), row.names = FALSE)
write.csv(sig, file.path(RES_DIR, "02_deg_sig.csv"), row.names = FALSE)

## ---- 5. 敏感性 1：校正脑区 ----
deg_reg <- tryCatch(
  run_limma(sub_expr, sub_pd, "~ age + sex + region + condition", coef_main),
  error = function(e) { message("region 模型失败：", conditionMessage(e)); NULL })
if (!is.null(deg_reg)) {
  write.csv(deg_reg, file.path(RES_DIR, "02_deg_sens_region.csv"), row.names = FALSE)
  cat("【校正脑区】显著 DEG：",
      nrow(subset(deg_reg, abs(logFC) >= THRESH$lfc & adj.P.Val <= THRESH$fdr)), "\n")
}

## ---- 6. 敏感性 2：成人子集 ----
ad <- sub_pd$age >= 18
if (sum(ad) > 10 && length(unique(sub_pd$condition[ad])) == 2) {
  deg_ad <- tryCatch(
    run_limma(sub_expr[, ad], droplevels(sub_pd[ad, ]), "~ age + sex + condition", coef_main),
    error = function(e) { message("成人子集失败：", conditionMessage(e)); NULL })
  if (!is.null(deg_ad)) {
    write.csv(deg_ad, file.path(RES_DIR, "02_deg_sens_adult.csv"), row.names = FALSE)
    cat("【成人子集 ", sum(ad), " 例】显著 DEG：",
        nrow(subset(deg_ad, abs(logFC) >= THRESH$lfc & adj.P.Val <= THRESH$fdr)), "\n")
  }
} else cat("成人子集样本不足，跳过\n")

## ---- 7. 火山图（重点基因 ggrepel 标注）----
deg$sig_flag <- ifelse(abs(deg$logFC) >= THRESH$lfc & deg$adj.P.Val <= THRESH$fdr,
                       ifelse(deg$logFC > 0, "up", "down"), "ns")
lab_genes <- intersect(c("TNF", "CSF1R", "CLEC7A", "SELPLG", "PLD4", "CD37",
                         "HK2", "GPD1", "PFKFB3", "NLRP3", "BDNF", "SIAH2"),
                       deg$gene)
p <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig_flag)) +
  geom_point(size = .8, alpha = .6) +
  scale_color_manual(values = c(up = "#C0392B", down = "#2471A3", ns = "grey75")) +
  geom_vline(xintercept = c(-THRESH$lfc, THRESH$lfc), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(THRESH$fdr), linetype = 2, color = "grey40") +
  ggrepel::geom_label_repel(data = subset(deg, gene %in% lab_genes),
                            aes(label = gene), size = 2.9, fill = "white",
                            label.size = 0.15, label.padding = 0.12,
                            box.padding = 0.6, point.padding = 0.25,
                            min.segment.length = 0, max.overlaps = Inf,
                            segment.size = 0.3, segment.colour = "grey40",
                            color = "black", show.legend = FALSE, seed = 2) +
  labs(x = "log2FC (TLE-HS vs Control)", y = expression(-log[10]~adjusted~italic(P)),
       color = "",
       title = paste0("GSE256068: TLE-HS vs Control (n=", ncol(sub_expr), ")")) +
  theme_bw(PLOT_THEME$base_size)
ggsave(file.path(FIG_DIR, "02_volcano.png"), p, width = 6, height = 5,
       dpi = PLOT_THEME$save_dpi)

## ---- 8. 热图 ----
gs <- read.csv(GENESET_CSV)
gs_in <- intersect(gs$gene_symbol, rownames(sub_expr))
plot_genes <- intersect(union(head(deg$gene, 50), gs_in), rownames(sub_expr))
if (length(plot_genes) > 2) {
  m2 <- t(scale(t(sub_expr[plot_genes, ]))); m2[is.na(m2)] <- 0
  ann <- data.frame(Group = sub_pd$condition, Region = sub_pd$region,
                    row.names = colnames(sub_expr))
  pheatmap(m2, annotation_col = ann, show_rownames = TRUE,
           show_colnames = FALSE, fontsize_row = 6,
           filename = file.path(FIG_DIR, "02_heatmap.png"), width = 9, height = 10)
}

saveRDS(list(deg = deg, sig = sig, sens_region = deg_reg,
             sub_expr = sub_expr, sub_pd = sub_pd),
        file.path(RES_DIR, "02_deg_main.rds"))
message("02 完成")
