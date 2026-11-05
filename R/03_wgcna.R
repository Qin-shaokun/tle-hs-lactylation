# ============================================================
# 03_wgcna.R — 加权共表达网络（主分析样本：TLE-HS vs 对照皮层）
# 【2026-11-04 修正版】对接 02_deg_main.rds（原 02_deg_hippocampus.rds 已废弃）
# 输入：results/02_deg_main.rds（sub_expr 基因×样本，sub_pd 含 group/age/sex/region）
# 输出：results/03_wgcna_modules.csv, 03_module_trait_cor.csv；
#       figures/03_soft_threshold.png, 03_module_traits.png, 03_dendro.png
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(WGCNA); library(ggplot2); library(pheatmap)})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 3)   # 8GB 内存机器，线程留余量

d  <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
expr <- d$sub_expr; pd <- d$sub_pd
cat("WGCNA 输入：", nrow(expr), "基因 ×", ncol(expr), "样本\n")

## 低表达/低方差过滤：取方差前 5000 基因建网（控制计算量与内存）
v <- apply(expr, 1, var)
keep <- names(sort(v, decreasing = TRUE))[1:min(5000, length(v))]
datExpr <- t(expr[keep, ])   # 样本×基因
cat("建网基因数：", ncol(datExpr), "\n")

## 软阈值选择
powers <- THRESH$wgcna_soft_range
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6
cat("softPower =", softPower,
    "（scale-free R^2 =", round(sft$fitIndices$SFT.R.sq[
      match(softPower, sft$fitIndices$Power)], 3), "）\n")

## 软阈值曲线图
png(file.path(FIG_DIR, "03_soft_threshold.png"), 2000, 1000, res = PLOT_THEME$save_dpi)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     type = "h", xlab = "Soft threshold power", ylab = "Scale-free topology R^2",
     main = "Soft-thresholding power selection", cex.main = 1.1)
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k., type = "h",
     xlab = "Soft threshold power", ylab = "Mean connectivity", main = "Mean connectivity")
dev.off()

## 建网（分块，控制内存）
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                        minModuleSize = THRESH$wgcna_min_module,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        maxBlockSize = 5000,   # 5000 基因单块即可，省内存
                        verbose = 0)
moduleColors <- labels2colors(net$colors)
cat("模块数：", length(unique(moduleColors)), "\n")
print(table(moduleColors))

meigs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

## 模块-性状关联
traits <- data.frame(
  TLE    = as.numeric(pd$condition == "TLE-HS"),
  Age    = pd$age,
  SexM   = as.numeric(pd$sex == unique(pd$sex)[1]) * 0 + as.numeric(pd$sex != levels(pd$sex)[1]),  # 男=1
  RegionTemp = as.numeric(pd$region == "Temporal"),
  row.names = rownames(pd)
)
## sex 列是 factor，编码保存在 02 日志；此处按"第二水平=1"更稳妥：
traits$SexM <- as.numeric(pd$sex == levels(pd$sex)[2])

mm  <- cor(meigs, traits, use = "p")
mmp <- corPvalueStudent(mm, nrow(traits))

## 输出模块归属 + 乳酸化基因所在模块
gs <- read.csv(GENESET_CSV)
modules <- data.frame(gene = colnames(datExpr),
                      module = moduleColors,
                      in_lactylation_geneset = colnames(datExpr) %in% gs$gene_symbol)
write.csv(modules, file.path(RES_DIR, "03_wgcna_modules.csv"), row.names = FALSE)
write.csv(cbind(module = sub("^ME", "", rownames(mm)), mm, mmp),
          file.path(RES_DIR, "03_module_trait_cor.csv"), row.names = FALSE)

## 摘要打印：与 TLE 最相关的模块 & 乳酸化基因分布
ord <- order(-abs(mm[, "TLE"]))
cat("\n== 与 TLE 状态相关性前 5 的模块 ==\n")
for (i in ord[1:5]) {
  mod <- sub("^ME", "", rownames(mm)[i])
  n_in <- sum(modules$module == mod & modules$in_lactylation_geneset)
  cat(sprintf("  %s: r=%.3f, p=%.2e, 含乳酸化基因 %d 个\n",
              mod, mm[i, "TLE"], mmp[i, "TLE"], n_in))
}
cat("\n乳酸化基因模块分布：\n")
print(table(modules$module[modules$in_lactylation_geneset]))

## 模块-性状热图（数值 + 红正蓝负，供投稿）
lab_mm <- matrix(sprintf("%.2f\n(%.1e)", mm, mmp), nrow = nrow(mm),
                 dimnames = dimnames(mm))
rownames(mm) <- sub("^ME", "", rownames(mm)); rownames(lab_mm) <- rownames(mm)
colnames(mm) <- c("TLE-HS", "Age", "Sex (male)", "Temporal cortex")[seq_len(ncol(mm))]
colnames(lab_mm) <- colnames(mm)
pheatmap(mm, display_numbers = lab_mm, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101), fontsize_number = 7,
         border_color = "white",
         main = "Module-trait associations (Pearson r, P)",
         filename = file.path(FIG_DIR, "03_module_traits.png"),
         width = 7.2, height = 4.6)

saveRDS(list(modules = modules, meigs = meigs, softPower = softPower,
             traits = traits, mm = mm, mmp = mmp),
        file.path(RES_DIR, "03_wgcna.rds"))
message("03 完成")
