# ============================================================
# 06_immune_infiltration.R — 免疫/细胞浸润分析（2026-11-04 对接 02_deg_main）
# 方法：ssGSEA（GSVA ≥2.0 API），脑组织细胞类型签名
# 输出：results/06_cell_fractions.csv, 06_cell_wilcox.csv,
#       06_corr_hub_cells.csv；figures/06_boxplot.png, 06_corr_hub.png
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(GSVA); library(GSEABase); library(ggplot2)})

d <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
expr <- d$sub_expr   # 基因×样本（TLE-HS 64 + 对照皮层 14）

## ---- 细胞类型签名（内置精简版，论文 Methods 注明来源）----
## 注意：签名中刻意排除 hub gene（TNFRSF13C），否则"hub 基因 vs 细胞评分"的相关为自相关
signatures <- list(
  Microglia    = c("CX3CR1", "P2RY12", "TMEM119", "AIF1", "C1QA", "C1QB", "CSF1R", "TYROBP"),
  Astrocyte    = c("GFAP", "AQP4", "ALDH1L1", "S100B", "SLC1A3", "SLC1A2", "AGT"),
  ExcitNeuron  = c("SLC17A7", "CAMK2A", "GRIN1", "SNAP25", "KCNB1"),
  InhibNeuron  = c("GAD1", "GAD2", "SLC6A1", "PVALB", "SST", "CALB1"),
  Oligodendro  = c("MBP", "PLP1", "MOG", "OLIG2", "MAG"),
  Endothelium  = c("PECAM1", "VWF", "CLDN5", "FLT1"),
  Tcell        = c("CD3D", "CD3E", "CD8A", "IL7R"),
  Bcell        = c("CD79A", "MS4A1", "BANK1", "CD19", "CD79B"),   # 已剔除 TNFRSF13C（hub gene）
  Monocyte     = c("LYZ", "S100A8", "S100A9", "FCGR3A")
)
signatures <- lapply(signatures, function(g) intersect(g, rownames(expr)))
cat("签名基因命中数：", sapply(signatures, length), "\n")
sig <- GeneSetCollection(mapply(function(n, g) GeneSet(g, setName = n),
                                names(signatures), signatures))

## ---- ssGSEA 评分 ----
gsva_par <- GSVA::ssgseaParam(exprData = as.matrix(expr), geneSets = sig)
score <- GSVA::gsva(gsva_par)   # 细胞类型 × 样本
frac <- as.data.frame(t(score))
frac$condition <- d$sub_pd$condition
frac$age       <- d$sub_pd$age
frac$sex       <- d$sub_pd$sex
frac$sample    <- rownames(frac)
write.csv(frac, file.path(RES_DIR, "06_cell_fractions.csv"), row.names = FALSE)

## ---- 组间 Wilcoxon 检验（含 age/sex 协变量的稳健性检查：残差法）----
cell_names <- names(signatures)
wx <- do.call(rbind, lapply(cell_names, function(ct) {
  x <- frac[[ct]]
  p_raw <- wilcox.test(x ~ frac$condition)$p.value
  # age/sex 残差化后再比（控制协变量）
  fit <- lm(x ~ age + sex, data = frac)
  res <- residuals(fit)
  p_adj <- wilcox.test(res ~ frac$condition)$p.value
  data.frame(celltype = ct,
             median_TLE = median(x[frac$condition == "TLE-HS"]),
             median_Ctrl = median(x[frac$condition == "ControlCortex"]),
             log2FC_med = log2(median(x[frac$condition == "TLE-HS"]) -
                               median(x[frac$condition == "ControlCortex"]) + 1),
             p_wilcox = p_raw, p_residual = p_adj)
}))
wx$p_bh <- p.adjust(wx$p_wilcox, "BH")
write.csv(wx, file.path(RES_DIR, "06_cell_wilcox.csv"), row.names = FALSE)
print(wx)

## ---- 组间箱线图 ----
long <- reshape(frac, direction = "long", varying = cell_names,
                v.names = "score", times = cell_names, idvar = "sample")
long$condition <- factor(long$condition, levels = c("ControlCortex", "TLE-HS"))
## 细胞类型缩写 → 标准英文名（图表统一）
cname_full <- c(Microglia = "Microglia", Astrocyte = "Astrocyte",
                ExcitNeuron = "Excitatory neuron", InhibNeuron = "Inhibitory neuron",
                Oligodendro = "Oligodendrocyte", Endothelium = "Endothelium",
                Tcell = "T cell", Bcell = "B cell", Monocyte = "Monocyte")
long$time <- factor(cname_full[as.character(long$time)],
                    levels = unname(cname_full[cell_names]))
## BH 后显著性星号（按 celltype 对应，名称同步映射）
sig_lab <- ifelse(wx$p_bh < 0.001, "***",
                  ifelse(wx$p_bh < 0.01, "**",
                         ifelse(wx$p_bh < 0.05, "*", "ns")))
names(sig_lab) <- unname(cname_full[as.character(wx$celltype)])
star_df <- data.frame(time = names(sig_lab), star = unname(sig_lab),
                      stringsAsFactors = FALSE)
star_df$y <- tapply(long$score, long$time, max)[star_df$time]
p <- ggplot(long, aes(condition, score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  facet_wrap(~time, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(ControlCortex = "#7FA6C9", `TLE-HS` = "#C94F4F"),
                    labels = c("Control cortex", "TLE-HS")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  geom_text(data = star_df, aes(x = 1.5, y = y, label = star),
            inherit.aes = FALSE, vjust = -0.35, size = 3.8) +
  scale_x_discrete(labels = c("Ctrl", "TLE-HS")) +
  labs(x = "", y = "ssGSEA score",
       title = "Cell-type infiltration: TLE-HS vs control cortex") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIG_DIR, "06_boxplot.png"), p, width = 10, height = 7, dpi = PLOT_THEME$save_dpi)

## ---- hub gene × 细胞评分相关 + 乳酸化基因集成员相关 ----
ml  <- readRDS(file.path(RES_DIR, "04_ml.rds"))
gs  <- read.csv(GENESET_CSV)
focus <- union(intersect(ml$hub, rownames(expr)),
               intersect(gs$gene_symbol, rownames(expr)))
cm <- cor(t(expr[focus, ]), t(score), method = "spearman")
write.csv(cm, file.path(RES_DIR, "06_corr_hub_cells.csv"))
long2 <- as.data.frame(as.table(cm)); names(long2) <- c("gene", "celltype", "rho")
## 缩写 → 标准英文名
lvl_map <- c(Microglia = "Microglia", Astrocyte = "Astrocyte",
             ExcitNeuron = "Excitatory neuron", InhibNeuron = "Inhibitory neuron",
             Oligodendro = "Oligodendrocyte", Endothelium = "Endothelium",
             Tcell = "T cell", Bcell = "B cell", Monocyte = "Monocyte")
new_lvl <- lvl_map[as.character(cell_names)]
long2$celltype <- factor(lvl_map[as.character(long2$celltype)], levels = unname(new_lvl))
long2$is_hub <- long2$gene %in% ml$hub
p2 <- ggplot(long2, aes(celltype, gene, fill = rho)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal(base_size = 10) +
  labs(title = "Spearman correlation: hub and lactylation-related genes vs cell-type scores") +
  theme(axis.text.y = element_text(face = ifelse(levels(long2$gene) %in% ml$hub,
                                                 "bold", "plain")),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "06_corr_hub.png"), p2, width = 7.5, height = 5.5,
       dpi = PLOT_THEME$save_dpi)
message("06 完成")
