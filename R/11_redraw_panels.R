## 11_redraw_panels.R
## 按印刷尺寸重画需要修正的面板（不重跑分析，全部从已存结果读取）
source("config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap); library(WGCNA); library(dplyr); library(tidyr)
})
OUT <- file.path(FIG_DIR, "print"); dir.create(OUT, showWarnings = FALSE)
BASE <- 8

## ============================================================
## 1) Fig 1B 热图：36 行（top 30 DEG + 6 个显著乳酸化基因），行名可读
## ============================================================
d    <- readRDS(file.path(RES_DIR, "02_deg_main.rds"))
deg  <- d$deg; expr <- d$sub_expr; pd <- d$sub_pd
gs   <- read.csv(GENESET_CSV)
lac6 <- c("GPD1", "HK2", "NLRP3", "PFKFB3", "BDNF", "BRDT")

top30 <- head(intersect(deg$gene, rownames(expr)), 30)
plot_genes <- unique(c(top30, intersect(lac6, rownames(expr))))
cat("热图行数:", length(plot_genes), "\n")

m2 <- t(scale(t(expr[plot_genes, ]))); m2[is.na(m2)] <- 0
ann <- data.frame(Group = pd$condition, Region = pd$region,
                  row.names = colnames(expr))
pheatmap(m2, annotation_col = ann, annotation_legend = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 4.6, fontsize_col = 6, fontsize = 6,
         treeheight_row = 12, treeheight_col = 10,
         border_color = NA,
         filename = file.path(OUT, "p1B.png"), width = 3.05, height = 3.55)
cat("p1B done\n")

## 完整 97 基因热图 → 补充图 S8
all_genes <- intersect(union(head(deg$gene, 50), intersect(gs$gene_symbol, rownames(expr))),
                       rownames(expr))
m3 <- t(scale(t(expr[all_genes, ]))); m3[is.na(m3)] <- 0
pheatmap(m3, annotation_col = ann, annotation_legend = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 5, fontsize_col = 6,
         treeheight_row = 14, treeheight_col = 10, border_color = NA,
         filename = file.path(OUT, "pS8.png"), width = 5.2, height = 8.2)
cat("pS8 done:", length(all_genes), "genes\n")

## ============================================================
## 2) Fig 2A 软阈值图（重算 sft，缩短标题）
## ============================================================
v <- apply(expr, 1, var)
keep <- names(sort(v, decreasing = TRUE))[1:min(5000, length(v))]
datExpr <- t(expr[keep, ])
sft <- pickSoftThreshold(datExpr, powerVector = THRESH$wgcna_soft_range, verbose = 0)
cat("softPower 复算 =", sft$powerEstimate, "\n")

png(file.path(OUT, "p2A.png"), width = 3.35, height = 1.75, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(3.4, 3.4, 1.5, 0.5), mgp = c(1.9, 0.5, 0), cex.axis = 0.58,
    cex.lab = 0.62, cex.main = 0.66, tcl = -0.25)
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     type = "h", xlab = "Soft threshold power", ylab = expression("Scale-free " * R^2),
     main = "Scale-free fit")
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k., type = "h",
     xlab = "Soft threshold power", ylab = "Mean connectivity",
     main = "Connectivity")
dev.off()
cat("p2A done\n")

## ============================================================
## 3) Fig 5A 细胞浸润箱线图（5 列 × 2 行，重新排版）
## ============================================================
frac <- read.csv(file.path(RES_DIR, "06_cell_fractions.csv"), stringsAsFactors = FALSE)
wx   <- read.csv(file.path(RES_DIR, "06_cell_wilcox.csv"), stringsAsFactors = FALSE)
cell_names <- c("Microglia", "Astrocyte", "ExcitNeuron", "InhibNeuron", "Oligodendro",
                "Endothelium", "Tcell", "Bcell", "Monocyte")
cname_full <- c(Microglia = "Microglia", Astrocyte = "Astrocyte",
                ExcitNeuron = "Excitatory neuron", InhibNeuron = "Inhibitory neuron",
                Oligodendro = "Oligodendrocyte", Endothelium = "Endothelium",
                Tcell = "T cell", Bcell = "B cell", Monocyte = "Monocyte")

long <- frac %>%
  select(all_of(cell_names), condition, sample) %>%
  pivot_longer(-c(condition, sample), names_to = "celltype", values_to = "score") %>%
  mutate(celltype = factor(cname_full[celltype], levels = unname(cname_full[cell_names])),
         condition = factor(condition, levels = c("ControlCortex", "TLE-HS")))

sig_lab <- ifelse(wx$p_bh < 0.001, "***",
                  ifelse(wx$p_bh < 0.01, "**", ifelse(wx$p_bh < 0.05, "*", "ns")))
names(sig_lab) <- cname_full[as.character(wx$celltype)]
star_df <- data.frame(celltype = factor(names(sig_lab), levels = unname(cname_full[cell_names])),
                      star = unname(sig_lab))
star_df$y <- tapply(long$score, long$celltype, max)[as.character(star_df$celltype)]

p5A <- ggplot(long, aes(condition, score, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.32, width = 0.62) +
  facet_wrap(~celltype, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(ControlCortex = "#7FA6C9", `TLE-HS` = "#C94F4F"),
                    labels = c("Control", "TLE-HS")) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  geom_text(data = star_df, aes(x = 1.5, y = y, label = star),
            inherit.aes = FALSE, vjust = -0.4, size = BASE * 0.32) +
  scale_x_discrete(labels = c("Ctrl", "TLE-HS")) +
  labs(x = NULL, y = "ssGSEA score", fill = NULL) +
  theme_bw(base_size = BASE) +
  theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = BASE - 1, face = "bold"),
        axis.text.x = element_text(size = BASE - 2, angle = 30, hjust = 1),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "p5A.png"), p5A, width = 6.45, height = 2.55, dpi = 300)
cat("p5A done\n")

## ============================================================
## 4) Fig 5B 相关性热图（精选 14 个基因）
## ============================================================
cm <- as.matrix(read.csv(file.path(RES_DIR, "06_corr_hub_cells.csv"),
                         row.names = 1, check.names = FALSE))
sel <- c("TNFRSF13C", "TMEM156", "HK2", "GPD1", "PFKFB3", "NLRP3",
         "BDNF", "BRDT", "LDHA", "SLC16A7", "SLC16A3", "DLG4", "EP300", "HDAC1")
sel <- intersect(sel, rownames(cm))
cm2 <- cm[sel, , drop = FALSE]
cn  <- colnames(cm2)
colnames(cm2) <- cname_full[cn]
hub <- c("TNFRSF13C", "TMEM156")

long2 <- as.data.frame(as.table(cm2)); names(long2) <- c("gene", "celltype", "rho")
long2$celltype <- factor(long2$celltype, levels = unname(cname_full[cn]))
long2$gene <- factor(long2$gene, levels = rev(sel))

p5B <- ggplot(long2, aes(celltype, gene, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho), colour = abs(rho) > 0.55),
            size = BASE * 0.30, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = BASE) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE - 1.5),
        axis.text.y = element_text(size = BASE - 1,
                                   face = ifelse(levels(long2$gene) %in% hub, "bold", "plain")),
        panel.grid = element_blank(),
        legend.key.size = unit(2.8, "mm"),
        legend.title = element_text(size = BASE - 1.5),
        legend.text = element_text(size = BASE - 2))
ggsave(file.path(OUT, "p5B.png"), p5B, width = 5.0, height = 2.95, dpi = 300)
cat("p5B done\n")

cat("REDRAW_DONE\n")
