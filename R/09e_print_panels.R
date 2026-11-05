## 09e_print_panels.R
## 按投稿规格重画 Fig 7 的五个面板并拼版
## 目标：成品宽 180 mm（7.09 in）@ 300 dpi = 2126 px，字号在印刷尺寸下可直接阅读
source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(png); library(grid)
})

obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
Idents(obj) <- "celltype"
OUT <- file.path(FIG_DIR, "print")
dir.create(OUT, showWarnings = FALSE)

BASE <- 7          # 面板基准字号
DPI  <- 300

## ---------- A / B：UMAP ----------
dim_theme <- theme_classic(base_size = BASE) +
  theme(legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = BASE - 0.5),
        legend.title = element_blank(),
        plot.title = element_text(size = BASE + 0.5, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = BASE - 1))

pA <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE,
              label.size = 2.2, repel = TRUE, pt.size = 0.02) +
  ggtitle("Cell types") + dim_theme + theme(legend.position = "right")
ggsave(file.path(OUT, "pA.png"), pA, width = 3.3, height = 2.65, dpi = DPI)

pB <- DimPlot(obj, reduction = "umap", group.by = "condition", pt.size = 0.02,
              cols = c(Control = "#E64B35", TLE = "#3C5488")) +
  ggtitle("Condition") + dim_theme + theme(legend.position = "right")
ggsave(file.path(OUT, "pB.png"), pB, width = 3.3, height = 2.65, dpi = DPI)

## ---------- C：每样本构成 ----------
comp <- read.csv(file.path(RES_DIR, "09b_composition_by_sample.csv"),
                 stringsAsFactors = FALSE) %>%
  pivot_longer(-c(sample_short, condition), names_to = "celltype", values_to = "frac")

ct_order <- c("ExNeuron", "Astrocyte", "Microglia", "Oligodendrocyte",
              "OPC", "InNeuron", "Endothelial")
comp$celltype <- factor(comp$celltype, levels = ct_order)
comp$condition <- factor(comp$condition, levels = c("Control", "TLE"))
comp$sample_short <- factor(comp$sample_short,
                            levels = c("Normal1", "Normal2", "Normal_3", "Normal_4",
                                       "TLE_1", "TLE_2", "TLE_3", "TLE_4"))
pal <- c(ExNeuron = "#D95F02", Astrocyte = "#1B9E77", Microglia = "#7570B3",
         Oligodendrocyte = "#E7298A", OPC = "#66A61E", InNeuron = "#E6AB02",
         Endothelial = "#666666")

pC <- ggplot(comp, aes(sample_short, frac, fill = celltype)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Nuclear fraction", fill = NULL) +
  theme_classic(base_size = BASE) +
  theme(legend.key.size = unit(2.4, "mm"),
        legend.text = element_text(size = BASE - 1),
        legend.title = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = BASE - 0.5, face = "bold"),
        axis.text.x = element_text(size = BASE - 1.8, angle = 40, hjust = 1),
        axis.text.y = element_text(size = BASE - 1))
ggsave(file.path(OUT, "pC.png"), pC, width = 3.3, height = 1.95, dpi = DPI)

## ---------- D：hub 基因表达比例 ----------
hub <- read.csv(file.path(RES_DIR, "09c_hub_localization.csv"), stringsAsFactors = FALSE) %>%
  filter(condition == "All", gene %in% c("TNFRSF13C", "TMEM156")) %>%
  mutate(celltype = factor(celltype, levels = rev(ct_order)),
         pct = pct_expr * 100,
         gene = factor(gene, levels = c("TNFRSF13C", "TMEM156")))

pD <- ggplot(hub, aes(pct, celltype, fill = gene)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
            position = position_dodge(width = 0.78), hjust = -0.15,
            size = BASE * 0.28) +
  scale_fill_manual(values = c(TNFRSF13C = "#2C7FB8", TMEM156 = "#D95F02")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Nuclei expressing the gene (%)", y = NULL, fill = NULL) +
  theme_classic(base_size = BASE) +
  theme(legend.position = "top",
        legend.key.size = unit(2.4, "mm"),
        legend.text = element_text(size = BASE - 0.5),
        legend.title = element_blank(),
        axis.text = element_text(size = BASE - 1),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.2))
ggsave(file.path(OUT, "pD.png"), pD, width = 3.3, height = 1.95, dpi = DPI)

## ---------- E：关键基因 dotplot ----------
genes <- c("TNFRSF13C", "TMEM156", "HK2", "PFKFB3", "NLRP3", "SLC16A7", "EP300", "BDNF")
pE <- DotPlot(obj, features = genes, group.by = "celltype",
              split.by = "condition", dot.scale = 4.2) +
  RotatedAxis() +
  labs(x = NULL, y = "Cell type / condition") +
  theme_classic(base_size = BASE) +
  theme(axis.text.x = element_text(size = BASE - 1.2, angle = 45, hjust = 1),
        axis.text.y = element_text(size = BASE - 1),
        legend.key.size = unit(2.4, "mm"),
        legend.text = element_text(size = BASE - 1),
        legend.title = element_text(size = BASE - 1))
ggsave(file.path(OUT, "pE.png"), pE, width = 4.5, height = 2.95, dpi = DPI)

## ---------- 拼版（英寸坐标） ----------
rd <- function(f) readPNG(file.path(OUT, f))
put_in <- function(img, x, y, w, label = NULL, lab_cex = 11) {
  ih <- dim(img)[1]; iw <- dim(img)[2]
  h <- w * ih / iw
  grid.raster(img, x = unit(x, "in"), y = unit(y, "in"),
              width = unit(w, "in"), height = unit(h, "in"),
              just = c("left", "bottom"), interpolate = TRUE)
  if (!is.null(label))
    grid.text(label, x = unit(x - 0.06, "in"), y = unit(y + h + 0.02, "in"),
              just = c("right", "top"),
              gp = gpar(fontface = "bold", fontsize = lab_cex))
  h
}

A <- rd("pA.png"); B <- rd("pB.png"); C <- rd("pC.png"); D <- rd("pD.png"); E <- rd("pE.png")

W <- 7.09; H <- 8.30
png(file.path(FIG_DIR, "assembled", "Fig7_print.png"),
    width = round(W * DPI), height = round(H * DPI), res = DPI, bg = "white")
grid.newpage()

gap <- 0.16
w2  <- 3.12
xL  <- 0.30; xR <- xL + w2 + 0.20

put_in(A, xL, H - 0.14 - 2.65, w2, "A")
put_in(B, xR, H - 0.14 - 2.65, w2, "B")

y2 <- H - 0.14 - 2.65 - gap
put_in(C, xL, y2 - 1.95, w2, "C")
put_in(D, xR, y2 - 1.95, w2, "D")

y3 <- y2 - 1.95 - 0.24
put_in(E, (W - 4.4) / 2, y3 - 2.90, 4.4, "E")

dev.off()
cat("PRINT_FIG7_DONE\n")
