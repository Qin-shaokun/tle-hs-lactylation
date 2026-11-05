## 09b3_composition_hub.R
## 重做 Fig 7C（按样本构成，修掉均值分母不一致导致的 >100%）与 Fig 7D（hub 基因表达比例条形图）
source("config.R")
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(tidyr)})

obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
md  <- obj@meta.data

ct_order <- c("ExNeuron", "Astrocyte", "Microglia", "Oligodendrocyte",
              "OPC", "InNeuron", "Endothelial")

## ---- C：每个样本的构成（8 根柱，按 condition 分面）----
comp <- md %>%
  count(sample, condition, celltype, name = "n") %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = ct_order),
         sample   = factor(sample, levels = c(
           "GSM5724021_Normal1", "GSM5724022_Normal2", "GSM5724023_Normal_3",
           "GSM6893724_Normal_4", "GSM5724024_TLE_1", "GSM5724025_TLE_2",
           "GSM5724026_TLE_3", "GSM6893725_TLE_4")),
         condition = factor(condition, levels = c("Control", "TLE")))

## 用样本短名做 x 轴标签
comp$sample_short <- factor(sub("^GSM[0-9]+_", "", as.character(comp$sample)),
                            levels = c("Normal1", "Normal2", "Normal_3", "Normal_4",
                                       "TLE_1", "TLE_2", "TLE_3", "TLE_4"))

pal <- c(ExNeuron = "#D95F02", Astrocyte = "#1B9E77", Microglia = "#7570B3",
         Oligodendrocyte = "#E7298A", OPC = "#66A61E", InNeuron = "#E6AB02",
         Endothelial = "#666666")

pC <- ggplot(comp, aes(sample_short, frac, fill = celltype)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.25) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Nuclear fraction", fill = NULL, title = "Cell-type composition per specimen") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "09b_composition_bar.png"), pC, width = 8.4, height = 4.8, dpi = 200)
cat("C done\n")

## ---- D：hub 基因表达比例（合并两组，条形图 + 核数标注）----
hub <- read.csv(file.path(RES_DIR, "09c_hub_localization.csv"), stringsAsFactors = FALSE) %>%
  filter(condition == "All", gene %in% c("TNFRSF13C", "TMEM156")) %>%
  mutate(celltype = factor(celltype, levels = rev(ct_order)),
         pct      = pct_expr * 100,
         gene     = factor(gene, levels = c("TNFRSF13C", "TMEM156")))

pD <- ggplot(hub, aes(pct, celltype, fill = gene)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%% (n=%s)", pct, format(n_cells, big.mark = ","))),
            position = position_dodge(width = 0.75), hjust = -0.08, size = 2.9) +
  scale_fill_manual(values = c(TNFRSF13C = "#2C7FB8", TMEM156 = "#D95F02")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +
  labs(x = "Nuclei expressing the gene (%)", y = NULL, fill = NULL,
       title = "Hub genes by cell type (pooled)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
ggsave(file.path(FIG_DIR, "09c_hub_barplot.png"), pD, width = 7.6, height = 4.4, dpi = 200)
cat("D done\n")

## 顺便导出 C 的样本级构成表与组均值（供正文与图注核对）
wid <- comp %>% select(sample_short, condition, celltype, frac) %>%
  pivot_wider(names_from = celltype, values_from = frac, values_fill = 0)
write.csv(wid, file.path(RES_DIR, "09b_composition_by_sample.csv"), row.names = FALSE)
grp <- wid %>% group_by(condition) %>%
  summarise(across(all_of(ct_order), ~ round(mean(.x) * 100, 2)), .groups = "drop")
print(as.data.frame(grp))
write.csv(grp, file.path(RES_DIR, "09b_composition_group_mean.csv"), row.names = FALSE)
cat("ALL_DONE\n")
