# 09b2_composition.R — 细胞类型构成堆叠条形图（Fig 7C）
source("config.R")
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
obj <- readRDS(file.path(RES_DIR, "09b_sc_annotated.rds"))
pt <- obj@meta.data %>%
  count(sample, condition, celltype, name = "n") %>%
  group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
  group_by(condition, celltype) %>% summarise(mean_frac = mean(frac), .groups = "drop")
ct_order <- pt %>% group_by(celltype) %>% summarise(s = sum(mean_frac)) %>%
  arrange(desc(s)) %>% pull(celltype)
pt$celltype <- factor(pt$celltype, levels = ct_order)
p <- ggplot(pt, aes(condition, mean_frac, fill = celltype)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Mean nuclear fraction", fill = NULL,
       title = "Cell-type composition, snRNA-seq (GSE190452)") +
  theme_classic(base_size = 13) + theme(legend.position = "right")
ggsave(file.path(FIG_DIR, "09b_composition_bar.png"), p, width = 7.5, height = 5, dpi = 200)
cat("done\n")
