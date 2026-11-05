# ============================================================
# run_all.R — 总控：按顺序执行，失败即停
# 注：07_single_cell.R 是早期占位方案（基于 GEO 处理矩阵），已被 09 系取代，不纳入流水线
# Usage: from the repository root, `Rscript run_all.R`
# ============================================================
source("config.R")

steps <- c("01_download_and_metadata.R", "02_deg_main.R", "02b_sensitivity.R",
           "03_wgcna.R", "04_ml_hub_selection.R", "04b_ml_venn_rf.R",
           "05_validation_roc.R", "05b_nomogram_firth.R", "05c_calibration.R",
           "06_immune_infiltration.R",
           "08_enrichment.R", "08b_kegg.R",
           "09_sc_qc.R", "09b_cluster.R", "09b3_composition_hub.R", "09c_localize.R",
           "09c_keygenes_dotplot.R", "09c2_markers_fast.R",
           "09d_assemble_figs.R", "09e_print_panels.R",
           "11_redraw_panels.R", "10_assemble_main.R")

## 依赖 R 包自检（首次运行前提示安装）
need <- c("GEOquery", "limma", "WGCNA", "glmnet", "e1071", "randomForest",
          "pROC", "rms", "logistf", "GSVA", "GSEABase", "Seurat", "dplyr",
          "clusterProfiler", "org.Hs.eg.db", "ggplot2", "plyr")
miss <- need[!need %in% rownames(installed.packages())]
if (length(miss)) {
  message("缺少以下包，请先安装：\n", paste(miss, collapse = ", "))
  message("可运行：install.packages(c('", paste(setdiff(miss, c('GEOquery','org.Hs.eg.db','clusterProfiler','GSVA','Seurat')),
                                                    collapse = "','"), "'))")
  message("Bioconductor 包：BiocManager::install(c('GEOquery','org.Hs.eg.db','clusterProfiler','GSVA'))")
  stop("安装后重跑")
}

for (s in steps) {
  message("\n====> 运行 ", s, " ...\n")
  source(file.path("R", s), echo = TRUE, max.deparse.length = 200)
}
message("\n全流程跑完。结果在 results/，图在 figures/。")
