# ============================================================
# config.R — global configuration: paths, dataset parameters, thresholds
# Run everything from the repository root (the directory containing this file).
# 项目：癫痫后认知障碍 × 蛋白乳酸化（纯生信）
# 修改本文件即可切换数据集/阈值，脚本本身不动
# ============================================================

## ---- 路径（相对项目根 03_分析/）----
ROOT      <- normalizePath(".", winslash = "/")
DATA_DIR  <- file.path(ROOT, "data")
RES_DIR   <- file.path(ROOT, "results")
FIG_DIR   <- file.path(ROOT, "figures")
DOC_DIR   <- file.path(ROOT, "docs")
for (d in c(DATA_DIR, RES_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 基因集 ----
GENESET_CSV <- file.path(DATA_DIR, "lactylation_geneset_v1.csv")
GENESET_VERSION <- "v1"   # 每次修订记得升版本号

## ---- 数据集定义 ----
# 主训练集：GSE256068（Nat Commun 2024）
# 【2026-11-03 metadata 逐例核实后的修正】
#   - 海马样本 13 例全部是对照（ControlHippocampus），**无 TLE 海马样本**
#   - TLE-HS 64 例全部为颞叶新皮层（tissue = Temporal）
#   - 对照皮层 14 例：额叶 8、颞叶 4、顶叶 2；其中 10 例 <18 岁
#   → 决策：主分析改为「TLE-HS 颞叶 vs 对照皮层」，年龄/性别作协变量，
#           另做成人子集（≥18 岁）与 region 校正两层敏感性分析
GSE256068 <- list(
  accession    = "GSE256068",
  platform     = "GPL24676",
  expr_source  = "processed_csv",   # 表达矩阵来自 suppl/processed_data.csv.gz
  group_col    = "disease state",   # metadata 实际列名（已核实）
  tissue_col   = "tissue",
  age_col      = "age",
  sex_col      = "Sex",
  tissue_keep  = NULL,              # NULL = 不按组织过滤（TLE-HS 全是 Temporal）
  groups_keep  = c("TLE-HS", "ControlCortex"),
  case_label   = "TLE-HS",
  ctrl_label   = "ControlCortex"
)

# 亚组/验证集
GSE71058  <- list(accession = "GSE71058",  platform = "GPL11154")  # HS vs noHS（LCM 齿状回颗粒细胞）
GSE205661 <- list(accession = "GSE205661", platform = NA)           # 外部验证 1
GSE190451 <- list(accession = "GSE190451", platform = "GPL24676")   # 外部验证 2

# 单细胞（皮层，snRNA-seq 4v4）
# 【2026-11-04 series matrix 核实】样本-条件映射与批次：
#   主批次（2021 提交）：GSM5724021-23 = Normal 1-3；GSM5724024-26 = TLE 1-3
#   补充批次（2022 补录）：GSM6893724 = Normal 4；GSM6893725 = TLE 4
#   → 共 4v4，但分两个批次；QC 与 09b 聚类时需检查批次效应
GSE190452 <- list(accession = "GSE190452", tissue = "temporal_cortex")
cfg_sc_condition_map <- list(
  GSM5724021 = "Control", GSM5724022 = "Control", GSM5724023 = "Control",
  GSM6893724 = "Control",
  GSM5724024 = "TLE",     GSM5724025 = "TLE",     GSM5724026 = "TLE",
  GSM6893725 = "TLE"
)
cfg_sc_batch_map <- list(
  GSM5724021 = "batch1", GSM5724022 = "batch1", GSM5724023 = "batch1",
  GSM5724024 = "batch1", GSM5724025 = "batch1", GSM5724026 = "batch1",
  GSM6893724 = "batch2", GSM6893725 = "batch2"
)

## ---- 统计阈值 ----
THRESH <- list(
  lfc        = 1,       # |log2FC|
  fdr        = 0.05,    # BH 校正
  wgcna_soft_range  = c(1:10, seq(12, 20, by = 2)),   # 1-10 全扫 + 12/14/16/18/20
  wgcna_min_module  = 30,
  lasso_lambda     = "lambda.min",
  rf_ntree         = 1000,
  roc_ci_method    = "delong"
)

## ---- 设计矩阵协变量 ----
COVARIATES <- c("age", "sex")   # 主分析设计：~ age + sex + condition

## ---- 通用绘图参数 ----
PLOT_THEME <- list(base_size = 12, save_dpi = 300)
