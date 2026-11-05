# ============================================================
# 09_sc_qc.R — GSE190452 单核测序（snRNA-seq）质控
# 数据：data/GSE190452_RAW.tar（8 样本 × 10x 三件套，4 TLE vs 4 对照颞叶皮层）
# 输出：
#   results/09_qc_per_sample.csv        每样本 QC 指标（过滤前后）
#   results/09_qc_cell_counts.csv       每样本过滤前后核数
#   results/09_sc_filtered.rds          过滤后合并对象（供 09b 聚类注释）
#   figures/09_qc_violin.png            过滤前 nFeature/nCount/mt 小提琴
#   figures/09_qc_after_violin.png      过滤后小提琴
#   figures/09_qc_scatter.png           nCount vs nFeature（按样本，标阈值线）
#   figures/09_mt_hist.png              percent.mt 分布直方图
# 日志：log_09qc.txt
# ============================================================
source("config.R")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(dplyr)
})
options(width = 120)

TAR     <- file.path(DATA_DIR, "GSE190452_RAW.tar")
RAWDIR  <- file.path(DATA_DIR, "GSE190452_raw")
dir.create(RAWDIR, showWarnings = FALSE, recursive = TRUE)

## ---- 0. 解包（幂等）----
if (length(list.files(RAWDIR, pattern = "\\.mtx\\.gz$")) == 0) {
  stopifnot(file.exists(TAR))
  utils::untar(TAR, exdir = RAWDIR)
}
allfiles <- list.files(RAWDIR, recursive = TRUE, full.names = TRUE)
mtx_f  <- allfiles[grepl("matrix\\.mtx\\.gz$", allfiles)]
cat("发现 matrix 文件:", length(mtx_f), "\n")
stopifnot(length(mtx_f) >= 8)

## 样本 id：支持两种布局
##  a) 平铺：GSMxxxx_title_matrix.mtx.gz        -> id = GSMxxxx_title
##  b) 目录：GSMxxxx_title/matrix.mtx.gz        -> id = GSMxxxx_title（目录名）
sample_of <- function(p) {
  b <- basename(p)
  if (grepl("_?matrix\\.mtx\\.gz$", b)) {
    sub("_?matrix\\.mtx\\.gz$", "", b)
  } else {
    basename(dirname(p))
  }
}

## ---- 1. 手工读取 10x 三件套 ----
read_one <- function(mtx_path) {
  s   <- sample_of(mtx_path)
  d   <- dirname(mtx_path)
  flat <- grepl("_?matrix\\.mtx\\.gz$", basename(mtx_path))   # 平铺布局
  bc  <- allfiles[grepl("barcodes\\.tsv\\.gz$", basename(allfiles)) &
                    (if (flat) grepl(paste0("^", s, "_?barcodes"), basename(allfiles)) else dirname(allfiles) == d)]
  ft  <- allfiles[grepl("features\\.tsv\\.gz$", basename(allfiles)) &
                    (if (flat) grepl(paste0("^", s, "_?features"), basename(allfiles)) else dirname(allfiles) == d)]
  stopifnot(length(bc) == 1, length(ft) == 1)
  bcv <- readLines(gzfile(bc)); bcv <- sub(" .*", "", bcv)
  ftt <- read.delim(gzfile(ft), header = FALSE, stringsAsFactors = FALSE)
  feat_name <- if (ncol(ftt) >= 2) ftt[[2]] else ftt[[1]]
  feat_name <- make.unique(as.character(feat_name), sep = "_")   # GEO features 有 10 个重复符号
  feat_type <- if (ncol(ftt) >= 3) ftt[[3]] else rep("Gene Expression", nrow(ftt))
  m <- readMM(gzfile(mtx_path))
  m <- as(m, "dgCMatrix")
  # 10x matrix 是 genes x cells
  if (nrow(m) != length(feat_name)) {
    stop(sprintf("维度不匹配: %s matrix %d 行 vs features %d", s, nrow(m), length(feat_name)))
  }
  rownames(m) <- feat_name; colnames(m) <- bcv
  # 只保留 Gene Expression（snRNA 包含 Antibody Capture 时）
  keep <- feat_type == "Gene Expression"
  list(mat = m[keep, , drop = FALSE], sample = s)
}

objs <- list()
for (f in mtx_f) {
  r <- read_one(f)
  s <- r$sample
  cat("读取:", s, " 基因:", nrow(r$mat), " 核:", ncol(r$mat), "\n")
  o <- CreateSeuratObject(r$mat, project = s, min.cells = 3, min.features = 100)
  o$sample <- s
  objs[[s]] <- o
}

# condition 映射：样本名含 GSM；TLE vs 对照来自 GEO metadata（4v4）
# 从 GEO series matrix 自动取 condition，取不到则按样本序（前 4 TLE 后 4 control 需人工确认）
meta_hint <- NULL
smx <- file.path(DATA_DIR, "GSE190452_series_matrix.txt.gz")
if (file.exists(smx)) {
  ln <- readLines(gzfile(smx), n = 3000)
  ch <- ln[grepl("!Sample_characteristics_ch1", ln)]
  if (length(ch) >= 1) meta_hint <- ch
  cat("series_matrix characteristics:\n"); print(ch)
}
# 人工映射（依据 GSE190452 描述：4 TLE vs 4 control；以 series matrix 输出为准修正）
cond_map <- cfg_sc_condition_map   # config.R：GSM号 -> "TLE"/"Control"
batch_map <- cfg_sc_batch_map
for (s in names(objs)) {
  gsm <- sub("^(GSM\\d+).*", "\\1", s)
  objs[[s]]$condition <- if (!is.null(cond_map[[gsm]])) cond_map[[gsm]] else NA
  objs[[s]]$batch     <- if (!is.null(batch_map[[gsm]])) batch_map[[gsm]] else NA
}
cat("condition 缺失样本数:", sum(is.na(sapply(objs, function(o) o$condition[1]))), "\n")

merged <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
merged <- PercentageFeatureSet(merged, pattern = "^MT-", col.name = "percent.mt")
merged <- PercentageFeatureSet(merged, pattern = "^RP[SL]-|^RPL|^RPS", col.name = "percent.rb")

## ---- 2. 过滤前 QC 表 ----
qc_tbl <- merged@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells_raw   = n(),
            median_nFeature = round(median(nFeature_RNA)),
            median_nCount   = round(median(nCount_RNA)),
            median_pct_mt   = round(median(percent.mt), 2),
            .groups = "drop")
print(qc_tbl)

p1 <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0, group.by = "sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(FIG_DIR, "09_qc_violin.png"), p1, width = 12, height = 4.5, dpi = 200)

p2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "sample", pt.size = 0.3) +
  geom_hline(yintercept = c(200, 7000), linetype = 2, col = "red2") +
  theme(axis.text.x = element_text(size = 8))
ggsave(file.path(FIG_DIR, "09_qc_scatter.png"), p2, width = 10, height = 6, dpi = 200)

p3 <- ggplot(merged@meta.data, aes(percent.mt, fill = condition)) +
  geom_histogram(bins = 80, position = "identity", alpha = 0.5) +
  geom_vline(xintercept = 10, linetype = 2, col = "red2") +
  labs(x = "percent mitochondrial (raw)", y = "nuclei")
ggsave(file.path(FIG_DIR, "09_mt_hist.png"), p3, width = 6, height = 4, dpi = 200)

## ---- 3. snRNA 常规过滤 ----
# nFeature 200–7000（snRNA 核范围）；percent.mt < 10；nCount > 500
merged <- subset(merged,
                 subset = nFeature_RNA >= 200 & nFeature_RNA <= 7000 &
                          percent.mt < 10 &
                          nCount_RNA  >= 500)

qc_tbl2 <- merged@meta.data %>%
  group_by(sample, condition) %>%
  summarise(n_cells_kept  = n(),
            median_nFeature = round(median(nFeature_RNA)),
            median_nCount   = round(median(nCount_RNA)),
            median_pct_mt   = round(median(percent.mt), 2),
            .groups = "drop")
cnt <- full_join(qc_tbl, qc_tbl2, by = c("sample", "condition")) %>%
  mutate(retained_pct = round(100 * n_cells_kept / n_cells_raw, 1))
write.csv(cnt, file.path(RES_DIR, "09_qc_cell_counts.csv"), row.names = FALSE)
cat("\n== 过滤前后核数 ==\n"); print(cnt)

p4 <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0, group.by = "sample", split.by = "condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(FIG_DIR, "09_qc_after_violin.png"), p4, width = 12, height = 4.5, dpi = 200)

saveRDS(merged, file.path(RES_DIR, "09_sc_filtered.rds"))
cat("过滤后总核数:", ncol(merged), " 已保存 09_sc_filtered.rds\n")
cat("QC 完成。下一步：09b 聚类 + 细胞类型注释 + hub/乳酸化基因定位。\n")
