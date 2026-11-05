# ============================================================
# 01_download_and_metadata.R — 下载数据 + 生成待人工核对的 metadata 表
# 输入：config.R
# 输出：data/GSE256068_expr.csv.gz, data/GSE256068_metadata_raw.csv
#       results/01_metadata_checklist.md  ← 人工核对清单
# ============================================================
source("config.R")
suppressPackageStartupMessages({library(GEOquery)})

dl <- function(acc, dest) {
  f <- file.path(DATA_DIR, paste0(acc, "_series_matrix.txt.gz"))
  if (!file.exists(f)) {
    message("下载 ", acc, " ...")
    g <- getGEO(acc, destdir = DATA_DIR, GSEMatrix = TRUE)
    saveRDS(g, file.path(DATA_DIR, paste0(acc, "_eset.rds")))
  } else {
    g <- readRDS(file.path(DATA_DIR, paste0(acc, "_eset.rds")))
  }
  g
}

eset <- dl(GSE256068$accession)
expr <- exprs(eset[[1]])
pdata <- pData(eset[[1]])

## 原始 metadata 全量导出（列名与取值必须人工过目！）
write.csv(pdata, file.path(DATA_DIR, "GSE256068_metadata_raw.csv"), row.names = FALSE)

## 探查关键列：把候选列打印出来供人工确认
cat("=== 可用的 sample characteristics 列 ===\n")
cols <- grep("characteristics|title|supplementary", colnames(pdata), value = TRUE, ignore.case = TRUE)
for (cl in cols) {
  cat("\n---", cl, "---\n")
  print(table(pdata[[cl]], useNA = "ifany"))
}

## 生成人工核对清单
checklist <- c(
  "# GSE256068 metadata 人工核对清单",
  "",
  paste0("- 总样本数：", nrow(pdata), "；基因/探针数：", nrow(expr)),
  paste0("- 平台：", unique(pdata$platform)), "",
  "## 需要人工确认的事项",
  "1. 哪一列是组织部位（hippocampus / cortex）？取值拼写是什么？",
  "2. 哪一列是病理分组（TLE-HS / FCD2b / TSC / Control）？对照的来源（尸检 or 手术）？",
  "3. 年龄、性别各在哪一列？单位？缺失比例多少？",
  "4. 对照海马到底 7 例还是 14 例？（文献说法不一，以本表为准）",
  "5. 确认后把结论写回 config.R 的列名，并在 02 脚本顶部的 group_map 中固化样本→分组映射。",
  "",
  "## 分组频数预览（原始 characteristics 列逐列）",
  capture.output(print(cols))
)
writeLines(checklist, file.path(RES_DIR, "01_metadata_checklist.md"))

saveRDS(list(expr = expr, pdata = pdata), file.path(DATA_DIR, "GSE256068_main.rds"))
message("完成。请打开 results/01_metadata_checklist.md 进行人工核对。")
