## 生成 docs/session_info.txt（R 版本与依赖包版本）
pk <- c("GEOquery", "limma", "WGCNA", "glmnet", "e1071", "randomForest",
        "pROC", "rms", "logistf", "GSVA", "GSEABase", "clusterProfiler",
        "org.Hs.eg.db", "Seurat", "SeuratObject", "ggplot2", "dplyr",
        "plyr", "Matrix", "patchwork", "png", "grid")

con <- file("docs/session_info.txt", "w", encoding = "UTF-8")
writeLines(c(
  paste("R version:", R.version.string),
  paste("Platform :", R.version$platform),
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "Package versions:"
), con)
for (p in pk) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  writeLines(sprintf("  %-20s %s", p, v), con)
}
close(con)
cat("written\n")
