## 快速自检：config 可加载、基因集可读、路径正确
source("config.R", encoding = "UTF-8")
cat("ROOT      :", ROOT, "\n")
cat("DATA_DIR  :", DATA_DIR, "\n")
cat("GENESET   :", GENESET_CSV, "exists =", file.exists(GENESET_CSV), "\n")
gs <- read.csv(GENESET_CSV, stringsAsFactors = FALSE)
cat("genes     :", nrow(gs), "| categories:", paste(sort(unique(gs$category)), collapse = "/"), "\n")
cat("main cohort:", GSE256068$accession, "| case/ctrl:", GSE256068$case_label, "/", GSE256068$ctrl_label, "\n")
cat("sc samples:", length(cfg_sc_condition_map), "| batches:", length(unique(unlist(cfg_sc_batch_map))), "\n")
