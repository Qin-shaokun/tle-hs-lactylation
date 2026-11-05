## 10_assemble_main.R
## 把 Fig 1–6 的单面板按 180 mm 宽（2126 px @300 dpi）拼成投稿整图
source("config.R")
suppressPackageStartupMessages({library(png); library(grid)})

OUT <- file.path(FIG_DIR, "assembled")
dir.create(OUT, showWarnings = FALSE)

DPI <- 300
W   <- 7.087          # 180 mm
MX  <- 0.28           # 左右边距
MT  <- 0.14; MB <- 0.14; GX <- 0.16; GY <- 0.18

## 每个 figure：rows -> 每行若干 panel（文件名, 目标宽度英寸）
plan <- list(
  Fig1 = list(
    list(c("02_volcano.png", 3.35), c("print/p1B.png", 2.95))
  ),
  Fig2 = list(
    list(c("print/p2A.png", 3.30), c("03_module_traits.png", 3.00))
  ),
  Fig3 = list(
    list(c("04_lasso.png", 4.30)),
    list(c("04_rf_importance.png", 3.15), c("04_venn.png", 3.15)),
    list(c("04_roc_train.png", 4.30))
  ),
  Fig4 = list(
    list(c("05_roc_valid.png", 2.60),  c("05_nomogram.png", 3.70)),
    list(c("05_calibration.png", 4.60))
  ),
  Fig5 = list(
    list(c("print/p5A.png", 6.45)),
    list(c("print/p5B.png", 5.00))
  ),
  Fig6 = list(
    list(c("08_enrich_GO.png", 3.15),  c("08_gsea_focus.png", 3.15)),
    list(c("08_enrich_KEGG.png", 3.15), c("08_gsea_KEGG.png", 3.15))
  )
)

lab <- LETTERS

for (fig in names(plan)) {
  rows <- plan[[fig]]
  imgs <- lapply(rows, function(r) lapply(r, function(p) {
    f <- file.path(FIG_DIR, p[1]); stopifnot(file.exists(f)); readPNG(f)
  }))
  ws <- lapply(rows, function(r) as.numeric(sapply(r, `[`, 2)))

  ## 行高
  hs <- numeric(length(rows))
  for (i in seq_along(rows)) {
    hs[i] <- max(mapply(function(im, w) w * dim(im)[1] / dim(im)[2], imgs[[i]], ws[[i]]))
  }
  H <- MT + sum(hs) + GY * (length(rows) - 1) + MB

  file <- file.path(OUT, sprintf("%s_print.png", fig))
  png(file, width = round(W * DPI), height = round(H * DPI), res = DPI, bg = "white")
  grid.newpage()

  k <- 0
  y <- H - MT
  for (i in seq_along(rows)) {
    ni <- length(rows[[i]])
    tot <- sum(ws[[i]]) + GX * (ni - 1)
    x <- MX + (W - 2 * MX - tot) / 2      # 行内水平居中
    for (j in seq_len(ni)) {
      k <- k + 1
      im <- imgs[[i]][[j]]; w <- ws[[i]][[j]]
      h <- w * dim(im)[1] / dim(im)[2]
      yy <- y - hs[i] + (hs[i] - h) / 2   # 行内垂直居中
      grid.raster(im, x = unit(x, "in"), y = unit(yy, "in"),
                  width = unit(w, "in"), height = unit(h, "in"),
                  just = c("left", "bottom"), interpolate = TRUE)
      grid.text(lab[k], x = unit(x - 0.05, "in"), y = unit(yy + h + 0.02, "in"),
                just = c("right", "top"),
                gp = gpar(fontface = "bold", fontsize = 11))
      x <- x + w + GX
    }
    y <- y - hs[i] - GY
  }
  dev.off()
  cat(sprintf("%s done: %.2f x %.2f in\n", fig, W, H))
}

cat("MAIN_FIGS_DONE\n")
