## 09d_assemble_figs.R
## 把单图的 PNG 面板拼成投稿用的整图：Fig 7（A-E）、Figure S6（QC 四联）、Figure S7
## 依赖 png + grid（不依赖 magick 的字体配置）
suppressPackageStartupMessages({library(png); library(grid)})

FIG <- "figures"
OUT <- file.path(FIG, "assembled")
dir.create(OUT, showWarnings = FALSE)

rd <- function(f) {
  p <- file.path(FIG, f)
  stopifnot(file.exists(p))
  readPNG(p)
}

## 画布：宽高像素 + 输出文件
new_canvas <- function(file, w, h) {
  png(file.path(OUT, file), width = w, height = h, res = 300, bg = "white")
  grid.newpage()
  c(w = w, h = h)
}

## 在 npc 坐标系放一张图：x,y = 左下角；w = npc 宽度
put <- function(img, cv, x, y, w, label = NULL, lab_size = 22) {
  iw <- dim(img)[2]; ih <- dim(img)[1]
  h_npc <- w * (ih / iw) * (cv[["w"]] / cv[["h"]])
  grid.raster(img, x = unit(x, "npc"), y = unit(y, "npc"),
              width = unit(w, "npc"), height = unit(h_npc, "npc"),
              just = c("left", "bottom"), interpolate = TRUE)
  if (!is.null(label)) {
    grid.text(label,
              x = unit(x, "npc") + unit(1.5, "mm"),
              y = unit(y + h_npc, "npc") - unit(1.0, "mm"),
              just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = lab_size))
  }
  invisible(h_npc)
}

## 找出最合适的一行高度：给一组 (图, 宽npc) 返回统一高 npc
row_h <- function(imgs, ws, cv) {
  max(mapply(function(im, w) w * (dim(im)[1] / dim(im)[2]) * (cv[["w"]] / cv[["h"]]),
             imgs, ws))
}

## 按统一行高把若干图垂直居中摆放
put_row <- function(imgs, labels, ws, y_bottom, h_row, cv, lab_size = 22) {
  x <- 0.03
  for (i in seq_along(imgs)) {
    iw <- dim(imgs[[i]])[2]; ih <- dim(imgs[[i]])[1]
    h_i <- ws[i] * (ih / iw) * (cv[["w"]] / cv[["h"]])
    y_i <- y_bottom + (h_row - h_i) / 2
    put(imgs[[i]], cv, x, y_i, ws[i], label = labels[i], lab_size = lab_size)
    x <- x + ws[i] + 0.02
  }
}

# ============================================================
# Figure 7 : A UMAP-celltype | B UMAP-condition
#            C composition   | D hub dotplot
#            E keygene dotplot（跨列居中）
# ============================================================
A <- rd("09b_umap_celltype.png")
B <- rd("09b_umap_condition.png")
C <- rd("09b_composition_bar.png")
D <- rd("09c_hub_barplot.png")
E <- rd("09c_keygene_dotplot.png")

W <- 3000; H <- 3900
cv <- new_canvas("Fig7_single_nucleus.png", W, H)

w2 <- 0.455
y1 <- 0.955
h1 <- row_h(list(A, B), c(w2, w2), cv)
put_row(list(A, B), c("A", "B"), c(w2, w2), y1 - h1, h1, cv)

y2 <- y1 - h1 - 0.035
h2 <- row_h(list(C, D), c(w2, w2), cv)
put_row(list(C, D), c("C", "D"), c(w2, w2), y2 - h2, h2, cv)

w3 <- 0.68
h3 <- w3 * (dim(E)[1] / dim(E)[2]) * (cv[["w"]] / cv[["h"]])
y3 <- y2 - h2 - 0.030
put(E, cv, x = (1 - w3) / 2, y = y3 - h3, w = w3, label = "E")

dev.off()
cat("Fig7 done\n")

# ============================================================
# Figure S6 : QC 四联（2x2）
# ============================================================
Q1 <- rd("09_qc_violin.png")
Q2 <- rd("09_qc_after_violin.png")
Q3 <- rd("09_qc_scatter.png")
Q4 <- rd("09_mt_hist.png")

cv2 <- new_canvas("FigS6_sc_qc.png", 2600, 2100)
wq <- 0.455
yq1 <- 0.960
hq1 <- row_h(list(Q1, Q2), c(wq, wq), cv2)
put_row(list(Q1, Q2), c("A", "B"), c(wq, wq), yq1 - hq1, hq1, cv2, lab_size = 20)

yq2 <- yq1 - hq1 - 0.040
hq2 <- row_h(list(Q3, Q4), c(wq, wq), cv2)
put_row(list(Q3, Q4), c("C", "D"), c(wq, wq), yq2 - hq2, hq2, cv2, lab_size = 20)
dev.off()
cat("S6 done\n")

# ============================================================
# Figure S7 : 49 基因 dotplot（单面板，加留白边）
# ============================================================
S7 <- rd("09c_lactylation_dotplot.png")
cv3 <- new_canvas("FigS7_lactylation_all.png", 1900, 2900)
put(S7, cv3, x = 0.02, y = 0.02, w = 0.96, label = NULL)
dev.off()
cat("S7 done\n")

cat("ALL_DONE\n")
