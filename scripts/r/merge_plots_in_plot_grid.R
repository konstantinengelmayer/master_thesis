# make sure you have these packages installed
# install.packages(c("cowplot","png","stringr"))

library(cowplot)
library(png)
library(stringr)

## 1. folders ------------------------------------------------------------------
vi_dir  <- "figures/VI_TRI_heatmaps_by_plot"   # VI-TRI heat-maps
spi_dir <- "figures/SPI_TRI_heatmaps_by_plot"  # SPI-TRI heat-maps
ts_dir  <- "figures/VI_TRI_by_plot"            # VI-TRI time-series

## 2. list PNG files -----------------------------------------------------------
vi_files  <- list.files(vi_dir,  pattern = "\\.png$", full.names = TRUE)
spi_files <- list.files(spi_dir, pattern = "\\.png$", full.names = TRUE)
ts_files  <- list.files(ts_dir,  pattern = "\\.png$", full.names = TRUE)

## 3. discover the plot IDs ----------------------------------------------------
all_files <- c(vi_files, spi_files, ts_files)
plots     <- unique(str_extract(basename(all_files), "(?:SF|PF)\\d+"))

## 4. output folder ------------------------------------------------------------
out_dir <- "figures/combined_heatmaps"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## 5. build the 2-row, 3-panel layout -----------------------------------------
for (pl in plots) {
  
  vf <- vi_files [str_detect(vi_files,  pl)]
  sf <- spi_files[str_detect(spi_files, pl)]
  tf <- ts_files [str_detect(ts_files,  pl)]
  
  if (length(vf)!=1 || length(sf)!=1 || length(tf)!=1) {
    warning(sprintf("Skipping %s: found %d VI, %d SPI, %d TS files",
                    pl, length(vf), length(sf), length(tf)))
    next
  }
  
  # read PNGs as raster grobs
  vi_plot  <- ggdraw() + draw_image(vf)  + theme(plot.margin = margin(t = 8))
  spi_plot <- ggdraw() + draw_image(sf)  + theme(plot.margin = margin(t = 8))
  ts_plot  <- ggdraw() + draw_image(tf)  + theme(plot.margin = margin(t = 8))
  
  
  # ── row 1 : VI heat-map ───────────────────────────────────────────────
  row1 <- plot_grid(
    vi_plot,
    labels       = "",
    label_size   = 12,
    label_fontface = "bold",
    label_y      = 1.05,   # ← shift label up ( >1 means above the panel)
    hjust        = 0,      # left-aligned
    vjust        = 1
  )
  
  # ── row 2 : SPI heat-map | VI time-series ─────────────────────────────
  row2 <- plot_grid(
    spi_plot, ts_plot,
    ncol        = 2,
    labels      = c("",
                    ""),
    label_size  = 12,
    label_fontface = "bold",
    label_y     = 1.05,    # ← same tweak for both labels
    hjust       = 0,
    vjust       = 1,
    rel_widths = c(1, 1.5)
  )
  
  
  ## combine rows -------------------------------------------------------------
  combined <- plot_grid(
    row1, row2,
    ncol        = 1,
    rel_heights = c(1.6, 1)   # give the heat-map row a bit more height
  )
  
  ## save ---------------------------------------------------------------------
  ggsave(
    filename = file.path(out_dir, sprintf("%s_combined_panel.png", pl)),
    plot     = combined,
    width    = 12,
    height   = 9,
    dpi      = 300,
    bg       = "white"
  )
  message("Saved combined panel for ", pl)
}
