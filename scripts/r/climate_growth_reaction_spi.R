################################################################################
# Climate vs SPI strength — combined current‐ & previous‐year heatmaps
################################################################################

## 0. PACKAGES & SETTINGS ------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(data.table)
  library(lubridate); library(SPEI); library(dplR)
  library(sf);    library(purrr);   library(writexl)
  library(ggplot2); library(RColorBrewer)
})

# — user paths  
xls_folder  <- "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS"
points_path <- "data/vector_data/trees_all_plots.gpkg"

ppt_tif   <- "data/climate/dwd_4_variables/precipitation_1991_2024_merged.tif"
other_vars<- c(
  tmax = "data/climate/dwd_4_variables/air_temp_max_1991_2024_merged.tif",
  pet  = "data/climate/dwd_4_variables/evapo_p_1991_2024_merged.tif",
  soi  = "data/climate/dwd_4_variables/soi_moist_1991_2024_merged.tif"
)
years_keep  <- 1991:2022
max_spi_k   <- 24   # SPI1…SPI9
out_xlsx    <- "plot_spi_vs_climate_strength_lagged.xlsx"

## assume `heat_all` is your VI–TRI table with columns:
## plot, VI, correlation, window, lag ("previous"/"current"), doy

vi_peaks <- heat_all %>%                 # or use best_combo if you have it
  dplyr::filter(correlation > 0) %>%            # keep positive r (or use abs())
  group_by(plot) %>%                     # ONE point per plot
  slice_max(correlation, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(
    month_vi   = month(as.Date(doy - 1, origin = "2000-01-01")),
    ext_month  = if_else(lag == "previous", month_vi, month_vi + 12L),
    win_len_vi = ceiling(window * 8 / 30.437)   # → months
  ) %>% 
  select(plot, ext_month, win_len_vi, VI, correlation)  # keep what we need


## 1. EXTRACT TREE‐RING CHRONOLOGIES -------------------------------------------
source("scripts/r/functions/extract_TRW.R")
TRW_df <- extract_TRW(xls_folder = xls_folder, points_path = points_path)

chron_by_plot <- TRW_df %>%
  group_by(plot) %>%
  group_modify(~{
    wide <- pivot_wider(.x, Year, names_from = tree_id, values_from = ring_width)
    rwi  <- detrend(as.rwl(wide[,-1]), method="Spline", nyrs=13, f=0.6)
    chron(as.data.frame(rwi), prefix = paste0("CH_", unique(.x$plot))) %>%
      tibble::rownames_to_column("year") %>%
      mutate(year = as.numeric(year) + 1649,
             plot = unique(.x$plot))
  }) %>%
  ungroup() %>%
  rename(TRI = starts_with("std"))

growth_dt <- as.data.table(
  chron_by_plot %>% select(plot, year, TRI) %>% dplyr::filter(!is.na(TRI))
)

## 2. EXTRACT MONTHLY CLIMATE & PREPARE ----------------------------------------
source("scripts/r/functions/extract_climate_variables.R")
cores <- sf::read_sf(points_path)
var_files <- c(ppt = ppt_tif, other_vars)

plot_clim <- extract_plot_climate(
  points_sf  = cores,
  var_files  = var_files,
  buffer_m   = 1,
  years      = years_keep,
  fun        = mean,
  na.rm      = TRUE,
  start_date = as.Date("1991-01-01"),
  plot_col   = "plot"
)

clim_dt <- as.data.table(plot_clim)
clim_dt[, date := as.IDate(sprintf("%04d-%02d-15", year, month))]
setorder(clim_dt, plot, date)

## 3. COMPUTE SPI‐k ------------------------------------------------------------
ppt_dt <- clim_dt[, .(plot, year, month, ppt)]

spi_by_scale <- function(x, scales=1:max_spi_k) {
  ts_ppt <- ts(x$ppt, start = c(x$year[1], x$month[1]), frequency = 12)
  out    <- lapply(scales, function(k) as.numeric(spi(ts_ppt, k)$fitted))
  setnames(as.data.table(out), paste0("SPI", scales))
}

ppt_spi <- ppt_dt[ , cbind(.SD, spi_by_scale(.SD)), by = plot ]

## 4. WIDE → LONG FOR SPI -----------------------------------------------------
spi_long <- melt(
  ppt_spi,
  id.vars        = c("plot","year","month"),
  measure.vars   = patterns("^SPI"),
  variable.name  = "k",
  value.name     = "spi",
  variable.factor= FALSE
)[, `:=`(
  window_length = as.integer(gsub("\\D+", "", k)),
  start_month   = month
)][, k := NULL][!is.na(spi)]

# keep only years with TRI
spi_long <- spi_long[year %in% growth_dt$year]

## 5. MERGE & CORRELATE — CURRENT YEAR (lag = 0) -------------------------------
setkey(spi_long,  plot, year)
setkey(growth_dt, plot, year)

spi_curr_dt     <- spi_long[growth_dt, nomatch=0L]
corr_spi_curr  <- spi_curr_dt[
  , .(r = cor(spi, TRI, use = "pairwise.complete.obs")),
  by = .(plot, start_month, window_length)
][, lag := "current"]

## 6. MERGE & CORRELATE — PREVIOUS YEAR (lag = -1) -----------------------------
# shift spi_long’s year forward so SPI from Y-1 pairs with TRI in Y
spi_prev       <- copy(spi_long)[, year := year + 1L]
spi_prev_dt <- inner_join(
  spi_prev,
  growth_dt,
  by = c("plot", "year")
)

corr_spi_prev  <- spi_prev_dt[
  , .(r = cor(spi, TRI, use = "pairwise.complete.obs")),
  by = .(plot, start_month, window_length)
][, lag := "previous"]

# 7. STACK & MAKE EXTENDED-MONTH INDEX numeric 1–24
corr_all <- rbindlist(list(corr_spi_curr, corr_spi_prev), use.names=TRUE)

corr_all[, ext_month := ifelse(
  lag == "previous",    # months  1–12 = SPI from Y−1
  start_month,
  start_month + 12      # months 13–24 = SPI from Y
)]


## 8. SAVE CORRELATION TABLES -------------------------------------------------
writexl::write_xlsx(corr_all, out_xlsx)

## 9. PLOTTING FUNCTION & LOOP ------------------------------------------------
draw_spi_heatmap_lagged <- function(data, plot_id,
                                    vi_peak_tbl,      # ← NEW
                                    r_break = 0.1) {
  
  df <- data[plot == plot_id & !is.na(r)]
  if (nrow(df) == 0) return(NULL)
  
  ## ------- bin colours (unchanged) ------------------------------
  r_min <- floor(min(df$r)/r_break)*r_break
  r_max <- ceiling(max(df$r)/r_break)*r_break
  edges <- seq(r_min, r_max, by = r_break)
  lbls  <- paste0(head(edges, -1), "–", tail(edges, -1))
  df[, r_bin := cut(r, breaks = edges,
                    include.lowest = TRUE, labels = lbls)]
  pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(length(lbls)))
  names(pal) <- lbls
  
  ## ------- point of strongest VI -------------------------------
  peak <- vi_peak_tbl %>% dplyr::filter(plot == plot_id)    # one row
  
  ggplot(df, aes(ext_month, factor(window_length), fill = r_bin)) +
    geom_tile(colour = "white") +
    
    ## the yellow diamond marking VI’s best window/month  ↑ NEW ↑
    geom_point(data = peak,
               aes(x = ext_month, y = factor(win_len_vi)),
               shape = 23, size = 1,
               fill  = "yellow", colour = "black", stroke = 0.8) +
    
    scale_fill_manual(values = pal,
                      drop  = FALSE,
                      name  = "Pearson r",
                      guide = guide_legend(reverse = TRUE)) +
    labs(
      title = sprintf("SPI–TRI correlation (lag −1 & 0) — %s", plot_id),
      x     = "Month (1–12 = SPI in Y−1; 13–24 = SPI in Y)",
      y     = "Window length (months)"
    ) +
    coord_equal() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x       = element_text(angle = 90, vjust = 0.5),
      panel.grid        = element_blank(),
      legend.key.height = unit(0.6, "cm")
    )
}


dir.create("figures/SPI_TRI_heatmaps_by_plot", showWarnings = FALSE)

plots <- unique(corr_all$plot)
for (p in plots) {
  g  <- draw_spi_heatmap_lagged(corr_all, p, vi_peaks)  # pass peaks
  if (is.null(g)) next
  fn <- sprintf("figures/SPI_TRI_heatmaps_by_plot/SPI_TRI_heatmap_%s.png", p)
  ggsave(fn, g, width = 6, height = 4, dpi = 300, bg = "white")
  message("saved: ", fn)
}

