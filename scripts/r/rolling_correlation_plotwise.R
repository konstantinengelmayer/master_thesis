################################################################################
# Rolling-Window Vegetation-Index (VI) vs. Tree-Ring Index (TRI) Analysis
# ─────────────────────────────────────────────────────────────────────────────
# This script
#   1. extracts Landsat surface-reflectance pixels (10 m buffer) for each plot
#   2. computes daily vegetation indices (VIs)
#   3. builds one detrended tree-ring chronology (TRI) per plot
#   4. rolls 1- to 24-day cumulative sums of every VI
#   5. finds, for every plot × VI, the (window, DOY) pair with the strongest
#      absolute Pearson correlation to TRI
#   6. draws per-plot time-series grids, global summary tables, and several
#      distribution / heat-map figures
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. LIBRARIES & HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)      # data wrangling
library(tidyr)      # pivot_* helpers
library(lubridate)  # date-time helpers (year(), yday(), …)
library(zoo)        # rollapply() for rolling sums
library(purrr)      # map_dfc(), walk()
library(ggplot2)    # all plotting
library(dplR)
library(terra)
library(sf)
library(zoo)        # for rolling means
library(signal)     # for Savitzky-Golay smoothing
library(readxl)
library(tibble)
library(RColorBrewer) 
library(data.table)

# custom extractors (paths relative to your repo)
source("scripts/r/functions/extract_TRW.R")          # tree-ring widths
source("scripts/r/functions/extract_SDC_buffered.R") # Landsat reflectance
source("scripts/r/functions/extract_spi_by_plot.R")
source("scripts/r/functions/extract_predictors.R")

# ─────────────────────────────────────────────────────────────────────────────
# 1. SATELLITE DATA  – daily VIs per plot
# ─────────────────────────────────────────────────────────────────────────────
## 1A  raw reflectance pixels --------------------------------------------------
SDC_df <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  buffer_m     = 1,
  method = "original",
)
names(SDC_df)[1] <- "plot"                          # first col → plot ID

#spi_summary <- extract_spi_by_plot(
#  points_path = "data/vector_data/trees_all_plots.gpkg",
#  spi_path    = "data/climate/SPI/SPI_1991_2024_12_month.tif",
#  buffer = 10,
#  years = 2000:2022,
#  month_of_interest = 12
#)

## 1B  add calendar fields -----------------------------------------------------
SDC_df <- SDC_df |>
  mutate(year = year(date),           # calendar year
         doy  = yday(date))           # day-of-year (1–365/366)

## 1C  daily mean vegetation indices ------------------------------------------
# if multiple Landsat observations fall on the same day, we average them
SDC_daily <- SDC_df |>
  group_by(plot, date) |>
  summarise(
    NDVI   = mean((nir-red)/(nir+red), na.rm = TRUE),
    EVI    = mean(2.5*(nir-red)/(nir+6*red-7.5*blue+1), na.rm = TRUE),
    NDMI   = mean((nir-swir1)/(nir+swir1), na.rm = TRUE),
    NDWI   = mean((green+nir)/(green-nir) * -nir, na.rm = TRUE),
    NBR    = mean((nir-swir2)/(nir+swir2), na.rm = TRUE),
    NBR2   = mean((swir1-swir2)/(swir1+swir2), na.rm = TRUE),
    OSAVI  = mean((nir-red)/(nir+red+0.16), na.rm = TRUE),
    NIRv   = mean(((nir-red)/(nir+red))*nir, na.rm = TRUE),
    NIRv   = mean(((nir-red)/(nir+red))*nir + 0.08, na.rm = TRUE),
    VMI    = mean((nir-swir1)*(nir-swir2)/(nir+swir1+swir2), na.rm = TRUE),
    NMDI   = mean((nir-(swir1-swir2))/(nir+(swir1-swir2)), na.rm = TRUE),
    MSI    = mean(nir/swir1, na.rm = TRUE),
    MSII   = mean(nir/swir2, na.rm = TRUE),
    MSAVI2 = mean((nir*swir1-swir2*red)/(nir*swir1+swir2*red+1e-3), na.rm = TRUE),
    VMSI   = mean((nir+swir1-swir2-red)/(nir+swir1+swir2+red+1e-3), na.rm = TRUE),
    EVDI   = mean((nir-red-swir2)/(nir+red+swir2), na.rm = TRUE),
    TBMI   = mean((swir1-swir2-red)/(swir1+swir2+red), na.rm = TRUE),
    year   = first(year),                       # keep calendar fields
    doy    = first(doy),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 2. TREE-RING CHRONOLOGIES  – detrended TRI per plot
# ─────────────────────────────────────────────────────────────────────────────
TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/trees_all_plots.gpkg"
)

chronologies_by_plot <- TRW_df |>
  group_by(plot) |>
  group_modify(~{
    wide  <- pivot_wider(.x, Year, names_from = tree_id, values_from = ring_width)
    rwi   <- detrend(as.rwl(wide[,-1]), method = "Spline", nyrs = 13, f = 0.6)         # detrend
    chron <- chron(as.data.frame(rwi), prefix = paste0("CH_", unique(.x$plot)))
    chron |>
      tibble::rownames_to_column("year") |>
      mutate(year = as.numeric(year) + 1649,     # calendar offset
             plot = unique(.x$plot))
  }) |>
  ungroup() |>
  rename(TRI = starts_with("std"))               # TRI column

# ─────────────────────────────────────────────────────────────────────────────
# 3. ROLLING-SUM VIs  (1–24-day windows, daily resolution)
# ─────────────────────────────────────────────────────────────────────────────
accum_windows <- 1:91                          # candidate window lengths

## 3A  long table: one VI per row ------------------------------------------------
SDC_long <- SDC_daily |>
  pivot_longer(
    -c(plot, date, year, doy),                  # keep date + calendar cols
    names_to  = "VI", values_to = "value"
  )

## 3B  rolling sums per window ---------------------------------------------------
SDC_cum_continuous <- SDC_long |>
  group_by(plot, VI) |>
  group_modify(~{
    bind_cols(
      .x,
      map_dfc(accum_windows, \(w)
              tibble(!!paste0("cum", w) :=
                       zoo::rollapply(.x$value, w, sum, align = "right", fill = NA)))
    )
  }) |>
  ungroup()

## 3C  tidy: one row = one window length -----------------------------------------
SDC_cum_long <- SDC_cum_continuous |>
  pivot_longer(starts_with("cum"),
               names_to  = "window", names_prefix = "cum",
               values_to = "cum_value") |>
  mutate(window = as.integer(window))

# ─────────────────────────────────────────────────────────────────────────────
# 4. JOIN TRI  +  tag forest type (coniferous / deciduous)
# ─────────────────────────────────────────────────────────────────────────────
cum_TRI_filtered <- SDC_cum_long |>
  left_join(chronologies_by_plot, by = c("plot","year")) |>
  dplyr::filter(year >= 2000) |>
  mutate(forest_type = case_when(
    plot %in% c("PF01","SF03","SF06","SF07","SF10", "SF11",
                "SF13","SF15","SF16", "SF21") ~ "coniferous",
    plot %in% c("PF02","PF03","SF01", "SF02","SF04","SF05","SF08",
                "SF09","SF12", "SF14","SF17", "SF18","SF19","SF20") ~ "deciduous",
    TRUE                                                  ~ NA_character_
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 5. CORRELATION HEAT-TABLE  &  best (window, DOY)
# ─────────────────────────────────────────────────────────────────────────────
heatmap_df <- cum_TRI_filtered |>
  group_by(plot, forest_type, VI, window, doy) |>
  summarise(
    correlation = if (sum(!is.na(cum_value) & !is.na(TRI)) > 1)
      cor(cum_value, TRI, use = "complete.obs")
    else NA_real_,
    .groups = "drop"
  )

heat_curr <- heatmap_df %>% 
  mutate(lag = 0)

# ——— 2. build the lag‐1 correlations ——————————————————————
# Start from your cum‐values + TRI joined table:
heat_prev <- cum_TRI_filtered %>%
  # pull out only what we need, and rename:
  select(
    plot, forest_type, VI, window, doy,
    cum_prev = cum_value,
    year_prev = year
  ) %>%
  # shift the year so that cum_prev comes from Y−1
  mutate(year = year_prev + 1) %>%
  # bring in *only* the current‐year TRI as TRI_current
  inner_join(
    chronologies_by_plot %>% 
      select(plot, year, TRI_current = TRI),
    by = c("plot", "year")
  ) %>%
  # now correlate the lagged cum_prev vs. TRI_current
  group_by(plot, forest_type, VI, window, doy) %>%
  summarise(
    correlation = if (n() > 1) 
      cor(cum_prev, TRI_current, use = "pairwise.complete.obs")
    else 
      NA_real_,
    .groups = "drop"
  ) %>%
  mutate(lag = 1)

# 3. Bind them and plot as before:
heat_all <- bind_rows(heat_curr, heat_prev)

saveRDS(heatmap_df, "corr_table_vi_tri.rds")

best_combo <- heatmap_df |>
  dplyr::filter(!is.na(correlation),          # keep non-missing
         correlation > 0) |>           # keep only positive r
  group_by(plot, VI) |>
  slice_max(correlation,               # largest positive r
            n = 1, with_ties = FALSE) |>
  ungroup()


# ─────────────────────────────────────────────────────────────────────────────
# 6. TIME-SERIES TABLE for plotting  (only best windows)
# ─────────────────────────────────────────────────────────────────────────────
plot_tbl <- SDC_cum_long |>
  inner_join(best_combo |>
               select(plot, VI, win_days = window, doy),
             by = c("plot","VI","window"="win_days","doy")) |>
  left_join(select(chronologies_by_plot, plot, year, TRI),
            by = c("plot","year")) |>
  rename(VI_val = cum_value) |>
  select(-value) |>
  pivot_longer(c(VI_val, TRI), names_to = "series", values_to = "value") |>
  group_by(plot, VI, series) |>
  mutate(value = scale(value)[,1]) |>
  ungroup()

plot_meta <- read_xlsx("data/vector_data/plot_metadata.xlsx")
plot_tbl <- left_join(plot_tbl, plot_meta, by = "plot")       # ▶︎ B

# ─────────────────────────────────────────────────────────────────────────────
# 7. PER-PLOT TIME-SERIES GRID  (one PNG per plot)
# ─────────────────────────────────────────────────────────────────────────────
dir.create("figures/VI_TRI_by_plot", showWarnings = FALSE)

walk(unique(plot_tbl$plot), \(pl){
  
  plt   <- dplyr::filter(plot_tbl, plot == pl)
  stats <- dplyr::filter(best_combo,   plot == pl)
  overall <- mean(abs(stats$correlation), na.rm = TRUE)
  
  ## ▲  meta for this plot
  species <- unique(plt$species)
  area    <- unique(plt$area)
  
  ## facet labels (unchanged)
  labels <- transmute(stats,
                      VI,
                      facet_lab = sprintf("%s\nr = %.2f | win = %dd | DOY = %d",
                                          VI, correlation, window, doy))
  plt <- left_join(plt, labels, by = "VI")
  
  ## plot
  p <- ggplot(plt, aes(year, value, colour = series)) +
    geom_line(linewidth = .4, na.rm = TRUE) +
    facet_wrap(~facet_lab, ncol = 4, scales = "free_y",
               strip.position = "top") +
    scale_colour_manual(values = c("VI_val" = "steelblue",
                                   "TRI"    = "black"),
                        breaks = c("VI_val","TRI"),
                        labels = c("VI","TRI")) +
    labs(title = sprintf("Best window & DOY — plot %s  (%s, %s)  |mean |r| = %.2f",
                         pl, species, area, overall),
         x = NULL, y = "z-score", colour = NULL) +
    theme(legend.position = "top",
          axis.text.x     = element_text(angle = 90,
                                         vjust = .5, hjust = 1),
          panel.spacing   = unit(3, "mm"),
          strip.placement = "outside")
  
  ## ▲  filename starts with species
  ggsave(file.path("figures/VI_TRI_by_plot",
                   sprintf("%s_%s_VI_TRI_best_windows.png", species, pl)),
         plot = p, width = 11, height = 8, dpi = 300, bg = "white")
})

dir.create("figures/VI_TRI_by_plot", showWarnings = FALSE)

walk(unique(plot_tbl$plot), \(pl){
  
  ## ── meta & helpers ─────────────────────────────────────────────
  stats_pl  <- dplyr::filter(best_combo, plot == pl)
  
  ## **1. pick the single “best” VI (|r| max) for this plot**
  best_stat <- slice_max(stats_pl, abs(correlation), n = 1, with_ties = FALSE)
  best_vi   <- best_stat$VI
  best_r    <- best_stat$correlation
  
  species <- unique(plot_tbl$species[plot_tbl$plot == pl])
  area    <- unique(plot_tbl$area   [plot_tbl$plot == pl])
  
  ## ── data for plotting ──────────────────────────────────────────
  ## **2. keep TRI + the best VI only**
  plt <- dplyr::filter(plot_tbl, plot == pl,
                series == "TRI" | VI == best_vi)
  
  ## single-facet label
  labels <- best_stat |> 
    transmute(VI,
              facet_lab = sprintf("%s\nr = %.2f | win = %dd | DOY = %d",
                                  VI, correlation, window, doy))
  plt <- left_join(plt, labels, by = "VI")
  
  ## ── plot ───────────────────────────────────────────────────────
  p <- ggplot(plt, aes(year, value, colour = series)) +
    geom_line(linewidth = .4, na.rm = TRUE) +
    scale_colour_manual(values = c("VI_val" = "steelblue",
                                   "TRI"    = "black"),
                        breaks = c("VI_val","TRI"),
                        labels = c(best_vi,"TRI")) +
    labs(title = sprintf("Best window %s & DOY %s — plot %s  (%s, %s)  | r = %.2f",
                         best_stat$window, best_stat$doy, pl, species, area, best_r),
         x = NULL, y = "z-score", colour = NULL) +
    theme(legend.position = "top",
          axis.text.x     = element_text(angle = 90, vjust = .5, hjust = 1),
          panel.spacing   = unit(3, "mm"),
          strip.placement = "outside")
  
  ## ── output ─────────────────────────────────────────────────────
  ## **3. include the best VI name in the file name**
  ggsave(file.path("figures/VI_TRI_by_plot",
                   sprintf("%s_%s_%s_TRI_best_window.png",
                           species, pl, best_vi)),
         plot = p, width = 11, height = 6, dpi = 300, bg = "white")
})


# ─────────────────────────────────────────────────────────────────────────────
# 8. SUMMARY TABLE  (mean |r|, window, DOY per VI × forest type)
# ─────────────────────────────────────────────────────────────────────────────
summary_tbl <- best_combo |>
  group_by(forest_type, VI) |>
  summarise(
    plots_used  = n(),
    mean_abs_r  = mean(abs(correlation)),
    mean_window = mean(window), sd_window = sd(window),
    mean_doy    = mean(doy),    sd_doy    = sd(doy),
    .groups     = "drop") |>
  arrange(forest_type, desc(mean_abs_r))
print(summary_tbl)

# ─────────────────────────────────────────────────────────────────────────────
# 9. DISTRIBUTION PLOTS  (window length & DOY)
# ─────────────────────────────────────────────────────────────────────────────
ggplot(best_combo, aes(window, fill = forest_type)) +
  geom_histogram(binwidth = 1, colour = "white", alpha = .7,
                 position = "identity") +
  facet_wrap(~VI, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(coniferous="forestgreen",deciduous="tan4")) +
  labs(title = "Best window length by VI & forest type",
       x = "window length (days)", y = "count", fill = NULL) +
  theme_minimal()

month_breaks <- cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31))
month_mids   <- month_breaks[-13] + diff(month_breaks)/2

ggplot(best_combo, aes(doy, fill = forest_type)) +
  geom_histogram(breaks = month_breaks, colour = "white",
                 alpha = .7, position = "identity", closed = "left") +
  facet_wrap(~VI, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = month_mids, labels = month.abb) +
  scale_fill_manual(values = c(coniferous="forestgreen",deciduous="tan4")) +
  labs(title = "Start-DOY of best window by VI & forest type",
       x = "month", y = "count", fill = NULL) +
  theme_minimal()

## --- 10B  output folder ------------------------------------------------------
dir.create("figures/VI_TRI_heatmaps_by_plot", showWarnings = FALSE)
heatmap_df <- dplyr::left_join(heatmap_df, plot_meta, by = "plot")

## --- 10C  loop over plots ----------------------------------------------------
for (pl in unique(heatmap_df$plot)) {
  
  df <- heatmap_df %>% 
    dplyr::filter(plot == pl, !is.na(correlation)) %>% 
    mutate(
      month_date   = as.Date(doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  if (nrow(df) == 0) next
  
  rng    <- range(df$correlation, na.rm = TRUE)          # data limits
  breaks <- seq(floor(rng[1] / 0.1) * 0.1,   # round down to nearest 0.1
                ceiling(rng[2] / 0.1) * 0.1, # round up   to nearest 0.1
                by = 0.1)                    # ← constant 0.1 gap
  
  # q_breaks <- quantile(df$correlation, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
  
  # how many bands (= how many colours do we need)?
  n_bands <- length(breaks) - 1         
  
  # build a long Spectral palette
  spectral_long <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "Spectral"))
  )(n_bands)
  
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(
      colour = "grey", size = 0.2, na.rm = TRUE, breaks = breaks
    ) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 1),
      expand = c(0, 0)
    ) +
    labs(
      title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)",
                      pl, unique(df$species), unique(df$area)),
      x = "Month", y = "Window length (months)"
    ) +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)   # discrete scale
  
  ggsave(
    filename = file.path("figures/VI_TRI_heatmaps_by_plot",
                         sprintf("%s_%s_VI_TRI_contour_%s.png",
                                 unique(df$species), unique(df$area), pl)),
    plot   = p,
    width  = 14, height = 12, dpi = 300, bg = "white"
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# Principal Component Analysis (PCA) of NIRv Correlation Patterns (PC1–PC6)
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)      # only needed for the optional scree & dendrogram
library(patchwork)    # idem
library(stats)        # princomp(), hclust(), dist()

# ------------------------------------------------------------------------------
# 1.  Keep rows that satisfy the temporal window you defined
#     (    ≥ 1 Mar   &   ≤ 31 Oct   )  for *all* VIs
# ------------------------------------------------------------------------------
valid_df <- heatmap_df %>% 
  dplyr::filter(!is.na(correlation)) %>% 
  mutate(start_doy = doy - window + 1) %>% 
  dplyr::filter(start_doy >= 60,        # ≥ 1 Mar
         doy       <= 304) %>%   # ≤ 31 Oct
  select(-start_doy)

# ------------------------------------------------------------------------------
# 2.  Build one huge “feature” for every combination VI × window × DOY,
#     then pivot wider so that each plot becomes a row and every feature a column
# ------------------------------------------------------------------------------
wide_df <- valid_df %>% 
  mutate(feature = paste(VI, sprintf("w%02d", window), sprintf("d%03d", doy),
                         sep = "_")) %>%           # e.g. NDVI_w03_d183
  select(plot, feature, correlation) %>% 
  pivot_wider(names_from = feature, values_from = correlation)

# ------------------------------------------------------------------------------
# 3.  Prepare a numeric matrix (plots × features)  -----------------------------
#     and compute the correlation matrix R across plots
# ------------------------------------------------------------------------------
X <- wide_df %>% select(-plot) %>% as.matrix()
rownames(X) <- wide_df$plot         # <-- 🔑 keep the plot IDs!


# Replace completely empty columns (all NAs) with 0 so cor() won’t choke
if (any(colSums(!is.na(X)) == 0)) {
  X[, which(colSums(!is.na(X)) == 0)] <- 0
}

R_all <- cor(t(X), use = "pairwise.complete.obs")    # plots × plots

# ------------------------------------------------------------------------------
# 4.  S-mode PCA on R_all, exactly the way you did for NIRv --------------------
# ------------------------------------------------------------------------------
pc_all <- princomp(covmat = R_all)

# Optional: variance explained
screeplot(pc_all, type = "lines", main = "S-mode PCA – all VIs")

# ------------------------------------------------------------------------------
# 5.  Classification = component with the strongest loading (PC1…PC6) ----------
# ------------------------------------------------------------------------------
load_mat <- unclass(loadings(pc_all))
load_df  <- as.data.frame(load_mat) %>% 
  rownames_to_column("plot")

classification <- load_df %>% 
  pivot_longer(cols      = starts_with("Comp."),
               names_to  = "PC",
               values_to = "loading") %>% 
  dplyr::filter(PC %in% paste0("Comp.", 1:3)) %>% 
  group_by(plot) %>% 
  slice_max(abs(loading), n = 1) %>% 
  ungroup()

print(classification)


# ─────────────────────────────────────────────────────────────────────────────
# Reconstruct and Visualize the Correlation Patterns for PC1–PC6
# ─────────────────────────────────────────────────────────────────────────────

# 10. Prepare loadings table restricted to PC1–PC6
pc_loadings <- load_df %>%
  select(plot, starts_with("Comp.")) %>%
  pivot_longer(
    cols = starts_with("Comp."),
    names_to = "PC",
    values_to = "loading"
  ) %>%
  dplyr::filter(PC %in% paste0("Comp.", 1:3))

# 11. Prepare list to store reconstructed maps for each PC
pc_maps_list <- list()

# 12. Loop over PC1 to PC6 to build normalized weighted correlation surfaces
for (pc_id in paste0("Comp.", 1:3)) {
  
  # Get loadings for the current PC
  loadings_this_pc <- pc_loadings %>%
    dplyr::filter(PC == pc_id) %>%
    select(plot, loading)
  
  # Calculate normalization factor to keep values bounded
  normalization_factor <- sum(abs(loadings_this_pc$loading), na.rm = TRUE)
  
  # Reconstruct weighted correlation surface for this PC
  pc_surface <- valid_df %>%
    inner_join(loadings_this_pc, by = "plot") %>%
    mutate(weighted_corr = correlation * loading) %>%
    group_by(window, doy) %>%
    summarise(
      correlation_estimate = sum(weighted_corr, na.rm = TRUE) / normalization_factor,
      .groups = "drop"
    ) %>%
    mutate(PC = pc_id)
  
  # Store the result
  pc_maps_list[[pc_id]] <- pc_surface
}

# 13. Combine all PC surfaces into one data frame
pc_maps <- bind_rows(pc_maps_list)

# ─────────────────────────────────────────────────────────────────────────────
# Plot Reconstructed PC1–PC6 Correlation Patterns as Contour Maps
# ─────────────────────────────────────────────────────────────────────────────
ggplot(pc_maps, aes(x = doy, y = window, z = correlation_estimate)) +
  geom_contour_filled(bins = 10) +
  geom_contour(bins = 10, colour = "grey30", size = 0.2) +
  facet_wrap(~PC, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 366, by = 30), name = "Day of Year") +
  scale_y_continuous(name = "Cumulative Window Length (days)") +
  labs(title = "Reconstructed Correlation Patterns for PC1–PC6") +
  theme_minimal()

overall_summary <- best_combo |>
  group_by(VI) |>
  summarise(
    plots_used  = n(),
    mean_abs_r  = mean(abs(correlation)), sd_abs_r = sd(abs(correlation)),
    mean_window = mean(window), sd_window = sd(window),
    mean_doy    = mean(doy),    sd_doy    = sd(doy),
    .groups     = "drop") |>
  arrange(desc(mean_abs_r))

#####################################################

pred_tbl <- extract_predictors(
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  raster_path  = "data/satellite_data/predictors/predictors.tif",
  buffer_width = 10
)

group_df_pca <- classification %>%              # already created above
  transmute(
    plot,
    group = factor(PC,                        # Comp.1, Comp.2, …
                   levels = paste0("Comp.", 1:6))  # keeps natural order
  )

# ──────────────────────────────────────────────────────────────
# B.  Merge with the predictor table  ─────────────────────────
# ──────────────────────────────────────────────────────────────
pred_tbl_grp <- pred_tbl %>%                  # your original table
  left_join(group_df_pca, by = "plot") %>%    # add the PCA group
  drop_na(group)                              # drop plots w/o a group

# ──────────────────────────────────────────────────────────────
# C.  Long form & facetted box-plots  ─────────────────────────
# ──────────────────────────────────────────────────────────────
pred_long <- pred_tbl_grp %>% 
  pivot_longer(
    cols      = c(slope, aspect, dem, ai),   # predictors to display
    names_to  = "variable",
    values_to = "value"
  )

ggplot(pred_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = 21, alpha = .6) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "PCA-based class",
       y = "Value",
       title = "Predictor distributions by PCA-derived group") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold")) +
  geom_point(position = position_jitter(width = .15), size = 1.2)








dir.create("figures/VI_TRI_heatmaps_by_plot", showWarnings = FALSE)
heat_all <- dplyr::left_join(heat_all, plot_meta, by = "plot")
for (pl in unique(heat_all$plot)) {
  
  df <- heat_all %>%
    dplyr::filter(plot == pl, !is.na(correlation)) %>%
    mutate(
      ext_doy      = if_else(lag == 1, doy, doy + 365),
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  if (nrow(df) == 0) next
  
  ## ── strongest positive r per VI ───────────────────────────────
  best_pts <- df %>%                               # ← NEW
    dplyr::filter(correlation > 0) %>%                    # keep positives only
    group_by(VI) %>%                               # one row per facet
    slice_max(correlation, n = 1, with_ties = FALSE) %>% 
    ungroup()                                      # ← NEW
  
  ## colour-band setup (unchanged) --------------------------------
  rng     <- range(df$correlation, na.rm = TRUE)
  breaks  <- seq(floor(rng[1]/0.1)*0.1,
                 ceiling(rng[2]/0.1)*0.1,
                 by = 0.1)
  n_bands <- length(breaks) - 1
  spectral_long <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "Spectral"))
  )(n_bands)
  
  ## ── plot ──────────────────────────────────────────────────────
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(colour = "grey", size = 0.2,
                        na.rm = TRUE, breaks = breaks) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "2 month",
                 date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 2),
                       expand = c(0, 0)) +
    
    ## — highlight the max-r point — -------------------------------
  geom_point(                                    # ← NEW
    data   = best_pts,
    aes(month_date, window_month),
    size   = 2,
    shape  = 21,
    stroke = 0.8,
    fill   = "yellow",
    colour = "black"
  ) +
    geom_text(                                     # ← optional NEW
      data  = best_pts,
      aes(month_date, window_month,
          label = sprintf("r = %.2f", correlation)),
      vjust = -0.5,
      size  = 3
    ) +
    
    labs(title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)",
                         pl, unique(df$species), unique(df$area)),
         x = "Month (prev year → current year)",
         y = "Window length (months)") +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)
  
  ggsave(
    filename = file.path(
      "figures/VI_TRI_heatmaps_by_plot",
      sprintf("%s_%s_VI_TRI_contour_%s.png",
              unique(df$species), unique(df$area), pl)),
    plot   = p,
    width  = 16, height = 9, dpi = 300, bg = "white"
  )
}


for (sp in unique(heat_all$species)) {
  
  ## ── prepare data ──────────────────────────────────────────────
  df <- heat_all %>%                                             # start with full table
    dplyr::filter(species == sp, !is.na(correlation)) %>%               # keep one species
    mutate(
      ext_doy      = if_else(lag == 1, doy, doy + 365),          # prev-year days →  366-730
      month_date   = as.Date(ext_doy - 1, origin = "2000-01-01"),# dummy calendar on year 2000
      window_month = window * 8 / 30.437                         # your old 8-day step in months
    ) %>%
    ## average the correlations across every plot that feeds this species
    group_by(VI, window_month, month_date) %>%
    summarise(correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
  
  if (nrow(df) == 0) next                                         # safety
  
  ## ── strongest positive r per VI (after averaging) ─────────────
  best_pts <- df %>%
    dplyr::filter(correlation > 0) %>%                                   # keep positives only
    group_by(VI) %>%
    slice_max(correlation, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ## ── colour-band setup (unchanged) ─────────────────────────────
  rng     <- range(df$correlation, na.rm = TRUE)
  breaks  <- seq(floor(rng[1] / 0.1) * 0.1,
                 ceiling(rng[2] / 0.1) * 0.1,
                 by = 0.1)
  n_bands <- length(breaks) - 1
  spectral_long <- colorRampPalette(
    rev(brewer.pal(11, "Spectral"))
  )(n_bands)
  
  ## ── plot ──────────────────────────────────────────────────────
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(colour = "grey", size = 0.2,
                        na.rm = TRUE, breaks = breaks) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "2 month",
                 date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(df$window_month, na.rm = TRUE)), 2),
                       expand = c(0, 0)) +
    geom_point(                                # highlight best point per VI
      data   = best_pts,
      aes(month_date, window_month),
      size   = 2,
      shape  = 21,
      stroke = 0.8,
      fill   = "yellow",
      colour = "black"
    ) +
    geom_text(
      data  = best_pts,
      aes(month_date, window_month,
          label = sprintf("r = %.2f", correlation)),
      vjust = -0.5,
      size  = 3
    ) +
    labs(title = sprintf("VI–TRI rolling correlation — species %s",
                         sp),
         subtitle = "Averaged across all plots where the species occurs",
         x = "Month (prev year → current year)",
         y = "Window length (months)") +
    theme(strip.text = element_text(face = "bold")) +
    scale_fill_manual(values = spectral_long)
  
  ggsave(
    filename = file.path(
      "figures/VI_TRI_heatmaps_by_species",
      sprintf("%s_VI_TRI_contour.png", sp)),
    plot   = p,
    width  = 16, height = 9, dpi = 300, bg = "white"
  )
}
