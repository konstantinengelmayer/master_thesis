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
  buffer_m     = 10,
  method = "sgolay"
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
    NDWI   = mean((green+nir)/(green-nir) * nir, na.rm = TRUE),
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
    rwi   <- detrend(as.rwl(wide[,-1]), method = "Friedman")         # detrend
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
accum_windows <- 1:46                          # candidate window lengths

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

best_combo <- heatmap_df |>
  dplyr::filter(!is.na(correlation)) |>
  group_by(plot, VI) |>
  slice_max(abs(correlation), n = 1, with_ties = FALSE) |>
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
                   sprintf("sgolay_%s_%s_VI_TRI_best_windows.png", species, pl)),
         plot = p, width = 11, height = 8, dpi = 300, bg = "white")
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

## 10A  prepare binning params -------------------------------------------------
heatmap_breaks  <- seq(-1, 1, by = 0.1)
heatmap_labels  <- sprintf("%.1f – %.1f",
                           head(heatmap_breaks, -1),
                           tail(heatmap_breaks, -1))
anchor_cols     <- c("#705276", "#584d77", "#434a74", "#26368e", "#346399", "#3d8a6a",
                     "#96ba58", "#ffe53c", "#ffb230", "#fe7c24", "#fd4415")

# Create a fixed color palette for all intervals
n_intervals <- length(heatmap_breaks) - 1
contour_pal <- colorRampPalette(anchor_cols)(n_intervals)


## 10B  output folder ----------------------------------------------------------
dir.create("figures/VI_TRI_heatmaps_by_plot", showWarnings = FALSE)

heatmap_df <- left_join(heatmap_df, plot_meta, by = "plot") 

## 10C  loop over plots --------------------------------------------------------
for (pl in unique(heatmap_df$plot)) {
  df <- heatmap_df %>%
    dplyr::filter(plot == pl, !is.na(correlation)) %>%
    mutate(
      month_date   = as.Date(doy - 1, origin = "2000-01-01"),
      window_month = window * 8 / 30.437
    )
  if (nrow(df)==0) next
  
  # rebuild same dynamic palette (pal) on the fly:
  corr_bins <- cut(df$correlation,
                   breaks = heatmap_breaks,
                   include.lowest = TRUE )
  pal <- if (nlevels(corr_bins) <= length(anchor_cols)) {
    anchor_cols[seq_len(nlevels(corr_bins))]
  } else {
    colorRampPalette(anchor_cols)(nlevels(corr_bins))
  }
  names(pal) <- levels(corr_bins)
  
  p <- ggplot(df, aes(month_date, window_month, z = correlation)) +
    geom_contour_filled(breaks = heatmap_breaks, na.rm = TRUE, colour = "grey", size = 0.2) +
    facet_wrap(~VI, ncol = 4) +
    scale_x_date(date_breaks = "1 month",
                 date_labels = "%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0,
                                    ceiling(max(df$window_month, na.rm = TRUE)),
                                    by = 1),
                       expand = c(0, 0)) +
    scale_fill_manual(values = contour_pal,
                      labels = heatmap_labels,
                      drop = FALSE,  # Show all breaks in legend
                      name = "Correlation") +
    labs(title = sprintf("VI–TRI rolling correlation — plot %s  (%s, %s)", df$pl, df$species, df$area),
         x = "Month", y = "Window length (months)") +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(file.path("figures/VI_TRI_heatmaps_by_plot",
                   sprintf("sgolay_%s_%s_VI_TRI_contour_%s.png",
                           unique(df$species),
                           unique(df$area),
                           pl)),
         plot   = p,
         width  = 14, height = 8, dpi = 300, bg = "white")
}

pred_tbl <- extract_predictors(
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  raster_path  = "data/satellite_data/predictors/predictors.tif",
  buffer_width = 10
)


#---------------------------------------------------------------
# 1. turn the groups vector into a data frame and join ----------
#---------------------------------------------------------------

# 6A  vectorise the heat-maps
mat <- heatmap_df |>
  dplyr::filter(!is.na(correlation)) |>
  pivot_wider(names_from = c(VI, window, doy),
              values_from = correlation) |>
  column_to_rownames("plot") |> as.matrix()

cl <- hclust(dist(mat), method = "ward.D2")
groups <- cutree(cl, k = 3)

group_df <- tibble(
  plot   = names(groups),          # row names from hclust object
  cluster = factor(groups)         # nicer on the x-axis
)

pred_tbl_grp <- pred_tbl %>% 
  left_join(group_df, by = "plot") %>%
  na.omit()

#---------------------------------------------------------------
# 2. long form + facetted box-plots -----------------------------
#---------------------------------------------------------------
pred_long <- pred_tbl_grp %>% 
  pivot_longer(
    cols = c(slope, aspect, dem, ai),   # predictors you want to show
    names_to  = "variable",
    values_to = "value"
  )

ggplot(pred_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(outlier.shape = 21, alpha = .6) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Cluster", y = "Value",
       title = "Predictor distributions by dendrogram cluster") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))+
  geom_point()
# ─────────────────────────────────────────────────────────────────────────────
# Principal Component Analysis (PCA) of NIRv Correlation Patterns (PC1–PC6)
# ─────────────────────────────────────────────────────────────────────────────

# 1. Extract only NIRv rows with valid correlation values
nirv_df <- heatmap_df %>%
  dplyr::filter(VI == "EVI", !is.na(correlation))

# 2. Pivot to wide format: one column per plot, one row per (window, DOY)
nirv_wide <- nirv_df %>%
  select(plot, window, doy, correlation) %>%
  tidyr::pivot_wider(
    names_from  = plot,
    values_from = correlation
  ) %>%
  arrange(window, doy)

# 3. Convert to matrix, dropping window and DOY columns
X <- nirv_wide %>% select(-window, -doy) %>% as.matrix()

# 4. Compute correlation matrix across plots (S-mode)
R <- cor(X, use = "pairwise.complete.obs")

# 5. Run PCA on this correlation matrix
pc_nirv <- princomp(covmat = R)

# 6. Visualize variance explained with a scree plot
screeplot(pc_nirv, type = "lines", main = "S-mode PCA (NIRv)")

# 7. Show summary of variance explained by components
summary(pc_nirv)

# 8. Extract plain numeric loadings matrix with plot IDs
load_mat <- unclass(loadings(pc_nirv))
load_df  <- as.data.frame(load_mat) %>%
  tibble::rownames_to_column("plot")

# 9. Classify each plot by its strongest loading among PC1–PC6
classification <- load_df %>%
  pivot_longer(
    cols      = starts_with("Comp."),
    names_to  = "PC",
    values_to = "loading"
  ) %>%
  dplyr::filter(PC %in% paste0("Comp.", 1:6)) %>%
  group_by(plot) %>%
  slice_max(abs(loading), n = 1) %>%
  ungroup()

# Display classification table
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
  dplyr::filter(PC %in% paste0("Comp.", 1:6))

# 11. Prepare list to store reconstructed maps for each PC
pc_maps_list <- list()

# 12. Loop over PC1 to PC6 to build normalized weighted correlation surfaces
for (pc_id in paste0("Comp.", 1:6)) {
  
  # Get loadings for the current PC
  loadings_this_pc <- pc_loadings %>%
    dplyr::filter(PC == pc_id) %>%
    select(plot, loading)
  
  # Calculate normalization factor to keep values bounded
  normalization_factor <- sum(abs(loadings_this_pc$loading), na.rm = TRUE)
  
  # Reconstruct weighted correlation surface for this PC
  pc_surface <- nirv_df %>%
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
                                                                                                                                                                              
