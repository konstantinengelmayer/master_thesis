# correlation_grid_search.R
# -----------------------------------------------------------------------------
# Purpose:  Grid‑search alternative **smoothing** methods (original, sgolay,
#           rolling, loess) **after** extracting Landsat reflectance (SDC)
#           exactly once.  This avoids calling the I/O‑heavy
#           `extract_SDC_buffered()` repeatedly and makes the search far faster.
# -----------------------------------------------------------------------------
#  Author:  <your-name>
#  Date:    2025-05-22
# -----------------------------------------------------------------------------

# ─────────────────────────────────────────────────────────────────────────────
# 0. PACKAGES & HELPERS  (multi‑core set‑up) ----------------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(terra)
library(future)
library(furrr)
library(progressr)
library(dplR)     # detrend(), chron()
library(zoo)      # rollapply
library(signal)   # sgolayfilt()
library(readxl)

# Custom extractors -----------------------------------------------------------
source("scripts/r/functions/extract_TRW.R")           # tree‑ring widths
source("scripts/r/functions/extract_SDC_buffered.R")  # Landsat reflectance

# ─────────────────────────────────────────────────────────────────────────────
# 1. STATIC DATA  -------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

## 1A.  Tree‑ring indices ------------------------------------------------------
TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/trees_all_plots.gpkg"
)

chronologies_by_plot <- TRW_df |>
  group_by(plot) |>
  group_modify(~{
    wide  <- pivot_wider(.x, Year, names_from = tree_id, values_from = ring_width)
    rwi   <- detrend(as.rwl(wide[,-1]), method = "Spline", nyrs = 13, f = 0.6)
    chron <- chron(as.data.frame(rwi), prefix = paste0("CH_", unique(.x$plot)))
    chron |>
      tibble::rownames_to_column("year") |>
      mutate(year = as.numeric(year) + 1649,
             plot = unique(.x$plot))
  }) |>
  ungroup() |>
  rename(TRI = starts_with("std"))

## 1B.  Raw (unsmoothed) SDC – extracted **once** -----------------------------
SDC_raw <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  buffer_m     = 1,
  method       = "original"          # <- always raw
)
setnames(SDC_raw, 1, "plot")

# ─────────────────────────────────────────────────────────────────────────────
# 2. HYPER‑PARAMETER GRID -----------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

fill_na_cols <- function(df) {
  needed <- c("sg_n", "sg_p", "rolling_n", "loess_span")
  for (n in needed) if (!n %in% names(df)) df[[n]] <- NA
  df |> relocate(all_of(needed), .after = method)
}

sgolay_grid <- expand_grid(method = "sgolay",
                           sg_n  = c(5, 7, 9, 11),
                           sg_p  = c(2, 3)) |> dplyr::filter(sg_p < sg_n)
rolling_grid <- tibble(method = "rolling", rolling_n = c(3,5,7,11))
loess_grid   <- tibble(method = "loess",   loess_span = c(0.10,0.25,0.50))
original_grid<- tibble(method = "original")

param_grid <- list(original_grid, sgolay_grid, rolling_grid, loess_grid) |>
  map(fill_na_cols) |>
  bind_rows() |>
  mutate(grid_id = row_number())

# ─────────────────────────────────────────────────────────────────────────────
# 3.   Smoothing helpers ------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

cols_to_smooth <- c("blue","green","red","nir","swir1","swir2")

apply_smoothing <- function(df, method,
                            sg_n=NULL, sg_p=NULL,
                            rolling_n=NULL,
                            loess_span=NULL) {
  # Work per‑plot to preserve temporal order ---------------------------------
  out <- split(df, df$plot) |>
    map_dfr(function(d) {
      d <- as.data.table(d)
      setorder(d, date)
      for (col in cols_to_smooth) {
        vec <- d[[col]]
        d[[col]] <- switch(method,
                           "original" = vec,
                           "sgolay"   = signal::sgolayfilt(vec, p = sg_p, n = sg_n),
                           "rolling"  = zoo::rollapply(vec, rolling_n, mean, align = "center", fill = NA),
                           "loess"    = predict(loess(vec ~ as.numeric(date), span = loess_span), d$date),
                           stop("Unknown method"))
      }
      d
    })
  return(out)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4.  GRID‑SEARCH (single run helper) ----------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

run_one <- function(params) {
  with(params, {
    message(sprintf("grid_id %s – %s", grid_id, method))
    message(sprintf("aplly smoothing"))
    # 4A.  Smooth the raw reflectance ----------------------------------------
    SDC_df <- apply_smoothing(SDC_raw, method,
                              sg_n = sg_n, sg_p = sg_p,
                              rolling_n = rolling_n,
                              loess_span = loess_span)
    message("calculate VIs")
    # 4B.  Calendar fields & daily means -------------------------------------
    SDC_daily <- SDC_df |>
      mutate(year = lubridate::year(date), doy = lubridate::yday(date)) |>
      group_by(plot, date) |>
      summarise(
        NDVI  = mean((nir-red)/(nir+red), na.rm = TRUE),
        EVI   = mean(2.5*(nir-red)/(nir+6*red-7.5*blue+1), na.rm = TRUE),
        NDMI  = mean((nir-swir1)/(nir+swir1), na.rm = TRUE),
        NDWI  = mean((green+nir)/(green-nir)*nir, na.rm = TRUE),
        NBR   = mean((nir-swir2)/(nir+swir2), na.rm = TRUE),
        NBR2  = mean((swir1-swir2)/(swir1+swir2), na.rm = TRUE),
        OSAVI = mean((nir-red)/(nir+red+0.16), na.rm = TRUE),
        NIRv  = mean(((nir-red)/(nir+red))*nir, na.rm = TRUE),
        VMI   = mean((nir-swir1)*(nir-swir2)/(nir+swir1+swir2), na.rm = TRUE),
        NMDI  = mean((nir-(swir1-swir2))/(nir+(swir1-swir2)), na.rm = TRUE),
        MSI   = mean(nir/swir1, na.rm = TRUE),
        MSII  = mean(nir/swir2, na.rm = TRUE),
        MSAVI2= mean((nir*swir1-swir2*red)/(nir*swir1+swir2*red+1e-3), na.rm = TRUE),
        VMSI  = mean((nir+swir1-swir2-red)/(nir+swir1+swir2+red+1e-3), na.rm = TRUE),
        EVDI  = mean((nir-red-swir2)/(nir+red+swir2), na.rm = TRUE),
        TBMI  = mean((swir1-swir2-red)/(swir1+swir2+red), na.rm = TRUE),
        year  = first(year),
        doy   = first(doy),
        .groups = "drop")
    
    message("perform rolling sum")
    
    # 4C.  Long VI & rolling cumulative sums (1–46 days) --------------------
    SDC_cum_long <- SDC_daily |>
      pivot_longer(-c(plot, date, year, doy), names_to = "VI", values_to = "value") |>
      group_by(plot, VI) |>
      group_modify(~bind_cols(.x, map_dfc(1:46, \(w) tibble(!!paste0("cum", w) :=
                                                              zoo::rollapply(.x$value, w, sum,
                                                                             align = "right", fill = NA))))) |>
      ungroup() |>
      pivot_longer(starts_with("cum"), names_to = "window", names_prefix = "cum",
                   values_to = "cum_value") |>
      mutate(window = as.integer(window))
    
    message("perform rolling correlation")
    # 4D.  Join TRI & correlations ------------------------------------------
    heatmap_df <- SDC_cum_long |>
      left_join(chronologies_by_plot, by = c("plot", "year")) |>
      dplyr::filter(year >= 2000) |>
      mutate(forest_type = case_when(
        plot %in% c("PF01","SF03","SF06","SF07","SF10","SF11","SF13","SF15","SF16","SF21") ~ "coniferous",
        plot %in% c("PF02","PF03","SF01","SF02","SF04","SF05","SF08","SF09","SF12","SF14","SF17","SF18","SF19","SF20") ~ "deciduous",
        TRUE ~ NA_character_)) |>
      group_by(plot, forest_type, VI, window, doy) |>
      summarise(correlation = if (sum(!is.na(cum_value) & !is.na(TRI)) > 1)
        cor(cum_value, TRI, use = "complete.obs") else NA_real_,
        .groups = "drop")
    
    message("perform summary")
    best_combo <- heatmap_df |>
      dplyr::filter(!is.na(correlation)) |>
      group_by(plot, VI) |>
      slice_max(abs(correlation), n = 1, with_ties = FALSE) |>
      ungroup()
    
    overall_summary <- best_combo |>
      group_by(VI) |>
      summarise(
        plots_used  = n(),
        mean_abs_r  = mean(abs(correlation)), sd_abs_r = sd(abs(correlation)),
        mean_window = mean(window), sd_window = sd(window),
        mean_doy    = mean(doy),    sd_doy    = sd(doy),
        .groups     = "drop") |>
      arrange(desc(mean_abs_r)) |>
      mutate(grid_id = grid_id, method = method,
             sg_n = sg_n, sg_p = sg_p, rolling_n = rolling_n, loess_span = loess_span)
    
    return(overall_summary)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# 5.  EXECUTE GRID  (sequential for-loop) -------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# Iterate **by row**, not by column! -----------------------------------------
param_rows <- split(param_grid, param_grid$grid_id)   # list of 1‑row data.frames
n_runs     <- length(param_rows)
results_list <- vector("list", n_runs)                # pre-allocate


for (i in seq_len(n_runs)) {
    message(sprintf("Running grid_id %s of %s", i, n_runs))
    results_list[[i]] <- run_one(param_rows[[i]])
}

results_df <- bind_rows(results_list)

best_per_vi <- results_df |>
  group_by(VI) |>
  slice_max(mean_abs_r, n = 1, with_ties = FALSE) |>
  select(VI, method, sg_n, sg_p, rolling_n, loess_span,
         mean_abs_r, plots_used, mean_window, mean_doy) |>
  arrange(desc(mean_abs_r))
         