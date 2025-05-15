# -----------------------------------------------------------------------------
# Fast Grid‑Search of 2‑ and 3‑Band Vegetation Indices versus Tree‑Ring Chronology
# (data.table + frollsum + furrr parallel version)
# -----------------------------------------------------------------------------
# Author   : <your name>
# Created  : 2025‑05‑07
# Purpose  :
#   • explore a large design‑space of simple, data‑driven vegetation indices
#   • score each candidate by the mean of the 10 strongest |Pearson r|
#   • drastically reduce run‑time compared with a tidyverse loop (≥ 10×)
# -----------------------------------------------------------------------------
#   0  Libraries & helpers
#   1  Load & pre‑aggregate satellite data (daily means per plot)
#   2  Load tree‑ring chronologies (detrended, one TRI per plot)
#   3  Define formula lists + band combinations
#   4  Core function `run_combo()` → returns a one‑row data.table with the score
#   5  Parallel grid search (2‑band then 3‑band)
#   6  Inspect / export winners
# -----------------------------------------------------------------------------

# ────────────────────────────────────────────────────────────────────────────
# 0 ▸  Libraries & helper functions
# ────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)     # fast aggregation & rolling windows
  library(future)         # parallel back‑end
  library(furrr)          # future‑aware purrr
  library(progressr)      # nice progress bars with furrr
  library(lubridate)      # date helpers (for year(), yday())
  library(dplR)           # detrend(), chron()
})

# For reproducibility when using multiple workers
set.seed(1)

# Helpers --------------------------------------------------------------
# Fast column‑wise correlation between a vector y and every column of X
colCor <- function(y, X) {
  # X must be matrix, y numeric vector same length as nrow(X)
  yc  <- y - mean(y, na.rm = TRUE)
  ys2 <- sum(yc^2, na.rm = TRUE)
  apply(X, 2, function(x) {
    keep <- !is.na(x) & !is.na(yc)
    if (sum(keep) > 1) {
      xc  <- x[keep] - mean(x[keep])
      sum(xc * yc[keep]) / sqrt(sum(xc^2) * ys2)
    } else NA_real_
  })
}

# A convenience wrapper around data.table::frollsum for multiple window sizes
rollsum_dt <- function(v, windows) {
  # returns a list() suitable for := assignment
  lapply(windows, function(k) frollsum(v, n = k, align = "right"))
}

# ────────────────────────────────────────────────────────────────────────────
# 1 ▸  Satellite data  (daily means per plot)  ------------------------------
# ────────────────────────────────────────────────────────────────────────────
source("scripts/r/functions/extract_SDC_buffered.R")

SDC_raw <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/tree_coordinates_SF_12_21_PF_02_03.gpkg",
  buffer_m     = 10
)
setDT(SDC_raw)
setnames(SDC_raw, 1L, "plot")      # first column = plot id

# daily mean reflectances ----------------------------------------------------
SDC_daily <- SDC_raw[ , .(
  blue  = mean(blue , na.rm = TRUE),
  green = mean(green, na.rm = TRUE),
  red   = mean(red  , na.rm = TRUE),
  nir   = mean(nir  , na.rm = TRUE),
  swir1 = mean(swir1, na.rm = TRUE),
  swir2 = mean(swir2, na.rm = TRUE)
), by = .(plot, date)][ , `:=`(
  year = year(date),
  doy  = yday(date)
)][order(plot, date)]

setkey(SDC_daily, plot, date)

# ────────────────────────────────────────────────────────────────────────────
# 2 ▸  Tree‑ring detrended chronologies  --------------------------------------
# ────────────────────────────────────────────────────────────────────────────
source("scripts/r/functions/extract_TRW.R")

TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/tree_coordinates_SF_12_21_PF_02_03.gpkg"
)
setDT(TRW_df)

chronologies_by_plot <- TRW_df[ , {
  wide <- dcast(.SD, Year ~ treeid, value.var = "ring_width")
  rwi  <- detrend(as.rwl(wide[ , -1]), method = "Friedman")
  ch   <- chron(as.data.frame(rwi), prefix = paste0("CH_", unique(plot)))
  out  <- as.data.table(ch, keep.rownames = "year")
  out[ , `:=`(
    year = as.numeric(year) + 1649,   # calendar offset
    plot = unique(plot)
  )]
}, by = plot]

setnames(chronologies_by_plot, grep("^std", names(chronologies_by_plot), value = TRUE), "TRI")
setkey(chronologies_by_plot, plot, year)

# ────────────────────────────────────────────────────────────────────────────
# 3 ▸  Index formulae & band combinations  ------------------------------------
# ────────────────────────────────────────────────────────────────────────────
all_bands <- c("blue", "green", "red", "nir", "swir1", "swir2")

# 3.1  Two‑band formulas (must be vectorised)
TwoBand <- list(
  SimpleRatio    = \(b1, b2)  b1 / b2,
  DiffNormalized = \(b1, b2) (b1 - b2) / (b1 + b2),
  DiffSimple     = \(b1, b2)  b1 - b2,
  LogRatio       = \(b1, b2)  log1p(b1) - log1p(b2)
)

pairs_dt   <- as.data.table(t(combn(all_bands, 2)))  # unique (b1 < b2)
setnames(pairs_dt, c("band1", "band2"))

# 3.2  Three‑band formulas
ThreeBand <- list(
  ThreeBandRatio = \(b1, b2, b3) (b1 - b2) * (b1 - b3) / (b1 + b2 + b3 + 1e-6)
)

triples_dt <- as.data.table(t(combn(all_bands, 3)))
setnames(triples_dt, c("band1", "band2", "band3"))

accum_windows <- 1:46         # candidate window lengths
TOP_N         <- 10           # how many top |r| values to average

# ────────────────────────────────────────────────────────────────────────────
# 4 ▸  Core worker function  ---------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
run_combo <- function(b1, b2, fun, fun_name, b3 = NULL) {
  # Compute custom index ********************************************************
  dt <- copy(SDC_daily)[ , custom_index := if (is.null(b3)) {
    fun(get(b1), get(b2))
  } else {
    fun(get(b1), get(b2), get(b3))
  }]
  
  # Rolling sums ***************************************************************
  dt[ , (paste0("cum", accum_windows)) := rollsum_dt(custom_index, accum_windows), by = plot]
  
  # Add TRI, keep complete cases ***********************************************
  dt[ , year := year(date)]
  dt <- chronologies_by_plot[dt, on = .(plot, year)]
  dt <- dt[!is.na(TRI)]
  if (nrow(dt) < 10)               # not enough data → NA score
    return(data.table(formula_name     = fun_name,
                      band1            = b1,
                      band2            = b2,
                      band3            = ifelse(is.null(b3), NA_character_, b3),
                      mean_top_10_corr = NA_real_))
  
  # Long melt once *************************************************************
  melt_dt <- melt(dt,
                  id.vars = c("plot", "doy", "TRI"),
                  measure.vars = patterns("^cum"),
                  variable.name = "window",
                  value.name   = "cum_value")
  melt_dt[ , window := as.integer(sub("cum", "", window))]
  
  # Correlation per plot × window × doy ****************************************
  corr_dt <- melt_dt[ , .(correlation = if (.N > 1)
    cor(cum_value, TRI, use = "complete.obs")
    else NA_real_),
    by = .(plot, window, doy)]
  
  # Score **********************************************************************
  cor_vals   <- abs(na.omit(corr_dt$correlation))
  score      <- if (length(cor_vals)) mean(head(sort(cor_vals, decreasing = TRUE), TOP_N)) else NA_real_
  
  data.table(formula_name     = fun_name,
             band1            = b1,
             band2            = b2,
             band3            = ifelse(is.null(b3), NA_character_, b3),
             mean_top_10_corr = score)
}

# ────────────────────────────────────────────────────────────────────────────
# 5 ▸  Parallel grid search  ---------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
plan(multisession, workers = max(1, future::availableCores() - 1))
handlers(global = TRUE)

with_progress({
  p <- progressor(steps = nrow(pairs_dt) * length(TwoBand) +
                    nrow(triples_dt) * length(ThreeBand))
  
  # 5A ▸  Two‑band ------------------------------------------------------------
  tasks_2 <- CJ(formula_name = names(TwoBand), pair_id = seq_len(nrow(pairs_dt)))
  
  results_2 <- future_map_dfr(seq_len(nrow(tasks_2)), function(j) {
    formlab  <- tasks_2$formula_name[j]
    pair_idx <- tasks_2$pair_id[j]
    fun      <- TwoBand[[formlab]]
    
    b1 <- pairs_dt$band1[pair_idx]
    b2 <- pairs_dt$band2[pair_idx]
    
    res <- run_combo(b1, b2, fun, formlab)
    p()
    res
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))
  
  # 5B ▸  Three‑band ----------------------------------------------------------
  tasks_3 <- CJ(formula_name = names(ThreeBand), triple_id = seq_len(nrow(triples_dt)))
  
  results_3 <- future_map_dfr(seq_len(nrow(tasks_3)), function(j) {
    formlab  <- tasks_3$formula_name[j]
    triple_i <- tasks_3$triple_id[j]
    fun      <- ThreeBand[[formlab]]
    
    b1 <- triples_dt$band1[triple_i]
    b2 <- triples_dt$band2[triple_i]
    b3 <- triples_dt$band3[triple_i]
    
    res <- run_combo(b1, b2, fun, formlab, b3)
    p()
    res
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))
})

all_results <- rbindlist(list(results_2, results_3))

# ────────────────────────────────────────────────────────────────────────────
# 6 ▸  Inspect winners  --------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
setorder(all_results, -mean_top_10_corr)
print(head(all_results, 20))

# Optionally write CSV ---------------------------------------------------------
# fwrite(all_results, "results/custom_index_gridsearch_fast.csv")

# -----------------------------------------------------------------------------
# END -------------------------------------------------------------------------
