################################################################################
# Script: Detrending and Vegetation Index Comparison
# Author: <Your Name>
# Date: YYYY-MM-DD
#
# This script processes tree ring width (TRW) data:
# 1. Loads required libraries and functions.
# 2. Loads and reshapes TRW data into a dplR rwl object.
# 3. Detrends TRW using several methods and computes ring-width indices (RWI).
# 4. Computes correlations between RWI and vegetation indices.
# 5. Plots comparisons of detrended RWI (z-scores) versus vegetation indices.
# 6. Plots comparisons between raw and detrended TRW values.
################################################################################

#### 1) LOAD LIBRARIES & FUNCTIONS ####
library(dplyr)      # Data manipulation
library(tidyr)      # Data reshaping (pivot_wider/pivot_longer)
library(dplR)       # Detrending and RWI calculation
library(ggplot2)    # Plotting
library(purrr)      # Functional programming (map_df)
library(tibble)     # For rownames_to_column

# Load custom functions (e.g., extract_TRW)
source("scripts/r/functions/extract_TRW.R")
source("scripts/r/functions/extract_SDC_buffered.R") # Landsat reflectance


# ─────────────────────────────────────────────────────────────────────────────
# 1. SATELLITE DATA  – daily VIs per plot
# ─────────────────────────────────────────────────────────────────────────────
## 1A  raw reflectance pixels --------------------------------------------------
SDC_df <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  buffer_m     = 1,
  method = "sgolay"
)
names(SDC_df)[1] <- "plot"                          # first col → plot ID

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

SDC_yearly <- SDC_daily %>%
  dplyr::filter(doy >= 90 & doy <= 300) %>%
  group_by(plot, year) %>%
  summarise(across(NDVI:TBMI, sum, na.rm = TRUE), .groups = "drop")



#### 2) LOAD & PREPARE TREE RING WIDTH (TRW) DATA ####
TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/trees_all_plots.gpkg"
)


# Reshape data from long to wide format:
# - Each column is one tree’s ring-width measurements.
# - 'Year' remains as a variable (later used as rownames).
TRW_df_wide <- TRW_df %>%
  select(Year, tree_id, ring_width) %>%
  pivot_wider(id_cols = Year,
              names_from = tree_id,
              values_from = ring_width) %>%
  arrange(Year)

# Convert to a dplR rwl object:
# - Set row names to be the years.
TRW_rwl <- as.data.frame(TRW_df_wide)
rownames(TRW_rwl) <- TRW_rwl$Year
TRW_rwl <- as.rwl(TRW_rwl[ , -1])  # Drop the Year column

# Print basic reports and statistics
rwl.report(TRW_rwl)
rwl.stats(TRW_rwl)

# Remove years where all tree measurements are NA
all_na_rows <- rowSums(is.na(TRW_rwl)) == ncol(TRW_rwl)
TRW_rwl <- TRW_rwl[!all_na_rows, ]


#### 3) DETRENDING & COMPUTING RWI WITH MULTIPLE METHODS ####
# Define the methods to apply
methods_available <- c("Ar", "Mean", "Spline", "Friedman",
                       "AgeDepSpline", "ModNegExp", "ModHugershoff")

# Initialize a list to store RWI results by method
rwi_list <- list()

# Loop through each method, detrend and store the result
for (m in methods_available) {
  cat("Detrending with method:", m, "\n")
  rwi_list[[m]] <- detrend(rwl = TRW_rwl, method = m)
}

plot.rwl(as.rwl(rwi_list$Spline), plot.type = "spag")
###############################################################################
## Hyper-parameter grid-search for dplR::detrend()                           ##
## - compares many detrending settings, computes RWI–vegetation correlations ##
###############################################################################
library(furrr)        # future + purrr syntax

## use all logical cores; adapt to your cluster if needed
future::plan(multisession)    # on Windows/macOS
# future::plan(multicore)     # on Linux

###############################################################################
## 1.  Build the parameter grid   (exactly as you already defined)           ##
###############################################################################
spline_grid <- crossing(
  method = "Spline",
  nyrs   = seq(1,60,1),
  f      = seq(0.5, 0.9, 0.05)
)



grid <- bind_rows(spline_grid) %>%
  mutate(method_id = sprintf("%s_%03d", method, row_number()))

###############################################################################
## 2.  A function that does *everything* for ONE grid row                    ##
###############################################################################
evaluate_detrend <- function(params_row,
                             trw  = TRW_df,
                             vix  = SDC_yearly) {
  
  ## ---- 2·A  collect arguments for detrend() -------------------------------
  detr_args <- list(
    method = params_row$method,
    nyrs   = params_row$nyrs,
    f      = params_row$f
  ) |>
    purrr::discard(is.na)
  
  ## ---- 2·B  build per-plot chronologies -----------------------------------
  chronologies_by_plot <- trw %>%
    group_by(plot) %>%
    group_modify(~{
      wide <- pivot_wider(.x,
                          id_cols     = Year,
                          names_from  = tree_id,
                          values_from = ring_width) %>%
        tibble::column_to_rownames("Year")
      
      rwi  <- do.call(detrend, c(list(as.rwl(wide)), detr_args))
      
      chron <- chron(as.data.frame(rwi),
                     prefix = paste0("CH_", unique(.x$plot)))
      
      chron |>
        tibble::rownames_to_column("year") |>
        mutate(year = as.numeric(year),
               plot = unique(.x$plot))
    }) |>
    ungroup() |>
    rename(TRI = starts_with("std"))
  
  ## ---- 2·C  merge with veg-indices & get correlations ----------------------
  df_cor <- vix |>
    left_join(chronologies_by_plot, by = c("plot", "year")) |>
    na.omit() |>
    group_by(plot) |>
    summarize(across(
      c(NDVI, EVI, NDMI, NDWI, NBR, NBR2, OSAVI, NIRv, VMI,
        NMDI, MSI, MSII, MSAVI2, VMSI, EVDI, TBMI),
      \(x) cor(x, TRI, use = "complete.obs"),
      .names = "cor_{.col}"
    ), .groups = "drop")
  
  ## ---- 2·D  collapse to a single row with mean over all indices -----------
  out <- df_cor |>
    summarize(across(starts_with("cor_"), mean, na.rm = TRUE)) |>
    mutate(mean_cor_all = rowMeans(across(starts_with("cor_")), na.rm = TRUE))
  
  ## attach parameter info so we know what produced the result
  bind_cols(params_row, out)
}

###############################################################################
## 3.  Run the function in parallel for every grid row                       ##
###############################################################################
## split the grid into a list of single-row data frames so the function
## receives one row at a time:
grid_list <- split(grid, seq_len(nrow(grid)))

df_method_summary <- future_map_dfr(
  grid_list,
  evaluate_detrend,
  .options = furrr_options(seed = TRUE),   # reproducible correlations
  .progress = TRUE                         # nice progress bar
) |>
  arrange(desc(mean_cor_all))

print(df_method_summary)

# Detrend again with Spline and Friedman; this time we assume the raw TRW is in long format.
# Ensure raw TRW data uses consistent naming (Year -> year, treeid remains)

# Detrend (and pivot) using Spline
rwi_spline_long <- detrend(TRW_rwl, method = "Spline") %>%
  tibble::rownames_to_column("year") %>%
  pivot_longer(-year, names_to = "tree_id", values_to = "RWI_Spline") %>%
  mutate(year = as.numeric(year),
         tree_id = as.numeric(tree_id))

# Detrend (and pivot) using Friedman
rwi_fried_long <- detrend(TRW_rwl, method = "Friedman") %>%
  tibble::rownames_to_column("year") %>%
  pivot_longer(-year, names_to = "tree_id", values_to = "RWI_Friedman") %>%
  mutate(year = as.numeric(year),
         tree_id = as.numeric(tree_id))

# Convert raw TRW data to long format (rename Year to year)
TRW_long <- TRW_df %>%
  rename(year = Year) %>%
  mutate(tree_id = as.numeric(tree_id))

# Merge raw TRW with detrendedtree_id# Merge raw TRW with detrended RWI data
df_compare <- TRW_long %>%
  left_join(rwi_spline_long, by = c("year", "tree_id")) %>%
  left_join(rwi_fried_long, by = c("year", "tree_id")) %>%
  na.omit()

# Compute z-scores for raw TRW and detrended RWI values (by treeid)
df_zscore <- df_compare %>%
  group_by(tree_id) %>%
  mutate(
    z_raw    = (ring_width - mean(ring_width, na.rm = TRUE)) / sd(ring_width, na.rm = TRUE),
    z_spline = (RWI_Spline - mean(RWI_Spline, na.rm = TRUE)) / sd(RWI_Spline, na.rm = TRUE),
    z_fried  = (RWI_Friedman - mean(RWI_Friedman, na.rm = TRUE)) / sd(RWI_Friedman, na.rm = TRUE)
  ) %>%
  ungroup()

# Reshape to long format for plotting
df_longplot <- df_zscore %>%
  select(year, tree_id, z_raw, z_spline, z_fried) %>%
  pivot_longer(cols = c(z_raw, z_spline, z_fried),
               names_to = "variable",
               values_to = "z_value") %>%
  dplyr::filter(year >= 1900 & year <= 2022)



# Get unique tree IDs
tree_ids <- unique(df_longplot$tree_id)

# Loop over each tree_id
for (tree in tree_ids) {
  
  # Filter for current tree
  df_tree <- df_longplot %>% dplyr::filter(tree_id == tree)
  
  # Create the plot
  p <- ggplot(df_tree, aes(x = year, y = z_value, color = variable)) +
    geom_line() +
    scale_color_manual(
      name = NULL,
      values = c("z_raw" = "black", "z_spline" = "blue", "z_fried" = "red"),
      labels = c("z_raw" = "Raw TRW", "z_spline" = "Spline RWI", "z_fried" = "Friedman RWI")
    ) +
    labs(
      title = paste("Z-Score Comparison for Tree ID:", tree),
      x = "Year",
      y = "Z-score"
    ) +
    theme_minimal()
  
  # Print or save the plot
  print(p)
  
  # Optional: Save to file, e.g., PNG
  # ggsave(paste0("Tree_", tree, "_Comparison.png"), plot = p, width = 8, height = 5)
}
