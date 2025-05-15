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


#### 2) LOAD & PREPARE TREE RING WIDTH (TRW) DATA ####
# Define paths
points_path <- "data/vector_data/Measured_Tree_coordinates.gpkg"
xls_folder  <- "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS"

# Extract TRW data using custom function
TRW_df <- extract_TRW(xls_folder, points_path, crop = FALSE)

# Reshape data from long to wide format:
# - Each column is one tree’s ring-width measurements.
# - 'Year' remains as a variable (later used as rownames).
TRW_df_wide <- TRW_df %>%
  select(Year, treeid, ring_width) %>%
  pivot_wider(id_cols = Year,
              names_from = treeid,
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
#### 4) CORRELATION ANALYSIS: RWI vs. VEGETATION INDICES ####
# (Assumes SDC_yearly exists and contains vegetation indices such as NDVI, EVI, etc.)

# Initialize list to store correlation results per method
method_cor_results <- list()

# Loop over each detrending method to compute correlations
for (m in methods_available) {
  cat("Computing correlations for method:", m, "\n")
  
  # a) Detrend data and convert to long format (IDs as 'point_id')
  rwi_long <- rwi_list[[m]] %>%
    tibble::rownames_to_column("year") %>%
    pivot_longer(cols = -year,
                 names_to = "point_id",
                 values_to = "RWI") %>%
    mutate(year = as.numeric(year),
           point_id = as.integer(point_id))
  
  # b) Merge RWI with vegetation indices (SDC_yearly) and drop NA rows
  df_merged <- SDC_yearly %>%
    left_join(rwi_long, by = c("point_id", "year")) %>%
    na.omit()
  
  # c) Compute per-tree correlations between RWI and each vegetation index
  df_cor <- df_merged %>%
    group_by(point_id) %>%
    summarize(
      cor_NDVI  = cor(NDVI,  RWI, use = "complete.obs"),
      cor_EVI   = cor(EVI,   RWI, use = "complete.obs"),
      cor_GNDVI = cor(GNDVI, RWI, use = "complete.obs"),
      cor_NDMI  = cor(NDMI,  RWI, use = "complete.obs"),
      cor_NDWI  = cor(NDWI,  RWI, use = "complete.obs"),
      cor_MNDWI = cor(MNDWI, RWI, use = "complete.obs"),
      cor_NBR   = cor(NBR,   RWI, use = "complete.obs"),
      cor_NBR2  = cor(NBR2,  RWI, use = "complete.obs"),
      cor_ARVI  = cor(ARVI,  RWI, use = "complete.obs"),
      cor_OSAVI = cor(OSAVI, RWI, use = "complete.obs"),
      cor_NIRv  = cor(NIRv,  RWI, use = "complete.obs")
    ) %>%
    ungroup()
  
  # d) Store the results
  method_cor_results[[m]] <- df_cor
}

# Summarize mean correlations by method using purrr::map_df
df_method_summary <- map_df(method_cor_results, function(df_cor) {
  df_cor %>%
    summarize(
      mean_cor_NDVI  = mean(cor_NDVI,  na.rm = TRUE),
      mean_cor_EVI   = mean(cor_EVI,   na.rm = TRUE),
      mean_cor_GNDVI = mean(cor_GNDVI, na.rm = TRUE),
      mean_cor_NDMI  = mean(cor_NDMI,  na.rm = TRUE),
      mean_cor_NDWI  = mean(cor_NDWI,  na.rm = TRUE),
      mean_cor_MNDWI = mean(cor_MNDWI, na.rm = TRUE),
      mean_cor_NBR   = mean(cor_NBR,   na.rm = TRUE),
      mean_cor_NBR2  = mean(cor_NBR2,  na.rm = TRUE),
      mean_cor_ARVI  = mean(cor_ARVI,  na.rm = TRUE),
      mean_cor_OSAVI = mean(cor_OSAVI, na.rm = TRUE),
      mean_cor_NIRv  = mean(cor_NIRv,  na.rm = TRUE)
    )
}, .id = "method")

# View summary correlations by detrending method
print(df_method_summary)


#### 5) PLOTTING: RWI (Friedman & Spline) vs. VEGETATION INDICES ####
# Detrend with Friedman and Spline and convert to long format

# Function to detrend and pivot to long format
get_rwi_long <- function(rwl, method, value_name) {
  detrend(rwl, method = method) %>%
    tibble::rownames_to_column("year") %>%
    pivot_longer(-year, names_to = "point_id", values_to = value_name) %>%
    mutate(year = as.numeric(year),
           point_id = as.integer(point_id))
}

rwi_friedman <- get_rwi_long(TRW_rwl, method = "Friedman", value_name = "RWI_Friedman")
rwi_spline   <- get_rwi_long(TRW_rwl, method = "Spline",   value_name = "RWI_Spline")




# Merge detrended data with SDC_yearly (vegetation indices)
df_rwi_compare <- rwi_friedman %>%
  left_join(rwi_spline, by = c("year", "point_id")) %>%
  left_join(SDC_yearly, by = c("year", "point_id")) %>%
  na.omit()

# Define vegetation indices to compare
vegetation_indices <- c("NDVI", "EVI", "GNDVI", "NDMI", 
                        "NDWI", "MNDWI", "NBR", "NBR2", 
                        "ARVI", "OSAVI", "NIRv")

# For each vegetation index, standardize (z-score) and plot comparisons
for (idx in vegetation_indices) {
  
  # A) Standardize by tree/point (z-score per point_id)
  df_normalized <- df_rwi_compare %>%
    group_by(point_id) %>%
    mutate(
      idx_value   = .data[[idx]],
      z_index     = (idx_value - mean(idx_value, na.rm = TRUE)) / sd(idx_value, na.rm = TRUE),
      z_fried     = (RWI_Friedman - mean(RWI_Friedman, na.rm = TRUE)) / sd(RWI_Friedman, na.rm = TRUE),
      z_spline    = (RWI_Spline - mean(RWI_Spline, na.rm = TRUE)) / sd(RWI_Spline, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # B) Pivot longer for plotting (one line per variable)
  df_longplot <- df_normalized %>%
    select(year, point_id, z_index, z_fried, z_spline) %>%
    pivot_longer(cols = c(z_index, z_fried, z_spline),
                 names_to = "variable",
                 values_to = "z_value")
  
  # C) Plot the z-scores (one panel per point)
  p <- ggplot(df_longplot, aes(x = year, y = z_value, color = variable)) +
    geom_line() +
    facet_wrap(~ point_id) +
    scale_color_manual(
      name = NULL,
      values = c("z_index"  = "red",
                 "z_fried"  = "blue",
                 "z_spline" = "green"),
      labels = c("z_index"  = idx,
                 "z_fried"  = "RWI (Friedman)",
                 "z_spline" = "RWI (Spline)")
    ) +
    labs(
      title = paste("Z-Score Comparison:", idx, "vs. RWI (Friedman & Spline)"),
      x = "Year",
      y = "Normalized Value (Z-score)"
    )
  
  print(p)
  # Optionally, save the plot:
  # ggsave(filename = paste0("zscore_comparison_", idx, ".png"), plot = p, width = 10, height = 6)
}


#### 6) PLOTTING: RAW TRW vs. DETRENDED RWI (Spline & Friedman) ####
# Detrend again with Spline and Friedman; this time we assume the raw TRW is in long format.
# Ensure raw TRW data uses consistent naming (Year -> year, treeid remains)

# Detrend (and pivot) using Spline
rwi_spline_long <- detrend(TRW_rwl, method = "Spline") %>%
  tibble::rownames_to_column("year") %>%
  pivot_longer(-year, names_to = "treeid", values_to = "RWI_Spline") %>%
  mutate(year = as.numeric(year),
         treeid = as.numeric(treeid))

# Detrend (and pivot) using Friedman
rwi_fried_long <- detrend(TRW_rwl, method = "Friedman") %>%
  tibble::rownames_to_column("year") %>%
  pivot_longer(-year, names_to = "treeid", values_to = "RWI_Friedman") %>%
  mutate(year = as.numeric(year),
         treeid = as.numeric(treeid))

# Convert raw TRW data to long format (rename Year to year)
TRW_long <- TRW_df %>%
  rename(year = Year)

# Merge raw TRW with detrended RWI data
df_compare <- TRW_long %>%
  left_join(rwi_spline_long, by = c("year", "treeid")) %>%
  left_join(rwi_fried_long, by = c("year", "treeid")) %>%
  na.omit()

# Compute z-scores for raw TRW and detrended RWI values (by treeid)
df_zscore <- df_compare %>%
  group_by(treeid) %>%
  mutate(
    z_raw    = (ring_width - mean(ring_width, na.rm = TRUE)) / sd(ring_width, na.rm = TRUE),
    z_spline = (RWI_Spline - mean(RWI_Spline, na.rm = TRUE)) / sd(RWI_Spline, na.rm = TRUE),
    z_fried  = (RWI_Friedman - mean(RWI_Friedman, na.rm = TRUE)) / sd(RWI_Friedman, na.rm = TRUE)
  ) %>%
  ungroup()

# Reshape to long format for plotting
df_longplot <- df_zscore %>%
  select(year, treeid, z_raw, z_spline, z_fried) %>%
  pivot_longer(cols = c(z_raw, z_spline, z_fried),
               names_to = "variable",
               values_to = "z_value") %>%
  filter(year >= 2000 & year <= 2022)

# Plot raw TRW and detrended RWI (Spline & Friedman) on the same chart
p <- ggplot(df_longplot, aes(x = year, y = z_value, color = variable)) +
  geom_line() +
  facet_wrap(~ treeid) +
  scale_color_manual(
    name = NULL,
    values = c("z_raw"    = "black",
               "z_spline" = "blue",
               "z_fried"  = "red"),
    labels = c("z_raw"    = "Raw TRW",
               "z_spline" = "Spline RWI",
               "z_fried"  = "Friedman RWI")
  ) +
  labs(
    title = "Z-Score Comparison: Non-detrended TRW vs. Spline & Friedman RWI",
    x = "Year",
    y = "Z-score (per tree)"
  )

print(p)
# Optionally, save the plot:
# ggsave("detrended_vs_raw.png", p, width = 10, height = 6)
