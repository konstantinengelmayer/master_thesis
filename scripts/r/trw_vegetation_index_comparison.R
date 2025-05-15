library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

source("scripts/r/functions/extract_SDC.R")
source("scripts/r/functions/extract_TRW.R")

# Define folder and list all TIFF files
folder_path <- "data/satellite_data/SDC/lindenberg_eifel_koenigsforst"
points_path <- "data/vector_data/Measured_Tree_coordinates.gpkg"
xls_folder  <- "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS"

SDC_df <- extract_SDC(folder_path, points_path)
TRW_df <- extract_TRW(xls_folder, points_path)

SDC_yearly <- SDC_df %>%
  mutate(
    year  = year(date),
    
    # 1) NDVI
    NDVI  = (nir - red) / (nir + red),
    # 2) EVI
    EVI   = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
    # 3) GNDVI
    GNDVI = (nir - green) / (nir + green),
    # 4) NDMI
    NDMI  = (nir - swir1) / (nir + swir1),
    # 5) NDWI
    NDWI  = (green - nir) / (green + nir),
    # 6) MNDWI
    MNDWI = (green - swir1) / (green + swir1),
    # 7) NBR
    NBR   = (nir - swir2) / (nir + swir2),
    # 8) NBR2
    NBR2  = (swir1 - swir2) / (swir1 + swir2),
    # 9) ARVI
    ARVI  = (nir - (2*red - blue)) / (nir + (2*red - blue)),
    # 10) OSAVI
    OSAVI = (nir - red) / (nir + red + 0.16)
  ) %>%
  # 11) NIRv (computed after NDVI)
  mutate(
    NIRv = NDVI * nir
  ) %>%
  group_by(point_id, year) %>%
  summarize(
    NDVI  = mean(NDVI,  na.rm = TRUE),
    EVI   = mean(EVI,   na.rm = TRUE),
    GNDVI = mean(GNDVI, na.rm = TRUE),
    NDMI  = mean(NDMI,  na.rm = TRUE),
    NDWI  = mean(NDWI,  na.rm = TRUE),
    MNDWI = mean(MNDWI, na.rm = TRUE),
    NBR   = mean(NBR,   na.rm = TRUE),
    NBR2  = mean(NBR2,  na.rm = TRUE),
    ARVI  = mean(ARVI,  na.rm = TRUE),
    OSAVI = mean(OSAVI, na.rm = TRUE),
    NIRv  = mean(NIRv,  na.rm = TRUE),
    .groups = "drop"
  )

# 4) Merge with Tree-Ring Width data -------------------------------------------
TRW_prepped <- TRW_df %>%
  rename(
    year     = Year,
    point_id = treeid
  )

df_merged <- SDC_yearly %>%
  left_join(TRW_prepped, by = c("point_id", "year"))

df_merged <- na.omit(df_merged)

# Check that df_merged has:
# point_id, year, NDVI, EVI, GNDVI, NDMI, NDWI, MNDWI, NBR, NBR2, ARVI, OSAVI, NIRv,
# ring_width, plot, (etc.)

# 5) Calculate correlation (per point) between each index and TRW --------------
df_correlations <- df_merged %>%
  group_by(point_id) %>%
  summarize(
    cor_NDVI  = cor(NDVI,  ring_width, use = "complete.obs"),
    cor_EVI   = cor(EVI,   ring_width, use = "complete.obs"),
    cor_GNDVI = cor(GNDVI, ring_width, use = "complete.obs"),
    cor_NDMI  = cor(NDMI,  ring_width, use = "complete.obs"),
    cor_NDWI  = cor(NDWI,  ring_width, use = "complete.obs"),
    cor_MNDWI = cor(MNDWI, ring_width, use = "complete.obs"),
    cor_NBR   = cor(NBR,   ring_width, use = "complete.obs"),
    cor_NBR2  = cor(NBR2,  ring_width, use = "complete.obs"),
    cor_ARVI  = cor(ARVI,  ring_width, use = "complete.obs"),
    cor_OSAVI = cor(OSAVI, ring_width, use = "complete.obs"),
    cor_NIRv  = cor(NIRv,  ring_width, use = "complete.obs")
  ) %>%
  ungroup()

# 6) Compute mean correlation for each vegetation index across all points ------
df_mean_cor <- df_correlations %>%
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

# Print mean correlations
print(df_mean_cor)

# 7) (Optional) Visual Comparison of Per-Point Correlations --------------------
df_cor_long <- df_correlations %>%
  pivot_longer(
    cols = starts_with("cor_"),
    names_to = "Index",
    values_to = "Correlation"
  )

ggplot(df_cor_long, aes(x = Index, y = Correlation)) +
  geom_boxplot() +
  labs(
    title = "Distribution of Per-Point Correlations with TRW",
    x     = "Vegetation Index",
    y     = "Correlation with TRW"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# End of Script ----------------------------------------------------------------

# Suppose df_merged has these columns:
#  point_id, year, ring_width, NDVI, EVI, GNDVI, NIRv, etc.

# 1) Define the vegetation indices to iterate over
vegetation_indices <- c("NDVI", "EVI", "GNDVI", "NIRv",
                        "NDMI", "NDWI", "MNDWI", "NBR", "NBR2", "ARVI", "OSAVI")

# 2) Loop over each vegetation index and create a faceted plot
for (idx in vegetation_indices) {
  
  # A) Standardize (z-score) the chosen vegetation index (idx) and TRW by point_id
  df_normalized <- df_merged %>%
    group_by(point_id) %>%
    mutate(
      # Evaluate the index column by its name
      idx_value   = .data[[idx]],
      norm_index  = (idx_value - mean(idx_value, na.rm = TRUE)) /
        sd(idx_value, na.rm = TRUE),
      norm_trw    = (ring_width - mean(ring_width, na.rm = TRUE)) /
        sd(ring_width, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # B) Plot the lines with ggplot
  p <- ggplot(df_normalized, aes(x = year)) +
    geom_line(aes(y = norm_index, color = idx)) +
    geom_line(aes(y = norm_trw, color = "TRW")) +
    facet_wrap(~ point_id) +
    scale_color_manual(
      name = NULL,
      values = c(idx = "red", "TRW" = "blue")
    ) +
    labs(
      title = paste("Z-score Comparison:", idx, "vs. TRW"),
      x = "Year",
      y = "Normalized Value (Z-score)"
    )
  
  # C) Print or save the plot (here we just print to the current device)
  print(p)
  
  # Optionally, save to file (uncomment if desired):
  # ggsave(filename = paste0("plot_", idx, "_vs_TRW.png"), plot = p, width = 8, height = 6)
}
