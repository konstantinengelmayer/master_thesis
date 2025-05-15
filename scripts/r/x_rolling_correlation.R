################################################################################
# Script: SDC (Plot-Level) & Tree-Ring Chronology Correlation (Daily Version)
# Author: <Your Name>
# Date: YYYY-MM-DD
#
# This script:
# 1) Extracts Landsat data for each plot using a 10m buffer, computes daily
#    vegetation indices (without monthly aggregation), and computes daily 
#    anomalies based on long-term day-of-year means.
# 2) Builds a tree-ring chronology (detrended, robust mean) for each plot.
# 3) Computes rolling (cumulative) sums of each vegetation index over 1-24-day
#    windows.
# 4) Merges the cumulative VIs with the tree-ring chronologies and calculates 
#    correlations, visualizing the results in a heatmap with the x-axis as day of year.
################################################################################

#-----------------------
# 1) LOAD LIBRARIES
#-----------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(zoo)      # For rollapply
library(purrr)
library(dplR)

#-----------------------
# 2) LOAD CUSTOM FUNCS
#-----------------------
source("scripts/r/functions/extract_TRW.R")
source("scripts/r/functions/extract_SDC_buffered.R")
# - extract_TRW:        Extracts tree-ring data from XLS & GPKG
# - extract_SDC_buffered: Extracts Landsat data in a 10m buffer around each plot

################################################################################
# 2) EXTRACT & PREPARE SATELLITE DATA (SDC)
################################################################################

# 2A) EXTRACT SATELLITE DATA AT PLOT LEVEL
SDC_df <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst","data/satellite_data/SDC/kellerwald_lahntal") ,
  points_path = "data/vector_data/tree_coordinates_SF_12_21_PF_02_03.gpkg",
  buffer_m = 10
)
names(SDC_df)[1] <- "plot"
# Now each row is: [plot, date, blue, green, red, nir, swir1, swir2]

# MODIFIED: Extract Year and Day-of-Year (DOY) instead of Month.
SDC_df <- SDC_df %>%
  mutate(
    year = year(date),
    doy  = yday(date)
  )

# 2B) COMPUTE DAILY MEAN VEGETATION INDICES
# (If there are multiple observations per day, average them.)
SDC_daily <- SDC_df %>%
  group_by(plot, date) %>%
  summarize(
    NDVI  = mean((nir - red) / (nir + red), na.rm = TRUE),
    EVI   = mean(2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1), na.rm = TRUE),
    GNDVI = mean((nir - green) / (nir + green), na.rm = TRUE),
    NDMI  = mean((nir - swir1) / (nir + swir1), na.rm = TRUE),
    NDWI  = mean((green - nir) / (green + nir), na.rm = TRUE),
    NBR   = mean((nir - swir2) / (nir + swir2), na.rm = TRUE),
    NBR2  = mean((swir1 - swir2) / (swir1 + swir2), na.rm = TRUE),
    ARVI  = mean((nir - (2 * red - blue)) / (nir + (2 * red - blue)), na.rm = TRUE),
    OSAVI = mean((nir - red) / (nir + red + 0.16), na.rm = TRUE),
    NIRv  = mean(((nir - red) / (nir + red)) * nir, na.rm = TRUE),
    VMI   = mean(((nir - swir1) * (nir - swir2)) / (nir + swir1 + swir2), na.rm = TRUE),
    NMDI  = mean( (nir - (swir1 - swir2)) / (nir + (swir1 - swir2)), na.rm = TRUE ),
    MSI   = mean((swir1/nir), na.rm = TRUE),
    MSAVI2 = mean(
      (nir * swir1 - swir2 * red) / (nir * swir1 + swir2 * red + 0.001), 
      na.rm = TRUE
    ),
    VMSI = mean(
      ((nir + swir1) - (swir2 + red)) / ((nir + swir1) + (swir2 + red) + 0.001), 
      na.rm = TRUE
    ),
    GSDI = mean(
      (green - swir2) / (green + swir2), 
      na.rm = TRUE
    ),
    EVDI = mean(
      (nir - red - swir2) / (nir + red + swir2), 
      na.rm = TRUE
    ),
    TBMI = mean(
      (swir1 - swir2 - red) / (swir1 + swir2 + red), 
      na.rm = TRUE
    ),
    year  = first(year),
    doy   = first(doy),
    .groups = "drop"
  )

################################################################################
# 3) EXTRACT & PREPARE TREE-RING CHRONOLOGIES
################################################################################

# 3A) EXTRACT TREE-RING DATA
TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/tree_coordinates_SF_12_21_PF_02_03.gpkg"
)



# 3B) BUILD DETRENDED CHRONOLOGY PER PLOT
chronologies_by_plot2 <- TRW_df %>%
  group_by(plot) %>%
  group_modify(~ {
    df_wide <- .x %>%
      select(Year, treeid, ring_width) %>%
      pivot_wider(names_from = treeid, values_from = ring_width) %>%
      arrange(Year)
    
    rownames(df_wide) <- as.character(df_wide$Year)
    rwl_obj <- as.rwl(df_wide[, -1])
    
    rwi_obj <- detrend(rwl_obj, method = "Friedman")
    
    chron_obj <- chron(rwi_obj, prefix = paste0("CH_", unique(.x$plot)))
    
    chron_df <- as.data.frame(chron_obj)
    chron_df$year <- as.numeric(rownames(chron_df))
    chron_df$plot <- unique(.x$plot)
    
    # Apply offset if needed (e.g., +1649)
    chron_df$year <- chron_df$year + 1649
    
    chron_df
  }) %>%
  ungroup()

# Rename the main chronology column to "TRI"
chronologies_by_plot2 <- chronologies_by_plot2 %>%
  rename(TRI = starts_with("std"))

################################################################################
# 4) COMPUTE CUMULATIVE VIs (ROLLING SUMS) FOR DIFFERENT WINDOWS
################################################################################

# 4A) PIVOT DAILY VIs TO LONG FORMAT
SDC_long2 <- SDC_daily %>%
  pivot_longer(
    cols = c(NDVI, EVI, GNDVI, NDMI, NDWI, NMDI, NBR, NBR2, ARVI, OSAVI, NIRv, VMI, MSAVI2, VMSI, GSDI, EVDI, TBMI, MSI),
    names_to = "VI",
    values_to = "value"
  ) %>%
  arrange(plot, VI, date)

# 4B) DEFINE ACCUMULATION WINDOWS (in days)
accum_windows <- 1:24  # These now represent rolling windows in days

# 4C) COMPUTE ROLLING SUMS (PER PLOT & VI, ACROSS TIME)
SDC_cum_continuous2 <- SDC_long2 %>%
  group_by(plot, VI) %>%
  group_modify(~ {
    df <- .x
    cum_df <- map_dfc(accum_windows, function(w) {
      roll_val <- if (nrow(df) >= w) {
        zoo::rollapply(df$value, width = w, FUN = sum, align = "right", fill = NA)
      } else {
        rep(NA_real_, nrow(df))
      }
      tibble(!!paste0("cum", w) := roll_val)
    })
    bind_cols(df, cum_df)
  }) %>%
  ungroup()

# 4D) RESHAPE CUMULATIVE COLUMNS TO LONG FORMAT
SDC_cum_long2 <- SDC_cum_continuous2 %>%
  pivot_longer(
    cols = starts_with("cum"),
    names_to = "window",
    names_prefix = "cum",
    values_to = "cum_value"
  ) %>%
  mutate(window = as.integer(window))

################################################################################
# 5) MERGE CUMULATIVE VIs WITH CHRONOLOGIES & COMPUTE CORRELATIONS
################################################################################

# 5A) MERGE ON "plot" & "year"
# Note: Since TRI is annual, each year's TRI will join with all daily records of that year.
cum_TRI2 <- SDC_cum_long2 %>%
  left_join(chronologies_by_plot, by = c("plot", "year"))

# 5B) FILTER OUT EARLY YEARS IF NECESSARY
cum_TRI_filtered2 <- cum_TRI2 %>%
  dplyr::filter(year >= 2000)  # Adjust the start year as needed

# NEW: ASSIGN FOREST TYPE BASED ON PLOT
# Classification based on provided information:
# critical plots Mischbestände: SF12, SF13, SF14, SF15, SF16, 
cum_TRI_filtered2 <- cum_TRI_filtered2 %>%
  mutate(forest_type = case_when(
    plot %in% c("PF01", "SF03", "SF06", "SF07", "SF10", "SF12", "SF13", "SF15", "SF16", "SF17", "SF21") ~ "coniferous",
    plot %in% c("PF02", "PF03" , "SF02", "SF04", "SF05", "SF08", "SF09",  "SF14", "SF18", "SF19", "SF20") ~ "deciduous",
    TRUE ~ NA_character_
  ))

# 5C) COMPUTE CORRELATION PER (VI, window, day-of-year, forest_type)
heatmap_df2 <- cum_TRI_filtered2 %>%
  group_by(VI, window, doy, forest_type) %>%
  summarize(
    correlation = if (sum(!is.na(cum_value) & !is.na(TRI)) > 1) {
      cor(cum_value, TRI, use = "complete.obs")
    } else {
      NA_real_
    },
    n = n(),
    .groups = "drop"
  )

################################################################################
# 6) PLOT HEATMAP (WITHOUT INTERPOLATION)
################################################################################
heatmap_df_coniferous <- heatmap_df[heatmap_df$forest_type == "coniferous",]
heatmap_df_deciduous <- heatmap_df[heatmap_df$forest_type == "deciduous",]
ggplot(heatmap_df_deciduous, aes(x = doy, y = window, fill = correlation)) +
  geom_tile() +
  facet_wrap(~ VI, ncol = 4) +
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +
  scale_y_continuous(breaks = seq(min(heatmap_df$window, na.rm = TRUE),
                                  max(heatmap_df$window, na.rm = TRUE), by = 5)) +
  scale_fill_gradientn(
    colors = c("navy", "blue", "cyan", "yellow", "red"),
    name = "Correlation"
  ) +
  labs(
    x = "Day of Year", 
    y = "Accumulation Window (days)",
    title = "Correlation between Cumulative Vegetation Indices and Detrended TRW"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggplot(heatmap_df, aes(x = doy, y = window, fill = correlation)) +
  geom_tile() +
  facet_grid(forest_type ~ VI) +
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +
  scale_y_continuous(breaks = seq(min(heatmap_df$window, na.rm = TRUE),
                                  max(heatmap_df$window, na.rm = TRUE), by = 5)) +
  scale_fill_gradientn(
    colors = c("navy", "blue", "cyan", "yellow", "red"),
    name = "Correlation"
  ) +
  labs(
    x = "Day of Year", 
    y = "Accumulation Window (days)",
    title = "Correlation between Cumulative Vegetation Indices and Detrended TRW by Forest Type"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

