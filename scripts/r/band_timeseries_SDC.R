library(terra)
library(sf)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(lubridate)
library(zoo)       # for rollapply
library(signal)    # for Savitzky-Golay filter (install if needed)
library(tidyr)

source("scripts/r/functions/extract_SDC.R")

# 1) Define folder and read data
df <- extract_SDC(folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                                   "data/satellite_data/SDC/kellerwald_lahntal"),
                  points_path <- "data/vector_data/trees_all_plots.gpkg")

# 2) Basic data prep
df$date <- as.Date(df$date)

df <- df %>%
  mutate(
    year   = year(date),
    period = ifelse(year <= 2012, "2000-2012", "2013+")
  )

bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
point_ids <- unique(df$point_id)

#-----------------------------
# Helper function to apply several smoothing methods
#-----------------------------
apply_smoothing_methods <- function(date_vec, value_vec) {
  # Make sure our input is in ascending date order
  # (rollapply etc. assumes chronological order)
  tmp_df <- data.frame(date = date_vec, value = value_vec)
  tmp_df <- tmp_df[order(tmp_df$date), ]
  
  # 1) Rolling mean (window=7)
  roll7 <- zoo::rollapply(tmp_df$value, width = 7, FUN = mean, align = "center", fill = NA)
  
  # 2) Rolling mean (window=21)
  roll21 <- zoo::rollapply(tmp_df$value, width = 21, FUN = mean, align = "center", fill = NA)
  
  # 3) LOESS smoothing
  #    For daily data over many years, you might want a smaller or larger span
  loess_fit <- loess(value ~ as.numeric(date), data = tmp_df, span = 0.05, na.action = na.exclude)
  loess_vals <- predict(loess_fit, newdata = tmp_df)
  
  # 4) Savitzky-Golay smoothing
  #    We pick a window length and polynomial order.
  #    Example: window = 11, poly = 2.
  #    The length must be odd, so 11 is typical for daily data.
  sg_window <- 11
  if (nrow(tmp_df) >= sg_window) {
    sg_vals <- signal::sgolayfilt(tmp_df$value, p = 2, n = sg_window)
  } else {
    sg_vals <- rep(NA_real_, nrow(tmp_df))
  }
  
  # Combine results
  out_df <- data.frame(
    date      = tmp_df$date,
    original  = tmp_df$value,
    roll7     = roll7,
    roll21    = roll21,
    loess     = loess_vals,
    savitzky  = sg_vals
  )
  return(out_df)
}

#-----------------------------
# Main loop: For each point -> For each band -> Smooth & Plot
#-----------------------------
for (point in point_ids) {
  
  # Filter data for the current point_id
  df_point <- df %>% dplyr::filter(point_id == point)
  
  # Create a list to store each band’s final ggplot
  plot_list <- list()
  
  for (band in bands) {
    # Prepare a small data frame with date, period, and raw band
    df_band <- df_point %>%
      select(date, period, all_of(band))
    
    colnames(df_band)[3] <- "value"  # rename the band column to a standard "value"
    
    # Apply smoothing methods
    smoothed_df <- apply_smoothing_methods(df_band$date, df_band$value)
    # Add 'period' info back in
    smoothed_df <- smoothed_df %>%
      left_join(df_band %>% select(date, period), by = "date")
    
    # Pivot to long format: columns = [date, period, method, value]
    long_df <- smoothed_df %>%
      tidyr::pivot_longer(
        cols      = c("original", "savitzky"),
        names_to  = "method",
        values_to = "smoothed_value"
      )
    
    # Plot: color by smoothing method, facet by period (optional)
    p <- ggplot(long_df, aes(x = date, y = smoothed_value, color = method)) +
      geom_line() +
      scale_color_manual(
        values = c("original" = "black", 
                   "roll7"    = "blue", 
                   "roll21"   = "red", 
                   "loess"    = "green",
                   "savitzky" = "purple")
      ) +
      ggtitle(paste("Point", point, "-", band, "Smoothing Comparison")) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom"
      )
    
    plot_list[[band]] <- p
  }
  
  # Arrange the six band plots in a grid (3x2) and save
  grid_obj <- grid.arrange(grobs = plot_list, ncol = 2, nrow = 3,
                           top = paste("Time Series Smoothing - Point", point))
  
  ggsave(
    filename = paste0("figures/timeseries_smoothing/point_nr_", point, ".png"),
    plot     = grid_obj,
    width    = 12,
    height   = 12,
    dpi      = 300
  )
}

source("scripts/r/functions/extract_SDC_buffered.R")
SDC_df <- extract_SDC_buffered(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  buffer_m     = 1,
  method = "original",
)
names(SDC_df)[1] <- "plot"   

SDC_df <- SDC_df %>%
  mutate(
    year   = year(date),
    period = ifelse(year <= 2012, "2000-2012", "2013+")
  )


plot_ids <- unique(SDC_df$plot)


# ── user inputs ────────────────────────────────────────────────────────────────
bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
outdir  <- "figures/timeseries"              # where to write the PNGs
# nice, high-contrast colours for the two periods
period_cols <- c("2000-2012" = "#1f78b4",   # blue
                 "2013+"     = "#e31a1c")   # red
library(patchwork) 
# ── iterate over every site (plot) ─────────────────────────────────────────────
for (pid in unique(SDC_df$plot)) {
  
  # 1. tidy data for this site --------------------------------------------------
  df_site <- SDC_df %>%
    dplyr::filter(plot == pid) %>%
    mutate(period = factor(period, levels = names(period_cols))) %>%
    select(date, period, all_of(bands)) %>%
    pivot_longer(cols   = all_of(bands),
                 names_to  = "band",
                 values_to = "value")
  
  # 2. build six individual plots (one per band) -------------------------------
  plot_list <- lapply(bands, function(b) {
    ggplot(df_site %>% dplyr::filter(band == b),
           aes(x = date, y = value,
               colour = period, group = period)) +
      geom_line() +
      scale_colour_manual(values = period_cols) +
      labs(title = b, x = NULL, y = "Reflectance (DN)") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")          # suppress legends inside panels
  })
  
  
  # 3. arrange them in a 3 × 2 grid with patchwork -----------------------------
  grid   <- wrap_plots(plotlist = plot_list, ncol = 2) +
    plot_annotation(
      title = paste("Original time series – site", pid),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  # 4. save to disk -------------------------------------------------------------
  ggsave(filename = file.path(outdir, paste0("orig_site_", pid, ".png")),
         plot     = grid,
         width    = 14, height = 8, dpi = 300)
}
