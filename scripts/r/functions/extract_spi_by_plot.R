extract_spi_by_plot <- function(points_path, spi_path,
                                  buffer = 10,
                                  years = 2000:2022,
                                  month_of_interest = 12) {

  
  # ---------------------------------------------------------------
  # 1. Load data --------------------------------------------------
  spi <- rast(spi_path)            # SpatRaster of SPI layers
  pts <- vect(points_path)         # SpatVector of tree points
  
  # ---------------------------------------------------------------
  # 3. Create buffer and dissolve by plot -------------------------
  buf      <- buffer(pts, width = buffer)
  plots    <- aggregate(buf, "plot")     # dissolve by plot
  
  # ---------------------------------------------------------------
  # 4. Select SPI layers for given month and years ---------------
  ras_dates <- time(spi)
  sel_idx   <- which(month(ras_dates) == month_of_interest &
                       year(ras_dates) %in% years)
  spi_sel   <- spi[[sel_idx]]
  names(spi_sel) <- year(ras_dates[sel_idx])
  
  # ---------------------------------------------------------------
  # 5. Extract mean SPI per plot polygon --------------------------
  spi_ext <- terra::extract(
    spi_sel,
    plots,
    fun   = mean,
    na.rm = TRUE,
    bind = TRUE
  )
  
  # ---------------------------------------------------------------
  # 6. Clean and return -------------------------------------------
  names(spi_ext)[1] <- "plot_id"
  spi_ext <- as.data.frame(spi_ext)
  names(spi_ext) <- gsub("^X", "", names(spi_ext))
  spi_ext <- spi_ext[, !(names(spi_ext) %in% c("agg_n", "area", "Species", "treeid"))]
  
  spi_plot_long <- spi_ext |>
    pivot_longer(
      cols      = -plot_id,
      names_to  = "year",
      values_to = "SPI12_Dec"
    ) |>
    mutate(year = as.integer(year)) |>
    arrange(plot_id, year)
  
  
  
  return(spi_plot_long)
}

