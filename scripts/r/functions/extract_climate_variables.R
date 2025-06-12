# ─────────────────────────────────────────────────────────────────────────────
# extract_plot_climate()
# ---------------------------------------------------------------------------
#  • points_sf : sf POINT layer with a column that identifies the plot
#  • var_files : named vector/list  c(tmax = "...tif", ppt = "...tif", …)
#  • buffer_m  : buffer radius in metres around each plot (default 30)
#  • years     : numeric vector of calendar years to keep (default 2000:2022)
#  • fun       : aggregation function across pixels (default mean)
#
# returns: one tibble  plot | year | month | <var1> | <var2> | …
# ---------------------------------------------------------------------------
extract_plot_climate <- function(points_sf,
                                 var_files,
                                 buffer_m   = 1,
                                 years      = 2000:2022,
                                 fun        = mean,
                                 na.rm      = TRUE,
                                 start_date = as.Date("1991-01-01"),
                                 plot_col   = "plot") {
  
  # ---- 0. libs (quietly) ---------------------------------------------------
  
  # ---- 1. polygon buffers per plot ----------------------------------------
  plots_sf <- points_sf |>
    group_by(.data[[plot_col]]) |>
    summarise(geometry = st_union(geom), .groups = "drop") |>
    st_buffer(buffer_m)
  
  # re-project polygons to match the first raster
  plots_v  <- terra::vect(plots_sf)                        # SpatVector
  
  # ---- 2. band → calendar date lookup -------------------------------------
  n_bands      <- terra::nlyr(terra::rast(var_files[[1]])) # should be 408
  band_dates   <- seq(start_date, by = "month", length.out = n_bands)
  
  keep_idx     <- which(lubridate::year(band_dates) %in% years)
  band_dates   <- band_dates[keep_idx]                     # 276 dates
  
  # ── 3. loop through each climate raster ─────────────────────────────────────
  climate_plot_month <- NULL   # will grow with each variable
  
  for (i in seq_along(var_files)) {
    
    varname <- names(var_files)[i]
    path    <- var_files[[i]]
    
    message(sprintf("Processing %s (%d of %d)", varname, i, length(var_files)))
    
    # 3A. read & subset raster to 2000-22
    r <- terra::rast(path)[[keep_idx]]
    names(r) <- as.character(seq_len(nlyr(r)))   # "1", "2", … "276"
    
    
    # 3B. extract mean inside each plot buffer
    raw <- terra::extract(r, plots_v, fun = mean, na.rm = TRUE, ID = FALSE)
    raw$plot <- unlist(as.vector(unname(plots_v[[plot_col]])))      # attach plot IDs
    
    # 3C. tidy: layers → long rows, add calendar fields
    tbl <- raw |>
      pivot_longer(-plot, names_to = "band", values_to = varname) |>
      mutate(
        band  = as.integer(band),     # now 1…276
        date  = band_dates[band],
        year  = year(date),
        month = month(date)
      ) |>
      dplyr::select(plot, year, month, all_of(varname))
    
    # 3D. merge with previously processed variables
    if (is.null(climate_plot_month)) {
      climate_plot_month <- tbl                     # first variable sets it up
    } else {
      climate_plot_month <- dplyr::full_join(       # subsequent vars are merged in
        climate_plot_month,
        tbl,
        by = c("plot", "year", "month")
      )
    }
  }
  climate_plot_month <- climate_plot_month %>%
    mutate(tmax = tmax / 10)
  return(climate_plot_month)
}
