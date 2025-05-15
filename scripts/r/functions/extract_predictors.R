#───────────────────────────────────────────────────────────────
# extract_predictors()
#   points_path  : vector file (GPKG/Shp) containing tree points.
#   raster_path  : single- or multi-band predictor raster.
#   buffer_width : radius around each point (same units as data;
#                  default = 10).
#───────────────────────────────────────────────────────────────
extract_predictors <- function(points_path,
                               raster_path,
                               buffer_width = 10) {
  
  #----- fixed internal settings -------------------------------------------
  plot_col    <- "plot"      # must exist in points file
  summary_fun <- mean        # pixel aggregation statistic
  keep_na     <- FALSE       # NA pixels discarded
  
  # 1 ─ read data -----------------------------------------------------------
  pts <- terra::vect(points_path)
  ras <- terra::rast(raster_path)
  
  # 2 ─ buffer (10 m by default) & dissolve --------------------------------
  buf10   <- terra::buffer(pts, width = buffer_width)
  bufPlot <- terra::aggregate(buf10, by = plot_col)
  
  # plain character vector of plot names (row-order = bufPlot)
  plot_names <- terra::values(bufPlot, df = TRUE)[[plot_col]]
  
  # 3 ─ extract raster values & summarise -----------------------------------
  plot_means <- terra::extract(
    ras, bufPlot,
    fun   = summary_fun,
    na.rm = !keep_na          # layer names become column names
  ) %>%
    mutate( !!plot_col := plot_names ) %>%
    relocate(!!plot_col) %>%          # make 'plot' the first column
    select(-any_of("ID"))             # drop the helper column
  
  return(plot_means)
}
