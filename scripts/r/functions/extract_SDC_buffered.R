# ---------------------------------------------------------------------------
#  Fast extract (+ optional smoothing) of multi-date, multi-band rasters
#  for plot-centre points • values are **averaged per plot & date**
#
#  Depends on: terra, sf, data.table, signal, zoo, parallel, dplyr
# ---------------------------------------------------------------------------

extract_SDC_points <- function(folder_paths,
                               points_path,
                               method       = c("original", "sgolay",
                                                "rolling7", "loess"),
                               sg_n         = 11,   # Savitzky-Golay window
                               sg_p         = 2,    # Savitzky-Golay poly order
                               roll_k       = 7,    # rolling-mean window
                               loess_span   = 0.1)  # LOESS span
{
  method <- match.arg(method)
  
  ## ── 0•packages & checks ────────────────────────────────────────────────
  pkgs <- c("terra", "sf", "data.table", "signal", "zoo", "parallel", "dplyr")
  for (p in pkgs)
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required but not installed.")
  
  if (!length(folder_paths))
    stop("folder_paths is empty.")
  if (!file.exists(points_path))
    stop("points_path does not exist: ", points_path)
  
  ## ── 1reference raster & points ────────────────────────────────────────
  first_tif <- list.files(folder_paths[1], "\\.tif$", full.names = TRUE)[1]
  if (is.na(first_tif)) stop("No .tif files in the first folder.")
  ref_rast  <- terra::rast(first_tif)
  
  pts <- terra::vect(points_path) |>      # read
    terra::project(ref_rast)         # re-project
  
  if (!all(c("area", "plot") %in% names(pts)))
    stop("Points layer must contain both 'area' and 'plot' columns.")
  
  ## helper: YYYYDDD → Date --------------------------------------------------
  fast_dates <- function(tifs) {
    as.Date(substr(basename(tifs),
                   nchar(basename(tifs)) - 10,
                   nchar(basename(tifs)) - 4),
            format = "%Y%j")
  }
  
  ## ── 2split points by area ─────────────────────────────────────────────
  pts_by_area <- split(pts, pts$area)
  
  ## ── 3loop over areas ─────────────────────────────────────────────────
  results <- vector("list", length(pts_by_area))
  names(results) <- names(pts_by_area)
  
  bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
  
  for (ar in names(pts_by_area)) {
    
    folder <- folder_paths[grepl(ar, folder_paths, ignore.case = TRUE)][1]
    if (is.na(folder)) next
    
    tif_paths <- list.files(folder, "\\.tif$", full.names = TRUE)
    if (!length(tif_paths)) next
    
    dates <- fast_dates(tif_paths)
    n_img <- length(tif_paths)
    
    ras <- terra::rast(tif_paths)                       # load all images
    names(ras) <- paste0(rep(bands, n_img), "_",
                         rep(format(dates, "%Y%j"), each = 6))
    
    ## exact pixel values at each point -------------------------------------
    vals <- terra::extract(ras,
                           pts_by_area[[ar]],
                           ID   = TRUE,    # row index of pts
                           bind = FALSE)   # lean output (matrix)
    
    vals_dt <- data.table::as.data.table(vals)
    vals_dt[, plot_id := pts_by_area[[ar]]$plot[ID]]
    vals_dt[, ID      := NULL]                       # drop helper
    
    ## tidy long → wide, **aggregate by plot/date (mean)** ------------------
    long <- data.table::melt(vals_dt,
                             id.vars        = "plot_id",
                             variable.name  = "band_date",
                             value.name     = "value",
                             variable.factor = FALSE)
    
    long[, c("band", "date_str") := tstrsplit(band_date, "_", fixed = TRUE)]
    long[, date := as.Date(date_str, format = "%Y%j")]
    
    wide <- data.table::dcast(long,
                              plot_id + date ~ band,
                              value.var      = "value",
                              fun.aggregate  = mean)  # average across points
    
    wide[, area := ar]
    results[[ar]] <- wide
  }
  
  df <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  
  ## ── 4scale reflectances ───────────────────────────────────────────────
  df[, (bands) := lapply(.SD, `*`, 0.0001), .SDcols = bands]
  
  ## ── 5optional temporal smoothing ─────────────────────────────────────
  if (method != "original") {
    
    sgolay <- function(x) if (length(x) < sg_n) NA_real_ else
      signal::sgolayfilt(x, p = sg_p, n = sg_n)
    roll7  <- function(x) zoo::rollmean(x, roll_k, fill = NA)
    loess1 <- function(x, d) {
      if (all(is.na(x))) return(x)
      predict(loess(x ~ d, span = loess_span, na.action = na.exclude), d)
    }
    
    data.table::setorder(df, plot_id, area, date)
    
    df[, (bands) :=
         lapply(.SD,
                function(v)
                  switch(method,
                         sgolay   = sgolay(v),
                         rolling7 = roll7(v),
                         loess    = loess1(v, as.numeric(date)))),
       by = .(plot_id, area),
       .SDcols = bands]
  }
  
  as.data.frame(df)
}
