extract_SDC_buffered <- function(folder_paths,
                                 points_path,
                                 buffer_m = 10,
                                 method = c("original", "sgolay", "rolling7", "loess")) {
  # Load required libraries
 
  
  # Ensure 'method' is one of the allowed values:
  method <- match.arg(method)
  
  # Read the points vector file and create an 'area' column if needed.
  # (Here we assume the points file already has an 'area' column.)
  points <- terra::vect(points_path)
  points_sf <- sf::st_as_sf(points)
  
  if (!"area" %in% names(points_sf)) {
    stop("Points file must have an 'area' column to match folder paths.")
  }
  
  # Use the first folder's first TIFF file as a reference for CRS projection.
  first_folder <- folder_paths[1]
  tif_files_sample <- list.files(path = first_folder, pattern = "\\.tif$", full.names = TRUE)
  if (length(tif_files_sample) == 0) {
    stop("No TIFF files found in the first folder provided.")
  }
  sample_rast <- terra::rast(tif_files_sample[1])
  
  # Reproject points to match the raster CRS
  points <- terra::project(points, sample_rast)
  points_sf <- sf::st_as_sf(points)
  
  # Prepare a list to store extracted data across areas
  data_records <- list()
  
  # Loop over each unique area in the points file
  unique_areas <- unique(points_sf$area)
  for (area in unique_areas) {
    # Find folder(s) that contain the area name (case-insensitive search)
    matched_folders <- folder_paths[grepl(area, folder_paths, ignore.case = TRUE)]
    if (length(matched_folders) == 0) {
      warning(paste("No folder path found that contains area", area, "- skipping this area."))
      next
    }
    # Use the first matching folder
    folder <- matched_folders[1]
    
    # List TIFF files in the area-specific folder
    tif_files <- list.files(path = folder, pattern = "\\.tif$", full.names = TRUE)
    if (length(tif_files) == 0) {
      warning(paste("No TIFF files found in folder for area", area, "- skipping this area."))
      next
    }
    
    # Filter points for the current area
    area_points_sf <- points_sf[points_sf$area == area, ]
    
    # Group by plot (assuming a column "plot") and buffer points
    plot_polygons <- area_points_sf %>%
      group_by(plot) %>%
      summarize(
        geometry = sf::st_union(sf::st_buffer(geometry, dist = buffer_m)),
        .groups  = "drop"
      ) %>%
      # Convert to SpatVector for terra::extract
      terra::vect()
    
    # Set up a progress bar for the current area
    pb <- txtProgressBar(min = 0, max = length(tif_files), style = 3)
    
    # Process each TIFF file for the current area
    for (i in seq_along(tif_files)) {
      file <- tif_files[i]
      filename <- basename(file)
      
      # Extract date from filename assuming pattern "YYYYDDD.tif" after last underscore
      parts <- strsplit(filename, "_")[[1]]
      date_str <- gsub("\\.tif$", "", tail(parts, 1))
      year <- as.integer(substr(date_str, 1, 4))
      day_of_year <- as.integer(substr(date_str, 5, nchar(date_str)))
      date_val <- as.Date(paste0(year, "-01-01")) + (day_of_year - 1)
      
      # Read the raster and assign band names
      r <- terra::rast(file)
      names(r) <- c("blue", "green", "red", "nir", "swir1", "swir2")
      
      # Extract mean values for each plot polygon
      vals <- terra::extract(r, plot_polygons, fun = mean, na.rm = TRUE)
      
      # Combine the extracted values with the plot ID, date, and area information
      extract_df <- cbind(
        data.frame(
          plot_id = plot_polygons$plot[vals$ID],
          date    = date_val,
          area    = area
        ),
        vals[, c("blue", "green", "red", "nir", "swir1", "swir2")]
      )
      
      # Store the result
      data_records[[length(data_records) + 1]] <- extract_df
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  # Combine all records into a single data frame
  df <- do.call(rbind, data_records)
  
  # Scale the band values by 0.0001
  numeric_cols <- c("blue", "green", "red", "nir", "swir1", "swir2")
  df[, numeric_cols] <- df[, numeric_cols] * 0.0001
  
  # Optionally apply a smoothing method per plot per band if requested
  if (method != "original") {
    smooth_sgolay <- function(x, n = 11, p = 2) {
      if (length(x) < n) rep(NA_real_, length(x)) else signal::sgolayfilt(x, p = p, n = n)
    }
    
    smooth_rolling <- function(x, k = 7) {
      zoo::rollapply(x, width = k, FUN = mean, fill = NA, align = "center")
    }
    
    smooth_loess <- function(x, date_vec, span = 0.1) {
      if (all(is.na(x))) return(x)
      fit <- loess(x ~ date_vec, span = span, na.action = na.exclude)
      predict(fit, newdata = date_vec)
    }
    
    df <- df %>%
      group_by(plot_id, area) %>%
      arrange(date) %>%
      mutate(
        dplyr::across(
          .cols = all_of(numeric_cols),
          .fns  = ~ {
            if (method == "sgolay") {
              smooth_sgolay(., n = 11, p = 2)
            } else if (method == "rolling7") {
              smooth_rolling(., k = 7)
            } else if (method == "loess") {
              numeric_date <- as.numeric(date)
              smooth_loess(., numeric_date, span = 0.1)
            } else {
              .
            }
          }
        )
      ) %>%
      ungroup()
  }
  
  return(df)
}
