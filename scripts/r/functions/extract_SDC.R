################################################################################
# Function: extract_point_values
# Author: [Your Name]
# Date:   [Date]
#
# Description:
#   This function reads all .tif files from a specified folder, extracts values 
#   at given point locations for each timestep, and combines the results into 
#   a single data frame. The function also applies a scale factor to the numeric 
#   columns. 
#
# Arguments:
#   folder_path (character): Path to the folder containing .tif files.
#   points_path (character): Path to the vector file (.gpkg, .shp, etc.) 
#                           containing the points of interest.
#
# Returns:
#   A data frame with the following columns:
#     - date:       The date extracted from the TIFF filename 
#                   (assumes 'yyyyddd' format after last underscore).
#     - point_id:   A unique integer ID for each point.
#     - red:        Extracted red band values (scaled).
#     - blue:       Extracted blue band values (scaled).
#     - green:      Extracted green band values (scaled).
#     - nir:        Extracted near-infrared band values (scaled).
#     - swir1:      Extracted shortwave-infrared band 1 values (scaled).
#     - swir2:      Extracted shortwave-infrared band 2 values (scaled).
#
# Dependencies:
#   - terra
#   - sf
#
# Example usage:
#   df_extracted <- extract_point_values(
#       folder_path = "path/to/tiff/folder",
#       points_path = "path/to/points.gpkg"
#   )
################################################################################

extract_SDC <- function(folder_paths, points_path) {
  #— load libraries
  library(terra)
  library(sf)
  
  #— read points as sf and check required columns
  points_sf <- sf::st_read(points_path, quiet=TRUE)
  if (!"area" %in% names(points_sf)) {
    stop("Points file must have an 'area' column to match folder paths.")
  }
  if (!"tree_id" %in% names(points_sf)) {
    stop("Points file must have a 'tree_id' column for unique IDs.")
  }
  
  #— get CRS from a sample raster in the first folder
  sample_tifs <- list.files(folder_paths[1], "\\.tif$", full.names=TRUE)
  if (length(sample_tifs)==0) {
    stop("No TIFF files found in the first folder provided.")
  }
  sample_rast <- terra::rast(sample_tifs[1])
  
  #— reproject points once to match the raster CRS
  points_proj_sf <- sf::st_transform(points_sf, crs(sample_rast))
  #— convert to terra vector for terra::extract
  points_proj_vec <- terra::vect(points_proj_sf)
  
  #— prepare container for all results
  all_records <- list()
  
  #— loop over each area, subset points, find matching folder, extract
  for (area in unique(points_proj_sf$area)) {
    # subset the points for this area
    pts_mask <- points_proj_sf$area == area
    if (!any(pts_mask)) next
    
    pts_area_vec <- points_proj_vec[pts_mask, ]
    
    # find the folder whose name contains the area (case‐insensitive)
    matched <- folder_paths[grepl(area, folder_paths, ignore.case=TRUE)]
    if (length(matched)==0) {
      warning("No folder matching area '", area, "' – skipping.")
      next
    }
    tif_files <- list.files(matched[1], "\\.tif$", full.names=TRUE)
    if (length(tif_files)==0) {
      warning("No TIFFs in folder for area '", area, "' – skipping.")
      next
    }
    
    
    # extract for each TIFF
    for (tif in tif_files) {
      fn <- basename(tif)
      # parse yyyyddd from filename (assumes last underscore + extension)
      date_str <- sub("\\.tif$", "", tail(strsplit(fn, "_")[[1]], 1))
      year <- as.integer(substr(date_str, 1, 4))
      doy  <- as.integer(substr(date_str, 5, nchar(date_str)))
      date_val <- as.Date(sprintf("%04d-01-01", year)) + (doy - 1)
      
      # read & name bands
      r <- terra::rast(tif)
      names(r) <- c("blue","green","red","nir","swir1","swir2")
      
      # extract and build a small data.frame
      vals <- terra::extract(r, pts_area_vec)
      df <- data.frame(
        date     = rep(date_val, nrow(vals)),
        point_id = pts_area_vec$tree_id,
        vals
      )
      
      all_records[[length(all_records)+1]] <- df
    }
  }
  
  #— combine, scale spectral columns, and return
  if (length(all_records)==0) {
    warning("No data extracted for any area.")
    return(NULL)
  }
  out <- do.call(rbind, all_records)
  
  # identify only the spectral columns (everything except date & point_id)
  spec_cols <- setdiff(names(out), c("date","point_id"))
  out[spec_cols] <- out[spec_cols] * 0.0001
  
  return(out)
}

