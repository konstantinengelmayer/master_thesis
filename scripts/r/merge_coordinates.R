library(terra)
library(dplyr)

# 1. Read in all the geopackage files
files <- list.files("data/vector_data/MW3-Koordinaten_final", full.names = TRUE)
list_vec <- lapply(files, vect)

v <- list_vec[[1]]

# --------------------------------------------------------------------
# 1. helper that cleans a single SpatVector --------------------------
# --------------------------------------------------------------------
clean_vec <- function(v, target_crs = "EPSG:25832") {
  
  # --- harmonise column names ---------------------------------------
  nms <- names(v)
  rename <- c(
    TreeID      = "tree_id",    # synonyms for ID ...
    TreeId      = "tree_id",
    treeid      = "tree_id",
    PT_NAME     = "tree_id",
    
    Plot        = "plot",       # synonyms for plot ...
    plot        = "plot",
    
    Treespecies = "species",    # synonyms for species ...
    Species     = "species"
  )
  present <- intersect(names(v), names(rename))          # which synonyms exist?
  for (old in present) {
    names(v)[names(v) == old] <- rename[old]             # overwrite with target
  }
  
  # --- keep only the three wanted attributes ------------------------
  v <- v[ , c("tree_id", "plot", "species"), drop = FALSE]
  
  # extract the attribute table
  vals <- values(v)
  
  # convert every field to character
  vals <- data.frame(lapply(vals, as.character),
                     stringsAsFactors = FALSE)
  
  # put it back
  values(v) <- vals
  
  # --- project to common CRS ----------------------------------------
  v <- project(v, target_crs)
  v
}

# --------------------------------------------------------------------
# 2. apply to every element in list_vec ------------------------------
# --------------------------------------------------------------------
cleaned <- lapply(list_vec, clean_vec)

# --------------------------------------------------------------------
# 3. merge and write -------------------------------------------------
# --------------------------------------------------------------------
all_plots <- do.call(rbind, cleaned)
writeVector(all_plots,
            filename = "trees_all_plots.gpkg",
            filetype = "GPKG",
            layer    = "trees_all_plots",
            overwrite = TRUE)
