################################################################################
# Function: extract_dendro_data
# Author: [Your Name]
# Date:   [Date]
#
# Description:
#   This function reads all .xlsx files from a specified folder, matches them 
#   to a vector of measured tree points, and extracts the tree ring width data 
#   for each tree within each plot. The function assumes that each Excel file 
#   contains a 'Year' column (in the first column) and columns named after 
#   each tree ID. 
#
# Arguments:
#   xls_folder (character): Path to the folder containing .xlsx files.
#   points_path  (character): Path to the vector file (.gpkg, .shp, etc.) 
#                            containing the measured tree coordinates. 
#                            Must contain at least 'plot' and 'tree_id' columns.
#
# Returns:
#   A data frame with the following columns:
#     - Year:   The year extracted from the Excel sheets (first column).
#     - ring_width: Tree ring width value for a given tree and year.
#     - plot:   The plot ID (e.g., "pf01").
#     - tree_id: The tree ID (numeric or integer).
#
# Dependencies:
#   - terra
#   - sf
#   - dplyr
#   - readxl
#
# Example usage:
#   df_ringwidth <- extract_dendro_data(
#       xls_folder = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
#       points_path  = "data/vector_data/Measured_Tree_coordinates.gpkg"
#   )
################################################################################

extract_TRW <- function(xls_folder, points_path, crop = FALSE) {
  
  # --- 1. Read your points file ---
  #     Each row in 'points' must have a 'plot' and a 'tree_id' column.
  points <- vect(points_path)
  df_points <- as.data.frame(points)
  
  # --- 2. List all Excel files in the folder ---
  excel_files <- list.files(xls_folder, pattern = "\\.xlsx$", full.names = TRUE)
  
  # --- 3. Get unique combinations of (plot, tree_id) from your points ---
  plot_tree_combos <- df_points %>%
    distinct(plot, tree_id)
  
  # --- 4. Initialize a data frame to store all the results ---
  final_data <- data.frame()
  
  # --- 5. Loop over each unique plot ---
  unique_plots <- unique(plot_tree_combos$plot)
  
  for (p in unique_plots) {
    
    # Convert the plot name to uppercase, e.g., "pf01" -> "PF01"
    # This helps match files named "PF01-xxx.xlsx" or containing "PF01" in the filename.
    p_upper <- toupper(p)
    
    # --- 5a. Find Excel files that match this plot. 
    #         We'll do a simple grepl search for the plot ID in the filename.
    candidate_file <- excel_files[grepl(p_upper, basename(excel_files))][1]
    
    if (is.na(candidate_file)) {
      message(paste0("No matching Excel file found for plot: ", p_upper))
      next
    }
    
    # --- 5b. Read the Excel file
    df_xl <- read_excel(candidate_file)
    
    # Ensure the first column is the Year (or rename accordingly).
    # If your file is indeed structured with the Year in the first column:
    colnames(df_xl)[1] <- "Year"
    
    # --- 5c. Find all tree IDs that belong to this plot
    trees_in_plot <- plot_tree_combos %>% 
      dplyr::filter(plot == p)
    
    # --- 5d. Loop over each tree in that plot
    for (t_id in trees_in_plot$tree_id) {
      
      # The Excel column name is assumed to match the tree ID (as a string).
      col_name <- as.character(t_id)
      
      # If the column does not exist, skip
      if (!col_name %in% colnames(df_xl)) {
        message(paste0("Tree ID ", col_name, " not found in file: ", basename(candidate_file)))
        next
      }
      
      # Extract the Year column + the tree column
      trw_vals <- df_xl %>%
        select(Year, all_of(col_name)) %>%
        rename(ring_width = all_of(col_name))
      
      # Add identifying columns
      trw_vals$plot   <- p
      trw_vals$tree_id <- t_id
      
      # Combine with the master data
      final_data <- bind_rows(final_data, trw_vals)
    }
    
    message(paste0("Processed plot: ", p))
  }
  
  if(crop == TRUE){
    # --- 6. Filter years if desired (e.g., 2000-2022 as in your example) ---
    final_data <- final_data[final_data$Year >= 2000 & final_data$Year <= 2022, ]
  }
 
  
  return(final_data)
}
