# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 08:14:21 2025

Description:
    - 
    -
    
Required User Input:
    -
    -
    
@author: konst
"""
import os
from glob import glob
from datetime import datetime, timedelta
from tqdm import tqdm
import geopandas as gpd
import pandas as pd
import rioxarray

# Define the folder containing the .tif files and get their paths
folder_path = "data/satellite_data/SDC/lindenberg_eifel_koenigsforst"
tif_files = glob(os.path.join(folder_path, "*.tif"))

# Define the band names (in the same order as they appear in the TIFF file)
band_names = ['red', 'blue', 'green', 'nir', 'swir1', 'swir2']

# Read the points from the geopackage
points = gpd.read_file("data/vector_data/Measured_Tree_coordinates.gpkg")

# Uncomment and set the CRS if you need to reproject the points to match the TIFF's CRS:
# points = points.to_crs("EPSG:your_raster_crs_here")

# Create an empty list to store extracted records
data_records = []

# Open one raster to check its CRS
sample_da = rioxarray.open_rasterio(tif_files[0])
print("Raster CRS:", sample_da.rio.crs)
print("Points CRS before reproject:", points.crs)

if points.crs != sample_da.rio.crs:
    points = points.to_crs(sample_da.rio.crs)
    print("Reprojected points CRS:", points.crs)

for file in tqdm(tif_files, desc="Processing files"):
    # Extract the filename (e.g., "CSDC30_32UMB_2022361.tif")
    filename = os.path.basename(file)
    # Extract the date part from the filename, e.g., "CSDC30_32UMB_2022361.tif" -> "2022361"
    date_str = filename.split('_')[-1].replace('.tif', '')
    year = int(date_str[:4])
    day_of_year = int(date_str[4:])
    # Convert year and day-of-year to a datetime object
    date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
    
    # Open the TIFF file with rioxarray; result is a DataArray with dims: (band, y, x)
    da = rioxarray.open_rasterio(file)
    
    # Assign meaningful band names to the 'band' coordinate
    da = da.assign_coords(band=band_names)
    
    # Loop through each point in your GeoDataFrame
    for idx, row in points.iterrows():
        # Get the point's x, y coordinates
        x, y = row.geometry.x, row.geometry.y
        
        # Use nearest neighbor selection to extract values at the point location
        sample = da.sel(x=x, y=y, method="nearest")
        # sample now contains a value for each band
        
        # Build a dictionary with the date and each band value.
        # Here, 'point_id' is either an "id" field in your GeoDataFrame or the row index.
        record = {
            "date": date,
            "point_id": row.get("id", idx)
        }
        # Add each band value to the record dictionary
        for band in band_names:
            # Ensure the value is a native Python type (e.g., float)
            record[band] = float(sample.sel(band=band).values)
        
        data_records.append(record)

# Create a DataFrame from the list of records
df = pd.DataFrame(data_records)
print(df.head())

scale_factor = 0.0001
for band in ['red', 'blue', 'green', 'nir', 'swir1', 'swir2']:
    df[band] = df[band] * scale_factor

print(df.head())


import pandas as pd
import matplotlib.pyplot as plt

# Assuming df is your DataFrame with the columns: date, point_id, red, blue, green, nir, swir1, swir2
# Make sure the date column is in datetime format
df['date'] = pd.to_datetime(df['date'])

# Calculate NDVI using the formula: NDVI = (nir - red) / (nir + red)
df['ndvi'] = (df['nir'] - df['red']) / (df['nir'] + df['red'])

# Get the unique point ids
point_ids = df['point_id'].unique()

# Loop through each point id and plot the NDVI time series
for pid in point_ids:
    point_data = df[df['point_id'] == pid].sort_values('date')
    
    plt.figure(figsize=(10, 5))
    plt.plot(point_data['date'], point_data['nir'], marker='o', linestyle='-')
    plt.title(f"NDVI Time Series for Point {pid}")
    plt.xlabel("Date")
    plt.ylabel("NIR")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    print(df.head())

