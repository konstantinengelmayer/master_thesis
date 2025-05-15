# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 10:32:04 2025

@author: konst
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from dask.distributed import Client

# Start a Dask client for parallel processing
client = Client()  # This will start a local Dask cluster
print(client)  # View the Dask dashboard link for monitoring progress

# Load datasets with Dask chunks
hlsl = xr.open_dataset("data/satellite_data/eifel/hls_landsat/hls_landsat_2013_2025.nc", chunks={'time': 100})
hlss = xr.open_dataset("data/satellite_data/eifel/hls_sentinel/hls_sentinel_2015_2024.nc", chunks={'time': 100})
landsat7 = xr.open_dataset("data/satellite_data/eifel/landsat_7_stricter/landsat_7_stricter_2000_2023.nc", chunks={'time': 100})

# Ensure the time coordinates are in the same format (e.g., datetime64)
hlsl['time'] = pd.to_datetime(hlsl['time'])
hlss['time'] = pd.to_datetime(hlss['time'])
landsat7['time'] = pd.to_datetime(landsat7['time'])

# Concatenate the datasets along the time dimension
combined = xr.concat([hlsl, hlss, landsat7], dim='time')

# Sort by time to ensure proper resampling
combined = combined.sortby('time')

# Resample to monthly frequency and calculate the median
monthly_median = combined.resample(time='1M').median(dim='time')

# Compute the result (this triggers the actual computation)
monthly_median = monthly_median.compute()
monthly_median = combined.resample(time='1M').median(dim='time')

# Load the shapefile containing the points
points = gpd.read_file("data/vector_data/Measured_Tree_coordinates.gpkg")  # Replace with the path to your shapefile

# Ensure the shapefile has the same CRS as the dataset (e.g., EPSG:4326 for WGS84)
points = points.to_crs(epsg=4326)  # Adjust the EPSG code if necessary

# Extract time series for each point
time_series = []
for idx, point in points.iterrows():
    # Extract the time series for the point using xarray's .sel() method
    x, y = point.geometry.x, point.geometry.y
    ts = monthly_median.sel(x=x, y=y, method="nearest")  # Use nearest neighbor interpolation
    time_series.append(ts)

# Plot the time series for the points
plt.figure(figsize=(12, 6))
for i, ts in enumerate(time_series):
    plt.plot(ts['time'], ts['NIR'], label=f'Point {i+1}')  # Replace 'your_variable_name' with the actual variable name

# Customize the plot
plt.title('Time Series of Monthly Median Mosaic for Points')
plt.xlabel('Time')
plt.ylabel('Value')  # Replace with the appropriate variable name
plt.legend()
plt.grid()
plt.show()
