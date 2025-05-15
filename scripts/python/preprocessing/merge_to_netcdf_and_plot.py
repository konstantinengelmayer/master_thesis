# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:10:48 2025

@author: konst
"""

import os
import glob
import numpy as np
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
from datetime import datetime

# 1) Adjust the working directory as appropriate
os.chdir("C:/Users/konst/Documents/Hiwi/mw3/master_thesis")

# folders = ["HLS_Sentinel_Kellerwald", "HLS_Landsat_Kellerwald", "Landsat_7_strictest_Kellerwald"]  # or other folders
folders = ["landsat_5"]
for folder in folders:
    print(f"Processing folder: {folder}")
    
    # 2) Find all .tif files
    tif_pattern = f"data/satellite_data/eifel/{folder}/*.tif"
    files = glob.glob(tif_pattern)
    if not files:
        print(f"No .tif files found for {folder}")
        continue
    
    # 3) Extract date strings from filenames
    date_strings = [os.path.basename(f).split("_")[-1].split(".")[0] for f in files]
    dates = [datetime.strptime(ds, "%Y%m%d").strftime("%Y-%m-%d") for ds in date_strings]

    # 4) Open each raster, store in a list
    rasters = [rxr.open_rasterio(f) for f in files]

    # 5) Concatenate along 'time' dimension
    # Shape becomes (time, band, y, x)
    stack = xr.concat(rasters, dim="time")
    stack = stack.assign_coords(time=("time", dates))

    # 6) Label the band dimension
    #    Adjust the band names to match your actual band order.
    band_names = ["Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"]
    stack = stack.assign_coords(band=band_names[: stack.sizes["band"]])

    # 7) Convert to a Dataset, splitting each band into its own variable
    ds = stack.to_dataset(dim="band")
    # Optionally rename them exactly as you like
    # ds = ds.rename({"Blue": "Blue", "Green": "Green", "Red": "Red", ...})

    # 8) Save as NetCDF
    nc_output = f"data/satellite_data/eifel/{folder}/{folder}_2015_2024.nc"
    ds.to_netcdf(nc_output)
    print(f"NetCDF saved: {nc_output}")
    
    # 9) Create figures
    fig_dir = f"figures/satellite_raw/eifel/{folder}/"
    os.makedirs(fig_dir, exist_ok=True)

    # 10) Compute global min and max as single scalars across all times + channels
    #     (You can pick your own percentile or do it band-by-band if you prefer.)
    #     We'll make one global min, one global max across Red,Green,Blue.
    rgb_da = ds[["Red", "Green", "Blue"]].to_array()  # shape: (band, time, y, x)
    # We flatten all dims to find a single 2% and 98% value
    low = float(rgb_da.quantile(0.02))
    high = float(rgb_da.quantile(0.98))

    def normalize_to_01(da, vmin, vmax):
        """
        Convert an xarray DataArray to [0,1] range using vmin and vmax as scalars.
        """
        arr = da.values  # convert to NumPy
        arr = (arr - vmin) / (vmax - vmin)
        arr = np.clip(arr, 0, 1)
        # Return as a DataArray so we keep dims & coords
        return xr.DataArray(arr, dims=da.dims, coords=da.coords, attrs=da.attrs)

    # 11) Plot each time step
    for i in range(ds.sizes["time"]):
        scene = ds.isel(time=i)
        date_str = str(scene.time.values)[:10]  # "YYYY-MM-DD"
        
        # Extract R,G,B as a single DataArray with dim "band"
        # shape => (band, y, x)
        rgb = xr.concat([scene["Red"], scene["Green"], scene["Blue"]], dim="band")
        
        # Normalize
        rgb = normalize_to_01(rgb, low, high)
        
        # Optional gamma
        gamma = 0.8
        rgb.values = rgb.values ** gamma

        # Prepare for plotting => (y, x, band)
        rgb_for_plot = rgb.transpose("y", "x", "band").values
        
        plt.figure(figsize=(10, 10 * (ds.sizes["y"] / ds.sizes["x"])))
        plt.imshow(rgb_for_plot)
        plt.axis("off")
        plt.title(f"{folder.capitalize()} - {date_str}")
        
        out_path = os.path.join(fig_dir, f"{folder}_{date_str}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        print(f"Saved figure: {out_path}")