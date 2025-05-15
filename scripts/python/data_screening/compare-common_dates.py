# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:46:21 2025

@author: konst
"""
import os
import xarray as xr
import rioxarray
import numpy as np
import matplotlib.pyplot as plt

# Set working directory
os.chdir("C:/Users/konst/Documents/Hiwi/mw3/master_thesis/")

# ---------------------------
# 1) Load datasets
# ---------------------------
hlsl = xr.load_dataset("data/satellite_data/eifel/hls_landsat/hls_landsat_2013_2025.nc")
hlss = xr.load_dataset("data/satellite_data/eifel/hls_sentinel/hls_sentinel_2015_2024.nc")
landsat7 = xr.load_dataset("data/satellite_data/eifel/landsat_7_stricter/landsat_7_stricter_2000_2023.nc")
landsat5 = xr.load_dataset("data/satellite_data/eifel/landsat_5/landsat_5_2015_2024.nc")

# ---------------------------
# 2) Add "satellite" coordinate to each dataset
# ---------------------------
hlsl = hlsl.assign_coords(satellite="HLS_Landsat")
hlss = hlss.assign_coords(satellite="HLS_Sentinel")
landsat7 = landsat7.assign_coords(satellite="Landsat-7")
landsat5 = landsat5.assign_coords(satellite="Landsat-5")

# ---------------------------
# 3) Set CRS for each dataset (replace "EPSG:4326" with the correct CRS if needed)
# ---------------------------
hlss.rio.write_crs("EPSG:4326", inplace=True)
hlsl.rio.write_crs("EPSG:4326", inplace=True)
landsat7.rio.write_crs("EPSG:4326", inplace=True)
landsat5.rio.write_crs("EPSG:4326", inplace=True)

# ---------------------------
# 4) (Optional) Load and apply a landcover mask
# ---------------------------
#landcover = rioxarray.open_rasterio("data/satellite_data/eifel/landcover/ESA_WorldCover_10m_2021_v200_N48E006_Map.tif")
#mask_reprojected = landcover.sel(band=1).rio.reproject_match(hlsl)
#mask_xr = mask_reprojected == 10
#hlsl = hlsl.where(mask_xr)
#hlss = hlss.where(mask_xr)
#landsat7 = landsat7.where(mask_xr)
#landsat5 = landsat5.where(mask_xr)

# ---------------------------
# 5) Exclude specific dates from each dataset
# ---------------------------
hlss_exclude = ["2018-02-14", "2018-03-01", "2024-01-09", "2024-01-19"]
hlsl_exclude = ["2014-01-31", "2017-01-14", "2018-02-18", "2018-02-27", "2019-01-20", "2022-04-03", "2024-01-11", "2024-01-12"]
landsat7_exclude = ["2000-01-24", "2000-03-21", "2000-07-27", "2000-10-22", "2001-01-19", "2001-02-27", "2001-08-31", "2001-11-03", 
                    "2001-11-10", "2002-04-12", "2003-01-09", "2003-01-25", "2003-02-10", "2003-02-26", "2003-03-30", "2004-04-01", 
                    "2004-04-17", "2004-07-22", "2005-02-06", "2006-08-04", "2007-04-17", "2007-06-13", "2008-05-30", "2008-05-30",
                    "2009-04-06", "2009-05-01", "2011-08-18", "2011-11-15", "2012-08-20", "2013-03-25", "2014-01-07", "2014-01-30",
                    "2014-05-22", "2014-06-07", "2014-10-29", "2024-11-23", "2015-04-16", "2015-05-18", "2015-12-12", "2016-04-02", 
                    "2016-11-28", "2016-12-30", "2017-01-22", "2018-01-18", "2018-03-07", "2018-03-30", "2018-09-22", "2018-10-17", 
                    "2019-01-21", "2019-05-20", "2020-08-19", "2020-11-23", "2021-02-11", "2021-05-09", "2021-07-12", "2022-02-21", 
                    "2022-03-02", "2022-04-03", "2022-08-28", "2023-09-20"]
# If needed, define exclusion dates for Landsat-5 (if none, leave as an empty list)
landsat5_exclude = [
    "2000-01-07", "2001-06-27", "2002-01-05", "2002-03-17", 
    "2003-09-28", "2004-02-28", "2004-12-03", "2005-01-04", 
    "2005-02-05", "2005-09-01", "2006-02-01", "2006-10-31", 
    "2009-07-19", "2009-10-30", "2010-12-04", "2011-01-30"
]  # e.g., ["2016-..."]

# Apply filtering using the exclusion lists
hlss = hlss.sel(time=~hlss.time.isin(hlss_exclude))
hlsl = hlsl.sel(time=~hlsl.time.isin(hlsl_exclude))
landsat7 = landsat7.sel(time=~landsat7.time.isin(landsat7_exclude))
landsat5 = landsat5.sel(time=~landsat5.time.isin(landsat5_exclude))

# ---------------------------
# 6) Extract time values and compute common dates between datasets
# ---------------------------
hlsl_dates = set(hlsl.time.values)
hlss_dates = set(hlss.time.values)
landsat7_dates = set(landsat7.time.values)
landsat5_dates = set(landsat5.time.values)

# Already existing common dates:
common_hlsl_hlss = sorted(str(date) for date in hlsl_dates.intersection(hlss_dates))
common_hlsl_landsat7 = sorted(str(date) for date in hlsl_dates.intersection(landsat7_dates))
common_hlss_landsat7 = sorted(str(date) for date in hlss_dates.intersection(landsat7_dates))

# New intersections including Landsat-5:
common_hlsl_landsat5 = sorted(str(date) for date in hlsl_dates.intersection(landsat5_dates))
common_hlss_landsat5 = sorted(str(date) for date in hlss_dates.intersection(landsat5_dates))
common_landsat7_landsat5 = sorted(str(date) for date in landsat7_dates.intersection(landsat5_dates))

# Print the results
print(f"Common dates between HLS-Landsat and HLS-Sentinel: {common_hlsl_hlss}")
print(f"Common dates between HLS-Landsat and Landsat-7: {common_hlsl_landsat7}")
print(f"Common dates between HLS-Sentinel and Landsat-7: {common_hlss_landsat7}")
print(f"Common dates between HLS-Landsat and Landsat-5: {common_hlsl_landsat5}")
print(f"Common dates between HLS-Sentinel and Landsat-5: {common_hlss_landsat5}")
print(f"Common dates between Landsat-7 and Landsat-5: {common_landsat7_landsat5}")

# ---------------------------
# 7) Function to create an RGB image for a given date
# ---------------------------
def get_rgb_image(ds, date, clip=(2, 98)):
    # Select data for the exact date
    ds_sel = ds.sel(time=date)  # no "nearest" here
    # Extract red, green, blue as numpy arrays
    r, g, b = [ds_sel[band].values for band in ["Red", "Green", "Blue"]]
    
    # Replace NaNs with 0
    r = np.nan_to_num(r, nan=0)
    g = np.nan_to_num(g, nan=0)
    b = np.nan_to_num(b, nan=0)
    
    # Stack bands into an RGB array
    rgb = np.dstack([r, g, b])
    
    # Percentile-based contrast stretching
    low, high = np.percentile(rgb, clip)
    print("Percentile low/high:", low, high)  # Debug print
    if high == low:  # Avoid division by zero
        high = low + 1e-10
    rgb_scaled = (rgb - low) / (high - low)
    rgb_scaled = np.clip(rgb_scaled, 0, 1)  # Ensure values are in [0, 1]
    
    return rgb_scaled

# ---------------------------
# 8) Function to plot two RGB images side by side and save
# ---------------------------
def plot_two_rgb_images(date, ds1, ds2, name1, name2, out_dir="figures"):
    os.makedirs(out_dir, exist_ok=True)
    rgb1 = get_rgb_image(ds1, date)
    rgb2 = get_rgb_image(ds2, date)
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    axs[0].imshow(rgb1)
    axs[0].axis("off")
    axs[0].set_title(f"{name1}\n{date}")
    axs[1].imshow(rgb2)
    axs[1].axis("off")
    axs[1].set_title(f"{name2}\n{date}")
    plt.tight_layout()
    fname = f"{name1}_vs_{name2}_{date}.png".replace(" ", "_")
    plt.savefig(os.path.join(out_dir, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)

# ---------------------------
# 9) Create side-by-side RGB plots for each pair of common dates
# ---------------------------
# Existing comparisons:
for date in common_hlsl_hlss:
    plot_two_rgb_images(date, hlsl, hlss, "HLS-Landsat", "HLS-Sentinel")

for date in common_hlsl_landsat7:
    plot_two_rgb_images(date, hlsl, landsat7, "HLS-Landsat", "Landsat-7")

for date in common_hlss_landsat7:
    plot_two_rgb_images(date, hlss, landsat7, "HLS-Sentinel", "Landsat-7")

# New comparisons including Landsat-5:
for date in common_hlsl_landsat5:
    plot_two_rgb_images(date, hlsl, landsat5, "HLS-Landsat", "Landsat-5")

for date in common_hlss_landsat5:
    plot_two_rgb_images(date, hlss, landsat5, "HLS-Sentinel", "Landsat-5")

for date in common_landsat7_landsat5:
    plot_two_rgb_images(date, landsat7, landsat5, "Landsat-7", "Landsat-5")

print("All side-by-side plots saved in the 'figures/' folder.")

# ---------------------------
# 10) Function to create scatter plots for a pair of datasets
# ---------------------------
from scipy.odr import ODR, Model, RealData  # Added for ODR
import matplotlib.cm as cm

coolwarm = cm.get_cmap("coolwarm")  # Get the colormap

def create_scatter_plots(
    ds_x, ds_y, common_dates, ds_x_name, ds_y_name, 
    band_names=None, num_pixels=20000
):
    print(f"Comparing {ds_x_name} vs. {ds_y_name}:")
    print(f"  - Number of common dates: {len(common_dates)}")

    if band_names is None:
        # Use only the bands common to both datasets
        band_names = list(set(ds_x.data_vars.keys()) & set(ds_y.data_vars.keys()))
        band_names = sorted(band_names)
    actual_bands = [b for b in band_names if b in ds_x.data_vars and b in ds_y.data_vars]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    axes = axes.flatten()

    for idx, band in enumerate(actual_bands[:6]):
        ax = axes[idx]
        x_vals_total, y_vals_total = [], []

        for date in common_dates:
            arr_x = ds_x[band].sel(time=date).values.ravel()
            arr_y = ds_y[band].sel(time=date).values.ravel()
            valid_mask = ~np.isnan(arr_x) & ~np.isnan(arr_y)
            arr_x, arr_y = arr_x[valid_mask], arr_y[valid_mask]
            
            if len(arr_x) >= num_pixels:
                rand_idx = np.random.choice(len(arr_x), size=num_pixels, replace=False)
                x_vals_total.append(arr_x[rand_idx])
                y_vals_total.append(arr_y[rand_idx])

        if not x_vals_total:
            ax.set_title(f"{band}\nNo valid samples")
            ax.grid(True, alpha=0.3)
            continue

        x_vals_total = np.concatenate(x_vals_total)
        y_vals_total = np.concatenate(y_vals_total)

        # Filter data to only include points between the 1st and 99th percentiles
        x_low, x_high = np.percentile(x_vals_total, [0, 100])
        y_low, y_high = np.percentile(y_vals_total, [0, 100])
        mask = (x_vals_total >= x_low) & (x_vals_total <= x_high) & \
               (y_vals_total >= y_low) & (y_vals_total <= y_high)
        x_vals_total = x_vals_total[mask]
        y_vals_total = y_vals_total[mask]

        hist, x_edges, y_edges = np.histogram2d(
            x_vals_total, y_vals_total, bins=50, 
            range=[[x_low, x_high], [y_low, y_high]]
        )
        x_bins = np.clip(np.digitize(x_vals_total, x_edges) - 1, 0, hist.shape[0]-1)
        y_bins = np.clip(np.digitize(y_vals_total, y_edges) - 1, 0, hist.shape[1]-1)
        density = hist[x_bins, y_bins]
        colors = coolwarm(density / density.max())

        def linear_model(beta, x):
            return beta[0] * x + beta[1]
        
        odr_model = Model(linear_model)
        data = RealData(x_vals_total, y_vals_total, 
                        sx=np.std(x_vals_total), 
                        sy=np.std(y_vals_total))
        odr = ODR(data, odr_model, beta0=[1.0, 0.0])
        odr_result = odr.run()
        slope, intercept = odr_result.beta
        sa, sb = odr_result.sd_beta

        # Calculate R² (using the ODR model)
        y_pred = slope * x_vals_total + intercept
        ss_res = np.sum((y_vals_total - y_pred) ** 2)
        ss_tot = np.sum((y_vals_total - np.mean(y_vals_total)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        # Calculate RMSE
        rmse = np.sqrt(np.mean((y_vals_total - y_pred) ** 2))

        # Plotting
        ax.scatter(x_vals_total, y_vals_total, c=colors, s=5, alpha=0.7, 
                   edgecolors='none', cmap='coolwarm')
        lim_min = min(x_low, y_low)
        lim_max = max(x_high, y_high)
        ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', alpha=0.6, lw=1)
        
        x_fit = np.linspace(lim_min, lim_max, 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, 'r-', lw=1.5, label='ODR Fit')

        # Annotation box with fit statistics
        textstr = '\n'.join((
            f'Slope: {slope:.3f} ± {sa:.3f}',
            f'Intercept: {intercept:.3f} ± {sb:.3f}',
            f'R²: {r_squared:.3f}',
            f'RMSE: {rmse:.3f}'
        ))
        
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))

        ax.set_title(f"{band}", fontsize=10)
        ax.set_xlabel(f"{band} ({ds_x_name})", fontsize=9)
        ax.set_ylabel(f"{band} ({ds_y_name})", fontsize=9)
        ax.set_aspect('equal')
        ax.grid(alpha=0.3)

    for j in range(len(actual_bands), 6):
        axes[j].axis('off')

    plt.tight_layout()
    plt.show()


# ---------------------------
# 11) Scatter plot comparisons between dataset pairs
# ---------------------------
# Existing comparisons:
create_scatter_plots(
    ds_x=hlsl,
    ds_y=hlss,
    common_dates=common_hlsl_hlss,
    ds_x_name="HLS-Landsat",
    ds_y_name="HLS-Sentinel"
)

create_scatter_plots(
    ds_x=hlsl,
    ds_y=landsat7,
    common_dates=common_hlsl_landsat7,
    ds_x_name="HLS-Landsat",
    ds_y_name="Landsat-7"
)

create_scatter_plots(
    ds_x=hlss,
    ds_y=landsat7,
    common_dates=common_hlss_landsat7,
    ds_x_name="HLS-Sentinel",
    ds_y_name="Landsat-7"
)

# New comparisons including Landsat-5:
create_scatter_plots(
    ds_x=hlsl,
    ds_y=landsat5,
    common_dates=common_hlsl_landsat5,
    ds_x_name="HLS-Landsat",
    ds_y_name="Landsat-5"
)

create_scatter_plots(
    ds_x=hlss,
    ds_y=landsat5,
    common_dates=common_hlss_landsat5,
    ds_x_name="HLS-Sentinel",
    ds_y_name="Landsat-5"
)

create_scatter_plots(
    ds_x=landsat7,
    ds_y=landsat5,
    common_dates=common_landsat7_landsat5,
    ds_x_name="Landsat-7",
    ds_y_name="Landsat-5"
)
