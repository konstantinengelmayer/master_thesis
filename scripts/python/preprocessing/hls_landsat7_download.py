# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 09:06:55 2025

@author: konst
"""

import ee
import geopandas as gpd
import json

# =============================================================================
# 0. Authenticate and initialize Earth Engine
# =============================================================================
# Uncomment the next line if you haven’t authenticated on this machine.
#ee.Authenticate()
ee.Initialize()

# =============================================================================
# 1. Load your AOI from a shapefile
# =============================================================================
shp_path = r"C:\Users\konst\Documents\Hiwi\mw3\drought_indicies\data\Untersuchungsgebiete\002a_Sen1-Subset_Eifel.shp"
gdf = gpd.read_file(shp_path)
geojson = gdf.to_json()
aoi = ee.Geometry(json.loads(geojson)['features'][0]['geometry'])

# =============================================================================
# 2. Set your date range
# =============================================================================
start_date = '2000-01-01'
end_date   = '2025-01-31'

# =============================================================================
# 3. Define helper functions
# =============================================================================

# 3a. Cloud & adjacent cloud masking for HLS (Sentinel and Landsat based)
def mask_clouds_adjacent(image):
    fmask = image.select('Fmask')
    # Extract the bits
    cloud_bit = fmask.rightShift(1).bitwiseAnd(1)      # Bit 1: Cloud
    adjacent_bit = fmask.rightShift(2).bitwiseAnd(1)     # Bit 2: Adjacent to Cloud/Shadow
    shadow_bit = fmask.rightShift(3).bitwiseAnd(1)       # Bit 3: Cloud Shadow

    #And(adjacent_bit.eq(0))
    # Create a mask: keep pixels clear of cloud, adjacent, and shadow flags.
    mask = cloud_bit.eq(0).And(shadow_bit.eq(0))
    
    return image.updateMask(mask)


# 3b. Cloud & adjacent cloud masking for Landsat 7
def mask_clouds_l7_adjacent(image):
    """
    For Landsat 7, use the QA_PIXEL band to detect and mask out:
      - Fill
      - Dilated clouds
      - Clouds
      - Cloud shadows
      - Pixels with high confidence for clouds
      - Pixels with high confidence for cirrus
    """
    qa = image.select('QA_PIXEL')
    
    # 1) Mask out fill, dilated cloud, cloud, cloud shadow:
    #    Fill (bit 0), Dilated Cloud (bit 1), Cloud (bit 3), Cloud Shadow (bit 4)
    #    We construct a bitmask with these bits set to 1, then keep only pixels
    #    where these bits are all zero (== 0).
    cloud_bits = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4)
    basic_mask = qa.bitwiseAnd(cloud_bits).eq(0)

    # 2) Exclude high-confidence clouds (bits 8–9 == 3):
    #    We shift right by 8 bits, then look at the 2-bit cloud confidence.
    #      0 = None, 1 = Low, 2 = Medium, 3 = High
    cloud_confidence = qa.rightShift(8).bitwiseAnd(3)
    cloud_conf_mask = cloud_confidence.lt(2)

    # 3) Exclude high-confidence cirrus (bits 14–15 == 3):
    #    Similarly, shift right by 14, then bitwiseAnd(3) to get cirrus confidence.
    #      0 = None, 1 = Low, 2 = Medium, 3 = High
    cirrus_confidence = qa.rightShift(14).bitwiseAnd(3)
    cirrus_conf_mask = cirrus_confidence.lt(1)
    
    # Combine all masks with logical AND
    final_mask = basic_mask.And(cloud_conf_mask).And(cirrus_conf_mask)
    
    # Update the image mask
    return image.updateMask(final_mask)


# 3c. Mask negative reflectance values by updating the mask.
def mask_negative_values(image, band_list):
    """
    For the given list of bands, compute the per-pixel minimum and mask out
    any pixels where that minimum is below 0.
    """
    valid_mask = image.select(band_list).reduce(ee.Reducer.min()).gte(0)
    return image.updateMask(valid_mask)

# 3d. Group images by day and mosaic them together.
def daily_mosaic(collection):
    """
    Creates daily mosaics from an image collection by grouping images by the
    day (using the system:time_start property) and mosaicking them.
    """
    # Get distinct dates in the collection as strings ("YYYY-MM-dd")
    dates = ee.List(
        collection.aggregate_array('system:time_start')
    ).map(lambda t: ee.Date(t).format('YYYY-MM-dd')).distinct()
    
    def mosaic_by_day(date_str):
        date = ee.Date(date_str)
        imgs = collection.filterDate(date, date.advance(1, 'day'))
        mosaic = imgs.mosaic().set('system:time_start', date.millis())
        return mosaic
    
    return ee.ImageCollection(dates.map(mosaic_by_day))

# 3e. Compute valid-data coverage using the Blue band’s mask as a proxy.
def add_coverage(image):
    """
    Adds a property 'coverage' to the image that is the mean of the Blue band’s
    mask over the AOI. (Since a mask value is 1 for valid data and 0 for masked.)
    """
    valid_mask = image.select('Blue').mask()
    stats = valid_mask.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=30,
        maxPixels=1e9
    )
    coverage = ee.Number(stats.get('Blue'))
    return image.set('coverage', coverage)

# 3f. Function to export a collection to a given Drive folder.
def export_collection(collection, folder, band_list):
    col_list = collection.toList(collection.size())
    count = col_list.size().getInfo()
    for i in range(count):
        image = ee.Image(col_list.get(i)).select(band_list).toFloat()  # select and cast to Float
        # Extract date from the image’s time property
        date_str = ee.Date(image.get('system:time_start')).format('YYYYMMdd').getInfo()
        task = ee.batch.Export.image.toDrive(
            image=image,
            description=f"{folder}_{date_str}",
            folder=folder,
            fileNamePrefix=f"{folder}_{date_str}",
            region=aoi,
            scale=30,
            maxPixels=1e13
        )
        task.start()
        print(f"Started export task for {folder}_{date_str}")


# =============================================================================
# 4. Load and process the HLS collections (Sentinel‐ and Landsat–based)
# =============================================================================

# HLS Sentinel–based
hlss30 = (ee.ImageCollection('NASA/HLS/HLSS30/v002')
          .filterBounds(aoi)
          .filterDate(start_date, end_date))
print("HLS Sentinel (raw) count:", hlss30.size().getInfo())

# HLS Landsat–based
hlsl30 = (ee.ImageCollection('NASA/HLS/HLSL30/v002')
          .filterBounds(aoi)
          .filterDate(start_date, end_date))
print("HLS Landsat (raw) count:", hlsl30.size().getInfo())

# Rename bands for consistency
hlss30_renamed = hlss30.map(lambda img: img.select(
    ['B2', 'B3', 'B4', 'B8A', 'B11', 'B12', 'B10', 'Fmask'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Cirrus', 'Fmask']
))
hlsl30_renamed = hlsl30.map(lambda img: img.select(
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B9', 'Fmask'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Cirrus', 'Fmask']
))

# =============================================================================
# 5. Load and process Landsat 7 and 5 data
# =============================================================================
landsat7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
            .filterBounds(aoi)
            .filterDate('2000-01-01', end_date))
print("Landsat 7 (raw) count:", landsat7.size().getInfo())

def apply_scale_factors(image):
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0)
    return image.addBands(opticalBands, None, True).addBands(thermalBand, None, True)

landsat7_scaled = landsat7.map(apply_scale_factors)
landsat7_renamed = landsat7_scaled.map(lambda img: img.select(
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL']
))

landsat5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
            .filterBounds(aoi)
            .filterDate('2000-01-01', end_date))
print("Landsat 5 (raw) count:", landsat5.size().getInfo())

landsat5_scaled = landsat5.map(apply_scale_factors)
landsat5_renamed = landsat5_scaled.map(lambda img: img.select(
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL']
)) 

# =============================================================================
# 6. Process each collection: apply adjacent cloud masking, mask negatives,
#    create daily mosaics, and filter mosaics with <20% valid data.
# =============================================================================

# For HLS Sentinel–based collection:
hlss30_processed = (hlss30_renamed
                    .map(mask_clouds_adjacent)
                    .map(lambda img: mask_negative_values(
                        img, ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Cirrus']
                    )))
hlss30_daily = (daily_mosaic(hlss30_processed)
                .map(add_coverage)
                .filter(ee.Filter.gte('coverage', 0.2)))

# For HLS Landsat–based collection:
hlsl30_processed = (hlsl30_renamed
                    .map(mask_clouds_adjacent)
                    .map(lambda img: mask_negative_values(
                        img, ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Cirrus']
                    )))
hlsl30_daily = (daily_mosaic(hlsl30_processed)
                .map(add_coverage)
                .filter(ee.Filter.gte('coverage', 0.2)))

# For Landsat 7 and 5 collection:
landsat7_processed = (landsat7_renamed
                      .map(mask_clouds_l7_adjacent)
                      .map(lambda img: mask_negative_values(
                          img, ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
                      )))
landsat7_daily = (daily_mosaic(landsat7_processed)
                  .map(add_coverage)
                  .filter(ee.Filter.gte('coverage', 0.2)))

landsat5_processed = (landsat5_renamed
                      .map(mask_clouds_l7_adjacent)
                      .map(lambda img: mask_negative_values(
                          img, ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
                      )))
landsat5_daily = (daily_mosaic(landsat5_processed)
                  .map(add_coverage)
                  .filter(ee.Filter.gte('coverage', 0.2)))
# -----------------------------------------------------------------------------
# Print out how many images there are for each daily mosaic collection:
print("HLS Sentinel daily mosaic count after filtering:", hlss30_daily.size().getInfo())
print("HLS Landsat daily mosaic count after filtering:", hlsl30_daily.size().getInfo())
print("Landsat 7 daily mosaic count after filtering:", landsat7_daily.size().getInfo())
print("Landsat 5 daily mosaic count after filtering:", landsat5_daily.size().getInfo())


# =============================================================================
# 7. Export each processed daily mosaic collection to separate folders
# =============================================================================
#export_collection(hlss30_daily, 'HLS_Sentinel_Kellerwald', ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
#export_collection(hlsl30_daily, 'HLS_Landsat_Kellerwald', ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
#export_collection(landsat7_daily, 'Landsat_7_strictest_Kellerwald', ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
export_collection(landsat5_daily, 'Landsat_5_strictest_Eifel', ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
