library(terra)

dem <- rast("data/satellite_data/predictors/dem.tif")

ai <- rast("data/satellite_data/predictors/ai.tif") %>%
  terra::project(y = dem) %>%
  terra::resample(y = dem) %>%
  terra::crop(y = dem)

predictors <- terrain(dem, v = c("slope", "aspect"))
predictors <- c(predictors, dem, ai)


writeRaster(predictors, "data/satellite_data/predictors/predictors.tif")

