
r_present <- raster::getData("worldclim", var = "bio", res = 10)
r_present  <- worldclim_global(var = "bio", res = 10,path = getwd())
r_present <- r_present[[c(1)]]

# Run necessary transformations on wordclim-provided temperature data
r_present$bio1 <- r_present$bio1 / 10
r_present$bio4 <- r_present$bio4 / 1000 # no need for further transformation

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland

north_america <- ne_countries(continent = "North America",type="map_units")

# north_america<- subset(north_america, admin != "Canada" & admin != "Greenland")

st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170, -55, 18, 72), "SpatialPolygons")

north_america_cropped <- st_crop(north_america, b)

plot(north_america_cropped)

# Crop to region
r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)

###
greatlakes <- rnaturalearth::ne_download(
  scale = 110, type = "lakes", category = "physical"
) %>%
  sf::st_as_sf(lakes110, crs = 4269) %>%
  dplyr::filter(name_en %in% c("Superior", "Michigan", "Huron", "Erie", "Ontario"))

clipOutPoly <- function(r, poly) {
  r_poly <- raster::mask(r, poly)
  r[which(!is.na(as.matrix(r_poly)))] <- NA
  r
}


r_present_northam <- clipOutPoly(r_present_northam, greatlakes) # the map based on climates so no need to add variables

## get the variables for each grid cell
# add the soil carbon data to the map
# Data downloaded from ISRIC https://files.isric.org/soilgrids/latest/data_aggregated/5000m/phh2o/
# and https://files.isric.org/soilgrids/latest/data_aggregated/5000m/soc/

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

# for soil pH (need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph) <- "ph"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized

# for soil organic carbon

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc) <- "organicCPercent"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)
r_soc_northam_resample <- raster::resample(r_soc_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_soc_northam_resample / 100) # need to know why 10 but all variables would be standardized

# add the quadratic term of the temperature and the map

mat_2 <- r_present_northam$bio1^2
map_2 <- r_present_northam$bio12^2
r_present_northam <- addLayer(r_present_northam, mat_2) 
r_present_northam <- addLayer(r_present_northam, map_2) 
names(r_present_northam) <- c("mat_celsius", "temp_seasonality", "map_mm", "soilInCaClpH", "organicCPercent", "mat_celsius_2", "map_mm_2")

# prepare the environmental variables for the rcp2-4.5 scenarios

future_data <- geodata::cmip6_world(model = "ACCESS-ESM1-5", ssp = "245", time = "2061-2080", var = "bioc", download = F, res = 10, path = "/Users/luowenqi/soil-sar/SDM")
future_data <- as(future_data, "Raster")
future_data <- stack(future_data)
# get the variables of our interest
r_future <- future_data[[c(1, 4, 12)]]

names(r_future) <- c("mat_celsius", "temp_seasonality", "map_mm")

# Run necessary transformations on world climate-provided temperature data

r_future$temp_seasonality <- r_future$temp_seasonality / 100

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <- ne_countries(continent = "North America")
library(sf)
st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170, -55, 18, 72), "SpatialPolygons")
north_america_cropped <- st_crop(north_america, b)

plot(north_america_cropped)

###

# Crop to region
r_future_northam <- raster::mask(raster::crop(r_future, north_america_cropped), north_america_cropped)

# Crop to remove Great Lakes
greatlakes <- rnaturalearth::ne_download(
  scale = 110, type = "lakes", category = "physical"
) %>%
  sf::st_as_sf(lakes110, crs = 4269) %>%
  dplyr::filter(name_en %in% c("Superior", "Michigan", "Huron", "Erie", "Ontario"))
clipOutPoly <- function(r, poly) {
  r_poly <- raster::mask(r, poly)
  r[which(!is.na(as.matrix(r_poly)))] <- NA
  r
}

r_future_northam <- clipOutPoly(r_future_northam, greatlakes)

names(r_future_northam) <- c("mat_celsius", "temp_seasonality", "map_mm")

# add the quadratic terms

mat_2 <- r_future_northam$mat_celsius^2
map_2 <- r_future_northam$map_mm^2
r_future_northam <- raster::addLayer(r_future_northam, mat_2) # need to know why 10 but all variables would be standardized
r_future_northam <- addLayer(r_future_northam, map_2) # need to know why 10 but all variables would be standardized


# add the soil variables

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

# for soil pH (need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph) <- "ph"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_future_northam)
r_future_northam <- addLayer(r_future_northam, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized

# for soil total organic carbon content

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc) <- "organicCPercent"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)

r_soc_northam_resample <- raster::resample(r_soc_northam, r_future_northam)
r_future_northam <- addLayer(r_future_northam, r_soc_northam_resample / 100) # need to know why 10 but all variables would be standardized

names(r_future_northam) <- c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")

setwd("/Users/luowenqi/soil-sar/SDM")
save(r_future_northam, file = "r_future_northam.RData")


# prepare the environmental variables for the rcp5-8.5 scenarios

future_data <- geodata::cmip6_world(model = "ACCESS-ESM1-5", ssp = "585", time = "2061-2080", var = "bioc", download = F, res = 10, path = "/Users/luowenqi/soil-sar/SDM")
future_data <- as(future_data, "Raster")
future_data <- stack(future_data)
# get the variables of our interest
r_future <- future_data[[c(1, 4, 12)]]

names(r_future) <- c("mat_celsius", "temp_seasonality", "map_mm")

# Run necessary transformations on wordclim-provided temperature data

r_future$temp_seasonality <- r_future$temp_seasonality / 100

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <- ne_countries(continent = "North America")
library(sf)
st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170, -55, 18, 72), "SpatialPolygons")
north_america_cropped <- st_crop(north_america, b)

plot(north_america_cropped)

###

# Crop to region
r_future_northam_rcp585 <- raster::mask(raster::crop(r_future, north_america_cropped), north_america_cropped)

# Crop to remove Great Lakes
greatlakes <- rnaturalearth::ne_download(
  scale = 110, type = "lakes", category = "physical"
) %>%
  sf::st_as_sf(lakes110, crs = 4269) %>%
  dplyr::filter(name_en %in% c("Superior", "Michigan", "Huron", "Erie", "Ontario"))
clipOutPoly <- function(r, poly) {
  r_poly <- raster::mask(r, poly)
  r[which(!is.na(as.matrix(r_poly)))] <- NA
  r
}

r_future_northam_rcp585 <- clipOutPoly(r_future_northam_rcp585, greatlakes)

names(r_future_northam_rcp585) <- c("mat_celsius", "temp_seasonality", "map_mm")

# add the quadratic terms

mat_2 <- r_future_northam_rcp585$mat_celsius^2
map_2 <- r_future_northam_rcp585$map_mm^2
r_future_northam_rcp585 <- raster::addLayer(r_future_northam_rcp585, mat_2) # need to know why 10 but all variables would be standardized
r_future_northam_rcp585 <- addLayer(r_future_northam_rcp585, map_2) # need to know why 10 but all variables would be standardized


# add the soil variables

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

# for soil pH (need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph) <- "ph"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_future_northam_rcp585)
r_future_northam_rcp585 <- addLayer(r_future_northam_rcp585, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized

# for soil total organic carbon content

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc) <- "organicCPercent"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)

r_soc_northam_resample <- raster::resample(r_soc_northam, r_future_northam_rcp585)
r_future_northam_rcp585 <- addLayer(r_future_northam_rcp585, r_soc_northam_resample / 100) # need to know why 10 but all variables would be standardized

names(r_future_northam_rcp585) <- c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")

setwd("/Users/luowenqi/soil-sar/SDM")
save(r_future_northam_rcp585, file = "r_future_northam_rcp585.RData")