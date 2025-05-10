#
###
#to get the sf object of the several states
# to crop part of the canada map based on the ranges

library(tigris)
library(rnaturalearth)
us_states <- states(cb = TRUE)
canadian_provinces <- ne_states(country = "Canada", returnclass = "sf")
cuba_provinces <- ne_states(country = "Cuba", returnclass = "sf")
mexico_states <- ne_states(country = "Mexico", returnclass = "sf")
rico_provinces <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Puerto Rico")
haiti_provinces <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Puerto Rico")

haiti_sf <- ne_countries(scale = "medium", country = "Haiti", returnclass = "sf")
dominican_republic <- ne_countries(scale = "medium", returnclass = "sf", country = "Dominican Republic")

baha_sf <- ne_countries(scale = "medium", country = "The Bahamas", returnclass = "sf")
jama_sf <- ne_countries(scale = "medium", country = "Jamaica", returnclass = "sf")


# to make projection for the sf objects

target_crs <- "EPSG:5070"
us_projected <- st_transform(us_states, crs = target_crs)
canadian_projected <- st_transform(canadian_provinces, crs = target_crs)
mexico_projected <- st_transform(mexico_states, crs = target_crs)
rico_projected <- st_transform(rico_provinces, crs = target_crs)
cuba_projected <- st_transform(cuba_provinces, crs = target_crs)

haiti_projected <- st_transform(haiti_sf, crs = target_crs)

dominican_projected <- st_transform(dominican_republic, crs = target_crs)

baha_projected <- st_transform(baha_sf, crs = target_crs)
jama_projected <- st_transform(jama_sf, crs = target_crs)

# to clip the canada map
bbox <- st_bbox(c(xmin = -2665004, ymin = 2157670 , xmax = 3188701, ymax = 5400000 ), crs = st_crs(canadian_projected))
bbox_sf <- st_as_sfc(bbox)
cropped_province <- st_crop(canadian_projected, bbox_sf)
canada_clipped=cropped_province



# the function to project the data based on a data.frame

my_function_project=function(data)
{
  if ("x"%in%colnames(data))
  {
    points <- vect(data, geom = c("x", "y"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    raster <- rasterize(points, raster_template, field = "value")
  }
  else{
    points <- vect(data, geom = c("lon", "lat"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    
    raster <- rasterize(points, raster_template, field = "group")
  }
  target_crs <- "EPSG:5070"
  raster_equal_area <- project(raster, target_crs,method="near")# there are options for the method used
  raster_df <- as.data.frame(raster_equal_area, xy = TRUE,)
  return(raster_df )
}

