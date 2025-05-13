#
library(tigris)
library(rnaturalearth)
library(sf)
library(plyr)
library(dplyr)
library(gstat)
library(raster)
library(ggplot2)
library(car)
library(classInt)
library(caret)
library(caretEnsemble)
library(doParallel)
library(gridExtra)
library(terra)
library(tidyr)
library(geodata)
library(reshape)
library(SSDM)
library(stars)
library(patchwork)
library(phyloseq)
library(dplyr)
library(iNEXT)

#to get the sf object of the several states
# to crop part of the Canada map based on the ranges


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

##add sf layers on the map

add_sf_layers <- function() {
  list(
    geom_sf(data = dominican_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = us_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = canada_clipped, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = rico_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = cuba_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = mexico_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = haiti_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = baha_projected, fill = NA, size = 0.01, color = "gray80"),
    geom_sf(data = jama_projected, fill = NA, size = 0.01, color = "gray80")
  )
}



#set different them for the maps
theme_map=theme(legend.spacing.y = unit(32, "pt"), 
                legend.position = c(0.15,0.35),
                legend.margin = margin(t = -30, r = 0, b = -1, l = 0),
                legend.text = element_text(size=8,angle=0),
                legend.box = "vertical",
                legend.justification = "center",
                legend.title = element_text(margin = margin(b = 4),size=10),
                text = element_text(size = 18),
                plot.title = element_text(size = 15, hjust = 0.5), 
                axis.text.y = element_blank(), 
                axis.text.x = element_blank(), 
                axis.title.y = element_text(size = 18), 
                axis.title.x = element_text(size = 18), 
                axis.ticks.x = element_blank(), 
                axis.ticks.y = element_blank(),
                plot.margin = unit(c(0.3, -5, -0.5, 0.5), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_blank())

theme_latitude=theme(legend.position = c(0.75,0.28),
                     legend.text = element_text(size=8),
                     legend.title  = element_text(size=10),
                     text = element_text(size = 18),
                     plot.title = element_text(size = 15, hjust = 0.5), 
                     axis.text.y = element_text(size=12), 
                     axis.text.x = element_text(size = 12), 
                     axis.title.y = element_text(size = 15), 
                     axis.title.x = element_text(size = 15), 
                     plot.margin = unit(c(0.3, 0.1, -.5, 0), "cm"),
                     panel.background = element_rect(fill = "NA"),
                     panel.border = element_rect(color = "black", size = 0.6, fill = NA))

##save the sf object

sf_layers <- list(
  dominican_projected, us_projected, canada_clipped, rico_projected,
  cuba_projected, mexico_projected, haiti_projected, baha_projected, jama_projected
)

# Create base plot




