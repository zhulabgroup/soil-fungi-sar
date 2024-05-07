
library(sf)
library(rnaturalearth)
newdata$lat=-1*newdata$lat
# add the variables of these sites

r_present <- raster::getData("worldclim",var="bio",res=10)

r_present <- r_present[[c(1,2,4,12,15,18)]]



# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$bio1/10
r_present$temp_seasonality <- r_present$bio4/1000# no need for further transformation

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland

north_america <- ne_countries(continent="North America")
st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)



b <- as(extent(-170,-55,18,72), "SpatialPolygons")
north_america_cropped  <- st_crop(north_america, b)
plot(north_america_cropped)


# Crop to region
r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)

###

greatlakes <- rnaturalearth::ne_download(
  scale = 110, type = 'lakes', category = 'physical'
) %>%
  sf::st_as_sf(lakes110, crs = 4269) %>%
  dplyr::filter(name_en %in% c("Superior", "Michigan", "Huron", "Erie", "Ontario"))
clipOutPoly <- function(r, poly) {
  r_poly <- raster::mask(r, poly)
  r[which(!is.na(as.matrix(r_poly)))] <- NA
  r
}
r_present_northam <- clipOutPoly(r_present_northam, greatlakes)# the map based on climates so no need to add variables

## get the variables for each of the grid
# add the soil carbon data to the map
# Data downloaded from ISRIC https://files.isric.org/soilgrids/latest/data_aggregated/5000m/phh2o/
# and https://files.isric.org/soilgrids/latest/data_aggregated/5000m/soc/



setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

# for the soil pH(need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph)="ph"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized

# for the soil soc

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc)="soc"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)
r_soc_northam_resample <- raster::resample(r_soc_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_soc_northam_resample / 10) # need to know why 10 but all variables would be standardized


r_nitrogen <- raster("nitrogen_5-15cm_mean_5000.tif")
names(r_nitrogen)="nitrogen"
r_nitrogen_reproj <- projectRaster(r_nitrogen, crs = crs(r_present_northam))
r_nitrogen_northam <- raster::mask(raster::crop(r_nitrogen_reproj, north_america_cropped), north_america_cropped)
r_nitrogen_northam_resample <- raster::resample(r_nitrogen_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_nitrogen_northam_resample / 10) # need to know why 10 but all variables would be standardized

###
#r_cec <- raster("cec_5-15cm_mean_5000.tif")

#r_cec_reproj <- projectRaster(r_cec, crs = crs(r_present_northam))
#r_cec_northam <- raster::mask(raster::crop(r_cec_reproj, north_america_cropped), north_america_cropped)
#r_cec_northam_resample <- raster::resample(r_cec_northam, r_present_northam)
#r_present_northam <- addLayer(r_present_northam, r_cec_northam_resample / 10) # need to know why 10 but all variables would be standardized


##
r_sand <- raster("sand_5-15cm_mean_5000.tif")
names(r_sand)="sand"
r_sand_reproj <- projectRaster(r_sand, crs = crs(r_present_northam))
r_sand_northam <- raster::mask(raster::crop(r_sand_reproj, north_america_cropped), north_america_cropped)
r_sand_northam_resample <- raster::resample(r_sand_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_sand_northam_resample / 10) # need to know why 10 but all variables would be standardized

names(r_present_northam) <- c("mat_celsius","MDR", "temp_seasonality","map_mm","pre_seas","pre_wq","ph", "soc",  "nitrogen", "sand")




# extract the six soil variables based on the coordinates of the plots
cl <- c("mat_celsius" ,  "MDR" ,  "temp_seasonality", "map_mm" ,    "pre_seas"  ,  "pre_wq","ph", "soc",  "nitrogen", "sand"  )

cl <- c("bio1",             "bio2"  ,           "bio4" ,            "bio12" ,           "bio15" ,          
        "bio18", "ph"  ,             "soc" , "nitrogen"    ,     "sand"  )


climate <- list()
for (i in 1:10) {
  climate[[i]] <- raster::extract(r_present_northam[[cl[i]]], newdata[, 1:2])
}

plot_loca_all_soil <- cbind(newdata[,1:2], climate[[1]], climate[[2]],climate[[3]], climate[[4]],climate[[5]], climate[[6]],climate[[7]], climate[[8]],climate[[9]], climate[[10]])

names(plot_loca_all_soil) <- c("lon", "lat" ,  cl)

newdata=plot_loca_all_soil
