library(sf)
library(rnaturalearth)
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


# (1). to determine the richness in each grid cell with SDM

# extract the climate and soil variables for each grid cell


r_present <- raster::getData("worldclim", var = "bio", res = 10)

r_present <- r_present[[c(1, 4, 12)]]

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

# add the quardratic term of the temperature and the map

mat_2 <- r_present_northam$bio1^2
map_2 <- r_present_northam$bio12^2
r_present_northam <- addLayer(r_present_northam, mat_2) # need to know why 10 but all variables would be standardized
r_present_northam <- addLayer(r_present_northam, map_2) # need to know why 10 but all variables would be standardized
names(r_present_northam) <- c("mat_celsius", "temp_seasonality", "map_mm", "soilInCaClpH", "organicCPercent", "mat_celsius_2", "map_mm_2")


# (2). construct SDMs based on current climate with 128 operational sites and then got the richness within each grid cell


agg_data <- readRDS("~/soil-sar/SDM/neon_dob_prevalent_v4.1.Rds")

env <- sample_data(agg_data) %>%
  data.frame() %>%
  dplyr::select(lon, lat, mat_celsius, temp_seasonality, soilInCaClpH, organicCPercent, mat_celsius_2, map_mm_2)

occ <- otu_table(agg_data) %>%
  data.frame() %>%
  bind_cols(env[, 1:2])

## need to convert the env data into roasters rather than data-frame
## transformed the occurrence data into desired format

occurrence_data <- melt(occ[, 1:8597])
df <- occ[, 8598:8599]
replicated_df <- do.call(rbind, replicate(8597, df, simplify = FALSE))
occurrence_data <- cbind(occurrence_data, replicated_df)
# check if all the species have occurrence data
d <- colSums(occ[, 1:8597])
sp <- colnames(occ[, 1:8597]) # model each species and save the output
cl <- makeCluster(4)
registerDoParallel(cl)

# each run involves modeling one species
sdm <- list()
for (i in 1:8597) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  sdm[[i]] <- modelling("RF", subset(occurrence_data, variable == sp[i]),
    r_present_northam,
    Xcol = "lon", Ycol = "lat", Pcol = "value", Spcol = "variable", verbose = FALSE
  )
  saveRDS(sdm[[i]], file = paste0("result_", i, ".rds"))
}

stopImplicitCluster()

# model all species with two different models (takes about 100 hours)

SSDM <- stack_modelling(c("GLM", "RF"), occurrence_data,
  r_present_northam,
  Xcol = "lon",
  Ycol = "lat",
  Pcol = "value",
  Spcol = "variable",
  rep = 1, ensemble.thresh = 0,
  verbose = FALSE
)

save(occurrence_data, file = "occurrence_data.RData")
save(r_present_northam, file = "r_present_northam.RData")


# (3). read in the data and quantify the total richness for each grid cell
# the coordinates of each cell was based on the land-use change data
# get the coordinates for each grid cell with the land-use change data in 2020

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
raster1 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2020.nc")
ext(raster1) <- c(-90, 90, -180, 180)
crs(raster1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent(-72, -18, -170, -55), "SpatialPolygons")

coarser_raster <- aggregate(raster1, fact = 3, fun = mean) # convert it to a coarse resolution

save(coarser_raster, file = "coarser_raster_2020.RData")


# (4). read in the species distribution data to get the total richness in each grid cell

setwd("/Users/luowenqi/soil-sar/SDM")

coords_present <- coords_present[, c(2, 1)]
coords_present[, 1] <- -1 * coords_present[, 1]

data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  df <- readRDS(paste0("result_", i, ".rds"))
  data[, i] <- raster::extract(df@binary, coords_present)
}
# did not use the two model

save(data, file = "presence_present.RData") # 0-1 data for 8597 species(can not find this data, perhaps, it was not saved)

# get the richness map for the current distribution

data <- data.frame(data)

d <- rowSums(data) # the total richness for each grid cell

richness_present <- cbind(coords_present, d) %>%
  data.frame() %>%
  rename(lon = x, lat = y, richness = d)
#this data should have been saved at great

# (5). model species distribution under climate change with SDM
# get the future climate variables
# https://bedatablog.netlify.app/post/download-and-illustrate-current-and-projected-climate-in-r/

# future_data <- raster::getData(name = "CMIP5", var = "bio", res = 10, rcp = 45, model = "IP", year = 70)

future_data <- geodata::cmip6_world(model = "ACCESS-ESM1-5", ssp = "245", time = "2061-2080", var = "bioc", download = F, res = 10, path = "/Users/luowenqi/soil-sar/SDM")
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



# (6). model each species' presence and get the total richness in each cell under future climate scenarios

sdm <- list()
for (i in 1:8597) {
  sdm[[i]] <- modelling("RF", subset(occurrence_data, variable == sp[i]),
    r_future_northam,
    Xcol = "lon", Ycol = "lat", Pcol = "value", Spcol = "variable", verbose = FALSE
  )
  saveRDS(sdm[[i]], file = paste0("result_future", i, ".rds"))
}

stopImplicitCluster()

# read in the data and get the total richness for each cell for future climate scenarios
# this is based on 
data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  df <- readRDS(paste0("result_future", i, ".rds"))
  data[, i] <- raster::extract(df@binary, coords_present)
}

save(data, file = "presence_future.RData") # 0-1 data for 8597 species

# get the richness map for the current distribution
data <- data.frame(data)
d <- rowSums(data) # the total richness for each grid cell
richness_future <- cbind(coords_present, d) %>%
  data.frame() %>%
  rename(lon = x, lat = y, richness = d)


species_change_climate <- cbind(richness_present, richness_future[, 3]) %>%
  rename_all(~ paste0(c("lon", "lat", "present_rich", "future_rich"))) %>%
  mutate(change = future_rich - present_rich, rate = change / present_rich)
save(species_change_climate, file = "species_change_climate.RData")

load("~/soil-sar/SDM/richness_present_RF.RData")# the stored data is richness_present

# just look at the cell with species loss, cells with increased richness were assigned to 0
species_change_climate1 <- cbind(richness_present, richness_future[, 3]) %>%
  rename_all(~ paste0(c("lon", "lat", "present_rich", "future_rich"))) %>%
  mutate(change = future_rich - present_rich, rate = change / present_rich) %>%
  mutate(rate = ifelse(rate > 0, 0, rate))

# just look at the cell with species increase, cells with decreased richness were assigned to 0
species_change_climate2 <- cbind(richness_present, richness_future[, 3]) %>%
  rename_all(~ paste0(c("lon", "lat", "present_rich", "future_rich"))) %>%
  mutate(change = future_rich - present_rich, rate = change / present_rich) %>%
  mutate(rate = ifelse(rate < 0, 0, rate))

# (7). look at land-use change-induced richness change

# for the year of 2020

PFT_2020 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster_2020 <- crop(flip(coarser_raster[[i]]), b)
  coords_present <- xyFromCell(cropped_raster_2020, cell = 1:ncell(cropped_raster_2020)) # get the coordinates
  cell_values <- raster::extract(cropped_raster_2020, coords_present) %>% as.matrix()
  PFT_2020[, i] <- cell_values
}

# total cover of non-crop land use types in each cell

PFT_2020 <- cbind(coords_present, PFT_2020) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster1))))

temp=PFT_2020[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2020%>%bind_cols(No_crop)->PFT_2020
###


# land use data in 2100 with global change scenario of SSP2 & RCP 4.5

raster2 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2100.nc")
ext(raster2) <- c(-90, 90, -180, 180)
crs(raster2) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent(-72, -18, -170, -55), "SpatialPolygons")

coarser_raster_2100 <- aggregate(raster2, fact = 3, fun = mean) # convert it to a coarse resolution

PFT_2100 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster_2100 <- crop(flip(coarser_raster_2100[[i]]), b)
  coords_present <- xyFromCell(cropped_raster_2100, cell = 1:ncell(cropped_raster_2100)) # get the coordinates
  cell_values <- raster::extract(cropped_raster_2100, coords_present) %>% as.matrix()
  PFT_2100[, i] <- cell_values
}
###


PFT_fine <- matrix(nrow = 2484000, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped<- crop(flip(raster1[[i]]), b)
  coords_present <- xyFromCell(cropped, cell = 1:ncell(cropped)) # get the coordinates
  cell_values <- extract(cropped, coords_present) %>% as.matrix()
  PFT_fine[, i] <- cell_values
}
## to see the total cover of non-crop data

PFT_fine%>%data.frame()%>%bind_cols(coords_present)%>%rename_all(~paste0(c(names(raster1),"lon","lat")))%>%mutate(nocrop = rowSums(across(PFT1:PFT14)))->PFT_fine

ggplot()+
  geom_point(data=PFT_fine,aes(x=lon,y=lat,color=nocrop/100))

PFT_2100 <- cbind(coords_present, PFT_2100) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster2)))) 

## get the total non-crop land cover value in each cell

temp=PFT_2100[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2100%>%bind_cols(No_crop)->PFT_2100
###

PFT_2100=PFT_2100%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14)

PFT_2020=PFT_2020%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14)

save(PFT_2100,file="PFT_2100_ssp245.RData")
save(PFT_2020,file="PFT_2020_ssp245.RData")

# (8). changes in different land-use types between the two time points (2020 vs 2100).

land_use_change <- matrix(nrow = dim(PFT_2020[1]), ncol = 36)
for (i in 1:dim(PFT_2020)[2])
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  land_use_change[, i] <- PFT_2100[, i] / PFT_2020[, i]
}

bind_cols(PFT_2020[, 1:2], land_use_change[, 3:36]) %>% dplyr::rename_all(~ paste0(c("lon", "lat", names(raster2), "nocrop"))) -> land_use_change

# if we have more explicit classification

land_use_change_explict <- matrix(nrow = dim(PFT_2020[1]), ncol = 39)
for (i in 1:dim(PFT_2020)[2])
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  land_use_change_explict[, i] <- PFT_2100[, i] / PFT_2020[, i]
}

bind_cols(PFT_2020[, 1:2], land_use_change_explict[, 3:39]) %>% dplyr::rename_all(~ paste0(c("lon", "lat", names(raster2), "nocrop","tree","shrub","grass"))) -> land_use_change_explict

# look at different changes

cbind(PFT_2020[,36:39],PFT_2100[,36:39],land_use_change_explict[,36:39])->temp

names(temp)=c("crop2020","tree2020","shrub2020","grass2020","crop2100","tree2100","shrup2100","grass2100","crop_ratio","tree_ratio","shrub_ratio","grass_ratio")


temp%>%filter(tree_ratio<1)%>%dim()/dim(temp)[1]
temp%>%filter(shrub_ratio<1)%>%dim()/dim(temp)[1]
temp%>%filter(grass_ratio<1)%>%dim()/dim(temp)[1]
# most change were caused by forest 

temp%>%filter(tree_ratio>1)%>%dim()/dim(temp)[1]#2%
temp%>%filter(shrub_ratio>1)%>%dim()/dim(temp)[1]
temp%>%filter(grass_ratio>1)%>%dim()/dim(temp)[1]#5%

# most change were caused by forest 

temp%>%filter(tree_ratio==1&shrub_ratio==1&grass_ratio<1)%>%dim()


## to look at changes in species richness

load("~/soil-sar/SDM/richness_present_RF.RData")

land_use_change <- land_use_change %>%dplyr::rename(lat = lon, lon = lat)

land_use_change$lat <- -1 * land_use_change$lat


# look at the changes in richness caused by land-use conversion
# we used the mean value of 0.3877299 for z across all cells
# k <- cbind(PFT_2020$nocrop, PFT_2100$nocrop, land_use_change$nocrop)

kk <- richness_present %>%
  mutate(future_rich = land_use_change$nocrop^0.3877299 * richness) %>%
  bind_cols(land_use_change[, c("nocrop")]) %>%
  rename_all(~ paste0(c("lon", "lat", "richness", "future_rich", "land_change_ratio"))) %>%
  mutate(species_change = future_rich - richness, rich_change_rate = species_change / richness)

land_use_change_rich <- kk

save(land_use_change_rich, file = "land_use_change_rich.RData")#save as ssp245

# for future projected richness, i made some conversions

# species richness was considered unchanged for the following scenarios
# 1. cells don't have land-use data (NA) but we have present richness data
# 2. cells don't have land-use data for both time points(0-0) but we have present richness data
# 3. cells that have land-use data increased from 0 to some value (Inf)
# 4. in these cases, future richness will be the same as it's current richness
# 5. some cells will have very high value of richness because of a high land cover change ratio
# these codes were run on the server and loaded the output
land_use_change_rich_nona <- subset(land_use_change_rich, richness > 0) # select cells with richness data

df <- land_use_change_rich_nona[, (3:7)]
dff <- matrix(ncol = 5, nrow = dim(land_use_change_rich_nona)[1])
for (i in 1:119607)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  if (is.nan(df[i, 3]) && !is.nan(df[i, 1])) # land-use type value in the cell is 0
    {
      dff[i, 1] <- df[i, 1]
      dff[i, 2] <- df[i, 1]
      dff[i, 3] <- NA
      dff[i, 4] <- 0
      dff[i, 5] <- 0 # we assume species does not change
    } else if (is.na(df[i, 3]) && !is.na(df[i, 1])) {
    dff[i, 1] <- df[i, 1]
    dff[i, 2] <- df[i, 1] # land-use type value in the cell is NA
    dff[i, 3] <- NA
    dff[i, 4] <- 0
    dff[i, 5] <- 0
  } else if (is.infinite(df[i, 3]) && !is.nan(df[i, 1])) {
    dff[i, 1] <- df[i, 1]
    dff[i, 2] <- df[i, 1] # land-use data changed from 0 to some value, in which case the equation indefinite
    dff[i, 3] <- NA
    dff[i, 4] <- NA
    dff[i, 5] <- NA # we do not evaluate such changes but richness should increase
  } else {
    dff[i, ] <- df[i, ] %>% as.numeric()
  }
}

dff %>%
  data.frame() %>%
  rename_all(~ paste0(c(colnames(df)))) %>%
  bind_cols(land_use_change_rich_nona[, 1:2]) -> land_use_change_rich_nona_conver 

# in some cases the data need to be transformed to make it more realistic

# if the estimated richness exceeds 8597, it is assigned to 8597, the total richness included in the analysis

land_use_change_rich_nona_conver%>%mutate(future_rich_adjust=ifelse(future_rich>8597,8597,future_rich),species_change_adjust=future_rich_adjust-richness,rich_change_rate_adjust=species_change_adjust/richness)->land_use_change_rich_nona_conver



# just look at cell with increased richness
land_use_change_rich_nona_conver %>% mutate(rich_change_rate_adjust = ifelse(rich_change_rate_adjust > 0, 0, rich_change_rate_adjust)) -> land_use_change_rich_nona_conver1

# just look at cells with decreased richness
land_use_change_rich_nona_conver %>% mutate(rich_change_rate_adjust = ifelse(rich_change_rate_adjust < 0, 0, rich_change_rate_adjust)) -> land_use_change_rich_nona_conver2


# (10) create maps for species distribution


# look at changes in land cover data among two time periods



ggplot(land_use_change_rich_nona_conver) +
  geom_point(data = land_use_change_rich_nona_conver, aes(x = lon, y = lat, color = land_change_ratio), size = 0.275)+

  scale_color_gradient2(expression("Ratio"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = 0.5, na.value = "white", limits = c(0, 1))+ 
theme(
  legend.position = c(0.2, 0.35),
  legend.key.size = unit(0.15, "inches"),
  guides(color = guide_legend(nrow = 2, byrow = TRUE)),
  legend.title = element_text(size = 8),
  text = element_text(size = 18),
  legend.text = element_text(size = 8),
  plot.title = element_text(size = 15, hjust = 0.5),
  axis.text.y = element_text(hjust = 0),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill = "NA"),
  panel.border = element_rect(color = "black", size = 1.5, fill = NA)
) +
  xlab("2100/2020") +
  ylab("") +
  ggtitle("Non-crop land cover %")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))




ggplot(PFT_2100) +
  geom_point(data = PFT_2100, aes(y = -1*lon, x = lat, color = No_crop), size = 0.275)+ 

  scale_color_gradient2(expression("Cover %"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = 50, na.value = "white", limits = c(0, 100))+ 
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Non-crop land cover in 2100") +
  ylab("") +
  ggtitle("")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))

#mapping richness in the region
custom_them=theme(
  legend.position = c(0.2, 0.35),
  legend.key.size = unit(0.15, "inches"),
  guides(color = guide_legend(nrow = 2, byrow = TRUE)),
  legend.title = element_text(size = 8),
  text = element_text(size = 18),
  legend.text = element_text(size = 8),
  plot.title = element_text(size = 15, hjust = 0.5),
  axis.text.y = element_text(hjust = 0),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill = "NA"),
  panel.border = element_rect(color = "black", size = 1.5, fill = NA)
) 

 p_land_future=ggplot(land_use_change_rich_nona_conver) +
  geom_point(data = land_use_change_rich_nona_conver, aes(x = lon, y = lat, color = future_rich), size = 0.275) +
  scale_color_gradient2(expression("Richness"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = 4000, na.value = "white", limits = c(0, 8000)) +
  custom_them
  xlab("Predicted richness in 2100") +
  ylab("") +
  ggtitle("RCP 4.5 & SSP 2")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))




p_land_loss=ggplot(land_use_change_rich_nona_conver1) +
  geom_point(data = land_use_change_rich_nona_conver1, pch = 15, 
  aes(x = lon, y = lat, color = rich_change_rate_adjust* 100), size = 0.275) +
  scale_color_gradient2(expression("%"), low = "red",mid="blue", high = "white", midpoint = -50,na.value = "white") +
  custom_them+
  xlab("Predicted species loss") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))




p_land_gain=ggplot(land_use_change_rich_nona_conver2) +
  geom_point(data = land_use_change_rich_nona_conver2, pch = 15, aes(x = lon, y = lat, color = rich_change_rate_adjust * 100), size = 0.275) +
  scale_color_gradient2(expression("%"), low = "white", mid="blue",high = "red",midpoint = 100, na.value = "blue",limits=c(0,200)) +
  custom_them+
  xlab("Predicted species gain") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))

wrap_plots(p1,p_land_future,p_land_loss,p_land_gain,ncol=2)+
  plot_annotation(tag_levels ="A",theme = theme(plot.tag = element_text(size = 5)))

###
#0.03% with ratio>1
#0.000259 with ratio>10
##0.000526725 with ratio>5
#0.7657411 with ratio <1
# suggest that most area will have few forest or grass land cover

ggplot(land_use_change_rich_nona_conver1) +
  geom_point(data = land_use_change_rich_nona_conver1, aes(x = lon, y = lat, color =  land_change_ratio), size = 0.275)+

  scale_color_gradient2(expression("Species change rate %"), low = "red", mid="blue",high = "seagreen1", midpoint = 300, na.value = "white") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species changes") +
  ylab("")


plot_grid(p1, p2, ncol = 2)

 p1=ggplot(species_change_climate) +
  geom_point(data = species_change_climate, aes(x = lon, y = lat, color = present_rich), size = 0.275) +
  scale_color_gradient2(expression("Richness"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = 4000, na.value = "white", limits = c(0, 8000)) +
   geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))+
   theme(
     plot.margin = unit(c(0, 0, 0, 0), "cm"),
     legend.position = c(0.2, 0.35),
     legend.key.size = unit(0.15, "inches"),
     guides(color = guide_legend(nrow = 2, byrow = TRUE)),
     legend.title = element_text(size = 8),
     text = element_text(size = 18),
     legend.text = element_text(size = 8),
     plot.title = element_text(size = 15, hjust = 0.5),
     axis.title.y = element_text(size = 18),
     axis.title.x = element_text(size = 18),
     axis.ticks.x= element_blank(),
     panel.background = element_rect(fill = "NA"),
     panel.border = element_rect(color = "black", size = 1.5, fill = NA)
   ) +
   xlab("Predicted richness in 2020") +
   ylab("") +
   ggtitle("Current")
 
 
p2=ggplot(species_change_climate) +
  geom_point(data = species_change_climate, aes(x = lon, y = lat, color = future_rich), size = 0.275) +
  scale_color_gradient2(expression("Richness"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = 4000, na.value = "white", limits = c(0, 8000)) +
  geom_sf(data=st_as_sf(north_america_cropped),
          size=0.1, col="black", fill=alpha("white", 0))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted richness in 2100") +
  ylab("") +
  ggtitle("RCP 4.5 & SSP 2")


p3=ggplot(species_change_climate1) +
  geom_point(data = species_change_climate1, pch = 15, aes(x = lon, y = lat, color = rate * 100), size = 0.275) +
  scale_color_gradient2(expression("%"), low = "blue", mid = "red", high = "white", midpoint = -50, na.value = "white", limits = c(-100, 0)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species loss") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  geom_sf(data=st_as_sf(north_america_cropped),
          size=0.1, col="black", fill=alpha("white", 0))


p4=ggplot(species_change_climate2) +
  geom_point(data = species_change_climate2, pch = 15, aes(x = lon, y = lat, color = rate * 100), size = 0.275) +
  scale_color_gradient2(expression("%"), low = "white", mid = "blue", high = "red", midpoint = 100, na.value = "white", limits = c(0, 200)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
  legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species gain") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  geom_sf(data=st_as_sf(north_america_cropped),
          size=0.1, col="black", fill=alpha("white", 0))

wrap_plots(p1,p2,p3,p4,ncol=2)+
  plot_annotation(
    tag_levels ="A",theme = theme(plot.title = element_text(size = 10))  # Labels plots with 'A', 'B', 'C', 'D'
  )
 

# the converted data was used to create the map


head(PFT_2020)

##changes in the land-use cover data

p1 <- ggplot(PFT_2020) +
  geom_point(data = PFT_2020, pch = 15, aes(y = -1 * lon, x = lat, color = nocrop), size = 0.275) +
  scale_color_gradient(expression("Non-crop cover"), low = "seagreen3", high = "purple", na.value = "white") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("2020") +
  ylab("")


p2 <- ggplot(PFT_2100) +
  geom_point(data = PFT_2100, pch = 15, aes(y = -1 * lon, x = lat, color = nocrop), size = 0.275) +
  scale_color_gradient(expression("Non-crop cover"), low = "seagreen3", high = "purple", na.value = "white") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("2100") +
  ylab("")

###the effects of land use change on richness




ggplot(dff) +
  geom_point(data = dff%>%filter(land_change_ratio <1), pch = 15, aes(x = lon, y = lat, color = land_change_ratio), size = 0.275) +
  scale_color_gradient(expression("Land-use change"), high = "gray", low = "purple",na.value = "red") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species changes") +
  ylab("")

### to see the spatial matching of richness loss caused by either by climate change and land use conversion
# cells with increased richness were assigned with 0 for better visualization
# threshold = -0.1; arbitrary, the loss of 10% of species was considered as high degree of species loss

threshold = -0.1

land_use_loss <- land_use_change_rich_nona_conver  %>% mutate(rich_change_rate_adjust = ifelse(rich_change_rate_adjust > 0, 0, rich_change_rate_adjust)) # range of 0-1

# assign different species loss degree to the data

land_use_loss %>% mutate(new_column = case_when(is.na(rich_change_rate_adjust) ~ NA_character_, rich_change_rate_adjust < threshold ~ "high", TRUE ~ "low")) -> land_loss_degree

# for the climate induced richness change

species_change_climate1 %>% mutate(new_column = case_when(is.na(rate) ~ NA_character_, rate < threshold ~ "high", TRUE ~ "low")) -> climate_loss_degree

## need to bind the two data.frames to show the impact of global change factors on richness

land_loss_degree %>% mutate(location = paste(lon, "*", lat)) -> land_loss_degree

names(land_loss_degree)[11] <- "land_degree"

climate_loss_degree %>% mutate(location = paste(lon, "*", lat)) -> climate_loss_degree

land_climate_loss_degree_rcp45_spp2 <- merge(climate_loss_degree, land_loss_degree[, c(1:5, 8:12)], by = "location") %>% mutate(com_degree = paste(new_column, "*", land_degree))

save(land_climate_loss_degree, file = "land_climate_loss_degree.RData")

save(land_climate_loss_degree_rcp45_spp2,file="land_climate_loss_degree_rcp45_spp2.RData")

# make some conversion so that some scenarios were grouped into one ca
#"NA * high" ="low * high" 
#"high * NA" ="high * low" 
#"NA * low"    
#"low * NA" 

land_climate_loss_degree%>%mutate(com_degree=gsub("NA-*-high","low-*-high",com_degree))->temp

land_climate_loss_degree$com_degree=gsub("NA-*-high","low-*-high",land_climate_loss_degree$com_degree)

land_climate_loss_degree$com_degree=as.character(land_climate_loss_degree$com_degree)


ggplot(land_climate_loss_degree_rcp45_spp2) +
  geom_point(data = land_climate_loss_degree_rcp45_spp2, pch = 15, aes(x = lon, y = lat, color = com_degree), size = 0.38)+

  scale_color_manual(expression("%"),
    breaks = unique(land_climate_loss_degree_rcp45_spp2$com_degree),
    labels = c("Land use high (23.9%)", "Both high (5.6%)",  "Both low (57.5%)", "Climate high (10.2%)","Uncertainty", "Uncertainty"),
    values = c("deepskyblue3", "black", "gray", "chocolate1", "gray", "deepskyblue3")
  ) +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 12),
    text = element_text(size = 18),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species loss") +
  ylab("") +
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0))
  



ggplot(dff) +
  geom_point(data = dff%>%filter(is.na(land_change_ratio)), pch = 15, aes(x = lon, y = lat, color = land_change_ratio), size = 0.275) +
  scale_color_gradient(expression("Land-use change"), high = "gray", low = "purple",na.value = "purple") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("Predicted species changes") +
  ylab("")

ggplot()+
  geom_point(data=PFT_fine,aes(y=-1*lon,x=lat,color=nocrop/100))+
  scale_color_gradient(expression("Land-use cover"), low = "gray", high = "purple",na.value = "white") +
  theme(
    legend.position = c(0.2, 0.35),
    legend.key.size = unit(0.15, "inches"),
    guides(color = guide_legend(nrow = 2, byrow = TRUE)),
    legend.title = element_text(size = 8),
    text = element_text(size = 18),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.y = element_text(hjust = 0),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  ) +
  xlab("") +
  ylab("")
  
## convert a data.frame into a raster

df=species_change_climate[,1:3]
coordinates(df) <- ~lon+lat

proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

extent_north_america <- extent(-170, -50, 5, 85)  # longitudes and latitudes

r1 <- raster(extent_north_america, ncol=1, nrow=275760)


r <- raster(df, r1,field="present_rich")

plot(r)

