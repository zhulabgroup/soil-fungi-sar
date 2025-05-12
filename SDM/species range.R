1.# determine species ranges at the continental scale with 8597 species as used in 

library(neonUtilities)
library(dplyr)
library(tidyr)
library(sp)
library(ggplot2)
library(rgeos)
library(sf)
library(rgdal)
library(raster)
library(rasterVis)
library(doParallel)
library(rnaturalearth)
library(phyloseq)

#functions

findThresholds <- function(physeq, models, terms = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality"),
                           s=c("min", "1se"), opt=c("eq","maxtss"), ncores=1) {
  message("findThresholds: Started at ", Sys.time())
  s <- match.arg(s, c("min", "1se"))
  opt <- match.arg(opt, c("eq","maxtss"))
  s <- switch(s, "min" = "lambda.min", "1se" = "lambda.1se")
  
  covars <- covarNamesFromLognets(models)
  pred_data <- as(sample_data(physeq)[,terms], "data.frame")
  pred_data_std <- scale(pred_data, center=TRUE, scale=TRUE)
  ind_cols <- match(covars, colnames(pred_data_std))
  pred_data_std <- pred_data_std[,ind_cols]
  
  registerDoParallel(ncores)
  thres <- foreach(i = seq_along(taxa_names(physeq)), .combine = c) %dopar% {
    if(identical(models[[i]], NA)) return(NA)
    sp <- taxa_names(physeq)[i]
    y01 <- as.vector(otu_table(physeq)[,sp]) > 0
    lognet_pred <- predict(models[[i]], newx = pred_data_std, type="response", s = s)
    t_grid <- seq(min(lognet_pred, na.rm=TRUE), max(lognet_pred, na.rm=TRUE),
                  by=(max(lognet_pred, na.rm=TRUE)-min(lognet_pred, na.rm=TRUE))/100)
    eval <- dismo::evaluate(p = lognet_pred[which(y01)], a = lognet_pred[which(!y01)], tr = t_grid)
    # plot(eval@TPR - eval@TNR)
    if(identical(opt, "eq")) {
      # Threshold to minimize diff b/w sensitivity and specificity
      t <- eval@t[which.min(abs(eval@TPR - eval@TNR))]
      return(t)
    }
    if(identical(opt, "maxtss")) {
      # threshold to maximize true skill statistic
      t <- eval@t[which.max(eval@TPR + eval@TNR)]
      return(t)
    }
  }
  stopImplicitCluster()
  message("findThresholds: Ended at ", Sys.time())
  return(thres)
}
###

covarNamesFromLognets <- function(models) {
  for(i in 1:length(models)) {
    if(identical(models[[i]], NA)) {
      if(i == length(models)) {
        warning("All models are NA! Returning all NAs.")
        return(rep(NA, length(models)))
      }
    } else {
      return(rownames(coef(models[[i]]))[-1])
    }
  }
}
## load the data


neon_dob <- readRDS("/Users/luowenqi/soil-sar/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)

#
xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 

# add a coordinates for each sample
sample_data(neon_dob) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon_dob)

###

cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  dplyr::filter(Project == "NEON") %>%
  group_by(Site) %>%
  summarise(n_coords = n_distinct(coordsSite)) %>%
  arrange(desc(n_coords)) # Number of operational sites that each NEON site gets split into
sites_to_project <- cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(coordsSite) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_agg <- merge_samples(neon_dob, group="coordsSite")
sample_data(neon_dob_agg) %>%
  as("data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_agg), sites_to_project$coordsSite)],
         lon = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[1])),
         lat = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[2])),
         coordsSite = rownames(.)) ->
  sample_data(neon_dob_agg)
# we have 128 sites now and the presence data is based on the sites within the 10-min degree



neon_dob_agg <- transform_sample_counts(neon_dob_agg, function(x) ifelse(x>0, 1, 0))

prevalence <- apply(otu_table(neon_dob_agg), 2, function(x) sum(x>0))
hist(prevalence)
sum(prevalence >= 20); mean(prevalence >= 20)
sum(prevalence >= 10); mean(prevalence >= 10)# the 8597 species included in the subsequent analysis
sum(prevalence >= 5); mean(prevalence >= 5)

neon_dob_prevalent <- prune_taxa(prevalence >= 10, neon_dob_agg)

# need to add this to the predicted occurrence

sample_data(neon_dob_prevalent)%>%data.frame()->temp

temp[,c("lon","lat")]->coordi_sample# the coordinates of the sampling points
# get the otu table based on the operational sites

otu_table(neon_dob_prevalent)%>%data.frame()->temp

# the absence data for the sample sites

sample_presence=cbind(coordi_sample,temp)

colnames(allpresence_absence)=colnames(sample_presence)
allpresence_absence=rbind(sample_presence,allpresence_absence)

colnames(allpresence_absence)=colnames(sample_presence)
allpresence_absence=rbind(sample_presence,allpresence_absence)

##

data=






## get soil data for neon sites

soil <- loadByProduct(dpID="DP1.10047.001")
soil$spc_biogeochem %>%
  # Get soil horizons that are completely less than 30 cm in depth
  dplyr::filter(horizonID %in% soil$spc_perhorizon$horizonID[which(soil$spc_perhorizon$horizonBottomDepth < 30)]) %>%
  left_join(dplyr::select(soil$spc_perplot, plotID, decimalLongitude, decimalLatitude), by="plotID") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(decimalLongitude, xgrid)], 2),
                            round(ygrid[findInterval(decimalLatitude, ygrid)], 2), sep="_")) %>%
  group_by(coordsSite) %>%
  summarise(across(c(phCacl2, phH2o, nitrogenTot, carbonTot, ctonRatio), .fns = ~ mean(.x, na.rm=FALSE))) %>%
  mutate(across(c(nitrogenTot, carbonTot), .fns = ~ .x/10)) %>%
  # NOTE: Although we call it "soilInCaClpH", we actually use pH base in a water solution.
  # Similarly, although we call it "organicCPercent," we actually use the percent TOTAL carbon.
  dplyr::select(coordsSite, soilInCaClpH = phH2o, nitrogenPercent = nitrogenTot, organicCPercent = carbonTot) ->
  soilchem_coordsSites

#
sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  dplyr::filter(Project=="DoB") %>%
  group_by(coordsSite) %>%
  summarise(across(c(soilInCaClpH, nitrogenPercent, organicCPercent), .fns = ~ mean(.x, na.rm=TRUE))) -> soilchem_coordsSites_dob

soilchem_coordsSites <- rbind(soilchem_coordsSites,
                              soilchem_coordsSites_dob[which(!soilchem_coordsSites_dob$coordsSite %in% soilchem_coordsSites$coordsSite),])


# Add soil data to neon_dob_prevalent

sample_data(neon_dob_prevalent) %>% as("data.frame") %>%
  mutate(soilInCaClpH = soilchem_coordsSites$soilInCaClpH[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         nitrogenPercent = soilchem_coordsSites$nitrogenPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         organicCPercent = soilchem_coordsSites$organicCPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)]) ->
  temp
# does not matche with Qin's result
apply(dplyr::select(temp, soilInCaClpH, nitrogenPercent, organicCPercent), 2, function(x) mean(!is.na(x)))

sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp

##traits




###


# Load present (actually 1970-2000 average) climate data




####

r_present <- getData("worldclim",var="bio",res=10)
r_present <- r_present[[c(1,4,12,15)]]
names(r_present) <- c("mat_celsius","temp_seasonality","map_mm","prec_seasonality")

# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$mat_celsius/10
r_present$temp_seasonality <- r_present$temp_seasonality/1000

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <- ne_countries(continent="North America")

st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
library(sf)
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

### to devide the map into grids and get the centroid of each

us_extent <- extent(-170,-55,18,72)

raster_layer <- raster(ext = us_extent, res = 1/6)

grid_polygons <- as(raster_layer, "SpatialPolygonsDataFrame")

grid_coordinates <- coordinates(grid_polygons)%>%data.frame()

## get the variables for each of the grid
# add the soil carbon data to the map
# Data downloaded from ISRIC https://files.isric.org/soilgrids/latest/data_aggregated/5000m/phh2o/
# and https://files.isric.org/soilgrids/latest/data_aggregated/5000m/soc/

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc) <- "organicCPercent"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)
r_soc_northam_resample <- raster::resample(r_soc_northam, r_present_northam)
plot(r_soc_northam_resample)
r_present_northam <- addLayer(r_present_northam, r_soc_northam_resample / 100)
# add the soil pH layer

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph) <- "soilInCaClpH"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10)



mat_celsius = raster::extract(r_present_northam[["mat_celsius"]], grid_coordinates)%>%data.frame()
map_mm = raster::extract(r_present_northam[["map_mm"]], grid_coordinates)%>%data.frame()
temp_seasonality = raster::extract(r_present_northam[["temp_seasonality"]], grid_coordinates)%>%data.frame()
prec_seasonality = raster::extract(r_present_northam[["prec_seasonality"]], grid_coordinates)%>%data.frame()
organicCPercent = raster::extract(r_present_northam[["organicCPercent"]], grid_coordinates)%>%data.frame()
soilInCaClpH = raster::extract(r_present_northam[["soilInCaClpH"]], grid_coordinates)%>%data.frame()


grid_coordinates_climate=cbind(grid_coordinates,mat_celsius,mat2=(mat_celsius)^2,map_mm,map2=(map_mm)^2, temp_seasonality,prec_seasonality,organicCPercent,soilInCaClpH)

names(grid_coordinates_climate)=c("lon","lat","mat_celsius","mat2", "map_mm", "map2","temp_seasonality","prec_seasonality","organicCPercent","soilInCaClpH")

# the variables of the sites are ready and now we can move forward to the prediction





### Get quadratic terms

sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2) ->
  temp
sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp

###

saveRDS(neon_dob_prevalent, "data/neon_dob_prevalent_v4.1.Rds")
 saveRDS(neon_dob_agg, "data/neon_dob_agg_v2.1.Rds")
 
 ##
 
 r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
 names(r_ph) <- "soilInCaClpH"
 r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
 r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
 r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
 r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10)
 
 rm(r_ph)
 rm(r_ph_reproj)
 rm(r_ph_northam)
 rm(r_ph_northam_resample)
 
 ## Soil organic carbon content, a proxy for total percent carbon
 
 r_soc <- raster("soc_5-15cm_mean_5000.tif")
 names(r_soc) <- "organicCPercent"
 r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
 r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)
 r_soc_northam_resample <- raster::resample(r_soc_northam, r_present_northam)
 plot(r_soc_northam_resample)
 r_present_northam <- addLayer(r_present_northam, r_soc_northam_resample / 100)
 
 rm(r_soc)
 rm(r_soc_reproj)
 rm(r_soc_northam)
 rm(r_soc_northam_resample)
 
 #
 
 sample_data(neon_dob_prevalent) %>%
   as("data.frame") %>%
   mutate(soilInCaClpH = if_else(is.na(soilInCaClpH), raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)), soilInCaClpH),
          organicCPercent = if_else(is.na(organicCPercent), raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent)) ->
   temp
 
 
 
 
 soil_comparison_df <- data.frame(
   variable = c(rep("pH", nrow(sample_data(neon_dob_prevalent))), rep("SOC", nrow(sample_data(neon_dob_prevalent)))),
   val_neon_dob = c(get_variable(neon_dob_prevalent, "soilInCaClpH"), get_variable(neon_dob_prevalent, "organicCPercent")),
   val_raster = c(
     raster::extract(r_present_northam[["soilInCaClpH"]], cbind(get_variable(neon_dob_prevalent, "lon"), get_variable(neon_dob_prevalent, "lat"))),
     raster::extract(r_present_northam[["organicCPercent"]], cbind(get_variable(neon_dob_prevalent, "lon"), get_variable(neon_dob_prevalent, "lat"))))
 )
 
 #
 
 ggplot(soil_comparison_df, aes(x = val_raster, y = val_neon_dob)) +
   facet_grid(variable~.) +
   geom_point() +
   geom_abline(slope=1, intercept=0)
 ##
 
 in.site.poly.inf <- function(sitedata, geodf, xvar, yvar, inflate = 0.25) {
   site.data = sitedata[, c(xvar, yvar)]
   x = site.data[,1]
   y = site.data[,2]
   Mx = mean(x)
   My = mean(y)
   CH=chull(site.data)
   BumpX = x[CH] + inflate*(x[CH]-Mx)
   BumpY = y[CH] + inflate*(y[CH]-My)
   geo.df = geodf[, c(xvar, yvar)]
   inpolygon(geo.df[,1], geo.df[,2], BumpX, BumpY, boundary = T)
   
   
   ###
   
   getMaskFromClimateHull <- function(sitedata, geo_r, inflate = 0.15) {
     geo_df <- as.data.frame(cbind(coordinates(geo_r), as.matrix(geo_r)))
     names(geo_df)[1:2] <- c("lon", "lat")
     geo_df <- geo_df[complete.cases(geo_df),]
     
     # matrix for logical test of variable combinations
     var_comb <- c("MAT.TMPSEA", "MAT.MAP", "MAT.PRSEA", "MAP.TMPSEA", "MAP.PRSEA", "TMPSEA.PRSEA")
     var_mat <- matrix(nrow=nrow(geo_df), ncol=6)
     colnames(var_mat) <- var_comb
     var_mat <- data.frame(var_mat)
     
     # logical tests
     var_mat$MAT.TMPSEA <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "temp_seasonality", inflate)
     var_mat$MAT.MAP <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "map_mm", inflate)
     var_mat$MAT.PRSEA <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "prec_seasonality", inflate)
     var_mat$MAP.TMPSEA <- in.site.poly.inf(sitedata, geo_df, "map_mm", "temp_seasonality", inflate)
     var_mat$MAP.PRSEA <- in.site.poly.inf(sitedata, geo_df, "map_mm", "prec_seasonality", inflate)
     var_mat$TMPSEA.PRSEA <- in.site.poly.inf(sitedata, geo_df, "temp_seasonality", "prec_seasonality", inflate)
     
     # rowSums to see total number of polygons in
     var_mat$in.or.out <- rowSums(var_mat)
     
     # rasterize the points
     r_within <- rasterFromXYZ(cbind(geo_df$lon, geo_df$lat, var_mat$in.or.out))
     
     r_within_01 <- r_within
     r_within_01[r_within==6] <- 1
     r_within_01[r_within!=6] <- 0
     
     # how much area inclusive
     within_area_df <- data.frame(rasterToPoints(raster::area(r_within_01, na.rm=T)))
     within_area <- sum(within_area_df$layer)
     
     # Define mask based on hull
     climatehull_mask <- r_within_01
     climatehull_mask[climatehull_mask==0] <- NA
     
     return(list(mask = climatehull_mask, area = within_area, map = r_within_01))
   }
 
###
   
getMaskFromEnvHull <- function(sitedata, geo_r,
                                  terms_nonquadratic = c("mat_celsius", "map_mm", "temp_seasonality", "soilInCaClpH", "organicCPercent"),
                                  inflate = 0.15) {
     geo_df <- as.data.frame(cbind(coordinates(geo_r), as.matrix(geo_r)))
     names(geo_df)[1:2] <- c("lon", "lat")
     geo_df <- geo_df[complete.cases(geo_df),]
     
     # matrix for logical test of variable combinations
     # var_comb <- c("MAT.TMPSEA", "MAT.MAP", "MAT.PRSEA", "MAP.TMPSEA", "MAP.PRSEA", "TMPSEA.PRSEA")
     var_mat <- matrix(nrow=nrow(geo_df), ncol=10)
     # colnames(var_mat) <- var_comb
     var_mat <- data.frame(var_mat)
     
     # logical tests
     inhull_list <- list()
     for(i in 1:(length(terms_nonquadratic)-1)) {
       for(j in (i+1):length(terms_nonquadratic)) {
         xvar <- terms_nonquadratic[i]
         yvar <- terms_nonquadratic[j]
         x <- sitedata[,xvar]
         y <- sitedata[,yvar]
         complete_ind <- which(complete.cases(cbind(x, y)))
         x <- x[complete_ind]
         y <- y[complete_ind]
         ch_ind <- chull(x, y)
         infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
         infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
         inhull_list <- c(inhull_list, list(inpolygon(geo_df[,xvar], geo_df[,yvar], infl_x, infl_y, boundary=TRUE)))
       }
     }
     
     # rowSums to see total number of polygons in
     inhull_bool <- rowSums(do.call(cbind, inhull_list)) == 10# there are 10 kinds of combinations
     
     # rasterize the points
     r_within_01 <- rasterFromXYZ(cbind(geo_df$lon, geo_df$lat, inhull_bool))
     
     # how much area inclusive
     within_area_df <- data.frame(rasterToPoints(raster::area(r_within_01, na.rm=T)))
     within_area <- sum(within_area_df$layer)
     
     # Define mask based on hull
     climatehull_mask <- r_within_01
     climatehull_mask[climatehull_mask==0] <- NA
     
     return(list(mask = climatehull_mask, area = within_area, map = r_within_01))
}

#
obsEnv <- as(sample_data(neon_dob_prevalent), "data.frame")


# filling the gaps soil variables

sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))) # n=8
sum(is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "nitrogenPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH")) |
      is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=10

sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%mutate(soilInCaClpH = if_else(is.na(soilInCaClpH), raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)), soilInCaClpH),
         organicCPercent = if_else(is.na(organicCPercent), raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent))->
  temp

sample_data(neon_dob_prevalent) <- temp

sample_data(neon_dob_prevalent)[which(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))),]

which(is.na(temp$organicCPercent))# which site that has NA for the data






na_ind <- which(is.na(temp$soilInCaClpH))
for(i in na_ind) {
  xy <- cbind(temp$lon[i], temp$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["soilInCaClpH"]], xy), is.na(r_present_northam[["soilInCaClpH"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  temp$soilInCaClpH[i] <- values["soilInCaClpH"]
}

na_ind <- which(is.na(temp$organicCPercent))
for(i in na_ind) {
  xy <- cbind(temp$lon[i], temp$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["organicCPercent"]], xy), is.na(r_present_northam[["organicCPercent"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  temp$organicCPercent[i] <- values["organicCPercent"]
}




sum(is.na(temp$soilInCaClpH))
sum(is.na(temp$organicCPercent))

sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp


# the data used for the model is not complete, we will not include soil nitrogen in the model

#3----construct the model

set.seed(10101)
train_set <- sample(c(TRUE, FALSE), size=length(sample_sums(neon_dob_prevalent)),
                    prob=c(0.70, 0.30), replace=TRUE)
test_set <- !train_set
neon_dob_prevalent_train <- prune_samples(train_set, neon_dob_prevalent)
neon_dob_prevalent_test <- prune_samples(test_set, neon_dob_prevalent)
n_taxa <- 8597
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent_train)), size=n_taxa))
neon_dob_prevalent_train_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_train)[spp_subset_ind],
                                              neon_dob_prevalent_train)
neon_dob_prevalent_test_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_test)[spp_subset_ind],
                                             neon_dob_prevalent_test)
##

# important_vars (a posteriori naming)

terms1 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")

library(phyloseq)

predictors_train1 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms1], "data.frame")
predictors_train_std1 <- scale(predictors_train1, center=TRUE, scale=TRUE)
complete_records1 <- which(complete.cases(predictors_train_std1))



library(glmnet)# do not forget to load the package
library(pbapply)
install.packages("foreach")
library(foreach)

registerDoParallel(5)
lognets_train1 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std1)[complete_records1,],
                     y01[complete_records1], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
  
} 
stopImplicitCluster()

###try to add an bar


library(doSNOW)

cl <- makeCluster(2)
registerDoSNOW(cl)

iterations <- 100
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

lognets_train1 <- foreach(i = 1:seq_along(spp_subset_ind),.combine = rbind, 
                          .options.snow = opts,.packages=c("phyloseq")) %dopar% 
  {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std1)[complete_records1,],
                     y01[complete_records1], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
  s <- summary(rnorm(1e6))[3]
  return(s)
} 
close(pb)
stopCluster(cl) 

stopImplicitCluster()




#
library(glmnet)

lognet_thresholds_train1 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train1,
  terms = terms1)


###the locations the model is going to predict

names(grid_coordinates_climate)[4]="mat_celsius_2"
names(grid_coordinates_climate)[6]="map_mm_2"

newdata1 <- as(grid_coordinates_climate, "data.frame")[,terms1]# variables in the test data

newdata_std1 <- scaleToReference(newdata1, predictors_train1)# variables in the train data


###3

registerDoParallel(5)
lognet_preds1 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train1[[i]], NA)) return(NA)
  prediction <- predict(lognets_train1[[i]],
                        newx = as.matrix(newdata_std1),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train1[[i]]
  as.numeric(classification)
}

# cbind all the results
allpresence_absence=matrix(ncol=n_taxa,nrow=dim(newdata1)[1])
for(i in 1:n_taxa){
  allpresence_absence[,i]=lognet_preds1[[i]]
}

## to plot the species range for one species
grid_coordinates_climate=cbind(code=rownames(grid_coordinates_climate),grid_coordinates_climate)
allpresence_absence=cbind(code=rownames(allpresence_absence),allpresence_absence)

d=cbind(grid_coordinates_climate[,1:4],allpresence_absence[,1:dim(allpresence_absence)[2]])
d=subset(d,mat_celsius!="NA")# leads to 98191 grids remained

# map the presence data on the map


#

length(unique(get_variable(neon_dob, "Site")))
sites_to_project <- cbind.data.frame(
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(Site) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_site <- merge_samples(neon_dob, group="Site")
sites_df <- as(sample_data(neon_dob_site), "data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_site), sites_to_project$Site)])


d=allpresence_absence[,1:10]
d=subset(d,1!="NA")

library(sf)

ggplot(sites_df) +
  geom_sf(data=st_transform(st_as_sf(north_america_cropped)),
          # geom_sf(data=st_as_sf(north_america_cropped),
          size=0.1, col="black", fill=alpha("white", 0)) +
  theme_custom_map+
  geom_point(data=subset(d,otu8>0),aes(x=lon,y=lat),color="blue",size=0.1)

names(d)=c("lon","lat","otu1","otu2","otu3","otu4","otu5","otu6","otu7","otu8")




stopImplicitCluster()

