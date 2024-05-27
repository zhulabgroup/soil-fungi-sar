

neon_dob  <- readRDS("~/soil-sar/phylo_V3.1.RDS")

neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))

neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)

# the extent of the sample data
#-156.64885  -66.82463 for the lontitude
#17.96382 71.28560 for the latitude

xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15

sample_data(neon_dob) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon_dob)

cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  dplyr::filter(Project == "NEON") %>%
  group_by(Site) %>%
  dplyr::summarise(n_coords = n_distinct(coordsSite)) %>%
  dplyr::arrange(desc(n_coords)) # Number of operational sites that each NEON site gets split into

sites_to_project <- cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(coordsSite) %>%
  dplyr::summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  dplyr::arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))

neon_dob_agg <- merge_samples(neon_dob, group="coordsSite")# all the elements would be grouped

sample_data(neon_dob_agg) %>%
  as("data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_agg), sites_to_project$coordsSite)],
         lon = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[1])),
         lat = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[2])),
         coordsSite = rownames(.)) ->
  sample_data(neon_dob_agg)

neon_dob_agg <- transform_sample_counts(neon_dob_agg, function(x) ifelse(x>0, 1, 0))

prevalence <- apply(otu_table(neon_dob_agg), 2, function(x) sum(x>0))

hist(prevalence)
sum(prevalence >= 20); mean(prevalence >= 20)
sum(prevalence >= 10); mean(prevalence >= 10)
sum(prevalence >= 5); mean(prevalence >= 5)

# Subset to OTUs present in at least 10 sites
neon_dob_prevalent <- prune_taxa(prevalence >= 10, neon_dob_agg)
# neon_dob_prevalent20 <- prune_taxa(prevalence >= 20, neon_dob_agg)




# get soil data

soil <- loadByProduct(dpID="DP1.10047.001")


soil$spc_biogeochem %>%
  # Get soil horizons that are completely less than 30 cm in depth
  dplyr::filter(horizonID %in% soil$spc_perhorizon$horizonID[which(soil$spc_perhorizon$horizonBottomDepth < 30)]) %>%
  left_join(dplyr::select(soil$spc_perplot, plotID, decimalLongitude, decimalLatitude), by="plotID") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(decimalLongitude, xgrid)], 2),
                            round(ygrid[findInterval(decimalLatitude, ygrid)], 2), sep="_")) %>%
  group_by(coordsSite) %>%dplyr::summarise(dplyr::across(c(phCacl2, phH2o, nitrogenTot, carbonTot, ctonRatio), .fns = ~ mean(.x, na.rm=FALSE))) %>%
  mutate(dplyr::across(c(nitrogenTot, carbonTot), .fns = ~ .x/10)) %>%
  # NOTE: Although we call it "soilInCaClpH", we actually use pH base in a water solution.
  # Similarly, although we call it "organicCPercent," we actually use the percent TOTAL carbon.
  dplyr::select(coordsSite, soilInCaClpH = phH2o, nitrogenPercent = nitrogenTot, organicCPercent = carbonTot) ->
  soilchem_coordsSites
# soil variables for the neon sites


sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  dplyr::filter(Project=="DoB") %>%
  group_by(coordsSite) %>%
  dplyr::summarise(dplyr::across(c(soilInCaClpH, nitrogenPercent, organicCPercent), .fns = ~ mean(.x, na.rm=TRUE))) -> soilchem_coordsSites_dob

soilchem_coordsSites <- rbind(soilchem_coordsSites,
                              soilchem_coordsSites_dob[which(!soilchem_coordsSites_dob$coordsSite %in% soilchem_coordsSites$coordsSite),])
# select the rows that do was not included in the neon data 




sample_data(neon_dob_prevalent) %>% as("data.frame") %>%
  mutate(soilInCaClpH = soilchem_coordsSites$soilInCaClpH[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         nitrogenPercent = soilchem_coordsSites$nitrogenPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         organicCPercent = soilchem_coordsSites$organicCPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)]) ->
  temp

apply(dplyr::select(temp, soilInCaClpH, nitrogenPercent, organicCPercent), 2, function(x) mean(!is.na(x)))

# soilInCaClpH nitrogenPercent organicCPercent# does not match with the initial data
#    0.9375000       0.9609375       0.9609375

sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp

## still not sure sure why the generated mean value does not match


r_present <- raster::getData("worldclim",var="bio",res=10)

r_present <- r_present[[c(1,4,12,15)]]
names(r_present) <- c("mat_celsius","temp_seasonality","map_mm","prec_seasonality")

# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$mat_celsius/10
r_present$temp_seasonality <- r_present$temp_seasonality/1000

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <-  ne_countries(continent="North America")
library(sf)
st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170,-55,18,72), "SpatialPolygons")
north_america_cropped  <- st_crop(north_america, b)
plot(north_america_cropped)

###

# Crop to region
r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)

# Crop to remove Great Lakes
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
r_present_northam <- clipOutPoly(r_present_northam, greatlakes)



# Load climate data into Phyloseq sample data
# Calculate decomposition coefficient (k) in sample data
addClimateToPhyseq <- function(physeq, r_climate) {
  data <- sample_data(physeq) %>% as("data.frame") %>%
    mutate(mat_celsius = raster::extract(r_climate[["mat_celsius"]], cbind(lon, lat)),
           map_mm = raster::extract(r_climate[["map_mm"]], cbind(lon, lat)),
           temp_seasonality = raster::extract(r_climate[["temp_seasonality"]], cbind(lon, lat)),
           prec_seasonality = raster::extract(r_climate[["prec_seasonality"]], cbind(lon, lat)))
  # If climate is NA for any sites, it's most likely because it falls just outside the raster cells.
  # Fill in with nearest cell.
  na_climate_ind <- which(is.na(data$mat_celsius))
  for(i in na_climate_ind) {
    xy <- cbind(data$lon[i], data$lat[i])
    nearest_ind <- which.min(replace(distanceFromPoints(r_climate[[1]], xy), is.na(r_climate[[1]]), NA))
    values <- r_climate@data@values[nearest_ind,]
    data$mat_celsius[i] <- values["mat_celsius"]
    data$map_mm[i] <- values["map_mm"]
    data$temp_seasonality[i] <- values["temp_seasonality"]
    data$prec_seasonality[i] <- values["prec_seasonality"]
  }
  # which(is.na(data$mat_celsius)) # All filled in!
  # data$k <- yasso_k(data$mat_celsius, data$map_mm)
  sample_data(physeq) <- data
  return(physeq)
}

neon_dob_prevalent <- addClimateToPhyseq(neon_dob_prevalent, r_present_northam)

# Get quadratic terms
sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2) ->
  temp
sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp
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

# Soil organic carbon content, a proxy for total percent carbon

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


r_nitrogen <- raster("nitrogen_5-15cm_mean_5000.tif")
names(r_nitrogen) <- "nitrogenPercent"

r_nitrogen_reproj <- projectRaster(r_nitrogen, crs = crs(r_present_northam))

r_nitrogen_northam <- raster::mask(raster::crop(r_nitrogen_reproj, north_america_cropped), north_america_cropped)

r_nitrogen_northam_resample <- raster::resample(r_nitrogen_northam, r_present_northam)
plot(r_nitrogen_northam_resample)
r_present_northam <- addLayer(r_present_northam, r_nitrogen_northam_resample / 100)

rm(r_soc)
rm(r_soc_reproj)
rm(r_soc_northam)
rm(r_soc_northam_resample)





###

sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))) # n=8
sum(is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "nitrogenPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH")) |
      is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=10


sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  mutate(soilInCaClpH = if_else(is.na(soilInCaClpH), raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)), soilInCaClpH),
         nitrogenPercent= if_else(is.na(nitrogenPercent), raster::extract(r_present_northam[["nitrogenPercent"]], cbind(lon, lat)), nitrogenPercent),
         organicCPercent = if_else(is.na(organicCPercent), raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent))->
  temp

sample_data(neon_dob_prevalent)=temp

sample_data(neon_dob_prevalent)[which(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))),]

which(is.na(temp$organicCPercent))# the 54th row has na



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

na_ind <- which(is.na(temp$nitrogenPercent))
for(i in na_ind) {
  xy <- cbind(temp$lon[i], temp$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["nitrogenPercent"]], xy), is.na(r_present_northam[["nitrogenPercent"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  temp$nitrogenPercent[i] <- values["nitrogenPercent"]
}



sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp


saveRDS(neon_dob_prevalent, "neon_dob_prevalent_v4.1.Rds")
saveRDS(neon_dob_agg, "neon_dob_agg_v2.1.Rds")
