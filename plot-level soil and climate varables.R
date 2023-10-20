# extract plot-level soil variables, climate variables and SPEI

neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)#
# the 
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")

plot_loca_all=sample_data(neon_dob)[,c("lon","lat","geneticSampleID","Project","Site","siteID")]
plotIDD=substr(plot_loca_all$geneticSampleID,1,8)
plot_loca_all=cbind(plotIDD,plot_loca_all)

# get the soil variables of all the plots but now the points are at the finer resolution as one site may composed miltiple points

#https://github.com/claraqin/fungal-climate-niche/blob/master/R/1-data.R# creat the map of the study region
# for soil bulk density
r_bdod <- raster("bdod_5-15cm_mean_5000.tif")
names(r_bdod) <- "bdod"

r_bdod_reproj <- projectRaster(r_bdod, crs = crs(r_present_northam))

r_bdod_northam <- raster::mask(raster::crop(r_bdod_reproj, north_america_cropped), north_america_cropped)
r_bdod_northam_resample <- raster::resample(r_bdod_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_bdod_northam_resample / 10)#need to know why 10 but all variables would be standardized
# add soil nitrogen

r_nitrogen <- raster("nitrogen_5-15cm_mean_5000.tif")
names(r_nitrogen) <- "nitrogen"

r_nitrogen_reproj <- projectRaster(r_nitrogen, crs = crs(r_present_northam))

r_nitrogen_northam <- raster::mask(raster::crop(r_nitrogen_reproj, north_america_cropped), north_america_cropped)
r_nitrogen_northam_resample <- raster::resample(r_nitrogen_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_nitrogen_northam_resample / 10)#need to know why 10 but all variables would be standardized
##add cec
r_cec <- raster("cec_5-15cm_mean_5000.tif")
names(r_cec) <- "cec"

r_cec_reproj <- projectRaster(r_cec, crs = crs(r_present_northam))

r_cec_northam <- raster::mask(raster::crop(r_cec_reproj, north_america_cropped), north_america_cropped)
r_cec_northam_resample <- raster::resample(r_cec_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_cec_northam_resample / 10)#need to know why 10 but all variables would be standardized
#add soil sand
r_sand <- raster("sand_5-15cm_mean_5000.tif")
names(r_sand) <- "sand"

r_sand_reproj <- projectRaster(r_sand, crs = crs(r_present_northam))

r_sand_northam <- raster::mask(raster::crop(r_sand_reproj, north_america_cropped), north_america_cropped)
r_sand_northam_resample <- raster::resample(r_sand_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_sand_northam_resample / 10)#need to know why 10 but all variables would be standardized

# extract the six soil variables based on the coordinates of the plots
cl=c("organicCPercent","ph","bdod","nitrogen","cec","sand")

climate=list()
for (i in 1:6){
  climate[[i]]=raster::extract(r_present_northam[[cl[i]]], plot_loca_all[,2:3]) 
}
plot_loca_all_soil=cbind(plot_loca_all,climate[[1]],climate[[2]],climate[[3]],climate[[4]],climate[[5]],climate[[6]])
names(plot_loca_all_soil)=c("plotIDD" , "lon" ,  "lat"  , "geneticSampleID", "Project"  ,  "Site" ,   "siteID" ,"organicCPercent","ph","bdod","nitrogen","cec","sand" )
write.csv(plot_loca_all_soil,"plot_loca_all_soil.csv")

# get the climate data of all the plots
r <- getData("worldclim",var="bio",res=10)
points <- SpatialPoints(plot_loca_all[,2:3], proj4string = r@crs)
values <- extract(r,points)
df <- cbind.data.frame(coordinates(points),values)
# select the seven climate variables that do not show colinearity
##

# select the seven climate variables
neon_dob_climate=df[,c("bio2","bio8","bio18","bio4","bio12","bio15","bio1")] 
# all plot level climate and soil variables
plot_loca_all_soil_climate=cbind(plot_loca_all_soil,neon_dob_climate)# For some reasons,there are 233 NAS (for one site)

write.csv(plot_loca_all_soil_climate,"plot_loca_all_soil_climate.csv")

#add the id for all the plots,924 dob cores and 6218 neon cores
iddob=plot_loca_all_soil_climate[,"Site"][1:924]
idneon=plot_loca_all_soil_climate[,"plotIDD"][925:7142]
id=data.frame(c(iddob,idneon))

plot_loca_all_soil_climate=cbind(id,plot_loca_all_soil_climate)
names(plot_loca_all_soil_climate)[1]="plotIDM"# the plot id used to get the mean values of the variales

plot_loca_all_soil_climate_mean=aggregate(plot_loca_all_soil_climate[,c(9:22)],by=list(plot_loca_all_soil_climate$plotIDM),mean,sort=FALSE)
names(plot_loca_all_soil_climate_mean)[1]="plotID"
write.csv(plot_loca_all_soil_climate_mean,"plot_loca_all_soil_climate_mean.csv")
## get the land aridity index

devtools::install_github('seschaub/getSpei') 
require(getSpei)
require(ncdf4)
require(chron)

plot_spei=plot_loca_all[,1:3]
names(plot_spei)=c("site","longitude","latitude")

site_spei <- spec_spei(spei_files = c("spei01"), start_y = 2010, end_y = 2018,locations=plot_spei)

## need to get the mean values of the variables



