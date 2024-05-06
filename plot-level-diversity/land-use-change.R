
# land use types change

# Define the extent of latitude and longitude coordinates
min_lat = -72  # Minimum latitude
max_lat = -18 # Maximum latitude
min_lon = -170 # Minimum longitude
max_lon =-55 # Maximum longitude


na_extent <- extent(-170, -55, 18, 72)

res <- 1/6  # 10 minutes in degrees
ncols <- diff(range(c(na_extent@xmin, na_extent@xmax))) / res
nrows <- diff(range(c(na_extent@ymin, na_extent@ymax))) / res
na_raster <- raster(extent(na_extent), ncols = ncols, nrows = nrows)


#b <- as(extent(-72,-18,-170,-55), "SpatialPolygons")

# Calculate the number of grid cells needed in latitude and longitude directions

# Create a grid of latitude and longitude points
lat_grid <- seq(from = min_lat, to = max_lat, by = 10 / 60)
lon_grid <- seq(from = min_lon, to = max_lon, by = 10 / 60)



setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
library(RColorBrewer)
library(raster)
library(terra)
library(dplyr)
library(ggplot2)

raster1 <- rast("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc")
ext(raster1)=c(-90,90,-180,180)

crs(raster1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

b <- as(extent(-72,-18,-170,-55), "SpatialPolygons")

na_extent <- extent(-168, -15, 15, 83)

res <- 1/6  # 10 minutes in degrees
ncols <- diff(range(c(na_extent@xmin, na_extent@xmax))) / res
nrows <- diff(range(c(na_extent@ymin, na_extent@ymax))) / res

na_raster <- raster(extent(na_extent), ncols = ncols, nrows = nrows)

set.seed(123)
sample_points <- data.frame(
  lon = runif(100, -168, -15),  # Sample points within the extent of North America
  lat = runif(100, 15, 83)
)

sample_points$grid_id <- cellFromXY(na_raster, cbind(sample_points$lon, sample_points$lat))

centroids <- rasterToPoints(na_raster, spatial = TRUE)


interval_size <- 1/6 # For example, 20 minutes in degrees

# Create intervals around central points
centroids$lon_interval <- cut(centroids[,1], breaks = seq(from = na_extent@xmin, to = na_extent@xmax, by = interval_size), labels = FALSE, include.lowest = TRUE)
centroids$lat_interval <- cut(centroids[,2], breaks = seq(from = na_extent@ymin, to = na_extent@ymax, by = interval_size), labels = FALSE, include.lowest = TRUE)


#b <- as(extent(-170,-55,18,72), "SpatialPolygons")

# get different layers

PFT=matrix(nrow=2484000,ncol=33)
for (i in 1:33){
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing

cropped_raster <- crop(flip(raster1[[i]]), b)

coords_present <- xyFromCell(cropped_raster, cell = 1:ncell(cropped_raster))# get the coordinates

cell_values <- extract(cropped_raster,coords_present)%>%as.matrix()
PFT[,i]=cell_values 
}
# each column indicates the proportion of each of the 32 types

# to check is the row sums equals 1

rowSums(PFT)%>%unique()# does not awalys equals 1

PFT=cbind(PFT,rowSums(PFT))%>%data.frame()
PFT=cbind(coords_present,PFT)

## aggregate the cells into a largr cell based on the resolution


# 
PFT$x=-1*PFT$x
names(PFT)[1]="lat"
names(PFT)[2]="lon"

data=PFT[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))%>%bind_cols(PFT[,3:35])

# add the grid to the land types data
data=cbind(data,PFT[,3:35])
# different land cover types and the central point of the grid

# get the mean value for the land cover types for each grid

data_mean=aggregate(data[,c(1,2,7:38)],by=list(data$grid),FUN=mean)%>%data.frame()# aggregate the value within a larger grid

data_mean=data_frame(data_mean)

# to see the change of land cover types for some grids




# differencce in land cover types between the year of 2020 and 2100

df=matrix(nrow=2484000,ncol=8)
for (i in 1:2484000){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  d1=PFT[i,][,3:10]
  d2=PFT_FUTURE[i,][,3:10]
  df[i,]=as.matrix(d2-d1)
}



PFT%>%filter(X33<100)

cell_values =cbind(coords_present,cell_values )%>%


plot(flip(raster1[[2]]))


##

# the present distribution of the land cover types


raster2 <- rast("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2020.nc")

ext(raster2)=c(-90,90,-180,180)

crs(raster2) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

b <- as(extent(-72,-18,-170,-55), "SpatialPolygons")

PFT_FUTURE=matrix(nrow=2484000,ncol=33)

for (i in 1:33){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
cropped_raster1 <- crop(flip(raster2[[i]]), b)

coords_futrue <- xyFromCell(cropped_raster1, cell = 1:ncell(cropped_raster1))# get the coordinates

cell_values_future <- extract(cropped_raster1,coords_futrue)%>%as.matrix()
PFT_FUTURE[,i]=cell_values_future 
}

PFT_FUTURE=cbind(PFT_FUTURE[,1:33],rowSums(PFT_FUTURE[,1:33]))%>%data.frame()

PFT_FUTURE=cbind(coords_futrue,PFT_FUTURE)


##
names(PFT_FUTURE)[1]="lat"
names(PFT_FUTURE)[2]="lon"

data_future=PFT_FUTURE[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))

# add the grid to the land types data
data_future=cbind(data_future,PFT_FUTURE[,3:35])

# different land cover types and the central point of the grid

# get the mean value for the land cover types for each grid

data_mean_future=aggregate(data_future[,c(1,2,7:38)],by=list(data_future$grid),FUN=mean)%>%data.frame()# aggregate the value within a larger grid



# to compare present and futurr land cover types

# for the original data
dkk=numeric()
for(i in c(2:25,27:35))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  dk=summary(lm(data_mean[,i]~data_mean_future[,i]))
  dkk[i]=dk$coefficients[2,1]
}
# to see each type of change that cause changes in the richness

# used the coordinates to estimate the richness and the z value with the regression kriging

#1. to look at the estimated richness in these sites

# get the soil and climate variables
env=matrix(nrow=515,ncol=13)
a1=unique(model_data$plotID)
for (i in 1:515)
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
 env1=subset(model_data,plotID==a1[i])[,c(5:17)]
 env[i,]=unique(env1)%>%as.numeric()
}
env=data.frame(env)

env=cbind(a1,env)%>%data.frame()
names(env)=names(model_data[,c(1,5:17)])

plot_diversity_env=merge(plot_diversity,env,by="plotID")

#data_plot_rich=subset(data_plot_rich,plotID!=c("GUAN_001", "GUAN_002", "GUAN_003", "GUAN_004", "GUAN_006","GUAN_007", "GUAN_042", "GUAN_043", "GUAN_048", "GUAN_049"))

save(plot_diversity_env,file="plot_diversity_env.RData")

## 


lat_grid <- seq(from = 17.96382 , to = 71.24142 , by = 10 / 60)
lon_grid <- seq(from = -156.50286 , to = -66.82463, by = 10 / 60)

variables=c("plot_rich", "lon","lat", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18"  )


data_plot_diversity=plot_diversity_env[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))

data_plot_diversity=cbind(data_plot_diversity,plot_diversity_env[,c(variables)])

# get the mean value of all the variables
#
data_plot_diversity_mean=aggregate(data_plot_diversity[,c(1,2,6,9:18)],by=list(data_plot_diversity$grid),FUN=mean,na.rm=TRUE)

data_plot_diversity_mean=subset(data_plot_diversity_mean,ph!="NaN")# data used for spatial modeling

## construct the model with the plot-level richness

powerTransform(data_plot_diversity_mean$plot_rich)

data_plot_diversity_mean$richbc<-bcPower(data_plot_diversity_mean$plot_rich, 0.2059074 )

mod_plot_rich=data_plot_diversity_mean

##
mod_plot_rich[,5:14]=apply(mod_plot_rich[,5:14],2,range01)

train_plot_rich=mod_plot_rich[,c("richbc","lon","lat")]

pre_plot_rich=mod_plot_rich[,c("richbc", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18")]

## the new data for with the land cover types


#newdata=data_mean[,c("lon","lat")]# data mean is based on the land cover data

newdata=newdata[complete.cases(newdata),]

newloca=newdata[,c(1,2)]
newpre=newdata[,c(3:12)]
newpre=newpre[complete.cases(newpre),]
newpre=apply(newpre,2,range01)

# define the response and explainatory variables

RESPONSE<-train_plot_rich$richbc

mod_predictors<-pre_plot_rich[2:11]

coordinates(train_plot_rich) = ~lon+lat
coordinates(newloca) = ~lon+lat

mc <- makeCluster(detectCores())
registerDoParallel(mc)

myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1808)

GLM_plot<-train(mod_predictors,
                RESPONSE,
                method = "glm",
                
                trControl=myControl,
                preProc=c('center', 'scale'))
print(GLM_plot)


train_plot_rich$residuals.glm<-resid(GLM_plot)


v.glm<-variogram(residuals.glm~ 1, data = train_plot_rich,cutoff=35, width=35/15)

v.glm<-variogram(residuals.glm~ 1, data = train_plot_rich,nint = )

plot(v.glm)

m.glm<-vgm(0.025,"Exp",11,0.012)# based on estimation of the semivariance

m.f.glm<-fit.variogram(v.glm, m.glm)
m.f.glm


plot(v.glm, pl=F, 
     model=m.f.glm,
     col="black", 
     cex=0.9, 
     lwd=0.5,
     lty=1,
     pch=19,
     main="Variogram and Fitted Model\n Residuals of GLM model",
     xlab="Distance (m)",
     ylab="Semivariance")

###
newpre=data.frame(newpre)

names(newpre)=c("bio1","bio2","bio4","bio12","bio15","bio18","ph","organicCPercent","nitrogen","sand")

newloca$GLM_plot <- predict(GLM_plot, as.matrix(newpre))# the names of the predictors should be exactly the same as the predictors

SK.GLM_plot<-krige(residuals.glm~ 1, 
                   loc=train_plot_rich,        # Data frame
                   newdata=newloca,     # Prediction location
                   model = m.f.glm,     # fitted varigram model
                   beta = 0)    


newloca$SK.GLM_plot<-SK.GLM_plot$var1.pred

# Add RF predicted + SK predicted residuals

newloca$RK.GLM.bc<-(newloca$GLM_plot+newloca$SK.GLM_plot)

k1<-1/0.2059074 

newloca$RK.GLM <-((newloca$RK.GLM.bc *0.9882336 +1)^k1)

summary(newloca)


attributes_df <- data.frame(newloca$RK.GLM)

coordinates_df <- data.frame(coordinates(newloca))

pred_plot=cbind(coordinates_df,attributes_df)
names(pred_plot)[3]="rich"



ggplot(pred_plot) +
  geom_point(data=pred_plot,pch=15,aes(x=lon, y=lat,color=rich), size=0.275)+
  scale_color_gradient(expression(Richness),low = "blue", high = "yellow")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted core-level richness")+
  ylab("")+
  xlim(-175,-42)








newdata$lat=-1*newdata$lat
# add the variables of these sites

r_present <- raster::getData("worldclim",var="bio",res=10)

r_present <- r_present[[c(1,2,4,12,15,18)]]

names(r_present) <- c("mat_celsius","MDR", "temp_seasonality","map_mm","pre_seas","pre_wq")

# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$mat_celsius/10
r_present$temp_seasonality <- r_present$temp_seasonality/1000# no need for further transformation

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
r_cec <- raster("cec_5-15cm_mean_5000.tif")

r_cec_reproj <- projectRaster(r_cec, crs = crs(r_present_northam))
r_cec_northam <- raster::mask(raster::crop(r_cec_reproj, north_america_cropped), north_america_cropped)
r_cec_northam_resample <- raster::resample(r_cec_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_cec_northam_resample / 10) # need to know why 10 but all variables would be standardized


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



climate <- list()
for (i in 1:10) {
  climate[[i]] <- raster::extract(r_present_northam[[cl[i]]], newdata[, 1:2])
}

plot_loca_all_soil <- cbind(newdata[,1:2], climate[[1]], climate[[2]],climate[[3]], climate[[4]],climate[[5]], climate[[6]],climate[[7]], climate[[8]],climate[[9]], climate[[10]])

names(plot_loca_all_soil) <- c("lon", "lat" ,  cl)

newdata=plot_loca_all_soil

##3








# to see if there are any columns with all values equal zero

colSums(data_mean[,2:35],na.rm = TRUE)

colSums(data_mean[,2:35],na.rm = TRUE)








## get the change of specific regions

# to see how different types change across times

# we just look at the change in forest, grassland and shrub lands
# the remaining are croplands

# to see change for each of the land cover types for all the land use types




cell_values_future <- cbind(coords_futrue,cell_values_future)

cell_values_comp=cbind(coords_present,cell_values[,"PFT1"],cell_values_future[,"PFT1"])%>%data.frame()%>%mutate(change=V3-V4)

my_colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)

ggplot()+
  geom_point(data=data_mean,aes(x=lon,y=lat),color="red",size=1)

  scale_color_gradient("Change (%)",low = "blue", high = "yellow",na.value="white")+
  


ggplot()+
  #scale_color_gradientn(colours = my_colormap,name="change(%)",na.value = "white")+
  scale_color_gradient("Change (%)",low = "blue", high = "yellow",na.value="white")+
  geom_point(data=data_mean,pch=21,aes(x=lon,y=-lat,color=X6/10))+
  
theme(legend.position =c(0.8,0.325), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted land cover change")+
  ylab("")
# for future secenarios
  
  ggplot()+
    geom_point(data=data_mean,aes(x=lon,y=lat),color="red",size=1)
  
  scale_color_gradient("Change (%)",low = "blue", high = "yellow",na.value="white")+
    
    
    
    ggplot()+
    #scale_color_gradientn(colours = my_colormap,name="change(%)",na.value = "white")+
    scale_color_gradient("Change (%)",low = "blue", high = "yellow",na.value="white")+
    geom_point(data=data_mean_future,pch=21,aes(x=lon,y=-lat,color=X6/10))+
    
    theme(legend.position =c(0.8,0.325), 
          legend.key.size = unit(0.15, "inches"),
          guides(color = guide_legend(nrow = 2, byrow = TRUE)),
          legend.title = element_text(size=8),
          text = element_text(size = 18), 
          legend.text = element_text(size=8),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(hjust = 0), 
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18),
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "NA"), 
          panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
    xlab("Predicted land cover change")+
    ylab("")




plot(flip(raster1[[2]]))

nc_file <- nc_open("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc")

pft4 <- ncvar_get(nc_file, "PFT4")# the fourth type of vegetation types

dim(pft4)# the dimension of the variable "PFT4"


raster1 <- rast(nc_file)





names(nc_file$var)

input_data="\\Users\\luowenqi\\soil-sar\\plot-sar-permutation\\GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc"

ra=stack(nc_file,varname="PFT4")


nc_data=nc_file

lat_var_name <- "latitude"
lon_var_name <- "longitude"

# Extract the latitude and longitude coordinate variables
latitude <- ncvar_get(nc_file, lat_var_name)
longitude <- ncvar_get(nc_file, lon_var_name)


bd=cbind(latitude,longitude)# the coordinates of the data

bd=data.frame(bd)

bd=melt(bd)


pft4 <- ncvar_get(nc_data, "PFT4")# the fourth type of vegetation types

# cbinding the coordinates and the variables within that

type_spa=cbind(bd,t(pft4))# what does the 3600 values indicate




fill.value <- ncatt_get(nc_data, "PFT4", "_FillValue")

pft4[pft4== fill.value$value] = NA

nc_close(nc_data)

r <- raster(t(pft4))

xmin(r) <- longitude[1]
xmax(r) <- longitude[length(longitude)]
ymin(r) <- latitude[length(latitude)]
ymax(r) <- latitude[1]
crs(r) <- "epsg:4326"

crs(r) <- "epsg:3857"

plot(r)

dset01_df <- as.data.frame(r,xy = TRUE)
head(dset01_df)

dset01_df  <-  dplyr::rename(dset01_df,long = x,lat=y)

dset01_df_nona <- dset01_df %>% filter(!is.na(layer))


ggplot()+
  geom_point(data=dset01_df ,aes(x=long,y=lat,fill=layer))
