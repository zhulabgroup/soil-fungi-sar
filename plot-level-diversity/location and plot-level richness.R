library(sf)
library(rnaturalearth)
library(plyr)
library(dplyr)
library(gstat)
library(raster)
library(ggplot2)
library(car)
library(classInt)
library(RStoolbox)
library(caret)
library(caretEnsemble)
library(doParallel)
library(gridExtra)
library(terra)

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
load("model_data.RData")# environmental variables
load("~/soil-sar/plot-level-diversity/plot-diversity.RData")
load("~/soil-sar/plot-sar-permutation/zvalue_all_env.RData")# 457 plot included

1# determine the location of the new grids, which was based on the land use change data


raster1 <- rast("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc")
ext(raster1)=c(-90,90,-180,180)
crs(raster1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent(-72,-18,-170,-55), "SpatialPolygons")

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

PFT=cbind(PFT,rowSums(PFT))%>%data.frame()
PFT=cbind(coords_present,PFT)

## aggregate the cells into a larger cell based on the resolution
# 

names(PFT)[1]="lat"
names(PFT)[2]="lon"

# define the range of the grids
min_lat = -72  # Minimum latitude
max_lat = -18 # Maximum latitude
min_lon = -170 # Minimum longitude
max_lon =-55 # Maximum longitude
# Create a grid of latitude and longitude points
lat_grid <- seq(from = min_lat, to = max_lat, by = 10 / 60)
lon_grid <- seq(from = min_lon, to = max_lon, by = 10 / 60)

data=PFT[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))%>%bind_cols(PFT[,3:35])

# add the grid to the land types data
data=cbind(data,PFT[,3:35])

# different land cover types and the central point of the grid

data=data[complete.cases(data),]

# get the mean value for the land cover types for each grid

#data_mean=aggregate(data[,c(1,2,7:38)],by=list(data$grid),FUN=mean)%>%data.frame()# aggregate the value within a larger grid

data_mean=aggregate(data[,c(1,2,7:38)],by=list(data$grid),FUN=mean,na.rm=TRUE)%>%data.frame()# aggregate the value within a larger grid

data_mean=data_frame(data_mean)

newdata=data_mean[,2:3]
newdata$lat=-1*newdata$lat
# add the variables of these sites

2# extract the climate and soil variables of the new sites


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

# for the soil pH (need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph)="ph"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized

# for the soil soc

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc)="organicCPercent"
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

cl <- c("bio1",   "bio2"  ,   "bio4" , "bio12" , "bio15" ,"bio18", "ph"  , "organicCPercent" , "nitrogen" ,     "sand"  )


climate <- list()
for (i in 1:10) {
  climate[[i]] <- raster::extract(r_present_northam[[cl[i]]], newdata[, 1:2])
}

plot_loca_all_soil <- cbind(newdata[,1:2], climate[[1]], climate[[2]],climate[[3]], climate[[4]],climate[[5]], climate[[6]],climate[[7]], climate[[8]],climate[[9]], climate[[10]])

names(plot_loca_all_soil) <- c("lon", "lat" ,  cl)

newdata=plot_loca_all_soil

3# to get the mean value of the z within each grid

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

## to aggregate the plots with estimated z at the 10-min resolution

lat_grid <- seq(from = 17.96382 , to = 71.24142 , by = 10 / 60)
lon_grid <- seq(from = -156.50286 , to = -66.82463, by = 10 / 60)

variables=c("plot_rich", "lon","lat", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18"  )

data_plot_diversity=plot_diversity_env[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))
data_plot_diversity=cbind(data_plot_diversity,plot_diversity_env[,c(variables)])

# get the mean value of all the variables
data_plot_diversity_mean=aggregate(data_plot_diversity[,c(1,2,6,9:18)],by=list(data_plot_diversity$grid),FUN=mean,na.rm=TRUE)

data_plot_diversity_mean=subset(data_plot_diversity_mean,ph!="NaN")# data used for spatial modeling

4.## construct the model with the plot-level richness

powerTransform(data_plot_diversity_mean$plot_rich)

data_plot_diversity_mean$richbc<-bcPower(data_plot_diversity_mean$plot_rich, 0.2059074 )

mod_plot_rich=data_plot_diversity_mean

##
mod_plot_rich[,5:14]=apply(mod_plot_rich[,5:14],2,range01)

train_plot_rich=mod_plot_rich[,c("richbc","lon","lat")]

pre_plot_rich=mod_plot_rich[,c("richbc", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18")]

## the new data for with the land cover types

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

v.glm<-variogram(residuals.glm~ 1, data = train_plot_rich)

plot(v.glm)

m.glm<-vgm(3,"Exp",11,5.8)# based on estimation of the semivariance

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

newloca$GLM_plot <- predict(GLM_plot, newpre)# the names of the predictors should be exactly the same as the predictors

SK.GLM_plot<-krige(residuals.glm~ 1, 
                   loc=train_plot_rich,        # Data frame
                   newdata=newloca,     # Prediction location
                   model = m.f.glm,     # fitted varigram model
                   beta = 0)    


newloca$SK.GLM_plot<-SK.GLM_plot$var1.pred

# Add RF predicted + SK predicted residuals

newloca$RK.GLM.bc<-(newloca$GLM_plot+newloca$SK.GLM_plot)

k1<-1/0.2059074 

newloca$RK.GLM <-((newloca$RK.GLM.bc *0.2059074  +1)^k1)

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

## 

6.# for the prediction of the z values 

d=sample_data(rare_all)%>%data.frame()
d=d[,c("plotIDM","lon","lat")]
names(d)[1]="plotID"
d=unique(d)

d=aggregate(d[,2:3],by=list(d$plotID),FUN=mean)
names(d)[1]="plotID"

zvalue_all_env=merge(zvalue_all_env,d,by="plotID")

# to aggregated the plots

lat_grid <- seq(from = 18.01780  , to = 68.66407 , by = 10 / 60)
lon_grid <- seq(from = -149.52853 , to = -67.06682, by = 10 / 60)

variables=c("zvalue", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18"  )


data_z_model=zvalue_all_env[,c("lon","lat")]%>%mutate(Grid_lon = cut(lon, breaks = lon_grid,dig.lab = 5), Grid_lat = cut(lat, breaks = lat_grid,dig.lab = 5))%>%mutate(grid=paste(Grid_lat,"*",Grid_lon))

data_z_model=cbind(data_z_model,zvalue_all_env[,variables])

# get the mean values

data_z_model_mean=aggregate(data_z_model[,c(1,2,6,7:16)],by=list(data_z_model$grid),FUN=mean,na.rm=TRUE)


# construe the model for the z value


### to try the own data

powerTransform(data_z_model_mean$zvalue)

data_z_model_mean$zbc<-bcPower(data_z_model_mean$zvalue, -1.81607  )


train.xy_z=data_z_model_mean[,c("zbc","lon","lat")]

train.df_z=data_z_model_mean[,c("zbc", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18")]


grid.xy_z=newdata[,c(1,2)]

grid.df_z=newdata[,c(3:12)]


RESPONSE<-train.df_z$zbc
train.x_z<-train.df_z[,2:11]

coordinates(train.xy_z) = ~lon+lat
coordinates(grid.xy_z) = ~lon+lat

mc <- makeCluster(detectCores())
registerDoParallel(mc)


myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1856)
GLM_z<-train(train.x_z,
             RESPONSE,
             method = "glm",
             trControl=myControl,
             preProc=c('center', 'scale'))
print(GLM_z)

train.xy_z$residuals.glm<-resid(GLM_z)

v.glm<-variogram(residuals.glm~ 1, data = train.xy_z,cutoff=35, width=35/15)
v.glm<-variogram(residuals.glm~ 1, data = train.xy_z)

plot(v.glm)

m.glm<-vgm(0.058,"Exp",10,0.09)

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

names(grid.df_z)[8]="organicCPercent"

grid.xy_z$GLM_z <- predict(GLM_z, grid.df_z)# the x and y disappear

SK.GLM_z<-krige(residuals.glm~ 1, 
                loc=train.xy_z,        # Data frame
                newdata=grid.xy_z,     # Prediction location
                model = m.f.glm,     # fitted varigram model
                beta = 0)    


grid.xy_z$SK.GLM_z<-SK.GLM_z$var1.pred
# Add RF predicted + SK preedicted residuals
grid.xy_z$RK.GLM.bc<-(grid.xy_z$GLM_z+grid.xy_z$SK.GLM_z)

k1<-1/-1.81607   

grid.xy_z$RK.GLM <-((grid.xy_z$RK.GLM.bc *-1.81607  +1)^k1)

summary(grid.xy_z)



attributes_df <- data.frame(grid.xy_z$RK.GLM)

coordinates_df <- data.frame(coordinates(grid.xy_z))

pred_zvalue=cbind(coordinates_df,attributes_df)
names(pred_zvalue)[3]="zvalue"



ggplot(pred_zvalue) +
  geom_point(data=pred_zvalue,pch=15,aes(x=lon, y=lat,color=zvalue), size=0.275)+
  scale_color_gradient(expression(italic(Z)*" value"),low = "blue", high = "yellow")+
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
