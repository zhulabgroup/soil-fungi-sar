# spatial 


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
####


train<-read.csv("train_data.csv", header= TRUE) 

state<-shapefile(paste0(dataFolder,"GP_STATE.shp"))

grid<-read.csv( "GP_prediction_grid_data.csv", header= TRUE) 

First, we will create a data.frame with SOC and continuous environmental data.

powerTransform(train$SOC)

train_data=
  state<-shapefile("GP_STATE.shp")

state<-shapefile(paste0("/Users/luowenqi/soil-sar/plot-sar-permutation/GP_STATE.shp"))

powerTransform(train$SOC)

train$SOC.bc<-bcPower(train$SOC, 0.2523339)

train.xy<-train[,c(1,24,8:9)]
train.df<-train[,c(1,24,11:21)]

grid.xy<-grid[,c(1,2:3)]
grid.df<-grid[,c(4:14)]

RESPONSE<-train.df$SOC.bc
train.x<-train.df[3:13]

coordinates(train.xy) = ~x+y
coordinates(grid.xy) = ~x+y


mc <- makeCluster(detectCores())
registerDoParallel(mc)

myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1856)
GLM<-train(train.x,
           RESPONSE,
           method = "glm",
           trControl=myControl,
           preProc=c('center', 'scale'))
print(GLM)



# Extract residuals
train.xy$residuals.glm<-resid(GLM)
# Variogram
v.glm<-variogram(residuals.glm~ 1, data = train.xy,cutoff=300000, width=300000/15)

# Intial parameter set by eye esitmation
m.glm<-vgm(0.15,"Exp",40000,0.05)
# least square fit
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

grid.xy$GLM <- predict(GLM, grid.df)


SK.GLM<-krige(residuals.glm~ 1, 
              loc=train.xy,        # Data frame
              newdata=grid.xy,     # Prediction location
              model = m.f.glm,     # fitted varigram model
              beta = 0)   


grid.xy$SK.GLM<-SK.GLM$var1.pred
# Add RF predicted + SK preedicted residuals
grid.xy$RK.GLM.bc<-(grid.xy$GLM+grid.xy$SK.GLM)


k1<-1/0.2523339                                   
grid.xy$RK.GLM <-((grid.xy$RK.GLM.bc *0.2523339+1)^k1)
summary(grid.xy)


GLM<-rasterFromXYZ(as.data.frame(grid.xy)[, c("x", "y", "GLM")])

SK.GLM<-rasterFromXYZ(as.data.frame(grid.xy)[, c("x", "y", "SK.GLM")])
RK.GLM.bc<-rasterFromXYZ(as.data.frame(grid.xy)[, c("x", "y", "RK.GLM.bc")])
RK.GLM.SOC<-rasterFromXYZ(as.data.frame(grid.xy)[, c("x", "y", "RK.GLM")])















glm1<-ggR(GLM, geom_raster = TRUE) +
  scale_fill_gradientn("", colours = c("orange", "yellow", "green",  "sky blue","blue"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle("GLM Predicted (BoxCox)")+
  theme(plot.title = element_text(hjust = 0.5))



### to try the own data

powerTransform(obs_aggre$meanz)

obs_aggre$zbc<-bcPower(obs_aggre$meanz, -1.847084  )


train.xy_z=obs_aggre[,c(1,7,2:3)]

train.df_z=obs_aggre[,c(1,7,5,6)]

grid.xy_z=newdata[,c(1,2)]

grid.df_z=newdata[,c(3,4)]

RESPONSE<-train.df_z$zbc
train.x_z<-train.df_z[3:4]

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

v.glm<-variogram(residuals.glm~ 1, data = train.xy_z,cutoff=40, width=40/15)

m.glm<-vgm(0.008,"Sph",20,0.05)

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


grid.xy_z$GLM_z <- predict(GLM_z, grid.df_z)# the x and y disappear

SK.GLM_z<-krige(residuals.glm~ 1, 
              loc=train.xy_z,        # Data frame
              newdata=grid.xy_z,     # Prediction location
              model = m.f.glm,     # fitted varigram model
              beta = 0)    


grid.xy_z$SK.GLM_z<-SK.GLM_z$var1.pred
# Add RF predicted + SK preedicted residuals
grid.xy_z$RK.GLM.bc<-(grid.xy_z$GLM_z+grid.xy_z$SK.GLM_z)

k1<-1/-1.847084   

grid.xy_z$RK.GLM <-((grid.xy_z$RK.GLM.bc *-1.847084 +1)^k1)

summary(grid.xy_z)

GLM0<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "RK.GLM ")])


##############

set.seed(66771)
train_set <- sample(c(TRUE, FALSE), size=dim(obs)[1],
                    prob=c(0.70, 0.30), replace=TRUE)

test_set=!train_set

train_set[!train_set]=NA
mm=train_set

which(!is.na(mm))

train_data=obs[which(!is.na(mm)),]

## for the test data
test_set[!test_set]=NA
mmt=test_set

which(!is.na(mmt))
test_data=obs[which(!is.na(mmt)),]


set.seed(56885)
dd=matrix(ncol=2,nrow=1000)
for (i in 1:1000){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  mod_cv <- cv.glmnet( as.matrix(train_data[,c(3:12)]), train_data[,"meanz"], family='gaussian', intercept = F, alpha=1)
  
  dd[i,1]=mod_cv$lambda.min
  dd[i,2]=mod_cv$lambda.1se
}
dd=data.frame(dd)

head(dd)

# will have less variables selected
# the best lambda is somewhat variable


lam=unique(dd$X1)
R=numeric()
term=list()
for (i in 1:length(lam)){
  
  la.eq <- glmnet( train_data[,c(3:12)], train_data[,"meanz"],lambda=lam[i], family='gaussian', intercept = F, alpha=1)
  matrix(coef(la.eq ))[2:11,]
  terms=matrix(nrow=10,ncol=2)
  terms[,1]=colnames(train_data[,c(3:12)])
  terms[,2]=coef(la.eq )[2:11,]
  terms=data.frame(terms)
  best_model <- glmnet( train_data[,c(subset(terms,X2!=0)[,"X1"])], train_data[,"meanz"], 
                        family='gaussian', intercept = F, alpha=1,lambda = lam[i],final.re=TRUE)
  pre_zvalue=predict(best_model,lambda = lam[i],newx=as.matrix(test_data[,c(subset(terms,X2!=0)[,"X1"])]))
  
  a4=summary(lm(pre_zvalue~test_data[,"meanz"]))
  R[i]=a4$adj.r.squared
  term[[i]]=terms
}
# the best lamda is  0.002165107






### to try the own data with all the plots

# with all the predictors included


head(zc_nest_mean_var)

# based on the selected variables, ten variables were included in the model



zc_nest_mean_var$meanz=as.numeric(zc_nest_mean_var$meanz)
zc_nest_mean_var$meanc=as.numeric(zc_nest_mean_var$meanc)

predictors=c("soc", "ph"  ,  "nitrogen"  ,   "sand"  ,  "mat_celsius","MDR", "temp_seasonality","map_mm","pre_seas","pre_wq"  )

obs=zc_nest_mean_var[,c("plotID","meanz",predictors)]

# add the coordinates of each plot

obs=merge(obs,coordinate,by="plotID")


obs=obs[complete.cases(obs),]

# based on the glm

powerTransform(obs$meanz)

obs$zbc<-bcPower(obs$meanz, -1.284961)# simply a way to make it more normal

#standardize the data
obs[,predictors]=apply(obs[,predictors],2,range01)# 



names(obs)=c("plotID","meanz",  "soc","ph" ,"nitrogen","sand", "mat_celsius","MDR", "temp_seasonality","map_mm","pre_seas","pre_wq","lon","lat","zbc")



train.xy_z=obs[,c("zbc","lon","lat")]

train.df_z=obs[,c("zbc",predictors)]

# add the 10 predictors for each of the new location


newdata=no_na_prediction[,1:2]

# coordinates for the new locations
# get soil variables 
# create the map


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



newdata=newdata[which(complete.cases(newdata)),]

newdata[,3:12]=apply(newdata[,3:12],2,range01)


# the new locations and the predictors

grid.xy_z=newdata[,c(1,2)]

grid.df_z=newdata[,c(3:12)]



# for the observed values
RESPONSE<-train.df_z$zbc
train.x_z<-train.df_z[2:11]

coordinates(train.xy_z) = ~lon+lat
coordinates(grid.xy_z) = ~lon+lat

mc <- makeCluster(detectCores())
registerDoParallel(mc)


myControl <- trainControl(method="repeatedcv", #"repeatedcv" for repeated cross-validation.
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1856)
GLM_z<-train(train.x_z,
             RESPONSE,
             method = "glm",
             trControl=myControl,
             preProc=c('center', 'scale'))
print(GLM_z)# the model residuals and then to look at how this varies among distances

#RMSE-0.13
#Rsquared 0.16

train.xy_z$residuals.glm<-resid(GLM_z)# look at the 

dk=variogram(zbc~1,data = train.xy_z)

dk=variogram(residuals.glm~1,data = train.xy_z)


v.glm<-variogram(residuals.glm~ 1, data = train.xy_z,cutoff=40, width=40/15)

m.glm<-vgm(0.015,"Exp",15,0.025)

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


grid.xy_z$GLM_z <- predict(GLM_z, grid.df_z)# the predicted values for new locations

SK.GLM_z<-krige(residuals.glm~ 1, 
                loc=train.xy_z,        # Data frame
                newdata=grid.xy_z,     # Prediction location
                model = m.f.glm,     # fitted varigram model
                beta = 0)    


grid.xy_z$SK.GLM_z<-SK.GLM_z$var1.pred

# Add RF predicted + SK predicted residuals
grid.xy_z$RK.GLM.bc<-(grid.xy_z$GLM_z+grid.xy_z$SK.GLM_z)

k1<-1/-1.284961   

grid.xy_z$RK.GLM <-((grid.xy_z$RK.GLM.bc *-1.284961  +1)^k1)

summary(grid.xy_z)

GLM0<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "RK.GLM")][1:224,])#some rows can not be converted

## extract the values and the coordinates

attributes_df <- data.frame(grid.xy_z$RK.GLM)
coordinates_df <- data.frame(coordinates(grid.xy_z))

predz=cbind(coordinates_df,attributes_df)
names(predz)[3]="zvalue"


ggplot(predz) +
  geom_point(data=predz,pch=15,aes(x=lon, y=lat,color=zvalue), size=0.275)+
  scale_color_gradient(expression(Z),low = "blue", high = "yellow")+
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
  xlab("")+
  ylab("")+
  xlim(-175,-42)










# for the random forest

set.seed(1856)
mtry <- sqrt(ncol(train.x_z))             # number of variables randomly sampled as candidates at each split.
tunegrid.rf <- expand.grid(.mtry=mtry)
RF<-train(train.x_z,
          RESPONSE,
          method = "rf",
          trControl=myControl,
          tuneGrid=tunegrid.rf,
          ntree= 100,
          preProc=c('center', 'scale'))
print(RF)# how this function works?

train.xy_z$residuals.rf<-resid(RF)

v.rf<-variogram(residuals.rf~ 1, data = train.xy_z)

m.rf<-vgm(0.5,"Exp",10,0.05)

m.f.rf<-fit.variogram(v.rf, m.rf)
m.f.rf


plot(v.rf, pl=F, 
     model=m.f.rf,
     col="black", 
     cex=0.9, 
     lwd=0.5,
     lty=1,
     pch=19,
     main="Variogram and Fitted Model\n Residuals of RF model",
     xlab="Distance (m)",
     ylab="Semivariance")



grid.xy_z$RF <- predict(RF, grid.df_z)


SK.RF<-krige(residuals.rf~ 1, 
             loc=train.xy_z,        # Data frame
             newdata=grid.xy_z,     # Prediction location
             model = m.f.rf,      # fitted varigram model
             beta = 0)            # residuals from a trend; expected value is 0     



grid.xy_z$SK.RF<-SK.RF$var1.pred
# Add RF predicted + SK preedicted residuals
grid.xy_z$RK.RF.bc<-(grid.xy_z$RF+grid.xy_z$SK.RF)


k1<-1/-1.284961                                   
grid.xy_z$RK.RF <-((grid.xy_z$RK.RF.bc *-1.284961 +1)^k1)
summary(grid.xy_z)


GLM<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "GLM_z")])

SK.GLM<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "SK.GLM")])
RK.GLM.bc<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "RK.GLM.bc")])
RK.GLM.SOC<-rasterFromXYZ(as.data.frame(grid.xy_z)[, c("lon", "lat", "RK.GLM")])


## assign all the plot into different grids at a larger spatial scale

# devide the grid at a resolution of 10 by 10

xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 

data=obs[,c("lon","lat")]
data <- mutate(data,Grid_Lat = cut(lon, breaks = xgrid), Grid_Lon = cut(lat, breaks = ygrid))

data=mutate(data,grid=paste(data$Grid_Lat,"*",data$Grid_Lon))

data=cbind(data,obs)


# get the mean value of the of the variables at each grid

data=aggregate(data[,c(1,2,7:17)],by=list(data$grid),FUN=mean)# aggregated at the 10 by 10 scale

# run the model when different plots were aggregated at one site

powerTransform(data$meanz)

data$zbc<-bcPower(data$meanz, -2.472824   )

train.xy_agg=data[,c("zbc","lon","lat")]

train.df_agg=data[,c("zbc",predictors)]

grid.xy_agg=newdata[,c(1,2)]

grid.df_agg=newdata[,c(3:12)]


RESPONSE<-train.df_agg$zbc
train.x_agg<-train.df_agg[2:11]

coordinates(train.xy_agg) = ~lon+lat
coordinates(grid.xy_agg) = ~lon+lat

mc <- makeCluster(detectCores())
registerDoParallel(mc)


myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1856)
GLM_agg<-train(train.x_agg,
             RESPONSE,
             method = "glm",
             trControl=myControl,
             preProc=c('center', 'scale'))
print(GLM_agg)

train.xy_agg$residuals.glm<-resid(GLM_agg)

dk=variogram(zbc~1,data = train.xy_agg)

v.glm<-variogram(residuals.glm~ 1, data = train.xy_agg,cutoff=35, width=40/15)

m.glm<-vgm(0.015,"Exp",8,0.032)# based on estimation of the semivariance

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


grid.xy_agg$GLM_agg <- predict(GLM_agg, grid.df_agg)# the x and y disappear

SK.GLM_agg<-krige(residuals.glm~ 1, 
                loc=train.xy_agg,        # Data frame
                newdata=grid.xy_agg,     # Prediction location
                model = m.f.glm,     # fitted varigram model
                beta = 0)    

###

grid.xy_agg$SK.GLM_agg<-SK.GLM_agg$var1.pred

# Add RF predicted + SK preedicted residuals
grid.xy_agg$RK.GLM.bc<-(grid.xy_agg$GLM_agg+grid.xy_agg$SK.GLM_agg)

k1<-1/-2.472824 

grid.xy_agg$RK.GLM <-((grid.xy_agg$RK.GLM.bc *-2.472824   +1)^k1)

summary(grid.xy_agg)



## extract the values and the coordinates

attributes_df <- data.frame(grid.xy_agg$RK.GLM)
coordinates_df <- data.frame(coordinates(grid.xy_agg))

predz_agg=cbind(coordinates_df,attributes_df)
names(predz_agg)[3]="zvalue"


ggplot(predz_agg) +
  geom_point(data=predz_agg,pch=15,aes(x=lon, y=lat,color=zvalue), size=0.275)+
  scale_color_gradient(expression(Z),low = "blue", high = "yellow")+
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
  xlab("")+
  ylab("")+
  xlim(-175,-42)

## based on the random forest approach

set.seed(1856)
mtry <- sqrt(ncol(train.x_agg))             # number of variables randomly sampled as candidates at each split.
tunegrid.rf <- expand.grid(.mtry=mtry)
RF<-train(train.x_agg,
          RESPONSE,
          method = "rf",
          trControl=myControl,
          tuneGrid=tunegrid.rf,
          ntree= 100,
          preProc=c('center', 'scale'))
print(RF)



train.xy_agg$residuals.rf<-resid(RF)
# Variogram
v.rf<-variogram(residuals.rf~ 1, data = train.xy_agg)
# Intial parameter set by eye esitmation
m.rf<-vgm(0.015,"Exp",8,0.032)
# least square fit
m.f.rf<-fit.variogram(v.rf, m.rf)
m.f.rf

plot(v.rf, pl=F, 
     model=m.f.rf,
     col="black", 
     cex=0.9, 
     lwd=0.5,
     lty=1,
     pch=19,
     main="Variogram and Fitted Model\n Residuals of RF model",
     xlab="Distance (m)",
     ylab="Semivariance")

grid.xy_agg$RF <- predict(RF, grid.df_agg)

SK.RF<-krige(residuals.rf~ 1, 
             loc=train.xy_agg,        # Data frame
             newdata=grid.xy_agg,     # Prediction location
             model = m.f.rf,      # fitted varigram model
             beta = 0)   

grid.xy_agg$SK.RF<-SK.RF$var1.pred
# Add RF predicted + SK predicted residuals
grid.xy_agg$RK.RF.bc<-(grid.xy_agg$RF+grid.xy_agg$SK.RF)

k1<-1/-2.472824 

grid.xy_agg$RK.RF <-((grid.xy_agg$RK.RF.bc *-2.472824 +1)^k1)

summary(grid.xy_agg)

####
attributes_df <- data.frame(grid.xy_agg$RK.RF)

coordinates_df <- data.frame(coordinates(grid.xy_agg))

predz_agg_rf=cbind(coordinates_df,attributes_df)
names(predz_agg_rf)[3]="zvalue"


ggplot(predz_agg_rf) +
  geom_point(data=predz_agg_rf,pch=15,aes(x=lon, y=lat,color=zvalue), size=0.275)+
  scale_color_gradient(expression(Z),low = "blue", high = "yellow")+
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
  xlab("")+
  ylab("")+
  xlim(-175,-42)











## look at the site based z values

siteid=matrix(ncol=1,nrow=483)
for (i in 1:483)
  {
  string_length <- nchar(obs$plotID[i])
  if(string_length<4)
  {
    siteid[i,]=substr(obs$plotID[i],1,3)
  }
  else{
    siteid[i,]=substr(obs$plotID[i],1,4) 
  }
  
}

obs=cbind(siteid,obs)

obs_aggre=aggregate(obs[,c(3,4,5,6,7)],by=list(obs$siteid),FUN=mean)



ggR <- function(img, layer = 1, maxpixels = 500000,  alpha = 1, hue = 1, sat = 0, stretch = "none", quantiles = c(0.02,0.98), 
                coord_equal = TRUE, ggLayer=FALSE, ggObj = TRUE, geom_raster = FALSE, forceCat = FALSE) {
  
  img   <- .toRaster(img)
  layer <- unlist(.numBand(img, layer))
  
  multLayers <- if (length(layer)>1) TRUE else FALSE
  if(multLayers & !geom_raster & ggObj) {
    warning("You asked for multiple layers but geom_raster is FALSE.",
            "\ngeom_raster will be reset to TRUE", 
            "\nHint: in case you're looking for a grayscale and facetted plot, use:",
            "\nggR(img, ..., geom_raster=TRUE)+scale_fill_gradientn(colors = grey.colors(100))",
            call. = FALSE)
    geom_raster <- TRUE
  }
  annotation <- !geom_raster
  
  xfort <- sampleRegular(img[[layer]], maxpixels, asRaster = TRUE)
  ex <- extent(xfort)
  dimImg <- dim(xfort)
  df <- lapply(names(xfort), function(layer) {
    df    <- data.frame(as.data.frame(xfort[[layer]],  xy = TRUE), 
                        layerName = factor(layer, levels = names(xfort)))
    colnames(df) <- c("x", "y", "value", "layerName") 
    df
  })
  df <- do.call(rbind, df)
  
  if(forceCat & !is.factor(df[,"value"])) df[,"value"] <- as.factor(df[,"value"])
  
  if(is.character(df[,"value"])) df[,"value"] <- factor(df[,"value"])
  fac <- is.factor(df[,"value"]) 
  if(fac & (annotation | !ggObj)) {
    .vMessage("img values are factors but annotation is TRUE. Converting factors as.numeric.")
    levelLUT   <- levels(df[,"value"])
    df[,"value"] <- as.numeric(df[,"value"])
  }
  
  if(!fac & stretch != "none")  {
    for(layer in levels(df$layerName)) {
      df[df$layerName==layer,"value"] <- .stretch(df[df$layerName==layer,"value"], method = stretch, quantiles = quantiles)    
    }
  }
  if(!(ggObj & !annotation)  ){
    ## Annotation processing
    ## Avoid rescaling "value"s with single values
    df <- lapply(levels(df$layerName), function(layer) {
      normVals    <- suppressWarnings(rescaleImage(df[df$layerName==layer,"value"], ymin = 0, ymax = 1))  ## suppress warnings due to single value bands. rescaleImage returns NA, which is fine.
      uval        <- unique(df[,"value"])
      if(sum(is.finite(uval)) == 1)   normVals[is.finite(df[,"value"])] <- 1
      nona         <- !is.na(normVals)
      df$fill      <- NA
      df[nona, "fill"] <- hsv(h = hue, s = sat, v = normVals[nona], alpha = alpha)
      df
    })
    df <- do.call(rbind, df)
  }
  x <- y <- value <- NULL  
  if(ggObj) {       
    if(annotation)  {        
      dmat <- matrix(df$fill, nrow=dimImg[1], ncol=dimImg[2], byrow = TRUE)  
      ggl  <- annotation_raster(raster = dmat, xmin = ex[1], xmax = ex[2], ymin = ex[3], ymax = ex[4], interpolate = FALSE)
    } else {
      ggl  <- geom_raster(data = df[,c("x","y","value","layerName")], aes(x = x, y = y, fill = value), alpha = alpha) 
    }
    if(multLayers) facetObj <- facet_wrap(~layerName)
    
    if(ggLayer) {
      if(multLayers) {
        return(list(ggl, facetObj))
      } else { 
        return(ggl) 
      }
    }
    if(annotation) {   
      dummy <- data.frame(x = ex[1:2],y = ex[3:4], layerName = rep(levels(df$layerName), each = 2) )       
      p <- ggplot()  + ggl + geom_blank(data = dummy, aes(x,y))
      if(coord_equal) p <- p + coord_equal()
      if(multLayers) p <- p + facet_wrap(~layerName)
      return(p)
    } else {
      p <- ggplot() + ggl
      if(coord_equal) p <- p + coord_equal()
      if(multLayers) p <- p + facetObj
      return(p)
    }
    
  } else {
    if(fac & (annotation | !ggObj)) df[,"value"] <- factor(levelLUT[df[, "value"]], levels=levelLUT)
    return(df)
  }
}



.toRaster <- function(x) {
  if (inherits(x, "SpatRaster")) {
    return(stack(x))
  } else {
    return(x)
  }
}

.numBand <- function(raster, ...){
  bands <- list(...)
  lapply(bands, function(band) if(is.character(band)) which(names(raster) == band) else band ) 
}
