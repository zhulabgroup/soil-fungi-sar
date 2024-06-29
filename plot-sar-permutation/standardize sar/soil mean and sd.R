## try to get the CV and the mean of the soil variables for each 

##
library(phyloseq)
library(dplyr)
library(tidyr)

d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site correspondes to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
rare_all<- merge_phyloseq(rare_all, plotIDM)# m


soil_data=sample_data(rare_all)%>%data.frame()%>%dplyr::select(Site,plotIDM,soilInCaClpH,nitrogenPercent, organicCPercent,soilMoisture)

# get the mean of the variables
a6=unique(soil_data$plotIDM)
soil_mean=matrix(ncol=5,nrow=length(a6))
soil_sd=matrix(ncol=5,nrow=length(a6))
#soil_sd=matrix(ncol=6,nrow=length(a6))
for(i in 1:length(a6))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  dk=soil_data%>%filter(plotIDM==a6[i])
 
  soil_mean[i,]=dk%>%group_by(plotIDM)%>%summarise(across(c(soilInCaClpH, nitrogenPercent, organicCPercent, soilMoisture), mean, na.rm = TRUE))%>%as.numeric()
  soil_sd[i,]=dk%>%group_by(plotIDM)%>%summarise(across(c(soilInCaClpH, nitrogenPercent, organicCPercent, soilMoisture), sd, na.rm = TRUE))%>%as.numeric()
  
  }


soil_mean=soil_mean%>%data.frame()%>%dplyr::select(-X1)%>%rename_all(~paste0(c("soilInCaClpH","nitrogenPercent","organicCPercent","soilMoisture")))%>%mutate(plotID=a6)



# get the coordinates of the plots

plot_coordinate=sample_data(rare_all)%>%data.frame()%>%dplyr::select(lon,lat,plotIDM)%>%distinct()

## load in the soil rasters


plot_coordinate=plot_coordinate%>%dplyr::select(plotIDM,lon,lat)%>%group_by(plotIDM)%>%summarize(across(c(lon, lat), mean, na.rm = TRUE))%>%rename(plotID=plotIDM)

soil_mean=left_join(plot_coordinate,soil_mean,by="plotID")

# fill the NAs with the soil raster



r_present <- raster::getData("worldclim",var="bio",res=10)

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


## get the variables for each of the grid
# add the soil carbon data to the map
# Data downloaded from ISRIC https://files.isric.org/soilgrids/latest/data_aggregated/5000m/phh2o/
# and https://files.isric.org/soilgrids/latest/data_aggregated/5000m/soc/
## ad the soil nitrogen to the raster

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

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


r_nitrogen <- raster("nitrogen_5-15cm_mean_5000.tif")
names(r_nitrogen) <- "nitrogenPercent"
r_nitrogen_reproj <- projectRaster(r_nitrogen, crs = crs(r_present_northam))
r_nitrogen_northam <- raster::mask(raster::crop(r_nitrogen_reproj, north_america_cropped), north_america_cropped)
r_nitrogen_northam_resample <- raster::resample(r_nitrogen_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_nitrogen_northam_resample / 1000)


##
r_cec <- raster("cec_5-15cm_mean_5000.tif")
names(r_cec) <- "cec"

r_cec_reproj <- projectRaster(r_cec, crs = crs(r_present_northam))

r_cec_northam <- raster::mask(raster::crop(r_cec_reproj, north_america_cropped), north_america_cropped)
r_cec_northam_resample <- raster::resample(r_cec_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_cec_northam_resample / 10) # need to know why 10 but all variables would be standardized

# add soil sand
r_sand <- raster("sand_5-15cm_mean_5000.tif")
names(r_sand) <- "sand"
r_sand_reproj <- projectRaster(r_sand, crs = crs(r_present_northam))
r_sand_northam <- raster::mask(raster::crop(r_sand_reproj, north_america_cropped), north_america_cropped)
r_sand_northam_resample <- raster::resample(r_sand_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_sand_northam_resample / 10) # need to know why 10 but all variables would be standardized

## extract the values from these rasters for the points without na
is.na(soil_mean$soilInCaClpH)%>%sum()
is.na(soil_mean$nitrogenPercent)%>%sum()
is.na(soil_mean$organicCPercent)%>%sum()
is.na(soil_mean$ soilMoisture)%>%sum()
# all the nas have been filled, soil moisture have some NAs

soil_mean=soil_mean%>%mutate(organicCPercent = if_else(is.na(organicCPercent), raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent),
                   nitrogenPercent = if_else(is.na(nitrogenPercent), raster::extract(r_present_northam[["nitrogenPercent"]], cbind(lon, lat)), nitrogenPercent),
                   soilInCaClpH = if_else(is.na(soilInCaClpH), raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)), soilInCaClpH))

# the unit for the soilN data for the NEON site is ?

soil_mean%>%mutate(cec =  raster::extract(r_present_northam[["cec"]], cbind(lon, lat)),
                   sand= raster::extract(r_present_northam[["sand"]], cbind(lon, lat)))->soil_mean

save(soil_mean,file="soil_mean.RData")

## the variation is low because of 

bind_cols(soil_mean[,1:3],soil_sd)%>%data.frame()%>%rename_all(~paste0(c("plotIDM","lon","lat","d","soilInCaClpH","nitrogenPercent","organicCPercent","kk","soilMoisture")))->temp

temp%>%dplyr::select(-d,-kk)->temp

# get the sd of the soil variables with kriging

soil_sd=soil_sd%>%data.frame()%>%select(X2,X3,X4,X5)%>%rename_all(~paste0(c("soilInCaClpH","nitrogenPercent","organicCPercent","soilMoisture")))%>%bind_cols(plot_coordinate[,1:3])




install.packages("gstat")
library(gstat)

library(gstat)
library(sp)

data_N=soil_sd%>%select(lon,lat,nitrogenPercent)%>%filter(!is.na(nitrogenPercent))

coordinates(data_N) <- ~lon+lat

variogram <- variogram(nitrogenPercent ~ 1, data_N)

model <- fit.variogram(variogram, model = vgm(1, "Sph", 1, 1))

## add the climate variables to the soil variables

plot_env=model_data%>%dplyr::select(plotID,siteIDD,bio1,bio2,bio4,bio8,bio12,bio15,bio18,richness)%>%distinct()


plot_soil_climate= left_join(soil_mean,plot_env, by="plotID")

model_data_SAR=data_zvalue%>%left_join(plot_soil_climate,by="plotID")

save(model_data_SAR,file="model_data_SAR.RData")

## some plots still don have sand and cec values


na_ind <- which(is.na(model_data_SAR$cec))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["cec"]], xy), is.na(r_present_northam[["cec"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$soilInCaClpH[i] <- values["cec"]
}

na_ind <- which(is.na(model_data_SAR$sand))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["sand"]], xy), is.na(r_present_northam[["sand"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$soilInCaClpH[i] <- values["sand"]
}





## for the full data
model_data_SAR$logc=2.71828^model_data_SAR$logc

model_data_SAR=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc, zvalue, soilMoisture ,soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)

model_data_SAR=model_data_SAR[complete.cases(model_data_SAR),]

data_all=model_data_SAR%>%filter(guild=="all")

  
ggcorrplot(cor(data_all%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#

##nitrogen and soil carbon were correlated 
ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[5]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[6]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[7]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[8]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[9]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
# bio1 and bio4 are correlated and i retained the former





# Find pairs with high correlation

data_all[,c(6:21)]=apply(data_all[,c(6:21)],2,range01)

mod <- lmer(zvalue ~ logc +nitrogenPercent +soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data_all)

mod <- lmer(zvalue ~ logc + soilInCaClpH +organicCPercent+soilMoisture +cec+ (1 |siteIDD), data = data_all)

mod <- glmmLasso(zvalue ~ logc + soilInCaClpH +organicCPercent+soilMoisture +cec+ (1|siteIDD), lambda = 5,data = data_all,family = gaussian(link = "identity"))

data_all$siteIDD=as.factor(data_all$siteIDD)



# when use a lambda of 5, several key climate predictors were selected, and when use 8, only 1-2 were selected
# when using 0.01 all were selected
#when using 0.1
# it looks increasing lambda does not lead to less variables selected by identities of the selected variables

mod <- cv.glmmLasso(fix=zvalue ~ logc + soilInCaClpH +organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                    rnd=list(siteIDD=~1),
                    dat = data[[1]],lambdas=10^seq(-2, 5, by = 0.1),family = gaussian(link = "identity"))#（0.01-10）




data_all$siteIDD=as.factor(data_all$siteIDD)

op=10^seq(-5, 2, by = 0.01)%>%round(digits = 3)

per=numeric()
for (i in 1:length(op)){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  model <- glmmLasso(fix=logc ~ soilInCaClpH +organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
   rnd=list(siteIDD=~1),
   dat = data_all,lambda=op[i],family = gaussian(link = "identity"),final.re = TRUE)#（0.01-10
per[i]=model$aic
  }

df=data.frame(per,op)

ggplot(data=df,aes(x=per,y=op))+
  geom_point(color="red")

plot(per~op)


model <- glmmLasso(fix=zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                   rnd=list(siteIDD=~1),
                   dat = data_all,lambda=7,family = gaussian(link = "identity"),final.re = TRUE)



# select the lambda for the linear mixed effect model

model_data_SAR$siteIDD=as.factor(model_data_SAR$siteIDD)
guild_select=unique(model_data_SAR$guild)

data=list()
for (i in 1:9)
{
  d=model_data_SAR%>%filter(guild==guild_select[i]) 
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
 data[[i]]=d
}


set.seed(235)
sel_vab=list()# the variables selected based on the optimal lambda
lambb=numeric()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  if(i<9)
    {
    mod1=cv.glmmLasso(fix=zvalue ~ logc + soilInCaClpH+organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                      rnd=list(siteIDD=~1), 
                      data = data[[i]],
                      family = gaussian(link = "identity"),
                      kfold = 5)#（0.01-10） 
    lambb[i]=mod1$lambda.min
    kk=mod1$glmmLasso.final$coefficients%>%data.frame()%>%filter(.!=0)
    sel_vab[[i]]=rownames(kk)[2:dim(kk)[1]]
    
  }
  else# here bio1 was excluded
    {
    mod1=cv.glmmLasso(fix=zvalue ~ logc + soilInCaClpH+ nitrogenPercent+organicCPercent+soilMoisture+ cec+sand+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                      rnd=list(siteIDD=~1), 
                      data = data[[i]],
                      family = gaussian(link = "identity"),
                      kfold = 5)#（0.01-10） 
    lambb[i]=mod1$lambda.min
    kk=mod1$glmmLasso.final$coefficients%>%data.frame()%>%filter(.!=0)
    sel_vab[[i]]=rownames(kk)[2:dim(kk)[1]]# the intercept was excluded
  }
  }
  
# based on the select lambda to select the variables with linear mixed effect model 
# extract the effect size of different variables

effect=list()
for (i in 1:9)
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  random_effect <- "(1|siteIDD)"
  response="zvalue"
  fixed_effects=sel_vab[[i]]
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")# the selected variables
  formula_string <- paste(response, "~", fixed_effects_formula, "+", random_effect)%>%as.formula()
  mod <- lmer(formula_string, data = data[[i]])
  mod2=summary(mod)
  effect[[i]]=mod2$coefficients[,c(1,2,5)]%>%data.frame()%>%mutate(var=rownames(mod2$coefficients))
}
## to plot the effect size
## to merge with the full variables
##get the modified effect size

modified_effect=list()
for (i in 1:9)
{
  modified_effect[[i]]=effect[[2]]%>%dplyr::select(var)%>%left_join(effect[[i]],by="var")%>%filter(!var%in%c("logc","(Intercept)"))%>%replace_na(list(Estimate=0))
}




ggplot()+
  geom_point(data=modified_effect[[1]],pch=21,color="black",aes(x= Estimate,y=1:dim(modified_effect[[1]])[1]),size=4,
  fill=rev(c("seagreen1", "royalblue","royalblue","royalblue","royalblue","royalblue", "royalblue","royalblue","peru","peru","peru", "peru","gray")))+
  geom_segment(data=modified_effect[[1]],size=.8,
  aes(x=modified_effect[[1]]$Estimate-1.96*modified_effect[[1]]$Std..Error,
  y=1:dim(modified_effect[[1]])[1],
  xend=modified_effect[[1]]$Estimate+1.96*modified_effect[[1]]$Std..Error,
  yend=1:dim(modified_effect[[1]])[1]),
  color=rev(c("seagreen1", "royalblue","royalblue","royalblue","royalblue","royalblue", "royalblue","royalblue","peru","peru","peru", "peru","peru")))+
  geom_vline(xintercept = 0,color="red",linetype="dashed",size=.8)+
  scale_y_continuous(breaks=1:13,labels=rev(c("Pla.rich", "Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.", "MDR","MAT","Sand","CEC","Moisture", "SoilC","pH")))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlab("Effect size")+
  ylab("")
  

### when only soil and climate variables were considered, more plots were included


effect_guild_specific=bind_rows(modified_effect[[2]][,c(1,2,4)],
          modified_effect[[3]][,c(1,2,4)],
          modified_effect[[4]][,c(1,2,4)],
          modified_effect[[5]][,c(1,2,4)],
          modified_effect[[6]][,c(1,2,4)],
          modified_effect[[7]][,c(1,2,4)],
          modified_effect[[8]][,c(1,2,4)],
          modified_effect[[9]][,c(1,2,4)])%>%mutate(guild=rep(c(guild_select[2:9]),each=13))%>%rename_all(~paste0(c("var","effect","pval","guild")))


effect_guild_specific$var=factor(effect_guild_specific$var,levels=c("soilInCaClpH", "organicCPercent", "soilMoisture" ,   "cec", "sand" , "bio1" ,    "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18" , "richness"))        
effect_guild_specific$guild=factor(effect_guild_specific$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           


effect_guild_specific$pval[effect_guild_specific$pval<0.01]="**"
effect_guild_specific$pval[effect_guild_specific$pval>0.01&effect_guild_specific$pval<0.05]="*"
effect_guild_specific$pval[effect_guild_specific$pval>0.05]=""



ggplot(effect_guild_specific, aes(x =var , y = guild, fill = effect)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  scale_x_discrete(breaks=as.character(unique(effect_guild_specific$var)),
                   labels = c("pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Tem.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" , "Pla.rich."))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15),
        plot.margin = margin(b=-0.5, unit="cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = effect_guild_specific$pval),size=6)

  save(effect_guild_specific,file="effect_guild_specific.RData")# need to be corrected





# use the more conventional way for model selection 
# select the data for each fungal guild

data=list()
for (i in 1:9)
{
  d=model_data_SAR%>%filter(guild==guild_select[i]) # this data set included both richness and climate
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}


sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
    {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }

  else{
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2  + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
      mod_sel=step(mod)
      kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
}
### if we don't do model selection and just focused on all the significant variables

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
   mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2  + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}


### we just focus on the soil and climate variables, load in the initial data

load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")

na_ind <- which(is.na(model_data_SAR$cec))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["cec"]], xy), is.na(r_present_northam[["cec"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$cec[i] <- values["cec"]
}

na_ind <- which(is.na(model_data_SAR$sand))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["sand"]], xy), is.na(r_present_northam[["sand"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$sand[i] <- values["sand"]
}


model_data_SAR$logc=2.71828^model_data_SAR$logc

# do not include plotID
# use this data for modeling

model_data_SAR_sub=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)

model_data_SAR_sub=model_data_SAR_sub[complete.cases(model_data_SAR_sub),]

data=list()
for (i in 1:9)
{
  d=model_data_SAR_sub%>%filter(guild==guild_select[i]) 
  
  d[,c(6:20)]=apply(d[,c(6:20)],2,range01)
  
  data[[i]]=d
}

# construct the model

sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+cec+ sand +bio1+ bio2  + bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])#bio4 removed
      mod_sel=step(mod)
      kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}



## when richness was included

model_data_SAR_rich=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18,richness)

model_data_SAR_rich=model_data_SAR_rich[complete.cases(model_data_SAR_rich),]

data=list()
for (i in 1:9)
{
  d=model_data_SAR_rich%>%filter(guild==guild_select[i]) 
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}

# construct the model

sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +soilMoisture+cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + richness+(1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+cec+ sand +bio1+ bio2  + bio12 + bio15 + bio18  + richness+(1 |siteIDD), data = data[[i]])#bio4 removed
      mod_sel=step(mod)
      kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}

                                      