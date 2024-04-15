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

# to see the fit for the neon sites

table(sar_neon_permutation$plotid)

point_number_neon=table(sar_neon_permutation$plotid)%>%data.frame()

point_number_five_neon=subset(point_number,Freq>4)
point_number_four_neon=subset(point_number_neon,Freq>3)

point_number_three_neon=subset(point_number_neon,Freq>2)

## to see the fits

for (i in 1:dim(point_number_five)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_five$Var1)[i])
  
  fit <-sar_power(data_subb[,3:2] )
  plot=plot(fit)
  
  # Print the plot
  print(plot)
  
}

for (i in 1:dim(point_number_five)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_five$Var1)[i])
  
  fit <-plot(log(area)~log(richness),data=data_subb[,3:2] )
 
  
  # Print the plot
  print(fit)
  
}


for (i in 1:dim(point_number_three_neon)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_three$Var1)[i])
  
  fit <-plot(log(area)~log(richness),data=data_subb[,3:2] )
  
  
  # Print the plot
  print(fit)
  
}

## estimate the z value



zvalue=numeric()
pvalue=numeric()
a5=unique(point_number_three_neon$Var1)
for (i in 1:length(a5))
{
  df1=subset(sar_neon_permutation,plotid==a5[i])
  if(dim(df1)[1]<3){
    zvalue[i]=NA
    pvalue[i]=NA
  }
  else{
    ft=lm(log(richness)~log(area),data=df1)%>%summary()
    zvalue[i]=ft$coefficients[2,1]
    pvalue[i]=ft$coefficients[1,4]
  }
  
}

df=cbind(data.frame(a5),zvalue)

df$zvalue=as.numeric(df$zvalue)
names(df)[1]="plotid"
# all the z values determined with more than three dots 

names(point_number)[1]="plotid"

df=merge(df,point_number,by="plotid")
z_neon=df
# when 4 and 5 points were included, there is no relation between core numbers and the estimated the z values

sar_dob_permutation$richness=as.numeric(sar_dob_permutation$richness)


table(sar_dob_permutation$plotid)

point_number_dob=table(sar_dob_permutation$plotid)%>%data.frame()

point_number_five=subset(point_number,Freq>=4)
point_number_four=subset(point_number,Freq==4)
point_number_three=subset(point_number,Freq>2)




for (i in 1:dim(point_number_five)[1])
{
  data_subb <- subset(sar_dob_permutation,plotid==unique(point_number_five$Var1)[i])
  
  fit <-plot(log(richness)~log(area),data=data_subb[,3:2] )
  
  
  # Print the plot
  print(fit)
  
}

# to estimate the z values of the dob sites

zvalue=numeric()
pvalue=numeric()
a6=unique(point_number_three$Var1)
for (i in 1:length(a6))
{
  df1=subset(sar_dob_permutation,plotid==a6[i])
  ft=lm(log(richness)~log(area),data=df1)%>%summary()
  zvalue[i]=ft$coefficients[2,1]
  pvalue[i]=ft$coefficients[1,4]
}
dob_z=cbind(data.frame(a6),zvalue)%>%data.frame()

names(dob_z)[1]="plotid"
names(point_number_dob)[1]="plotid"
df=merge(dob_z,point_number_dob,by="plotid")

# combind the two datasets for modeling

zvalue_all=rbind(z_neon,dob_z)
names(zvalue_all)[1]="plotID"

# get the environmental variables

# plot-level variables

plot_env=model_var[,c(1:16)]
plot_env=unique(plot_env)

zvalue_all_env=merge(zvalue_all,plot_env,by="plotID")

zvalue_all_env=zvalue_all_env[,c(1,3,2,4:17)]

zvalue_all_env[,3:17]=apply(zvalue_all_env[,3:17],2,range01)
# to see the coorelation among variables

ggcorrplot(cor(zvalue_all_env[,c(3:17)]), hc.order = TRUE, type = "lower", lab = TRUE)


mod=lmer(zvalue~organicCPercent  +  ph  +  nitrogen   +   cec +    sand+ bio1+ bio2+bio4 +bio8+ bio12+ bio15+ bio18  +    spei+(1|siteid),data= zvalue_all_env)

mod_best=lmer(zvalue ~ organicCPercent + ph + bio2 + bio15 + spei + (1 | siteid),data=zvalue_all_env)

la.eq <- glmnet( zvalue_all_env[,5:18], zvalue_all_env[,"zvalue"],lambda=0.006607138, family='gaussian', intercept = F, alpha=1)

mod_cv <- cv.glmnet( as.matrix(zvalue_all_env[,5:18]), zvalue_all_env[,"zvalue"], family='gaussian', intercept = F, alpha=1)

plot(mod_cv)

best_lambda=mod_cv$lambda.min

best_lambda=mod_cv$lambda.1se# will have less variables selected


best_model <- glmnet( zvalue_all_env[,c("ph","nitrogen","cec","sand","bio1","bio2","bio4","bio12")], zvalue_all_env[,"zvalue"], family='gaussian', intercept = F, alpha=1,lambda = 0.006607138,final.re=TRUE)

coef(best_model )

## map the current z value based on the best-selected model
# the variables include ph, nitrogen, cec, sand, bio1,bio2,bio4,bio12

# 
newdata=no_na_prediction[,1:2]
# coordinates for the new locations
# get soil variables 
# create the map


r_present <- getData("worldclim",var="bio",res=10)
r_present <- r_present[[c(1,2,4,12)]]
names(r_present) <- c("mat_celsius","MDR", "temp_seasonality","map_mm")

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



setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")

# for the soil pH(need to add and open the image in the work dic)

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10) # need to know why 10 but all variables would be standardized





r_nitrogen <- raster("nitrogen_5-15cm_mean_5000.tif")

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

r_sand_reproj <- projectRaster(r_sand, crs = crs(r_present_northam))
r_sand_northam <- raster::mask(raster::crop(r_sand_reproj, north_america_cropped), north_america_cropped)
r_sand_northam_resample <- raster::resample(r_sand_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_sand_northam_resample / 10) # need to know why 10 but all variables would be standardized

names(r_present_northam) <- c("mat_celsius","MDR", "temp_seasonality","map_mm","nitrogen","cec","ph","sand")




# extract the six soil variables based on the coordinates of the plots
cl <- c("nitrogen", "cec", "ph", "sand")

climate <- list()
for (i in 1:4) {
  climate[[i]] <- raster::extract(r_present_northam[[cl[i]]], newdata[, 1:2])
}
plot_loca_all_soil <- cbind(newdata, climate[[1]], climate[[2]], climate[[3]], climate[[4]])
names(plot_loca_all_soil) <- c("lon", "lat",  "nitrogen", "cec", "ph",  "sand")


# get the climate data of all the plots
r <- getData("worldclim", var = "bio", res = 10)
points <- SpatialPoints(newdata[, 1:2], proj4string = r@crs)
values <- extract(r, points)
df <- cbind.data.frame(coordinates(points), values)

df=df[,c("lon","lat","bio1","bio2","bio4","bio12")]


newdata=cbind(plot_loca_all_soil,df)

newdata=newdata[,-c(7,8)]

newdata=newdata[which(complete.cases(newdata)),]

newdata[,3:10]=apply(newdata[,3:10],2,range01)
newdata=newdata[,]
# make predictions


pre_value=predict(best_model,newx=as.matrix(newdata[,3:10]))

pre_value=cbind(newdata[,1:2],pre_value)

ggplot(pre_value) +
  geom_point(data=pre_value,pch=15,aes(x=lon, y=lat,color=s0), size=0.275)+
  scale_color_gradient(expression(Present *italic( z)* value),low = "blue", high = "yellow")+
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
  ylab("")
  


##

ggplot(richness_predic) +
  geom_point(data=richness_predic,pch=15,aes(x=lon, y=lat,color=richness), size=0.275)+
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
  xlab("")+
  ylab("")+
  xlim(-175,-42)


  
## the predicted species richness
richness

richness=rowSums(no_na_prediction[1:96928,3:8468])
 
       

#to create a map to show their current distributions

richness_predic=cbind(richness,no_na_prediction[,1:2])






# select the seven climate variables that do not show colinearity
##







new = matrix(rnorm(20), nrow=1, ncol=14) 
predict(best_model, s = best_lambda, newx = new)

# with another approach
zvalue_all_env$siteIDD=as.factor(zvalue_all_env$siteIDD)


model <- glmmLasso(zvalue~organicCPercent  +  ph  +  nitrogen  +cec +sand+ bio1+ bio2+bio4 +bio8+ bio12+ bio15+ bio18  +spei, rnd = list(siteid=~1),family=gaussian(link = "identity"),final.re=TRUE,data = zvalue_all_env,lambda = 1.068678)





mod1 <- cv.glmmLasso(fix = zvalue~organicCPercent  +  ph  +  nitrogen  +cec +sand+ bio1+ bio2+bio4 +bio8+ bio12+ bio15+ bio18  +spei, rnd = list(siteid=~1),data = zvalue_all_env, 
                     family = gaussian(link = "identity"), kfold = 10, lambda.final = 'lambda.1se')

## the function


cv.glmmLasso= function(fix, rnd, data, 
           family = stats::gaussian(link = "identity"), 
           kfold = 5, lambdas = NULL, nlambdas = 100, 
           lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001), 
           loss,
           lambda.final=c('lambda.1se', 'lambda.min'),
           ...)
{
  lambda.final <- match.arg(lambda.final)
  
  if(missing(loss))
  {
    # switch allows us to do take the family arg as assign the appropriate 
    # loss function 
    loss <- switch(family$family, 
                   'gaussian' = calc_mse,
                   'binomial' = calc_logloss,
                   'multinomial' = calc_multilogloss,
                   'poisson' = calc_deviance)
  }
  
  x <- useful::build.x(fix, data)
  nobs <- nrow(x)
  nvars <- ncol(x)
  
  # if lambda isn't specified by user, build the lambdas vector, this is 
  # static for all k folds
  if (is.null(lambdas))
  {
    # building the lambda vector
    lambdas <- buildLambdas(fix = fix,
                            rnd = rnd,
                            data = data, 
                            nlambdas = nlambdas, 
                            lambda.min.ratio= lambda.min.ratio)   
  }
  rowDF <- tibble::tibble(
    row = seq(nobs),
    group = sample(rep(seq(kfold), length.out=nobs), replace = FALSE)
  )
  
  # sorting by group 
  rowDF <-  dplyr::arrange(rowDF, .data$group)
  
  
  #instantiating list to hold loss and models for each fold
  lossVecList <- vector(mode = 'list', length = kfold)
  modList_foldk <- vector(mode = 'list', length = kfold)
  
  for(k in 1:kfold)
  {
    testIndices <- dplyr::filter(rowDF, .data$group == k) %>% dplyr::pull(row)
    trainIndices <- rowDF$row[-testIndices]
    modList_foldk[[k]] <- glmmLasso_MultLambdas(fix = fix,
                                                rnd = rnd,
                                                data = data %>% dplyr::slice(trainIndices),
                                                family = family,
                                                lambdas = lambdas,
                                                nlambdas = nlambdas,
                                                lambda.min.ratio = lambda.min.ratio,
                                                ...)
    
    
    
    # hacky way of getting the response variable out of the         
    response_var <- fix[[2]] %>% as.character()
    
    # pulling out actual data
    actualDataVector <- data %>% dplyr::slice(testIndices) %>% 
      dplyr::pull(response_var)
    predictionMatrix <- predict.glmmLasso_MultLambdas(
      object = modList_foldk[[k]],
      newdata = data %>% dplyr::slice(testIndices)
    )
    
    # employing the loss function in form loss(actual,predicted)
    # using loss function, calculating a list of loss values for each vector 
    # of prediction
    # which comes from a glmmLasso model with a specific lambda 
    # storing loss values for each fold
    
    # TODO: think an error is thrown here 
    lossVecList[[k]] <- loss(actual = actualDataVector, predicted = predictionMatrix)
    # each element of this list should be 1 x nlambdas
  }
  
  #building matrix (k by nlambdas) to help calculate cross-validated mean error
  cvLossMatrix <- do.call(what = rbind, args = lossVecList)
  
  cvm = colMeans(cvLossMatrix)
  cvsd <- apply(cvLossMatrix, 2, stats::sd, na.rm = TRUE)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  
  
  # finding the minimum cvm value in order pull out the lambda.min out of 
  # list of lambda
  minIndex <- which.min(cvm)    
  lambda.min <- lambdas[minIndex]
  
  # finding 1se index by doing vectorized comparision such that cvm <= cvup 
  # of minIndex
  my1seIndex <- min(which(cvm <= cvup[minIndex]))
  lambda.1se <- lambdas[my1seIndex]
  
  # chosing lambda.final to use by checking lambda.final option
  # note that first element lambda.final default value will return true for
  # lambda.1se 
  chosenLambda <- if(lambda.final == 'lambda.1se')
  {
    lambda.1se
  }else if(lambda.final == 'lambda.min')
  {
    lambda.min
  }
  
  
  
  glmmLasso.final <- glmmLasso::glmmLasso(fix = fix,
                                          rnd = rnd,
                                          data = data,
                                          family = family,
                                          lambda = chosenLambda)
  
  # add control list argument to this to make converge faster form one that 
  # create lambda.1se
  # TODO: (maybe) For final model fit, supply control list from the model that led to     either lambda.1se or lambda.min
  
  # mimicking cv.glmnet return objects
  return_List <- list(lambdas=lambdas,
                      cvm=cvm,
                      cvsd=cvsd,cvup=cvup,
                      cvlo=cvlo,
                      glmmLasso.final=glmmLasso.final,
                      lambda.min=lambda.min,
                      lambda.1se=lambda.1se)
  
  
  class(return_List) <- 'cv.glmmLasso'
  
  
  return(return_List)
  
}

## to make predictions

# add the ID to each plot
dm=zvalue_all_env$siteIDD

siteid=matrix(ncol=1,nrow=dim(zvalue_all_env)[1])
for (i in 1:dim(zvalue_all_env)[1])
{
  dk=dm[i]
  if(nchar(dk)<4)
    {
    siteid[i,1]=substr(dk,1,2)
  }
  else{
    siteid[i,1]=substr(dk,1,4) 
  }
  
}

zvalue_all_env=cbind(siteid,zvalue_all_env)
