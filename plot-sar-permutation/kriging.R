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






# for the random forest

set.seed(1856)
mtry <- sqrt(ncol(train.x))             # number of variables randomly sampled as candidates at each split.
tunegrid.rf <- expand.grid(.mtry=mtry)
RF<-train(train.x,
          RESPONSE,
          method = "rf",
          trControl=myControl,
          tuneGrid=tunegrid.rf,
          ntree= 100,
          preProc=c('center', 'scale'))
print(RF)

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
