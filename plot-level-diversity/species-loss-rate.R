# changes in the land cover types
library(ncdf4)

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
1# get the changes in land use types between 2020 and 2100 year
# for the present land use types

nc_file <- nc_open("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc")

raster1 <-  rast("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2020.nc")
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
data=data[complete.cases(data),]
data_mean=aggregate(data[,c(1,2,7:38)],by=list(data$grid),FUN=mean,na.rm=TRUE)%>%data.frame()# aggregate the value within a larger grid
data_mean=data_frame(data_mean)



2.# the future distribution of the land cover types

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
raster2 <- rast("GCAM_Demeter_LU_ssp1_rcp26_hadgem_2100.nc")

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

data_future=cbind(data_future,PFT_FUTURE[,3:35])

data_future=data_future[complete.cases(data_future),]

data_mean_future=aggregate(data_future[,c(1,2,7:38)],by=list(data_future$grid),FUN=mean)%>%data.frame()# aggregate the value within a larger grid

save(data_mean,file="data_mean_present.RData")
save(data_mean_future,file="data_mean_future.RData")

names(data_mean)=c("group","lon","lat",names(raster1)[2:33])

names(data_mean_future)=c("group","lon","lat",names(raster1)[2:33])

# compute changes in the land use, the ratio in land use cover between the two time periods
data_mean=data.frame(data_mean)

ratio=matrix(ncol=35,nrow=103567)
for (i in 4:35){
  ratio[,i]=data_mean_future[,i]/data_mean[,i]
  
}%>%data.frame()

names(ratio)=c("group","lon","lat",names(raster1)[2:33])

# can we just look at part of the habitats

data_mean=mutate(data_mean,loca=paste(data_mean$lon,"_",-1*data_mean$lat))

data_mean_future=mutate(data_mean_future,loca=paste(data_mean_future$lon,"_",-1*data_mean_future$lat))

richness_zva=cbind(pred_plot,pred_zvalue[,"zvalue"])%>%mutate(loca=paste(lon,"_",lat))

names(richness_zva)[4]="zvalue"

# to see some change in the land cover types for some grids
#
ratio=cbind(data_mean$loca,ratio[,9])%>%data.frame()# to see this type of the richness

ratio$X2=as.numeric(ratio$X2)
names(ratio)=c("loca","change")

species_change=merge(richness_zva,ratio,by="loca")

# the ratio between the future and present species richness

species_change=species_change%>%mutate(ratio_rich=change^zvalue)%>%mutate(future_rich=ratio_rich*rich)%>%mutate(change_rich=future_rich-rich)%>%mutate(change_rate=change_rich/rich)

# convert the NA values

save(species_change,file="species_change.RData")

df=species_change[,c("rich","change","ratio_rich","future_rich","change_rich","change_rate")]

dff=matrix(ncol=6,nrow=96693)
for (i in 1:dim(df)[1])
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  if(is.nan(df[i,2]))
    {
    dk=df[i,1]
    dff[i,1]=dk
    dff[i,2]=0
    dff[i,3]=0
    dff[i,4]=dk
    dff[i,5]=0
    dff[i,6]=0
    }
   
  else {
       dff[i,]=df[i,]%>%as.numeric()
     }
}


save(pred_plot,file="predict_plot_rich.RData")

save(pred_zvalue,file="predict_zvalue.RData")




# changes in the land cover of interest

p1=ggplot(data_mean) +
  geom_point(data=data_mean,pch=15,aes(x=lon, y=-lat,color=PFT14), size=0.275)+
  scale_color_gradient(expression("Cover %"),low = "gray", high = "purple")+
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
  xlab("PFT14 in 2020")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

###

p2=ggplot(data_mean_future) +
  geom_point(data=data_mean_future,pch=15,aes(x=lon, y=-lat,color=PFT14), size=0.275)+
  scale_color_gradient(expression("Cover %"),low = "gray", high = "purple",limit=c(0,100))+
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
  xlab("PFT14 in 2100")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

# changes in the richness
ggplot(species_change1) +
  geom_point(data=species_change,pch=15,aes(x=lon, y=lat,color=change_rate1), size=0.275)+
  scale_color_gradient(expression("Species change %"),low = "seagreen1", high = "purple",na.value = "yellow")+
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
  xlab("Predicted species changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

# mapping the richness

ggplot(species_change) +
  geom_point(data=species_change,pch=15,aes(x=lon, y=lat,color=change1), size=0.275)+
  scale_color_gradient(expression("Species change %"),low = "blue", high = "purple")+
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
  xlab("Predicted species changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

ggplot(species_change1) +
  geom_point(data=species_change1,pch=15,aes(x=lon, y=lat,color=change1), size=0.275)+
  scale_color_gradient(expression("Species change %"),low = "blue", high = "purple")+
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
  xlab("Predicted species changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

dd=cbind(species_change[,2:3],dff[,2])

ggplot(dd) +
  geom_point(data=dd,pch=15,aes(x=lon, y=lat,color=change), size=0.275)+
  scale_color_gradient(expression(""),low = "gray", high = "purple")+
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
  xlab("Predicted land cover changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations




