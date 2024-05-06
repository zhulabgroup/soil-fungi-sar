
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

# used the coordinates to estimate the richnes and the z value with the regression kriging

#1. to look at the estimated richness in these sites





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
