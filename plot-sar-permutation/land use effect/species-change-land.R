# get the land use data for each grid cell

nc_data <-nc_open ("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2015.nc")
longitude <- ncvar_get(nc_data , "longitude")
latitude<- ncvar_get(nc_data , "latitude")
data_var <- ncvar_get(nc_data, "PFT2")
r <- raster(matrix(data_var, ncol = length(longitude), nrow = length(latitude)))
extent(r) <- extent(range(longitude), range(latitude))
crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent ( -170, -55,18,72), "SpatialPolygons")
r_cropped=crop(r,b)
# change the resolution
coarser_raster <- aggregate(r_cropped, fact = 3, fun = mean) # convert it to a coarse resolution

##get the coordinates for each cell
coords_present <- xyFromCell(coarser_raster, cell = 1:ncell(coarser_raster)) 
cell_values <- raster::extract(coarser_raster, coords_present) %>% as.matrix()



equal_area_crs <- "EPSG:5070"

raster_equal<- terra::project(coarser_raster, equal_area_crs)

biomes1=biomes

st_crs(biomes1)=crs(r_cropped)

polygons_sf <- st_transform(biomes, crs(r_cropped))

values <- extract(r_cropped, polygons_sf)


## transform the coordinates

points <- st_as_sf(coords_present, coords = c("lon", "lat"), crs = st_crs(biomes))

biomes.t <- terra::vect(biomes)
points.t <- terra::vect(points)

values <- terra::extract(biomes.t, points.t)


# get the biomes type for each grid

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

biomes <- st_read("wwf_terr_ecos.shp")

biomes <- group_by(biomes, BIOME) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup() %>%
  mutate(BIOME = as.character(BIOME))

biome_labels <- data.frame(
  BIOME = as.character(c(seq(1, 14), 98, 99)),
  label = c(
    "Tropical & Subtropical Moist Broadleaf Forests",
    "Tropical & Subtropical Dry Broadleaf Forests",
    "Tropical & Subtropical Coniferous Forests",
    "Temperate Broadleaf & Mixed Forests",
    "Temperate Conifer Forests",
    "Boreal Forests/Taiga",
    "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    "Temperate Grasslands, Savannas & Shrublands",
    "Flooded Grasslands & Savannas",
    "Montane Grasslands & Shrublands",
    "Tundra",
    "Mediterranean Forests, Woodlands & Scrub",
    "Deserts & Xeric Shrublands",
    "Mangroves",
    "Undefined",
    "Undefined2"
  ),
  stringsAsFactors = FALSE
)

biomes$LABEL <- biome_labels$label[match(biomes$BIOME, biome_labels$BIOME)]
cbind(biomes$BIOME, biomes$LABEL)

biomes <- st_crop(biomes, c(xmin=-170,xmax=-55,ymin=17,ymax=72))



# just create an empty raster to save the values.
r <- rast(ext(biomes),resolution = res(coarser_raster),   # the resolution of your targeted raster
           crs = "EPSG:4326")
#r <- rasterize(biomes, r, field = "BIOME")  # 'field' can be a column name or a constant value

r <- rasterize(biomes, r, field = "LABEL")  # 'field' can be a column name or a constant value

df <- as.data.frame(r, xy = TRUE)

writeRaster(r,"/Users/luowenqi/soil-sar/plot-sar-permutation/biomes.tif",overwrite=T)

# get the c and z values for different biomes
# add the coordinates for each plot for the determination of the plot-level biomes

sample_data(rare_all_assign)%>%data.frame()%>%dplyr::select(plotIDM,lon,lat)%>%group_by(plotIDM)%>%
  summarise(lon=mean(lon),lat=mean(lat))%>%rename(plotID=plotIDM)->plot_coordinates

full_parameter_data%>%left_join(plot_coordinates,by="plotID")->full_parameter_data

plot_biomes=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct()

extract_biomes=terra::extract(r,plot_biomes[,c("lon","lat")])

plot_biomes%>%bind_cols(extract_biomes%>%dplyr::select(LABEL))%>%dplyr::select(plotID,LABEL)->plot_biomes

# the NAs should be filled 

plot_biomes%>%mutate(LABEL = replace_na(LABEL,"Temperate Broadleaf & Mixed Forests"), LABEL)->plot_biomes

full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%mutate(cvalue=2.71828^logc)%>%
  left_join(plot_biomes,by="plotID")%>%
  group_by(LABEL,guild)%>%summarise(mean_cvalue=mean(cvalue,na.rm=TRUE),mean_zvalue=mean(zvalue,na.rm=TRUE))->parameter_CZ

# just to check which plots were assigned as mangroves
full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%left_join(plot_biomes,by="plotID")%>%filter(LABEL=="Mangroves")%>%
  dplyr::select(plotID)%>%distinct()->plot_mangrove
# LAJA_001
# LAJA_002
# LAJA_003
# LAJA_004
# LAJA_005
# LAJA_015
# LAJA_042
# LAJA_044
# LAJA_046
# LAJA_051

# just to check which plots have not been assigned to specific biomes
full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%left_join(plot_biomes,by="plotID")%>%filter(is.na(LABEL))%>%dplyr::select(plotID)%>%
  distinct()->plot_na

# none plots

#to check the land cover type of these plots

#plot_diversity_env_land %>%filter(plotID%in%plot_na$plotID)%>%dplyr::select(lon,lat)->tem_coord

# these points were assigned as the Temperate Broadleaf & Mixed Forests

#plot_diversity_env_land%>%filter(plotID%in%plot_mangrove$plotID)


## add a column to specify the type based on the land use data
## but this will lead to the loss of some biomes

plot_diversity_env_land%>%dplyr::select(plotID,type)%>%filter(!is.na(type))->k

#some NA plots will be converted into natural communities
# to convert
#model_data_SAR%>%dplyr::select(logc,zvalue,plotID,guild)%>%left_join(k,by="plotID")%>%
#mutate(type=if_else(type%in% c("cultivatedCrops","pastureHay"),"modified","Natural"))%>%
#left_join(plot_biomes,by="plotID")%>%mutate(cvalue=2.71828^logc)%>%
#group_by(LABEL,guild,type)%>%
#summarise(mean_cvalue=mean(cvalue,na.rm=TRUE),mean_zvalue=mean(zvalue,na.rm=TRUE))->parameter_CZ_with_landuse


# to get the biome and guild specific c and z values

full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%
  left_join(plot_biomes,by="plotID")%>%mutate(cvalue=2.71828^logc)%>%
  group_by(LABEL,guild)%>%summarise(mean_cvalue=mean(cvalue,na.rm=TRUE),mean_zvalue=mean(zvalue,na.rm=TRUE))->parameter_CZ_no_landuse
  



# the richness among seven different guilds

mean_richness_guild=bind_rows(land_rich_AM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_EM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_soilsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_littersap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_woodsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_plapat_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_para_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_epiphy_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_all_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)))%>%
  mutate(guild=rep(c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all"),each=7))%>%data.frame()

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

# sensitivity isdefined as the richness ratio between the natural and human-modified landscapes

sensitivity=numeric()
for (i in 1:9)
{
  df=mean_richness_guild%>%filter(guild==guild_type[i])
  sensitivity[i]=df[1,2]/df[2:7,2] %>%mean()
}


sensitivity%>%data.frame()%>%bind_cols(guild_type)%>%data.frame()%>%rename_all(~paste0(c("sensitivity","guild")))%>%
  left_join(parameter_CZ_no_landuse,by="guild")%>%mutate(affinity=sensitivity^(1/mean_zvalue))->affinity_no_landuse



## to see how the plots are distributed

ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  #geom_point(data=model_data_SAR%>%dplyr::select(lon,lat)%>%distinct(),aes(x=lon,y=lat))+
  guides(position="bottom")+
  geom_point(data=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct(),aes(x=lon,y=lat),size=1)+
  theme(legend.position = "right",
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(0.3, "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("")
  


ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  #geom_point(data=plot_diversity_env_land%>%dplyr::select(lon,lat)%>%distinct(),aes(x=lon,y=lat))+
  guides(position="bottom")+
  geom_text(data=plot_diversity_env_land%>%dplyr::select(lon,lat,plotID)%>%distinct(),aes(x=lon,y=lat,label=plotID),size=3)

## read the file and get the proportional data for each grid cell
setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

raster1 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2015.nc")# use the 2015 data as the baseline

ext(raster1) <- c(-90, 90, -180, 180)
crs(raster1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent (-72, -18, -170, -55), "SpatialPolygons")

coarser_raster <- aggregate(raster1, fact = 3, fun = mean) # convert it to a coarse resolution
cropped_raster<- crop(coarser_raster, b)

cropped_raster=flip(t(cropped_raster),direction = "horizontal")# transpose the raster to a right position

new_raster =list()
for (i in 1:33){
temp<- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
                   xmin=-169.95, xmax=-55.05, 
                   ymin=18, ymax=72, 
                   crs=crs(cropped_raster))# change the ranges for the raster
values(temp) <- values(cropped_raster[[i]]) 
new_raster[[i]]=temp#each returns a transformed raster
}

# make equal area projection for each raster and extract the values
# need to convert the matrix into a data.frame for the coordinates

coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")))->coords_present

land_use_data=matrix(ncol=33,nrow=dim(coords_present)[1])
for (i in 1:33)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing

  r_projected <- terra::project(new_raster[[i]], "EPSG:5070")
# transformed the coordinates
  
#points <- vect(coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")), crs=crs(r_projected )))
points_sf <- st_as_sf(coords_present, coords = c("lon", "lat"), crs = 4326) # CRS 4326 is WGS84
coordinates_equal_area <- st_transform(points_sf, "EPSG:5070")
extracted_values <- terra::extract(r_projected, coordinates_equal_area)%>%as.matrix()
land_use_data[,i]=extracted_values[,2]
}

# to see if the sum of the propositional data is approximately 100% for each grid
# the range is 99.85-100.05
apply(land_use_data,1,sum)%>%data.frame()%>%rename_all(~paste0("cover"))->total_area


# get the natural and human-dominated land use type cover for each grid

PFT_2015=coords_present%>%bind_cols(land_use_data%>%data.frame()%>%rename_all(~paste0(names(raster1))))

PFT_2015=PFT_2015%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2015[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"

PFT_2015%>%bind_cols(No_crop)%>%mutate(area=11637.87^2)->PFT_2015

# the area is fixed for each grid 11637.87 m^2 at the resolution of 10 min


#get the richness within each grid based on the SAR
#S=cA^z
#get the biomes for each gird

grid_level_biomes=terra::extract(r,PFT_2015[,c("lon","lat")])

PFT_2015%>%bind_cols(biomes=grid_level_biomes%>%dplyr::select(LABEL))->PFT_2015

# based on the SAR parameters to determine the richness for each grid

# if we do not consider the habitat affinity and focused on the total richness

# the best way is not to differentiate the natural and modified habitats



guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")
# just based on the total area to determine the richness
richness_2015_no_landuse=matrix(ncol=9,nrow=dim(PFT_2015)[1])
for (i in 1:9)
{
  richness_temp=PFT_2015%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(affinity_no_landuse%>%filter(guild==guild_type[i]),by="LABEL")%>%
    mutate(richness=mean_cvalue*(area)^mean_zvalue)%>%as.matrix()
  richness_2015_no_landuse[,i]=richness_temp[,"richness"]
}



# consider the land use composition to estimated the richness
# note, the z should power on the total area rather than individual land area
guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

richness_2015=matrix(ncol=9,nrow=dim(PFT_2015)[1])
for (i in 1:9)
{
  richness_temp=PFT_2015%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(affinity_no_landuse%>%filter(guild==guild_type[i]),by="LABEL")%>%
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2015[,i]=richness_temp[,"richness"]
}


# for future richness in 2100 for the scenario of rcp245

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

raster2 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2100.nc")# use the 2015 data as the baseline

ext(raster2) <- c(-90, 90, -180, 180)
crs(raster2) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent (-72, -18, -170, -55), "SpatialPolygons")
coarser_raster <- aggregate(raster2, fact = 3, fun = mean) # convert it to a coarse resolution
cropped_raster<- crop(coarser_raster, b)
cropped_raster=flip(t(cropped_raster),direction = "horizontal")# transpose the raster to a right position
new_raster =list()
for (i in 1:33){
  temp<- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
              xmin=-169.95, xmax=-55.05, 
              ymin=18, ymax=72, 
              crs=crs(cropped_raster))# change the ranges for the raster
  values(temp) <- values(cropped_raster[[i]]) 
  new_raster[[i]]=temp#each returns a transformed raster
}

# make equal-area projection for each raster and extract the land-use cover data

land_use_data_2100=matrix(ncol=33,nrow=dim(coords_present)[1])
for (i in 1:33)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  r_projected <- terra::project(new_raster[[i]], "EPSG:5070")
  # transformed the coordinates
  
  #points <- vect(coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")), crs=crs(r_projected )))
  points_sf <- st_as_sf(coords_present, coords = c("lon", "lat"), crs = 4326) # CRS 4326 is WGS84
  coordinates_equal_area <- st_transform(points_sf, "EPSG:5070")
  extracted_values <- terra::extract(r_projected, coordinates_equal_area)%>%as.matrix()
  land_use_data_2100[,i]=extracted_values[,2]
}

# check if the total land use cover is about 100%
#apply(land_use_data,1,sum)%>%data.frame()%>%rename_all(~paste0("cover"))->oop



PFT_2100=coords_present%>%bind_cols(land_use_data_2100%>%data.frame()%>%rename_all(~paste0(names(raster2))))


PFT_2100=PFT_2100%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2100[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"

PFT_2100%>%bind_cols(No_crop)%>%mutate(area=11637.87^2)->PFT_2100

# the area is fixed as 11637.87^2 for each grid at the resolution of 10 min


#get the richness within each grid
#s=cA^z
# assign the biome type for each gird

grid_level_biomes=terra::extract(r,PFT_2100[,c("lon","lat")])

PFT_2100%>%bind_cols(biomes=grid_level_biomes%>%dplyr::select(LABEL))->PFT_2100



# get the future richness for each guild and grid cell

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

richness_2100=matrix(ncol=9,nrow=dim(PFT_2100)[1])
for (i in 1:9)
  {
  richness_temp=PFT_2100%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(affinity_no_landuse%>%filter(guild==guild_type[i]),by="LABEL")%>%
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2100[,i]=richness_temp[,10]
}

# difference in the richness among the two time points

richness_2100=richness_2100%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

richness_2015=richness_2015%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

richness_2015_no_landuse=richness_2015_no_landuse%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

colnames(richness_2100)=guild_type

colnames(richness_2015)=guild_type

colnames(richness_2015_no_landuse)=guild_type


species_change_land_rcp245=(richness_2100-richness_2015)/richness_2015

species_change_land_rcp245_no_landuse=(richness_2100-richness_2015_no_landuse )/richness_2015_no_landuse


# differences in the richness for each guild

present_future_richness0=list()
for (i in 1:9)
{
  present_future_richness0[[i]]=bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))
}

# bind all the data.frames

present_future_richness=present_future_richness0[[1]]
for (i in 2:9)
{
  present_future_richness=rbind(present_future_richness,present_future_richness0[[i]])
}

# to see the richness among time points and the change rates

present_future_richness%>%bind_cols(species_change_land_rcp245%>%melt())->species_change_land_rcp245

## there are 162638 cells with NAs for both time periods
d=list()
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&is.na(x2))%>%dim()
  
}

# there are 0 cases where there are species presently but NA in the future

d=list()
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(!is.na(x1)&is.na(x2))%>%dim()
  
}

d=list()# there were 0 cases where there are NA species presently but some species in future
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&!is.na(x2))%>%dim()
  
}

d=list()# there were cases where some species present while 0 species in the future (AM-158,epiphy-667 grids))
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1!=0&x2==0)%>%dim()
  
}


d=list()# for two guilds, we 0 cases where  0 species but some in the future (AM-158,ephiphy-667 grids)
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2!=0)%>%dim()
  
}

d=list()#  we have 0 cases where 0 species present for both time points
for (i in 1: 9)
{
  d[[i]]= bind_cols(richness_2015[,i],richness_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2==0)%>%dim()
  
}

# no need make some changes for the species change ratio



#present_future_richness%>%bind_cols(df_rcp585%>%melt())->species_change_temp

species_change_land_rcp245%>%melt()->species_change_land_rcp245

# for all the guilds were combined
species_change_land_rcp245%>%filter(variable=="all")%>%
  bind_cols(coords_present) ->change_richness_rcp245



species_change_land_rcp245_no_landuse%>%melt()->species_change_land_rcp245_no_landuse
species_change_land_rcp245_no_landuse%>%filter(variable=="all")%>%bind_cols(coords_present) ->change_richness_rcp245_no_landuse



#species_change_land_rcp245%>%mutate(value=if_else(x1==0&x2==0,0,value))->species_change_land_rcp245


species_change_land_rcp245%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
  summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->tem_df_rcp245_land
#the net effect of land use change on the richness

species_change_land_rcp245%>%group_by(variable)%>%
  summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%
  data.frame()->overall_change

tem_df_rcp245_land%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")->tem_df_rcp245_land




# test the significance of each group

T_test_result=matrix(ncol=3,nrow = 9)
for (i in 1:9)
{
  species_change_land_rcp245%>%filter(variable==guild_type[i])->guild_mean
  meam_compare=t.test(guild_mean$value,mu=0)
  T_test_result[i,2:3]=meam_compare$conf.int[1:2]
  T_test_result[i,1]=meam_compare$estimate%>%as.numeric()
  #T_test_result[i,4]=meam_compare$p.value
  
}

T_test_result%>%data.frame()%>%rename_all(~paste0(c("mean","low","up")))%>%mutate(variable=guild_type)->T_test_result


tem_df_rcp245_land%>%left_join(T_test_result,by="variable")->tem_df_rcp245_land


# to add a polygon of the north america


r_present <- raster::getData("worldclim", var = "bio", res = 10)

r_present <- r_present[[c(1, 4, 12)]]

# Run necessary transformations on wordclim-provided temperature data
r_present$bio1 <- r_present$bio1 / 10
r_present$bio4 <- r_present$bio4 / 1000 # no need for further transformation

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland

north_america <- ne_countries(continent = "North America",type="map_units")

# north_america<- subset(north_america, admin != "Canada" & admin != "Greenland")

st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170, -55, 18, 72), "SpatialPolygons")
north_america_cropped <- st_crop(north_america, b)
plot(north_america_cropped)

# Crop to region
r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)


p1=ggplot(change_richness_rcp245) +
  geom_point(data = change_richness_rcp245, pch=21,aes(x = lon, y = lat, color = value), size = 0.275) +
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid="yellow",high = "purple", na.value = "white")+ 
  xlab("Predicted species loss") +
  ylab("")+
  theme(legend.position = "bottom",
        legend.margin = margin(t = -15, r = -5, b = 5, l = 0),
        legend.text = element_text(size=8,angle=90),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.3, -0.5, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
  xlab("")+
  ylab("Land-use impact")+
  #ggtitle("RCP4.5 & SSP2")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("gray80", 0.2),linetype = "solid")+
  scale_size(range = c(0.5, 2))

  
  
summary_data <- change_richness_rcp245 %>% group_by(lat) %>%
   summarise(  mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE))


# the richness map based on sar

present_richness=richness_2015%>%bind_cols(coords_present)%>%dplyr::select(lon,lat,all)

rcp245_richness=richness_2100%>%bind_cols(coords_present)%>%dplyr::select(lon,lat,all)


present_richness$all_normalized <- rescale(present_richness$all, to = c(0, 1), na.rm = TRUE)

rcp245_richness$all_normalized <- rescale(rcp245_richness$all, to = c(0, 1), na.rm = TRUE)


summary_data_richness_sar_present <- present_richness%>% group_by(lat) %>%
  summarise(  mean_value = mean(all_normalized,na.rm=TRUE),sd_value = sd(all_normalized,na.rm=TRUE))

summary_data_richness_sar_rcp245 <- rcp245_richness%>% group_by(lat) %>%
  summarise(  mean_value = mean(all_normalized,na.rm=TRUE),sd_value = sd(all_normalized,na.rm=TRUE))


richness_change_curve=bind_rows(summary_data_richness_sar_present,summary_data_richness_sar_rcp245)%>%
  mutate(type=rep(c("present","rcp245"),each=360))

ggplot()+
  geom_line(data = richness_change_curve, aes(x = lat, y = mean_value,color=type), size = 1) +  # Mean trend line
  geom_ribbon(data =richness_change_curve, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value+ sd_value,fill=type), alpha = 0.2)+
  coord_flip()+
  scale_color_manual("",breaks = c("present","rcp245"), labels = c("present","rcp245"),  values=c("blue","red"))+
  scale_fill_manual("",breaks = c("present","rcp245"), labels = c("present","rcp245"),  values=c("blue","red"))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("Latitude")+
  ylab("Relative richness")+
  ggtitle("")+
  coord_flip()


ggplot()+
  geom_point(data=present_richness,aes(x=lon,y=lat,color=all_normalized),size=0.275)+
  scale_color_gradient2("Relative\n richness", na.value = "white" )+
  theme(legend.position = "right",
        legend.margin = margin(t = -15, r = -5, b = 5, l = 0),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0.2),linetype = "solid")




ggplot()+geom_point(data=richness_2100%>%bind_cols(coords_present),aes(x=lon,y=lat,color=all),size=0.275)+
  scale_color_gradient2(low = "white", mid = "gray", high = "purple", na.value = "white" )





  
p2=ggplot()+
  geom_line(data = summary_data, aes(x = lat, y = mean_value), color = "blue", size = 0.5) +  # Mean trend line
  geom_ribbon(data = summary_data, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value +sd_value), fill = "blue", alpha = 0.2)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("Latitude")+
  ylab("")+
  ggtitle("RCP4.5 & SSP2")+
  coord_flip()+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylim(-0.08,0.08)




tem_df_rcp245_land$variable=factor(tem_df_rcp245_land$variable,levels=guild_type)

p3=ggplot(data=tem_df_rcp245_land,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5,color="black")+
  geom_errorbar(data=tem_df_rcp245_land, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
  scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#8fd1e1","#fedc5e"))+
theme(legend.position = c(0.8,0.87),
      legend.text = element_text(size=8),
      legend.title  = element_text(size=10),
      text = element_text(size = 18),
      plot.title = element_text(size = 15, hjust = 0.5), 
      axis.text.y = element_text(hjust = 0), 
      axis.text.x = element_text(hjust = 1), 
      axis.title.y = element_text(size = 18), 
      axis.title.x = element_text(size = 18), 
      axis.ticks.x = element_blank(), 
      legend.key.size = unit(0.3, "cm"),
      panel.background = element_rect(fill = "NA"),
      panel.border = element_rect(color = "black", size = 1, fill = NA))+
  geom_vline(xintercept =0,color="gray",linetype="dashed")+
  ylab("")+
  xlab("")+
  scale_y_discrete(breaks=guild_type,position="right",labels=c("AM","EM","Soil sapro.","Litter sapro.","Wood sapro.","Plant patho.","Parasite","Epiphyte","All"))+
  geom_segment(data=tem_df_rcp245_land,size=0.35,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable))+
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1",alpha=0.5)+
  geom_hline(yintercept = 8.5,color="red",size=1,alpha=0.3,linetype="dotted")+
  ggtitle("RCP4.5 & SSP2")+
  xlim(-0.5,0.5)+
  geom_text(data=tem_df_rcp245_land,size=6,color="black",
            aes(x=rep(c(0.120640385390, -0.1320765014266, -0.1420753516064,  0.1320636994127, -0.11321197594057,  0.127649490386, -0.12421426806,
                     -0.132754493784, -0.10314108882),each=2),y=variable),label="***")


p1=ggplotGrob(p1)

p2=ggplotGrob(p2)
p3=ggplotGrob(p3)


p2$heights=p3$heights

p7=plot_grid(p1,p2,p3,ncol=3,rel_heights = c(1,0.6,0.6),rel_widths  = c(1,0.6,0.8))



###

# for future richness in the year of 2100

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

raster2 <- rast("GCAM_Demeter_LU_ssp5_rcp85_hadgem_2100.nc")# use the 2015 data as the baseline

ext(raster2) <- c(-90, 90, -180, 180)
crs(raster2) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent (-72, -18, -170, -55), "SpatialPolygons")

coarser_raster <- aggregate(raster2, fact = 3, fun = mean) # convert it to a coarse resolution


cropped_raster<- crop(coarser_raster, b)

cropped_raster=flip(t(cropped_raster),direction = "horizontal")# transpose the raster to a right position


new_raster =list()
for (i in 1:33){
  temp<- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
              xmin=-169.95, xmax=-55.05, 
              ymin=18, ymax=72, 
              crs=crs(cropped_raster))# change the ranges for the raster
  values(temp) <- values(cropped_raster[[i]]) 
  new_raster[[i]]=temp#each returns a transformed raster
}

# make equal area projection for each raster and extract the values

land_use_data_2100=matrix(ncol=33,nrow=dim(coords_present)[1])
for (i in 1:33)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  r_projected <- terra::project(new_raster[[i]], "EPSG:5070")
  # transformed the coordinates
  
  #points <- vect(coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")), crs=crs(r_projected )))
  
  points_sf <- st_as_sf(coords_present, coords = c("lon", "lat"), crs = 4326) # CRS 4326 is WGS84
  
  coordinates_equal_area <- st_transform(points_sf, "EPSG:5070")
  
  
  extracted_values <- terra::extract(r_projected, coordinates_equal_area)%>%as.matrix()
  
  land_use_data_2100[,i]=extracted_values[,2]
}

#apply(land_use_data,1,sum)%>%data.frame()%>%rename_all(~paste0("cover"))->oop




PFT_2100_rcp585=coords_present%>%bind_cols(land_use_data_2100%>%data.frame()%>%rename_all(~paste0(names(raster2))))


PFT_2100_rcp585=PFT_2100_rcp585%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2100_rcp585[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"

PFT_2100_rcp585%>%bind_cols(No_crop)%>%mutate(area=11637.87^2)->PFT_2100_rcp585

# the area is fixed 11637.87^2 with a resolution of 10 min


#get the richness within each grid
#s=cA^z
# get the biomes for each gird

grid_level_biomes=terra::extract(r,PFT_2100_rcp585[,c("lon","lat")])

PFT_2100_rcp585%>%bind_cols(biomes=grid_level_biomes%>%dplyr::select(LABEL))->PFT_2100_rcp585



# get the richness for each guild and each grid

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

richness_2100_rcp585=matrix(ncol=9,nrow=dim(PFT_2100_rcp585)[1])
for (i in 1:9)
{
  richness_temp=PFT_2100_rcp585%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(affinity_no_landuse%>%filter(guild==guild_type[i]),by="LABEL")%>%
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2100_rcp585[,i]=richness_temp[,10]
}

# difference in the richness among the two time points

richness_2100_rcp585=richness_2100_rcp585%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

richness_2015=richness_2015%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

colnames(richness_2100_rcp585)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

colnames(richness_2015)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")


species_change_land_rcp585=(richness_2100_rcp585-richness_2015)/richness_2015
species_change_land_rcp585%>%melt()->species_change_land_rcp585
species_change_land_rcp585%>%filter(variable=="all")%>%bind_cols(coords_present) ->change_richness_rcp585


#species_change_land_rcp245%>%mutate(value=if_else(x1==0&x2==0,0,value))->species_change_land_rcp245



species_change_land_rcp585%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
  summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->tem_df_rcp585_land

species_change_land_rcp585%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
  overall_change

tem_df_rcp585_land%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")->tem_df_rcp585_land

# add the test results
T_test_result_rcp585=matrix(ncol=4,nrow = 9)
for (i in 1:9)
{
  species_change_land_rcp585%>%filter(variable==guild_type[i])->guild_mean
  meam_compare=t.test(guild_mean$value,mu=0)
  T_test_result_rcp585[i,2:3]=meam_compare$conf.int[1:2]
  T_test_result_rcp585[i,1]=meam_compare$estimate%>%as.numeric()
  T_test_result_rcp585[i,4]=meam_compare$p.value
}

T_test_result_rcp585%>%data.frame()%>%rename_all(~paste0(c("mean","low","up","pva")))%>%mutate(variable=guild_type)->T_test_result_rcp585

tem_df_rcp585_land%>%left_join(T_test_result_rcp585,by="variable")->tem_df_rcp585_land



p4=ggplot(change_richness_rcp585) +
  geom_point(data = change_richness_rcp585, pch = 15, aes(x = lon, y = lat, color = 100*value), size = 0.275) +
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid = "yellow", high = "purple", na.value = "white")+ 
  xlab("Predicted species loss") +
  ylab("")+
  theme(legend.position ="bottom",
        legend.margin = margin(t = -15, r = -5, b = 5, l = 0),
        legend.text = element_text(size=8,angle=90),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        plot.margin = unit(c(0.3, -0.5, 0.5, 0.5), "cm"),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
  xlab("")+
  ylab("Land-use impact")+
  #ggtitle("RCP8.5 & SSP5")+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("gray80", 0.2),linetype = "solid")

  
summary_data_rcp585 <- change_richness_rcp585 %>% group_by(lat) %>%
  summarise(
    mean_value = mean(value,na.rm=TRUE),
    sd_value = sd(value,na.rm=TRUE))

###


p5=ggplot()+
  geom_line(data = summary_data_rcp585, aes(x = lat, y = mean_value), color = "blue", size = 0.5) +  # Mean trend line
  geom_ribbon(data = summary_data_rcp585, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value + sd_value), fill = "blue", alpha = 0.2)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("Latitude")+
  ylab("")+
  ggtitle("RCP8.5 & SSP5")+
  coord_flip()+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylim(-0.1,0.1)





tem_df_rcp585_land$variable=factor(tem_df_rcp585_land$variable,levels=guild_type)

p6=ggplot(data=tem_df_rcp585_land,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5,color="black")+
  #geom_errorbar(data=tem_df_rcp585_land, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
  scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#8fd1e1","#fedc5e"))+
  theme(legend.position = c(0.8,0.87),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(0.3, "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  geom_vline(xintercept =0,color="gray",linetype="dashed")+
  ylab("")+

  geom_segment(data=tem_df_rcp585_land,size=0.35,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable)
  )+
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1",alpha=0.5)+
  xlab("")+
  scale_y_discrete(breaks=guild_type,position="right",labels=c("AM","EM","Soil sapro.","Litter sapro.","Wood sapro.","Plant patho.","Parasite","Epiphyte","All"))+
  ggtitle("RCP8.5 & SSP5")+
  geom_hline(yintercept = 8.5,color="red",size=1,alpha=0.3,linetype="dotted")+
  xlim(-0.15,0.15)+
  geom_text(data=tem_df_rcp585_land,size=6,color="black",
            aes(x=rep(c(0.032182824253641, -0.0315022376859, -0.03021121617,  0.03022191932, -0.03004074219,  0.03030850445, -0.03010492210,
                         -0.03021995706, -0.03007214357),each=2),y=variable),label="***")

  
ggplotGrob(p4)

p5=ggplotGrob(p5)
p6=ggplotGrob(p6)
p5$heights=p6$heights



p8=plot_grid(p4,p5,p6,ncol=3,rel_heights = c(1,0.6,0.6),rel_widths  = c(1,0.6,0.8))

plot_grid(p7,p8,ncol=1)

plot_grid(p4,p5,p6,ncol=3,rel_heights = c(1,1,1),rel_widths  = c(1,0.6,0.8))
  
  
#d=p4+p5+p6+plot_layout(widths = c(2.2,0.7, 1),heights = c(2.5,1, 1))

  d&theme(plot.margin = unit(c(0.5, 0.2, -2, -1), "cm"))

  d1&theme(plot.margin = unit(c(0.5, 2, -2, 0), "cm"))
 
  
# changes in the crop and non-crop land
  op=data.frame(df=PFT_2100$No_crop-PFT_2015$No_crop,coords_present)
  op_crop=data.frame(df=PFT_2100$crop-PFT_2015$crop,coords_present)
  
head(op)

ggplot()+
  geom_point(data = op, pch = 15, aes(x = lon, y = lat, color = df), size = 0.275) +
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid = "yellow", high = "purple", midpoint = 0, na.value = "white")+ 
  xlab("Predicted species loss") 

ggplot()+
  geom_point(data = op_crop, pch = 15, aes(x = lon, y = lat, color = df), size = 0.275) +
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid = "yellow", high = "purple",midpoint = 0, na.value = "white")+ 
  xlab("changes in modified land") 



