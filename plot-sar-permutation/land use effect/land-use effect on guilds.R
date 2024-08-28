# for individual guilds
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

load("~/soil-sar/plot-sar-permutation/rare_all_assign.RData")
full_parameter_data <- readRDS("~/full_parameter_data.rds")

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

load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_all_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_AM_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_EM_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_epiphy_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_littersap_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_para_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_plapat_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_plapt_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_soilsap_updated.RData")
load("~/soil-sar/plot-sar-permutation/land use effect/land_rich_woodsap_updated.RData")


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

# sensitivity was defined as the richness ratio between the natural and human-modified landscapes

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
  #geom_point(data=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct(),aes(x=lon,y=lat),size=1)+
  theme(legend.position = "right",
        legend.margin = margin(t = -15, r = -5, b = 0, l = 0),
        legend.text = element_text(size=10),
        legend.title  = element_text(size=10),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.3,0, 0.2, 0), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
  xlab("")+
  ylab("")+
  geom_point(data=d,aes(x=lon,y=lat),color="black")




ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  #geom_point(data=plot_diversity_env_land%>%dplyr::select(lon,lat)%>%distinct(),aes(x=lon,y=lat))+
  guides(position="bottom")+
  geom_text(data=plot_diversity_env_land%>%dplyr::select(lon,lat,plotID)%>%distinct(),aes(x=lon,y=lat,label=plotID),size=3)

# to create a map show the habitat affinity

affinity_no_landuse%>%filter(guild=="soilsap")%>%dplyr::select(LABEL, mean_cvalue, mean_zvalue , affinity)%>%melt->habitat_affinity_map


ggplot(habitat_affinity_map%>%filter(variable!="mean_cvalue"), aes(y =variable , x = LABEL, fill = value)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("",low = "blue", mid = "#FFFFCC", high = "red",midpoint = 0.65)+
  theme(axis.text.x = element_text(angle=90,size=12),
        axis.text.y = element_text(size = 12))+
  scale_y_discrete(breaks=c("affinity","mean_zvalue"),labels=c("Habitat affinity",expression("Mean "*italic(Z))))+
  theme(legend.position = "right",
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1,angle=90), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(0.3, "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))


# create a graph to show the difference in the habitat affinity among different guilds 


ggplot(habitat_affinity_map%>%filter(variable=="mean_cvalue"), aes(y =variable , x = LABEL, fill = value)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("",low = "blue", mid = "#FFFFCC", high = "red",midpoint = 0.7)+
  theme(axis.text.x = element_text(angle=90,size=12),
        axis.text.y = element_text(size = 12))+
  scale_y_discrete(breaks=c("affinity","mean_zvalue"),labels=c("Habitat affinity",expression("Mean "*italic(Z))))



sensitivity%>%data.frame()%>%bind_cols(guild_type)%>%data.frame()%>%rename_all(~paste0(c("sensitivity","guild")))->sensitivity_guild

sensitivity_guild$guild=factor(sensitivity_guild$guild,levels=c("plapat","AM","littersap","woodsap", "all","para", "EM","epiphy","soilsap"))

ggplot(data=sensitivity_guild,aes(x=guild,y=sensitivity))+
  geom_bar(stat = "identity",width=0.5)+
  scale_x_discrete(breaks=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"),
                   labels=c("All","AM","EM","Epiphyte","Litter sprotroph","Parasite","Plant pathogen","Soil saprotroph","Wood saprotroph"))+
  theme(legend.position = c(0.8,0.87),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1,angle=90), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        legend.key.size = unit(0.3, "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  geom_vline(xintercept =0,color="gray",linetype="dashed")+
  geom_hline(yintercept = 1,color="red",linetype="dashed",size=1.1)+
  ylab("Response ratio")+
  xlab("")





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

load("~/soil-sar/plot-sar-permutation/land use effect/coords_present_new.RData")

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


# get the coverage of natural and human-dominated land use type for each grid

PFT_2015=coords_present%>%bind_cols(land_use_data%>%data.frame()%>%rename_all(~paste0(names(raster1))))

PFT_2015=PFT_2015%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2015[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"

PFT_2015%>%bind_cols(No_crop)%>%mutate(area=11637.87^2)->PFT_2015

# the area is fixed for each grid 11637.87 m^2 at the resolution of 10 min of a degree


#get the richness within each grid based on the SAR
#S=cA^z
#get the biomes for each gird

grid_level_biomes=terra::extract(r,PFT_2015[,c("lon","lat")])

PFT_2015%>%bind_cols(biomes=grid_level_biomes%>%dplyr::select(LABEL))->PFT_2015

# based on the SAR parameters to determine the total fungal richness for each grid

# if we do not consider the habitat affinity and focused on the total richness

# the best way is not to differentiate the natural and modified habitats


# do not consider the land use composition to estimated the richness
guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")
just based on the total area to determine the richness， this was not included in in the analysis
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

richness_2015=matrix(ncol=9,nrow=dim(PFT_2015)[1])
for (i in 1:9)
{
  richness_temp=PFT_2015%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(affinity_no_landuse%>%filter(guild==guild_type[i]),by="LABEL")%>%
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2015[,i]=richness_temp[,"richness"]
}


# for future richness in 2100 in the scenario of rcp245

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

# check if the total land use cover approximates 100%
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

# differences in the richness among the two time points in the 

richness_2100=richness_2100%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

richness_2015=richness_2015%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

richness_2015_no_landuse=richness_2015_no_landuse%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))

colnames(richness_2100)=guild_type

colnames(richness_2015)=guild_type

colnames(richness_2015_no_landuse)=guild_type


species_change_land_rcp245=(richness_2100-richness_2015)/richness_2015

#species_change_land_rcp245_no_landuse=(richness_2100-richness_2015_no_landuse )/richness_2015_no_landuse


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

# for individual guilds, rename the guilds names so that they are in the same of that in the climate change effect

guild_type=c("soilsap","plapat","EM","littersap","woodsap","AM","epiphy","para","all")

data_subset_land=list()
for(i in 1:9)
{
  species_change_land_rcp245%>%filter(variable==guild_type[i])%>%
    bind_cols(coords_present) ->data_subset_land[[i]]
}


pp=list()
for (i in 1:9)
  {
pp[[i]]=ggplot(data_subset[[i]]) +
  geom_point(data = data_subset[[i]], pch=21,aes(x = lon, y = lat, color = value), size = 0.275) +
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
  ylab("Land-use impact")
}

plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],ncol=3)


# for the impact of climate change on individual guilds




setwd("/Volumes/seas-zhukai/proj-soil-fungi/sdm-richness-data")

load("present_richness_all.RData")
load("present_richness_guild.RData")
load("future_richness_rcp245_guild.RData")
load("future_richness_rcp245_all.RData")
load("future_richness_rcp585_guild.RData")
load("future_richness_rcp585_all.RData")


guild_model=c("soil_saprotroph","plant_pathogen","ectomycorrhizal","litter_saprotroph","wood_saprotroph","arbuscular_mycorrhizal","epiphyte","para")

current_richness_climate=present_richness_all%>%dplyr::select(-x,-y)%>%bind_cols(present_richness_guild)%>%data.frame()%>%
  rename_all(~paste0(c("all",guild_model)))


richness_climate_rcp245=future_richness_rcp245_all%>%bind_cols(future_richness_rcp245_guild)%>%data.frame()%>%dplyr::select(-x,-y)%>%
  rename_all(~paste0(c("all",guild_model)))

richness_climate_rcp585=future_richness_rcp585_all%>%bind_cols(future_richness_rcp585_guild)%>%data.frame()%>%dplyr::select(-x,-y)%>%
  rename_all(~paste0(c("all",guild_model)))
## mapping current richness 

present_richness_all%>%rename_all(~paste0(c("lon","lat","richness")))->present_richness_all



### current richness change trend
richness_climate_rcp245%>%bind_cols(coords_present)%>%
  dplyr::select(lon,lat,all)->richness_climate_rcp245_all

richness_climate_rcp585%>%bind_cols(coords_present)%>%
  dplyr::select(lon,lat,all)->richness_climate_rcp585_all



## change in richness in the scenario of rcp245

climate_induced_species_change_rcp245= (richness_climate_rcp245-current_richness_climate)/current_richness_climate

climate_induced_species_change_rcp245%>%melt()->df_rcp245_climate


present_future_richness0=list()
for (i in 1:9)
{
  present_future_richness0[[i]]=bind_cols(current_richness_climate[,i],richness_climate_rcp245[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))
}

# bind all the data.frames

present_future_richness=present_future_richness0[[1]]
for (i in 2:9)
{
  present_future_richness=rbind(present_future_richness,present_future_richness0[[i]])
}

present_future_richness%>%bind_cols(df_rcp245_climate)->species_change_climate_rcp245

# to see different scenarios in terms of changes in richness among the two time points

#x1=NA&x2=NA(no need to make changes)
#X1!=NA&x2=NA(3376), leads to NAs
#x1=NA&x2!=NA(none)
#x1==0&x2!=0(92),species change rate was assigned to 1
#X1!=0&x2=0(336)#no need to change the ratio
#x1==0&x2==0(242)# species change rate was assigned to 0
#(x1==0&is.na(x2))(none)
#(is.na(x1)&x2==0)(none)

species_change_climate_rcp245%>%mutate(value=if_else(x1==0&x2!=0,1,value))%>%mutate(value=if_else(x1==0&x2==0,0,value))->
  species_change_climate_rcp245

guild_model=c("soil_saprotroph","plant_pathogen","ectomycorrhizal","litter_saprotroph","wood_saprotroph","arbuscular_mycorrhizal","epiphyte","para","all")


data_subset_climate=list()
for(i in 1:9)
{
  species_change_climate_rcp245%>%filter(variable==guild_model[i])%>%dplyr::select(variable, value)%>%
    bind_cols(coords_present) ->data_subset_climate[[i]]
}


pp=list()
for (i in 1:9)
{
  pp[[i]]=ggplot(data_subset[[i]]) +
    geom_point(data = data_subset[[i]], pch=21,aes(x = lon, y = lat, color = value), size = 0.275) +
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
    ylab("Climate impact")
}

plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],ncol=3)

# to map the relative impact of climate change and land use conversion on the fungal diversity

data=list()
for(i in 1:9){

bind_cols(data_subset[[i]]%>%rename_all(~paste0(c("guild","climate_effect","lon","lat"))),
          data_subset_land[[i]]%>%dplyr::select(variable, value )%>%rename_all(~paste0(c("variable","land_effect"))))->data[[i]]
}



my_function_guild_loss=function(data){
data%>%filter(climate_effect<0)->data_loss_climate
data%>%filter(climate_effect>0)->data_gain_climate
summary(data_loss_climate$climate_effect,na.rm=TRUE)->loss_percentile_climate
summary(data_gain_climate$climate_effect,na.rm=TRUE)->gain_percentile_climate
# to assign a category for each cell based on the magnitude
data%>%mutate(type_C=case_when(climate_effect<=loss_percentile_climate[2]%>%as.numeric()~"C-loss-high",
                               climate_effect>loss_percentile_climate[5]%>%as.numeric()&climate_effect<0~"C-loss-low",
                               climate_effect>=gain_percentile_climate[5]%>%as.numeric()~"C-gain-high",
                               climate_effect<gain_percentile_climate[2]%>%as.numeric()&climate_effect>0~"C-gain-low"))->df2
# for the impact of land use change 
data%>%filter(land_effect<0)->data_loss_land
data%>%filter(land_effect>0)->data_gain_land
summary(data_loss_land$land_effect,na.rm=TRUE)->loss_percentile
summary(data_gain_land$land_effect,na.rm=TRUE)->gain_percentile
data%>%mutate(type_L=case_when(land_effect<=loss_percentile[2]%>%as.numeric()~"L-loss-high",
                               land_effect>loss_percentile[5]%>%as.numeric()&land_effect<0~"L-loss-low",
                               land_effect>=gain_percentile[5]%>%as.numeric()~"L-gain-high",
                               land_effect<gain_percentile[2]%>%as.numeric()&land_effect>0~"L-gain-low"))->df1
# combing both effects of land use conversion and climate change
bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3
df3%>%mutate(joint_effect=paste(df3$type_L,"_", df3$type_C))%>%
  mutate(group=case_when(joint_effect=="L-loss-high _ C-loss-high"~"both-high-loss",
                         joint_effect=="L-loss-low _ C-loss-low"~"both-low-loss",
                         joint_effect%in%c("L-loss-high _ C-loss-low",
                                           "L-loss-high _ NA","L-loss-high _ C-gain-high",
                                           "L-loss-high _ C-gain-low")~"L-high-loss",
                         joint_effect%in%c("L-loss-low _ C-loss-high", 
                                           "NA _ C-loss-high",
                                           "L-gain-high _ C-loss-high",
                                           "L-gain-low _ C-loss-high" )~"C-high-loss",
                         TRUE~"other"))->df4
df4=my_function(df4)
return(df4)
}




df5=my_function_guild_loss(data[[1]])

# the function to mask the great lakes
my_function=function(data)
{
  data$group<- as.numeric(as.factor(data$group))
  r <- rasterFromXYZ(data[, c("lon", "lat", "group")])
  # to mask the object based on the NA map
  # to match the extent of the raster
  #r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)
  #r_present_northam <- clipOutPoly(r_present_northam, greatlakes) # the map based on climates so no need to add variables
  #r_present_northam=raster(r_present_northam)
  r <- crop(r, extent(r_present_northam)) 
  crs(r)="+proj=longlat +datum=WGS84 +no_defs"
  r=projectRaster(r, crs = projection(r_present_northam))
  r <- resample(r, r_present_northam, method = "bilinear") # or "ngb"
  masked_raster_sf <- terra::mask(r, r_present_northam)
  data <- as.data.frame(masked_raster_sf , xy = TRUE)
  data$group=round(data$group,0)
  data$group=as.character(data$group) 
  data%>%mutate(group = replace_na(group, "other"))->data
  return(data)
}





ggplot()+
  geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
  
  scale_color_manual("RCP4.5 & SSP2",breaks=c("1","2","3","4","5","other"),
                     labels=c("Both high(0.80%)","Both low(0.70%)","Climate high(3.8%)","Land-use high(8.3%)","Other",""),
                     values=c("violetred1","black","seagreen","tan1","darkseagreen2","white"))+
  main_theme+
  geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
  ggtitle("Species loss rate")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  ggtitle("Species loss rate")+
  xlab("")+
  ylab("")







