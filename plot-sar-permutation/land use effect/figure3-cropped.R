#estimate the richness
setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")

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


# get the land use data

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

###

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

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")


full_guild_affinity_data=readRDS("full_guild_affinity_data.rds")

full_guild_affinity_data%>%rename(LABEL=biome)->full_guild_affinity_data


richness_2015=matrix(ncol=9,nrow=dim(PFT_2015)[1])
for (i in 1:9)
{
  richness_temp=PFT_2015%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(full_guild_affinity_data%>%dplyr::select(guild, mean_zvalue, mean_cvalue,LABEL,affinity)%>%
                filter(guild==guild_type[i]),by="LABEL")%>%
    
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2015[,i]=richness_temp[,"richness"]
}

richness_2015%>%data.frame()%>%mutate(across(everything(), as.numeric))->richness_2015


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




richness_2100=matrix(ncol=9,nrow=dim(PFT_2100)[1])
for (i in 1:9)
{
  richness_temp=PFT_2100%>%dplyr::select(crop,No_crop,LABEL,area)%>%
    left_join(full_guild_affinity_data%>%dplyr::select(guild, mean_zvalue, mean_cvalue,LABEL,affinity)%>%
                filter(guild==guild_type[i]),by="LABEL")%>%
    
    mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
  richness_2100[,i]=richness_temp[,"richness"]
}

richness_2100%>%data.frame()%>%mutate(across(everything(), as.numeric))->richness_2100

species_change_land_rcp245=(richness_2100-richness_2015)/richness_2015

colnames(species_change_land_rcp245)=guild_type

species_change_land_rcp245=melt(species_change_land_rcp245)

species_change_land_rcp245=bind_cols(species_change_land_rcp245,coords_present)

richness_2015=bind_cols(richness_2015,coords_present)

# for the scenario of rcp585 for land use change effect

bbox <- st_bbox(c(xmin = -2665004, ymin = 2157670 , xmax = 3188701, ymax = 5400000 ), crs = st_crs(canadian_projected))
bbox_sf <- st_as_sfc(bbox)
cropped_province <- st_crop(canadian_projected, bbox_sf)
canada_clipped=cropped_province
bbox <- st_bbox(sf_object)#to see the range of the sf object


# project the results

###############the effect of land use change######################

# the spatial patterns 
species_change_land_rcp245_all=species_change_land_rcp245%>%
  filter(variable=="all")%>%bind_cols(coords_present)%>%rename(x=lon,y=lat)

species_change_land_rcp585_all=species_change_land_rcp585%>%
  filter(variable=="all")%>%bind_cols(coords_present)%>%rename(x=lon,y=lat)

species_change_land_rcp245_all=my_function_project(species_change_land_rcp245_all)
species_change_land_rcp585_all=my_function_project(species_change_land_rcp585_all)



# get the latitude pattern
summary_data=my_function_latitude_patterns(species_change_land_rcp245%>%filter(variable=="all"))

summary_data_land_rcp585=my_function_latitude_patterns(species_change_land_rcp585%>%filter(variable=="all"))

# get the overall effect of different guilds

effect_land_overall_245=my_function_overall_effect(species_change_land_rcp245)
effect_land_overall_585=my_function_overall_effect(species_change_land_rcp585)






###############the effect of climate change######################


species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")
species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")


# to crop the study region based on the land use data
#defin the region to be retained 

land_effect=bind_cols(coords_present,species_change_land_rcp245%>%filter(variable=="all"))

crs(r_climate) <- CRS("+proj=longlat +datum=WGS84")

r_land <- rasterFromXYZ(land_effect[, c("lon", "lat", "value")])

r_tem=as.data.frame(r_land, xy = TRUE, na.rm = TRUE)

# see the overall effects of climate change
#need to maks the data and then bind different guilds

crop_data_rcp585_climate=list()
for (i in 1:9)
{
  climate_effect=bind_cols(coords_present,species_change_climate_rcp585%>%filter(variable==guild_type[i]))
  r_climate <- rasterFromXYZ(climate_effect[, c("lon", "lat", "value")])
  crs(r_climate) <- CRS("+proj=longlat +datum=WGS84")
  masked_raster <- mask(r_climate, r_land)
  raster_df <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE)
  colnames(raster_df) <- c("x", "y", "value")
  crop_data_rcp585_climate[[i]]=raster_df%>%mutate(variable=rep(guild_type[i],dim(raster_df)[1]))
}


crop_data=list()
for (i in 1:9)
{
  climate_effect=bind_cols(coords_present,species_change_climate_rcp245%>%filter(variable==guild_type[i]))
  r_climate <- rasterFromXYZ(climate_effect[, c("lon", "lat", "value")])
  crs(r_climate) <- CRS("+proj=longlat +datum=WGS84")
  masked_raster <- mask(r_climate, r_land)
  raster_df <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE)
  colnames(raster_df) <- c("x", "y", "value")
  crop_data[[i]]=raster_df%>%mutate(variable=rep(guild_type[i],40650))
}

#to bind all the data for climate effect
combined_data=bind_rows(crop_data[[1]],crop_data[[2]],crop_data[[3]],
                        crop_data[[4]],crop_data[[5]],crop_data[[6]],
                        crop_data[[7]],crop_data[[8]],crop_data[[9]])

combined_data_climate_rcp585=bind_rows(crop_data_rcp585_climate[[1]],crop_data_rcp585_climate[[2]],crop_data_rcp585_climate[[3]],
                                       crop_data_rcp585_climate[[4]],crop_data_rcp585_climate[[5]],crop_data_rcp585_climate[[6]],
                                       crop_data_rcp585_climate[[7]],crop_data_rcp585_climate[[8]],crop_data_rcp585_climate[[9]])



#for the latitude trend

my_function_latitude_patterns_crop=function(data)
{
  data%>%rename(lon=x,lat=y)->temp_data
  latitude_patterns <- temp_data%>% group_by(lat) %>%
    summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE))
  return(latitude_patterns )
}

#for the climate change impact on overal fungal diversity

# for the latitude patterns for climate change impacts



climate_induced_change_richness_rcp245=my_function_project(crop_data[[9]])

climate_induced_change_richness_rcp585=my_function_project(crop_data_rcp585_climate[[9]])


summary_data_climate_rcp245=my_function_latitude_patterns_crop(crop_data[[9]])

summary_data_climate_rcp585=my_function_latitude_patterns_crop(crop_data_rcp585_climate[[9]])

# for the overall effect for climate change
climate_245_overall_effect=my_function_overall_effect(combined_data)

climate_245_overall_effect$variable=factor(climate_245_overall_effect$variable,
                                           levels=c("AM","EM","soilsap","littersap" ,
                                                    "woodsap","plapat" ,"para" ,"epiphy","all"))


climate_585_overall_effect=my_function_overall_effect(combined_data_climate_rcp585)

climate_585_overall_effect$variable=factor(climate_585_overall_effect$variable,
                                           levels=c("AM","EM","soilsap","littersap" ,
                                                    "woodsap","plapat" ,"para" ,"epiphy","all"))
###



effect_land_overall_585=my_function_overall_effect(species_change_land_rcp585)
effect_land_overall_585$variable=factor(effect_land_overall_585$variable,levels=guild_type)


# the effect of land use change on the overall fungal diversity

species_change_land_rcp585%>%filter(variable=="all")->species_change_land_rcp585_all
species_change_land_rcp585_all%>%bind_cols(coords_present)%>%rename(x=lon,y=lat)->temp
change_richness_rcp585=my_function_project(temp)

species_change_land_rcp585_all=my_function_project(species_change_land_rcp585_all)

summary_data_land_rcp585=my_function_latitude_patterns(species_change_land_rcp585_all)

# for the climate change impact


#create the plots
# for the overall climate effect

p_land_245=ggplotGrob(ggplot(species_change_land_rcp245_all) +
                        geom_point(data = species_change_land_rcp245_all, pch=15,aes(x = x, y = y, color = last*100), size = 0.175) +
                        scale_color_gradient2(expression("Change %"), high = "#1173EE", mid="white",low = "#ee8c11", na.value = "white")+
                        geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
                        geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = baha_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = jama_projected, fill = NA, size=0.01,color = "black")+
                        coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                        theme(legend.position = c(0.2,0.35),
                              legend.margin = margin(t = -30, r = -4, b = -1, l = 0),
                              legend.text = element_text(size=8,angle=0),
                              legend.title  = element_text(size=10),
                              text = element_text(size = 18),
                              plot.title = element_text(size = 15, hjust = 0.5), 
                              axis.text.y = element_blank(), 
                              axis.text.x = element_blank(), 
                              axis.title.y = element_text(size = 18), 
                              axis.title.x = element_text(size = 18), 
                              axis.ticks.x = element_blank(), 
                              axis.ticks.y = element_blank(),
                              plot.margin = unit(c(0.3, -5, -.5, 0.5), "cm"),
                              panel.background = element_rect(fill = "NA"),
                              panel.border = element_blank())+
                        xlab("")+
                        ylab("Land-use impact")+
                        ggtitle("RCP4.5-SSP2"))

## for the latidude patterns

p_land_latitude_245=ggplotGrob(ggplot()+
                                 geom_line(data = summary_data, aes(x = lat, y = 100*mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
                                 geom_ribbon(data = summary_data, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value +100*sd_value), fill = "#1173EE", alpha = 0.2)+
                                 theme(legend.position = c(0.75,0.28),
                                       legend.text = element_text(size=8),
                                       legend.title  = element_text(size=10),
                                       text = element_text(size = 18),
                                       plot.title = element_text(size = 15, hjust = 0.5), 
                                       axis.text.y = element_text(size=12), 
                                       axis.text.x = element_text(size = 12), 
                                       axis.title.y = element_text(size = 15), 
                                       axis.title.x = element_text(size = 15), 
                                       plot.margin = unit(c(0.3, 0.1, -.5, 0), "cm"),
                                       panel.background = element_rect(fill = "NA"),
                                       panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                 xlab("Latitude")+
                                 ylab("")+
                                 ggtitle("RCP4.5-SSP2")+
                                 coord_flip()+
                                 geom_hline(yintercept = 0,linetype="dashed")+
                                 xlim(15,65)+
                                 ylim(-20,20))
# for the overall effects


p_land_overall_245=ggplotGrob(ggplot(data=effect_land_overall_245,aes(fill=type,y=variable ,x=100*mean_value))+
                                geom_col(width = 0.5)+
                                geom_errorbar(data=effect_land_overall_245, aes(xmin = 100*mean_value- 100*sd_value/sqrt(count), xmax = 100*mean_value +100*sd_value/sqrt(count)),width=0.2)+
                                scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#1173EE","#ee8c11"))+
                                theme(legend.position = c(0.8,0.75),
                                      legend.text = element_text(size=8),
                                      legend.title  = element_text(size=10),
                                      text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size = 12), 
                                      axis.text.x = element_text(size = 12), 
                                      axis.title.y = element_text(size = 18), 
                                      axis.title.x = element_text(size = 18), 
                                      legend.key.size = unit(0.3, "cm"),
                                      plot.margin = unit(c(0.3, 0.5, -0.5, 0.1), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                ylab("")+
                                xlab("")+
                                xlim(-20,20)+
                                scale_y_discrete(breaks=guild_type,position="right",labels=c("AM (+9.4%)","EM (-6.7%)","Soil saprotroph (-4.7%)","Litter saprotroph (-0.3%)","Wood saprotroph (-1.9%)","Plant pathogen (+0.8%)","Parasite (-3.8%)","Epiphyte (-6.3%)","All (-2.4%)"))+
                                geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                ggtitle("RCP4.5-SSP2")+
                                geom_point(data=effect_land_overall_245,aes(y=variable,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
                                geom_errorbar(data=effect_land_overall_245,size=0.35,width=0,color="black",aes(xmin=100*low,xmax=100*up)))

#geom_text(data=effect_land_overall_245,size=5,color="black",
#aes(x=rep(c(0.25020640385390, -0.2011087320765014266, -0.1312420753516064,  0.1287910320636994127, -0.086540611321197594057,  0.1027649490386, -0.128586038102421426806,
#-0.13132754493784, -0.10245710314108882),each=2),y=variable),label="***"))


#


p_climate_245=ggplotGrob(ggplot(climate_induced_change_richness_rcp245) +
                           geom_point(data = climate_induced_change_richness_rcp245, pch=15,aes(x = x, y = y, color = last*100), size = 0.175) +
                           scale_color_gradient2(expression("Change %"), high = "#1173EE", mid="white",low = "#ee8c11", na.value = "white",limits=c(-50,50))+
                           geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = cropped_province, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = baha_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = jama_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
                           coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                           theme(legend.position = c(0.2,0.35),
                                 legend.margin = margin(t = -30, r = -4, b = -1, l = 0),
                                 legend.text = element_text(size=8,angle=0),
                                 legend.title  = element_text(size=10),
                                 text = element_text(size = 18),
                                 plot.title = element_text(size = 15, hjust = 0.5), 
                                 axis.text.y = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.title.y = element_text(size = 18), 
                                 axis.title.x = element_text(size = 18), 
                                 axis.ticks.x = element_blank(), 
                                 axis.ticks.y = element_blank(),
                                 plot.margin = unit(c(0.3, -5, -.5, 0.5), "cm"),
                                 panel.background = element_rect(fill = "NA"),
                                 panel.border = element_blank())+
                           xlab("")+
                           ylab("Climate impact")+
                           ggtitle("RCP4.5-SSP2"))


p_climate_latitude_245=ggplotGrob(ggplot()+
                                    geom_line(data = summary_data_climate_rcp245, aes(x = lat, y = 100*mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
                                    geom_ribbon(data = summary_data_climate_rcp245, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "#1173EE", alpha = 0.2)+
                                    theme(
                                      legend.text = element_text(size=8),
                                      legend.margin = margin(t = -30, r = -5, b = -1, l = 0),
                                      legend.title  = element_text(size=10),
                                      text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size=12), 
                                      axis.text.x = element_text(size = 12), 
                                      axis.title.y = element_text(size = 15), 
                                      axis.title.x = element_text(size = 15), 
                                      
                                      plot.margin = unit(c(0.3, 0.1, -0.5, 0), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                    xlab("Latitude")+
                                    ylab("")+
                                    ggtitle("")+
                                    coord_flip()+
                                    geom_hline(yintercept = 0,linetype="dashed")+
                                    ggtitle("RCP4.5-SSP2")+
                                    xlim(15,65)+
                                    ylim(-40,40))


p_climate_overall_245=ggplotGrob(ggplot(data=climate_245_overall_effect,aes(fill=type,y=variable ,x=100*mean_value))+
                                   geom_col(width = 0.5)+
                                   geom_errorbar(data=climate_245_overall_effect, aes(xmin = 100*mean_value- 100*sd_value/sqrt(count), xmax = 100*mean_value +100*sd_value/sqrt(count)),width=0.2)+
                                   scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#1173EE","#ee8c11"))+
                                   theme(legend.position ="none",
                                         legend.text = element_text(size=8),
                                         legend.title  = element_text(size=10),
                                         legend.margin = margin(t = -30, r = -5, b = -1, l = 0),
                                         text = element_text(size = 15),
                                         plot.title = element_text(size = 15, hjust = 0.5), 
                                         axis.text.y = element_text(size=12), 
                                         axis.text.x = element_text(size=12), 
                                         axis.title.y = element_text(size = 15), 
                                         axis.title.x = element_text(size = 15), 
                                         legend.key.size = unit(0.3, "cm"),
                                         plot.margin = unit(c(0.3, 0.5, -0.5, 0.1), "cm"),
                                         panel.background = element_rect(fill = "NA"),
                                         panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                   geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                   xlab("")+
                                   ylab("")+
                                   geom_point(aes(y=variable,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
                                   geom_errorbar(data=climate_245_overall_effect,size=0.35,width=0,color="black",aes(xmin=100*low,xmax =100*up)
                                   )+
                                   scale_y_discrete(breaks=c("AM","EM","soilsap","littersap" ,
                                                             "woodsap","plapat" ,"para" ,"epiphy","all"),position="right",
                                                    labels=c("AM (+1.9%)","EM (+3.2%)","Soil saprotroph (+3.5%)","Litter saprotroph (+4.4%)","Wood saprotroph (+2.2%)","Plant pathogen (+2.4%)","Parasite (+2.8%)","Epiphyte (+1.8%)","All (+3.4%)"))+
                                   ggtitle("RCP4.5-SSP2")+
                                   geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                   xlim(-20,20))

#geom_segment(data=tem_df_rcp245_climate,size=0.35,color="black",aes(x=overal_mean-overal_sd,xend=overal_mean+overal_sd,y=variable,yend=variable))+

#geom_text(data=climate_245_overall_effect,size=5,color="black",
# aes(x=rep(c(0.230632182824253641, 0.2330315022376859, 0.23403021121617,  0.2531245325303022191932, 0.24201853004074219,  0.253231703030850445, 0.3025423286213010492210,
# 0.283202963021995706, 0.264321503007214357),each=2),y=variable),label="***")



p_land_585=ggplotGrob(ggplot(change_richness_rcp585) +
                        geom_point(data = change_richness_rcp585, pch=15,aes(x = x, y = y, color = last*100), size = 0.175) +
                        scale_color_gradient2(expression("Change %"), high = "#1173EE", mid="white",low = "#ee8c11", na.value = "white")+
                        geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
                        geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = baha_projected, fill = NA, size=0.01,color = "black")+
                        geom_sf(data = jama_projected, fill = NA, size=0.01,color = "black")+
                        coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                        theme(legend.position = c(0.2,0.35),
                              legend.margin = margin(t = -30, r = -5, b = -1, l = 0),
                              legend.text = element_text(size=8,angle=0),
                              legend.title  = element_text(size=10),
                              text = element_text(size = 18),
                              plot.title = element_text(size = 15, hjust = 0.5), 
                              axis.text.y = element_blank(), 
                              axis.text.x = element_blank(), 
                              axis.title.y = element_text(size = 18), 
                              axis.title.x = element_text(size = 18), 
                              axis.ticks.x = element_blank(), 
                              axis.ticks.y = element_blank(),
                              plot.margin = unit(c(0.3, -5, -0.5, 0.5), "cm"),
                              panel.background = element_rect(fill = "NA"),
                              panel.border = element_blank())+
                        xlab("")+
                        ylab("Land-use impact")+
                        ggtitle("RCP8.5-SSP5"))



p_land_latitude_585=ggplotGrob(ggplot()+
                                 geom_line(data = summary_data_land_rcp585, aes(x = lat, y = 100*mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
                                 geom_ribbon(data = summary_data_land_rcp585, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "#1173EE", alpha = 0.2)+
                                 theme(
                                   legend.text = element_text(size=8),
                                   legend.title  = element_text(size=10),
                                   text = element_text(size = 18),
                                   plot.title = element_text(size = 15, hjust = 0.5), 
                                   axis.text.y = element_text(size = 12), 
                                   axis.text.x = element_text(size = 12), 
                                   axis.title.y = element_text(size = 15), 
                                   axis.title.x = element_text(size = 15), 
                                   plot.margin = unit(c(0.3, 0.1, -0.5, 0), "cm"),
                                   panel.background = element_rect(fill = "NA"),
                                   panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                 xlab("Latitude")+
                                 ylab("")+
                                 ggtitle("RCP8.5-SSP5")+
                                 coord_flip()+
                                 geom_hline(yintercept = 0,linetype="dashed")+
                                 ylim(-20,20)+
                                 xlim(15,65))


p_land_overall_585=ggplotGrob(ggplot(data=effect_land_overall_585,aes(fill=type,y=variable ,x=100*mean_value))+
                                geom_col(width = 0.5)+
                                geom_errorbar(data=effect_land_overall_585, aes(xmin = 100*mean_value- 100*sd_value/sqrt(count), xmax = 100*mean_value +100*sd_value/sqrt(count)),width=0.2)+
                                scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#1173EE","#ee8c11"))+
                                theme(legend.position = "none",
                                      legend.text = element_text(size=8),
                                      legend.title  = element_text(size=10),
                                      text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size = 12), 
                                      axis.text.x = element_text(size= 12), 
                                      axis.title.y = element_text(size = 15), 
                                      axis.title.x = element_text(size = 15), 
                                      legend.key.size = unit(0.3, "cm"),
                                      plot.margin = unit(c(0.3, 0.5, -0.5, 0.1), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                ylab("")+
                                xlab("")+
                                geom_point(aes(y=variable,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
                                geom_errorbar(data=effect_land_overall_585,size=0.3,width=0,color="black",aes(xmin=100*low,xmax=100*up)
                                )+
                                scale_y_discrete(breaks=guild_type,position="right",labels=c("AM (+0.2%)","EM (+0.7%)","Soil saprotroph (+0.3%)","Litter saprotroph (+0.5%)","Wood saprotroph (+0.4%)","Plant pathogen (+0.1%)","Parasite (+0.3%)","Epiphyte (-0.05%)","All (+0.3%)"))+
                                ggtitle("RCP8.5-SSP5")+
                                geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                xlim(-20,20))

#geom_text(data=effect_land_overall_585,size=5,color="black",
#aes(x=rep(c(0.05432182824253641, 0.065315022376859, 0.063021121617,  0.053022191932, 0.053004074219,  0.03030850445, 0.053010492210,
#0.053021995706, 0.053007214357),each=2),y=variable),label="***")



#for the climate effect in the scenario of the rcp585

p_climate_585=ggplotGrob(ggplot(climate_induced_change_richness_rcp585) +
                           geom_point(data = climate_induced_change_richness_rcp585, pch=15,aes(x = x, y = y, color = last*100), size = 0.175) +
                           scale_color_gradient2(expression("Change %"), high = "#1173EE", mid="white",low = "#ee8c11", na.value = "white")+
                           geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
                           geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = baha_projected, fill = NA, size=0.01,color = "black")+
                           geom_sf(data = jama_projected, fill = NA, size=0.01,color = "black")+
                           
                           coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                           theme(legend.position = c(0.2,0.35),
                                 legend.margin = margin(t = -30, r = -5, b = -1, l = 0),
                                 legend.text = element_text(size=8,angle=0),
                                 legend.title  = element_text(size=10),
                                 text = element_text(size = 18),
                                 plot.title = element_text(size = 15, hjust = 0.5), 
                                 axis.text.y = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.title.y = element_text(size = 18), 
                                 axis.title.x = element_text(size = 18), 
                                 axis.ticks.x = element_blank(), 
                                 axis.ticks.y = element_blank(),
                                 plot.margin = unit(c(0.3, -5, -0.5, 0.5), "cm"),
                                 panel.background = element_rect(fill = "NA"),
                                 panel.border = element_blank())+
                           xlab("")+
                           ylab("Climate impact")+
                           ggtitle("RCP8.5-SSP5"))

# latitude patterns for the rcp585

p_climate_latitude_585=ggplotGrob(ggplot()+
                                    geom_line(data = summary_data_climate_rcp585, aes(x = lat, y = 100*mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
                                    geom_ribbon(data = summary_data_climate_rcp585, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "#1173EE", alpha = 0.2)+
                                    theme(
                                      legend.text = element_text(size=8),
                                      legend.title  = element_text(size=10),
                                      text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size = 12), 
                                      axis.text.x = element_text(size = 12), 
                                      axis.title.y = element_text(size = 15), 
                                      axis.title.x = element_text(size = 15), 
                                      plot.margin = unit(c(0.3, 0.1, 0.5, 0), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                    xlab("Latitude")+
                                    ylab("")+
                                    ggtitle("RCP8.5-SSP5")+
                                    coord_flip()+
                                    geom_hline(yintercept = 0,linetype="dashed")+
                                    xlim(15,65)+
                                    ylim(-60,60))



p_climate_overall_585=ggplotGrob(ggplot(data=climate_585_overall_effect,aes(fill=type,y=variable ,x=100*mean_value))+
                                   geom_col(width = 0.5)+
                                   geom_errorbar(data=climate_585_overall_effect, aes(xmin = 100*mean_value- 100*sd_value/sqrt(count), xmax = 100*mean_value +100*sd_value/sqrt(count)),width=0.2)+
                                   scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#1173EE","#ee8c11"))+
                                   theme(legend.position ="none",
                                         legend.text = element_text(size=8),
                                         legend.title  = element_text(size=10),
                                         text = element_text(size = 18),
                                         plot.title = element_text(size = 15, hjust = 0.5), 
                                         axis.text.y = element_text(size = 12), 
                                         axis.text.x = element_text(size= 12), 
                                         axis.title.y = element_text(size = 15), 
                                         axis.title.x = element_text(size = 15), 
                                         legend.key.size = unit(0.3, "cm"),
                                         plot.margin = unit(c(0.3, 0.5, 0.5, 0.1), "cm"),
                                         panel.background = element_rect(fill = "NA"),
                                         panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                   geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                   ylab("")+
                                   #geom_errorbar(data=tem_df_rcp585_climate,size=0.25,color="black",aes(xmin=mean_value-overal_sd,xmax=mean_value+overal_sd,width=0.2))+
                                   geom_point(aes(y=variable,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
                                   xlab("")+
                                   scale_y_discrete(breaks=guild_type,position="right",
                                                    labels=c("AM (0.4%)","EM (+1.2%)","Soil saprotroph (+1.4%)","Litter saprotroph (+3.5%)","Wood saprotroph (+2.0%)","Plant pathogen (+2.1%)","Parasite (+1.9%)","Epiphyte (-0.8%)","All (+2.2%)"))+
                                   ggtitle("RCP 8.5-SSP5")+
                                   geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                   xlim(-20,20)+
                                   geom_errorbar(data=climate_585_overall_effect,width=0,size=0.35,color="black",aes(xmin=100*low,xmax=100*up)))


# geom_text(data=climate_585_overall_effect,size=5,color="black",
#aes(x=rep(c(0.3024032182824253641, 0.235315022376859, 0.23021121617,  0.26713022191932, 0.243004074219,  0.25253030850445, 0.28543010492210,
#0.28653021995706, 0.263007214357),each=2),y=variable),label="***")




p_land_latitude_245$heights=p_land_overall_245$heights
p_land_overall_245$widths=p_climate_overall_245$widths

p1=plot_grid(p_land_245,p_land_latitude_245,p_land_overall_245,ncol=3,rel_heights = c(1,0.4,0.4),rel_widths  = c(1,0.4,0.8))


p_climate_latitude_245$heights=p_climate_overall_245$heights


p2=plot_grid(p_climate_245,p_climate_latitude_245,p_climate_overall_245,ncol=3,rel_heights = c(1,0.4,0.4),rel_widths  = c(1,0.4,0.8))


p_land_latitude_585$heights=p_land_overall_585$heights

p_land_overall_585$widths=p_climate_overall_585$widths

p_land_overall_585$widths=p_land_overall_245$widths



p3=plot_grid(p_land_585,p_land_latitude_585,p_land_overall_585,ncol=3,rel_heights = c(1,0.4,0.4),
             rel_widths  = c(1,0.4,0.8))

p_climate_latitude_585$heights=p_climate_overall_585$heights

p_climate_latitude_585$heights=p_climate_overall_245$heights

p_climate_overall_585$heights=p_land_overall_245$heights
p_climate_overall_585$widths=p_land_overall_245$widths

p_land_overall_585$widths=p_land_overall_245$widths

p4=plot_grid(p_climate_585,p_climate_latitude_585,p_climate_overall_585,ncol=3,rel_heights = c(1,0.4,0.4),
             rel_widths  = c(1,0.4,0.8))


plot_grid(p1,p2,p3,p4,ncol=1)






