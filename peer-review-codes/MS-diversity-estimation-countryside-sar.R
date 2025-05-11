
#estimate fungal diversity for each 10-min grid cell
# assign the biome type for each cell

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

#coordinates for each grid cell

load("coords_present_new.RData")
coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")))->coords_present


# based on the SAR parameters to determine the total fungal diversity for each grid
#initial analysis included eight fungal guilds

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")


#land use data for different time points

raster1 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2015.nc")#
raster2 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2100.nc")#
raster3 <- rast("GCAM_Demeter_LU_ssp5_rcp85_hadgem_2100.nc")# 



#diversity within a grid cell was estimated based on the total area weighted by species relative affinity

habitat_affinity_with_land_history_rarefy_consider_nature_history=readRDS("habitat_affinity_with_land_history_rarefy_consider_nature_history.rds")

# the function to estimate grid-cell-level fungal diversity
my_function_raster=function(data)
{
  ext(data) <- c(-90, 90, -180, 180)
  crs(data) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  b <- as(extent (-72, -18, -170, -55), "SpatialPolygons")
  coarser_raster <- aggregate(data, fact = 3, fun = mean) # convert it to a coarse resolution
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
  df_temp=matrix(ncol=33,nrow=dim(coords_present)[1])
  for (i in 1:33)
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    r_projected <- terra::project(new_raster[[i]], "EPSG:5070")
    # transformed the coordinates
    #points <- vect(coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")), crs=crs(r_projected )))
    points_sf <- st_as_sf(coords_present, coords = c("lon", "lat"), crs = 4326) # CRS 4326 is WGS84
    coordinates_equal_area <- st_transform(points_sf, "EPSG:5070")
    extracted_values <- terra::extract(r_projected, coordinates_equal_area)%>%as.matrix()
    df_temp[,i]=extracted_values[,2]
  }
  #apply(land_use_data,1,sum)%>%data.frame()%>%rename_all(~paste0("cover"))->oop
  df_temp1=coords_present%>%bind_cols(df_temp%>%data.frame()%>%rename_all(~paste0(names(data))))
  df_temp1=df_temp1%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,
                             grass=PFT12+PFT13+PFT14,
                             crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)
  temp=df_temp1[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]
  No_crop=apply(temp, 1, sum)%>%data.frame()
  names(No_crop)="No_crop"
  df_temp1%>%bind_cols(No_crop)%>%mutate(area=11637.87^2)->df_temp1
  
  #the area is fixed 11637.87^2 with a resolution of 10-min of degree
  #get the diversity within each grid
  #s=cA^z
  # get the biomes for each gird
  grid_level_biomes=terra::extract(r,df_temp1[,c("lon","lat")])
  df_temp1%>%bind_cols(biomes=grid_level_biomes%>%dplyr::select(LABEL))->df_temp1
  
  # get the richness for each guild and each grid
  guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")
  richness_derived=matrix(ncol=9,nrow=dim(df_temp1)[1])
  for (i in 1:9)
  {
    richness_temp=df_temp1%>%dplyr::select(crop,No_crop,LABEL,area)%>%
      left_join(habitat_affinity_with_land_history_rarefy_consider_nature_history%>%
                  dplyr::select(guild, mean_zvalue, mean_cvalue,LABEL,affinity)%>%
                  filter(guild==guild_type[i]),by="LABEL")%>%
      mutate(richness=mean_cvalue*(affinity*crop/100*area+1*No_crop/100*area)^mean_zvalue)%>%as.matrix()
    richness_derived[,i]=richness_temp[,"richness"]#
  }
  richness_derived=richness_derived%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))
  #colnames(richness_derived)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")
  #richness_derived=melt(richness_derived)
  return(richness_derived)
}



# estimated grid-cell-level fungal diversity for different scenarios

richness_2015=my_function_raster(raster1)

richness_2100_rcp245=my_function_raster(raster2)
richness_2100_rcp585=my_function_raster(raster3)

colnames(richness_2100_rcp585)=guild_type

colnames(richness_2100_rcp245)=guild_type

colnames(richness_2015)=guild_type

# fungal diversity losses and gains among different scenarios

species_change_land_rcp245=(richness_2100_rcp245-richness_2015)/richness_2015
colnames(species_change_land_rcp245)=guild_type
species_change_land_rcp245=melt(species_change_land_rcp245)


species_change_land_rcp585=(richness_2100_rcp585-richness_2015)/richness_2015
species_change_land_rcp585%>%melt()->species_change_land_rcp585



saveRDS(species_change_land_rcp245,file="species_change_land_rcp245.rds")#save the data with the most updated affinity
saveRDS(species_change_land_rcp585,file="species_change_land_rcp585.rds")


