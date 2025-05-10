
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


load("coords_present_new.RData")
coords_present%>%data.frame()%>%rename_all(~paste0(c("lon","lat")))->coords_present


# based on the SAR parameters to determine the total fungal richness for each grid
# if we do not consider the habitat affinity and focused on the total richness
# the best way is not to differentiate the natural and modified habitats

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")



# for the scenario of rcp585 for land use change effect

# the effect of land use change in the scenario of rcp585

raster1 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2015.nc")# use the 2015 data as the baseline
raster2 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2100.nc")# use the 2015 data as the baseline
raster3 <- rast("GCAM_Demeter_LU_ssp5_rcp85_hadgem_2100.nc")# use the 2015 data as the baseline


load("coords_present_new.RData")

#diversity within a grid cell was estimated based on the total area weigthted by species relative affinity

habitat_affinity_with_land_history_rarefy_consider_nature_history=readRDS("habitat_affinity_with_land_history_rarefy_consider_nature_history.rds")

# the function to estimate fungal diversity for each grid-cell
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
  
  #the area is fixed 11637.87^2 with a resolution of 10 min
  #get the richness within each grid
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
      mutate(richness=mean_cvalue*(affinity*crop/100*area+No_crop/100*area)^mean_zvalue)%>%as.matrix()
    richness_derived[,i]=richness_temp[,"richness"]
  }
  richness_derived=richness_derived%>%data.frame%>%mutate(across(everything(), ~ as.numeric(as.character(.))))
  #colnames(richness_derived)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")
  #richness_derived=melt(richness_derived)
  return(richness_derived)
}



# get the richness data for different scenarios

richness_2015=my_function_raster(raster1)

richness_2100_rcp245=my_function_raster(raster2)
richness_2100_rcp585=my_function_raster(raster3)


# difference in the pixel-level richness over time


colnames(richness_2100_rcp585)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

colnames(richness_2100_rcp245)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

colnames(richness_2015)=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

species_change_land_rcp245=(richness_2100_rcp245-richness_2015)/richness_2015
colnames(species_change_land_rcp245)=guild_type
species_change_land_rcp245=melt(species_change_land_rcp245)


species_change_land_rcp585=(richness_2100_rcp585-richness_2015)/richness_2015
species_change_land_rcp585%>%melt()->species_change_land_rcp585



saveRDS(species_change_land_rcp245,file="species_change_land_rcp245.rds")#save the data with the most updated affinity

saveRDS(species_change_land_rcp585,file="species_change_land_rcp585.rds")

species_change_land_rcp585=readRDS("species_change_land_rcp585.rds")

# project the results
# the function to project the data based on a data.frame

my_function_project=function(data)
{
  if ("x"%in%colnames(data))
  {
    points <- vect(data, geom = c("x", "y"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    raster <- rasterize(points, raster_template, field = "value")
  }
  else{
    points <- vect(data, geom = c("lon", "lat"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    
    raster <- rasterize(points, raster_template, field = "group")
  }
  target_crs <- "EPSG:5070"
  raster_equal_area <- project(raster, target_crs,method="near")# there are options for the method used
  raster_df <- as.data.frame(raster_equal_area, xy = TRUE,)
  return(raster_df )
}
