

#(1) get the biomes type for each grid cell where we model species richness
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
# the resolution is different 
r <- rast(ext(biomes),resolution = 0.15,   # the resolution of your targeted raster
          crs = "EPSG:4326")
#r <- rasterize(biomes, r, field = "BIOME")  # 'field' can be a column name or a constant value

r <- rasterize(biomes, r, field = "LABEL")  # 'field' can be a column name or a constant value



#(2) get the plot-level c and z values for each biome

full_parameter_data <- readRDS("full_parameter_data.rds")
plot_coordinates=readRDS("plot_coordinates.rds")
full_parameter_data%>%left_join(plot_coordinates,by="plotID")->full_parameter_data
plot_biomes=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct()
extract_biomes=terra::extract(r,plot_biomes[,c("lon","lat")])
plot_biomes%>%bind_cols(extract_biomes%>%dplyr::select(LABEL))%>%dplyr::select(plotID,LABEL)->plot_biomes
# the NAs should be filled 
plot_biomes%>%mutate(LABEL = replace_na(LABEL,"Temperate Broadleaf & Mixed Forests"), LABEL)->plot_biomes
full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%
  left_join(plot_biomes,by="plotID")%>%mutate(cvalue=2.71828^logc)%>%
  group_by(LABEL,guild)%>%summarise(mean_cvalue=mean(cvalue,na.rm=TRUE),mean_zvalue=mean(zvalue,na.rm=TRUE))->parameter_CZ_no_landuse

#(3) get the richness response ratio for different guilds

biome_site_level_richness_ratio_consider_nature_history=readRDS("biome_site_level_richness_ratio_consider_nature_history.rds")

# to merge with the response ratio

parameter_CZ_no_landuse%>%mutate(guild_biome=paste(guild,"_",LABEL))->df1

# mismatch of the biomes because the choic of the resolution
biome_site_level_richness_ratio_consider_nature_history$guild_biome=gsub("Moist","Dry",biome_site_level_richness_ratio_consider_nature_history$guild_biome)

biome_site_level_richness_ratio_consider_nature_history%>%left_join(df1,by="guild_biome")%>%
  filter(!is.na( mean_ratio))%>%mutate(affinity=mean_ratio^(1/mean_zvalue))%>%
  rename(guild=guild.x)->habitat_affinity_with_land_history_rarefy_consider_nature_history

saveRDS(habitat_affinity_with_land_history_rarefy_consider_nature_history,file="habitat_affinity_with_land_history_rarefy_consider_nature_history.rds")

