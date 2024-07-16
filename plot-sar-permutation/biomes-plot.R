
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

ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  geom_point(data=tem_coord,aes(x=lon,y=lat))


ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  geom_point(data=(plot_diversity_env_land%>%filter(type=="mixedForest" )),aes(x=lon,y=lat))

