

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
    "Lake",
    "Rock and ice"
  ),
  stringsAsFactors = FALSE
)

biomes$LABEL <- biome_labels$label[match(biomes$BIOME, biome_labels$BIOME)]
cbind(biomes$BIOME, biomes$LABEL)


biomes <- st_crop(biomes, c(xmin=-170,xmax=-55,ymin=17,ymax=72))




biomes <- st_set_crs(biomes, 4326) 


neon_dob_prevalent_v4.1=readRDS("neon_dob_prevalent_v4.1.rds")


sample_data(neon_dob_prevalent_v4.1)%>%data.frame()%>%dplyr::select(lon,lat)->study_site

my_multipolygon_equal_area <- sf::st_transform(biomes,  5070)

df_sf <- st_as_sf(study_site, coords = c("lon", "lat"), crs = 4326)
df_equal_area <- sf::st_transform(df_sf, 5070)

df_sf_crop=df_sf <- st_as_sf(df4, coords = c("lon", "lat"), crs = 4326)
df_equal_area_crop <- sf::st_transform(df_sf_crop, 5070)


ggplot() +
  geom_sf(data=my_multipolygon_equal_area, aes(fill=LABEL))+
  geom_sf(data = df_equal_area, color = "black", size = 2)+
  geom_sf(data = df_equal_area_crop , pch=15,color = "yellow", alpha=1,size = 2)+
  theme(legend.position = "none",
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
        plot.margin = unit(c(0.3, 0.3,0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
  geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
  geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = baha_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = jama_projected, fill = NA, size=0.01,color = "black")+
  
  coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))
  
  
## to see how many soil samples are collected

site=c("BLAN", "SCBI", "STER", "SERC", "DSNY", "JERC", "KONA", "ORNL", "LAJA")
a=numeric()
for (i in 1:9){

subset_samples(get(data[m]),plotIDM%in%df4$plotIDM&Site==site[i])%>%nsamples()->a[i]
}
# add the polygons for different states


subset_sf <- my_multipolygon_equal_area %>% 
  filter(LABEL%in%c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Dry Broadleaf Forests"  ))


ggplot() +
  geom_sf(data=subset_sf, aes(fill=LABEL))+
  theme(legend.position = c(0.15,0.5),
        
        legend.key.height = unit(0.8, "cm"),
        legend.margin = margin(t = -30, r = -4, b = -1, l = 0),
        legend.text = element_text(size=25,angle=0),
        
       
        legend.title  = element_text(size=25),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.3, 0.3,0.5, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())+
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
  scale_fill_manual("Biomes",breaks = subset_sf$LABEL,
                    values=colors_2,
                    labels=c("Dry forests", 
                             "Temperate forests" ,        
                             "Conifer forests",
                             "Grasslands" 
                              ))
# subset the data

colors_1 <- c("#008080", "#FF6F61", "#FFB81C", "#333333")

# Color combination 2: Navy Blue, Soft Pink, Mint Green, and Gold
colors_2 <- c("#003366", "#F4C2C2", "#98FF98", "#FFD700")

# Color combination 3: Turquoise, Peach, Lavender, and Cream
colors_3 <- c("#40E0D0", "#FFDAB9", "#E6E6FA", "#FFF5E1")

# Color combination 4: Olive Green, Burgundy, Light Grey, and Ivory
colors_4 <- c("#556B2F", "#800020", "#D3D3D3", "#FFFFF0")

colors_2=c("#B983FF",  "#00C0AF","#00BCD8","#00B0F6")


[1] "#FD61D1" "#B983FF" "#619CFF" "#00C0AF" "#00BCD8" "#F8766D" "#E76BF3"
[8] "#00B0F6" "#C99800" "#FF67A4" "#00BA38" "#E58700" "#6BB100" "#A3A500"
[15] "#00BF7D"
