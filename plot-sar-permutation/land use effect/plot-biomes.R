# distribution of croplands
# get the updated plot level information

#https://data.neonscience.org/documents/-/document_library_display/kV4WWrbEEM2s/view_file/3411446?_110_INSTANCE_kV4WWrbEEM2s_redirect=https%3A%2F%2Fdata.neonscience.org%2Fdocuments%2F-%2Fdocument_library_display%2FkV4WWrbEEM2s%2Fview%2F2233450%3F_110_INSTANCE_kV4WWrbEEM2s_redirect%3Dhttps%253A%252F%252Fdata.neonscience.org%252Fdocuments%253Fp_p_id%253D110_INSTANCE_kV4WWrbEEM2s%2526p_p_lifecycle%253D0%2526p_p_state%253Dnormal%2526p_p_mode%253Dview%2526p_p_col_id%253Dcolumn-1%2526p_p_col_count%253D1

crop=read.csv("All_NEON_TOS_Plot_Centroids_V11.csv",sep=",")
crop%>%select(plotID,latitude,longitude,nlcdClass)%>%
  filter(nlcdClass%in%c("pastureHay","cultivatedCrops"))->crop_select


# select the plots with soil samples

plot_id_with_samples=sample_data(rare_all_assign)%>%data.frame()%>%
  select(plotIDM)%>%rename(plotID=plotIDM)
# merge with the updated data

load("~/soil-sar/plot-sar-permutation/plot_diversity_env_land.RData")
plot_diversity_env_land%>%select(plotID,type)%>%left_join(crop_select,by="plotID")%>%filter(type%in%c("pastureHay","cultivatedCrops"))->crop_plots



# get the biomes type for each grid

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

# create the plots

ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))+
  #geom_point(data=model_data_SAR%>%dplyr::select(lon,lat)%>%distinct(),aes(x=lon,y=lat))+
  guides(position="bottom")+
  #geom_point(data=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct(),aes(x=lon,y=lat),size=1)+
  theme(legend.position = "bottom",
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
  guides(fill = guide_legend(nrow = 5, byrow = TRUE))+
  xlab("")+
  ylab("")+
  geom_point(data=df4,size=1,alpha=0.8,aes(y=lat,x= lon),pch=15,fill="yellow",color="yellow")

