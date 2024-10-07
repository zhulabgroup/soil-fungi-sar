# figure4 in the main text
# to map the joint effect of climate change and land use change
# load the data stored in Turbo

setwd("/Volumes/seas-zhukai/proj-soil-fungi/sdm-richness-data")
species_change_land_rcp585=readRDS("species_change_land_rcp585.rds")
species_change_land_rcp245=readRDS("species_change_land_rcp245.rds")
species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")
species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")

############
species_change_climate_rcp245 %>%
  mutate(variable = case_when(
    variable == "soil_saprotroph" ~ "soilsap",
    variable == "plant_pathogen" ~ "plapat",
    variable == "ectomycorrhizal" ~ "EM",
    variable == "litter_saprotroph" ~ "littersap",
    variable == "wood_saprotroph" ~ "woodsap",
    variable == "arbuscular_mycorrhizal" ~ "AM",
    variable == "epiphyte" ~ "epiphy",
    TRUE ~ variable  # Keep other values unchanged
  ))->species_change_climate_rcp245

species_change_climate_rcp585 %>%
  mutate(variable = case_when(
    variable == "soil_saprotroph" ~ "soilsap",
    variable == "plant_pathogen" ~ "plapat",
    variable == "ectomycorrhizal" ~ "EM",
    variable == "litter_saprotroph" ~ "littersap",
    variable == "wood_saprotroph" ~ "woodsap",
    variable == "arbuscular_mycorrhizal" ~ "AM",
    variable == "epiphyte" ~ "epiphy",
    TRUE ~ variable  # Keep other values unchanged
  ))->species_change_climate_rcp585

####select the data set for each guild

data_subset_land_rcp245=list()
for(i in 1:9)
{
  species_change_land_rcp245%>%filter(variable==guild_type[i])%>%dplyr::select(variable, value)%>%
    bind_cols(coords_present) ->data_subset_land_rcp245[[i]]
}

data_subset_climate_rcp245=list()
for(i in 1:9)
{
  species_change_climate_rcp245%>%filter(variable==guild_type[i])%>%dplyr::select(variable, value)%>%
    bind_cols(coords_present) ->data_subset_climate_rcp245[[i]]
}

data_subset_land_rcp585=list()
for(i in 1:9)
{
  species_change_land_rcp585%>%filter(variable==guild_type[i])%>%dplyr::select(variable, value)%>%
    bind_cols(coords_present) ->data_subset_land_rcp585[[i]]
}

data_subset_climate_rcp585=list()
for(i in 1:9)
{
  species_change_climate_rcp585%>%filter(variable==guild_type[i])%>%dplyr::select(variable, value)%>%
    bind_cols(coords_present) ->data_subset_climate_rcp585[[i]]
}

# combined the two effects for the two scenarios

data_both_effect_rcp245=list()
for(i in 1:9){
  
  bind_cols(data_subset_climate_rcp245[[i]]%>%rename_all(~paste0(c("guild","climate_effect","lon","lat"))),
            data_subset_land_rcp245[[i]]%>%dplyr::select(variable, value )%>%rename_all(~paste0(c("variable","land_effect"))))->data_both_effect_rcp245[[i]]
}


data_both_effect_rcp585=list()
for(i in 1:9){
  
  bind_cols(data_subset_climate_rcp585[[i]]%>%rename_all(~paste0(c("guild","climate_effect","lon","lat"))),
            data_subset_land_rcp585[[i]]%>%dplyr::select(variable, value )%>%rename_all(~paste0(c("variable","land_effect"))))->data_both_effect_rcp585[[i]]
}



load("~/soil-sar/SDM/r_present_northam.RData")

r_present_northam=r_present_northam$mat_celsius

# the function to mask the great lakes
my_function=function(data)
{
  north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
  lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")
  north_america <- st_transform(north_america, "+proj=longlat +datum=WGS84 +no_defs")
  lakes <- st_transform(lakes, crs(r_present_northam))
  north_america_raster <- rasterize(north_america, r_present_northam, field = 1)
  lakes_raster <- rasterize(lakes, r_present_northam, field = 1)
  r_masked <- mask(r_present_northam, north_america_raster)
  # Mask the raster to exclude lake regions
  r_exclude_lakes <- mask(r_masked, lakes_raster, maskvalue = 1, updatevalue = NA)
  data$group<- as.numeric(as.factor(data$group))
  r <- rasterFromXYZ(data[, c("lon", "lat", "group")])
  r <- crop(r, extent(r_exclude_lakes ))
  crs(r)="+proj=longlat +datum=WGS84 +no_defs"
  r=projectRaster(r, crs = projection(r_exclude_lakes ))
  r <- resample(r, r_exclude_lakes  , method = "ngb") 
  masked_raster_sf <- terra::mask(r, r_exclude_lakes )
  data <- as.data.frame(masked_raster_sf , xy = TRUE)
  data$group=round(data$group,0)
  data$group=as.character(data$group) 
  data%>%mutate(group = replace_na(group, "other"))->data
  return(data)
}

# the function to project the raster
my_function_project=function(data)
{
  points <- vect(data, geom = c("x", "y"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
  raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
  raster <- rasterize(points, raster_template, field = "group")
  target_crs <- "EPSG:5070"
  raster_equal_area <- project(raster, target_crs,method="near")# there are options for the method used
  raster_df <- as.data.frame(raster_equal_area, xy = TRUE,)
  return(raster_df )
}


# function to determine the species loss rate

my_function_guild_loss=function(data){
  data%>%filter(climate_effect<0)->data_loss_climate 
  data%>%filter(climate_effect>0)->data_gain_climate
  summary(data_loss_climate$climate_effect,na.rm=TRUE)->loss_percentile_climate 
  summary(data_gain_climate$climate_effect,na.rm=TRUE)->gain_percentile_climate# to assign a category for each cell based on the magnitude
  data%>%mutate(type_C=case_when(climate_effect<=loss_percentile_climate[2]%>%as.numeric()~"C-loss-high", 
                                 climate_effect>loss_percentile_climate[5]%>%as.numeric()&climate_effect<0~"C-loss-low", 
                                 climate_effect>=gain_percentile_climate[5]%>%as.numeric()~"C-gain-high",
                                 climate_effect<gain_percentile_climate[2]%>%as.numeric()&climate_effect>0~"C-gain-low"))->df2# for the impact of land use change
  data%>%filter(land_effect<0)->data_loss_land 
  data%>%filter(land_effect>0)->data_gain_land
  summary(data_loss_land$land_effect,na.rm=TRUE)->loss_percentile 
  summary(data_gain_land$land_effect,na.rm=TRUE)->gain_percentile
  data%>%mutate(type_L=case_when(land_effect<=loss_percentile[2]%>%as.numeric()~"L-loss-high", land_effect>loss_percentile[5]%>%as.numeric()&land_effect<0~"L-loss-low", 
                                 land_effect>=gain_percentile[5]%>%as.numeric()~"L-gain-high",
                                 land_effect<gain_percentile[2]%>%as.numeric()&land_effect>0~"L-gain-low"))->df1# combing both effects of land use conversion and climate change
  bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3 

  df3%>%mutate(joint_effect=paste(df3$type_L,"_",df3$type_C))%>%
  mutate(group=case_when(joint_effect=="L-loss-high _ C-loss-high"~"both-high-loss", 
                         joint_effect=="L-loss-low _ C-loss-low"~"both-low-loss", 
                         joint_effect%in%c("L-loss-high _ C-loss-low","L-loss-high _ NA","L-loss-high _ C-gain-high","L-loss-high _ C-gain-low")~"L-high-loss", 
                         joint_effect%in%c("L-loss-low _ C-loss-high","NA _ C-loss-high",
                         "L-gain-high _ C-loss-high","L-gain-low _ C-loss-high")~"C-high-loss",
                         TRUE~"other"))->df4
  df4=my_function(df4)
  return(df4)
}
        
        
# function to determine the species gain rate

my_function_guild_gain=function(data){
  data%>%filter(climate_effect<0)->data_loss_climate 
  data%>%filter(climate_effect>0)->data_gain_climate
  summary(data_loss_climate$climate_effect,na.rm=TRUE)->loss_percentile_climate 
  summary(data_gain_climate$climate_effect,na.rm=TRUE)->gain_percentile_climate# to assign a category for each cell based on the magnitude
  data%>%mutate(type_C=case_when(climate_effect<=loss_percentile_climate[2]%>%as.numeric()~"C-loss-high", 
                                 climate_effect>loss_percentile_climate[5]%>%as.numeric()&climate_effect<0~"C-loss-low", 
                                 climate_effect>=gain_percentile_climate[5]%>%as.numeric()~"C-gain-high",
                                 climate_effect<gain_percentile_climate[2]%>%as.numeric()&climate_effect>0~"C-gain-low"))->df2# for the impact of land use change
  data%>%filter(land_effect<0)->data_loss_land 
  data%>%filter(land_effect>0)->data_gain_land
  summary(data_loss_land$land_effect,na.rm=TRUE)->loss_percentile 
  summary(data_gain_land$land_effect,na.rm=TRUE)->gain_percentile
  data%>%mutate(type_L=case_when(land_effect<=loss_percentile[2]%>%as.numeric()~"L-loss-high", land_effect>loss_percentile[5]%>%as.numeric()&land_effect<0~"L-loss-low", 
                                 land_effect>=gain_percentile[5]%>%as.numeric()~"L-gain-high",
                                 land_effect<gain_percentile[2]%>%as.numeric()&land_effect>0~"L-gain-low"))->df1# combing both effects of land use conversion and climate change
  bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3 
  df3%>%mutate(joint_effect=paste(df3$type_L,"_",df3$type_C))%>%
    mutate(group=case_when(joint_effect=="L-gain-high _ C-gain-high"~"both-high-gain", 
                           joint_effect=="L-gain-low _ C-gain-low"~"both-low-gain", 
                           joint_effect%in%c("L-gain-high _ C-loss-high","L-gain-high _ C-gain-low","L-gain-high _ C-loss-low","L-gain-high _ NA")~"L-high-gain", 
                           joint_effect%in%c("L-loss-low _ C-gain-high","L-gain-low _ C-gain-high","NA _ C-gain-high"  ,"L-loss-high _ C-gain-high")~"C-high-gain",
                           TRUE~"other"))->df4
  
  df4=my_function(df4)
  return(df4)
}

# the common theme for the plots

common_theme=theme(legend.position = c(0.18,0.60),
        legend.margin = margin(t = -5, r = -5, b = -5, l = 0),
        legend.text = element_text(size=12),
        legend.title  = element_text(size=12),
        text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0, 0, 0), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())

  
map_data=c("data_both_effect_rcp245","data_both_effect_rcp585")
map_loss=list()
  for (i in 1:2)
    {
data=get(map_data[i])[[9]]    
df5=my_function_guild_loss(data) 
df5=my_function_project(df5)
df5$group=as.factor(df5$group)
data_percent=table(df5$group)[1:5]/sum(table(df5$group)[1:5])->d
d%>%data.frame()%>%
  mutate(group=c("Both high","Both low","Climate high","Land-use high","Other"))->data_percent
if(i<2){
map_loss[[i]]=ggplot()+
  geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
  scale_color_manual("RCP4.5-SSP2",breaks=c("0","1","2","3","4","5"),
                     labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                     values=c("#E92316","#DA5725","#32CD92","#CD326D","gray88","white"))+

  #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
  ggtitle("Species loss rate")+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  xlab("")+
  ylab("")+
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
  common_theme
  
}
else{
  map_loss[[i]]=ggplot()+
    geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
    scale_color_manual("RCP5.8-SSP5",breaks=c("0","1","2","3","4","5"),
                       labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                       values=c("#E92316","yellow","#32CD92","#CD326D","gray88","white"))+
    
    #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
    ggtitle("Species loss rate")+
    guides(color = guide_legend(override.aes = list(size = 2)))+
    xlab("")+
    ylab("")+
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
    common_theme
}
}
  
# for species gain

map_gain=list()
for (i in 1:2)
{
  data=get(map_data[i])[[9]]    
  df5=my_function_guild_gain(data)  
  df5=my_function_project(df5)
  df5$group=as.factor(df5$group)
  data_percent=table(df5$group)[1:5]/sum(table(df5$group)[1:5])->d
  d%>%data.frame()%>%
    mutate(group=c("Both high","Both low","Climate high","Land-use high","Other"))->data_percent
  if(i<2){
    map_gain[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP4.5-SSP2",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","yellow","#32CD92","#CD326D","gray88","white"))+
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species gain rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
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
      common_theme
  }
  else{
    map_gain[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      
      scale_color_manual("RCP5.8-SSP5",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","yellow","#32CD92","#CD326D","gray88","white"))+
      
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species gain rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
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
      common_theme
  }
}

plot_grid(map_loss[[1]],map_gain[[1]],map_loss[[2]],map_gain[[2]],ncol=2,

          labels = c("(a)","(b)","(c)","(d)"),label_y = 1,label_x = 0.2,label_size = 20)


##if we just look at the regions with both species loss

my_function_guild_loss_overall=function(data){
  
  data%>%mutate(type_C=case_when(climate_effect<0~"C-loss", 
                                 climate_effect>0~"C-gain",
                                 is.na(climate_effect)~"C-NA",
                                 climate_effect==0~"C-no"))->df2
  
  
  data%>%mutate(type_L=case_when(land_effect<0~"L-loss", 
                                 land_effect>0~"L-gain",
                                 is.na(land_effect)~"L-NA",
                                  land_effect==0~"L-none")
                                 )->df1# combing both effects of land use conversion and climate change
  bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3 
  
  df3%>%mutate(joint_effect=paste(df3$type_L,"_",df3$type_C))%>%
    mutate(group=case_when(joint_effect=="L-loss _ C-loss"~"both-loss", 
                           #joint_effect=="L-gain _ C-gain"~"both-gain", 
                           joint_effect=="L-NA _ C-NA"~"both-NA", 
                           joint_effect%in%c("L-loss _ NA",  "L-loss _ C-gain", "L-loss _ C-none")~"L-loss", 
                           joint_effect%in%c( "L-none _ C-loss", "NA _ C-loss", "L-gain _ C-loss")~"C-loss",
                           TRUE~"other"))->df4
  
  #df4=my_function(df4)
  return(df4)
}

## for the species gain rate 

my_function_guild_gain_overall=function(data){
  
  data%>%mutate(type_C=case_when(climate_effect<0~"C-loss", 
                                 climate_effect>0~"C-gain",
                                 is.na(climate_effect)~"C-NA",
                                 climate_effect==0~"C-no"))->df2
  
  
  data%>%mutate(type_L=case_when(land_effect<0~"L-loss", 
                                 land_effect>0~"L-gain",
                                 is.na(land_effect)~"L-NA",
                                 land_effect==0~"L-none")
  )->df1# combing both effects of land use conversion and climate change
  bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3 
  
  df3%>%mutate(joint_effect=paste(df3$type_L,"_",df3$type_C))%>%
    mutate(group=case_when(
                           joint_effect=="L-gain _ C-gain"~"both-gain", 
                           joint_effect=="L-NA _ C-NA"~"both-NA", 
                           joint_effect%in%c("L-gain _ C-loss", "L-gain _ C-NA", "L-gain _ C-no" )~"L-gain", 
                           joint_effect%in%c( "L-none _ C-gain","L-NA _ C-gain","L-loss _ C-gain")~"C-gain",
                           TRUE~"other"))->df4
  
  #df4=my_function(df4)
  return(df4)
}
## for species loss rate for the two factors



map_loss=list()
for (i in 1:2)
{
  data=get(map_data[i])[[9]]    
  df5=my_function_guild_loss_overall(data)  
  df5=my_function_project(df5)
  df5$group=factor(df5$group,levels=c("0","2","3","4","1"))
  
  data_percent=table(df5$group)[c(1:4)]/sum(table(df5$group)[c(1:4)])->d
  d%>%data.frame()%>%
    mutate(group=c("Both loss","Climate loss","Land loss","Other"))->data_percent
  if(i<2){
    map_loss[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP4.5-SSP2",breaks=c("0","2","3","4","1"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#CD326D","#E3A72F","#32CD92","gray88","white"))+
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species loss rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
      geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
      geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
      coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
      common_theme
  }
  else{
    map_loss[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      
      scale_color_manual("RCP5.8-SSP5",breaks=c("0","2","3","4","1"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#CD326D","#E3A72F","#32CD92","gray88","white"))+
      
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species loss rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
      geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
      geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
      coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
      common_theme
  }
}


map_gain=list()
for (i in 1:2)
{
  data=get(map_data[i])[[9]]    
  df5=my_function_guild_gain_overall(data)  
  df5=my_function_project(df5)
  df5$group=factor(df5$group,levels=c("0","2","3","4","1"))
  
  data_percent=table(df5$group)[c(1:4)]/sum(table(df5$group)[c(1:4)])->d
  d%>%data.frame()%>%
    mutate(group=c("Both gain","Climate gain","Land gain","Other"))->data_percent
  if(i<2){
    map_gain[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP4.5-SSP2",breaks=c("0","2","3","4","1"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#CD326D", "#E3A72F","#32CD92","gray88","white"))+
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species gain rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
      geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
      geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
      coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
      common_theme
  }
  else{
    map_gain[[i]]=ggplot()+
      geom_point(data=df5,aes(x=x,y=y,color=group),size=0.01)+
      
      scale_color_manual("RCP5.8-SSP5",breaks=c("0","2","3","4","1"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#CD326D", "#E3A72F","#32CD92","gray88","white"))+
      
      #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
      ggtitle("Species gain rate")+
      guides(color = guide_legend(override.aes = list(size = 2)))+
      xlab("")+
      ylab("")+
      geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
      geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
      geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
      coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
      common_theme
  }
}


plot_grid(map_loss[[1]],map_gain[[1]],map_loss[[2]],map_gain[[2]],ncol=2,
          
          labels = c("(a)","(b)","(c)","(d)"),label_y = 0.98,label_x = 0.2,label_size = 20)



                                                                                                                                                                                                                                                                                                 