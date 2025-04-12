library(ggplot2)
library(tigris)
library(rnaturalearth)
library(sf)

#figures in the main text
#load the data required to create the figures
# for the impact of land use change
# when all the guilds were combined
##functions for projection

my_function_project=function(data)
{
  points <- vect(data, geom = c("x", "y"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
  raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
  raster <- rasterize(points, raster_template, field = "value")
  target_crs <- "EPSG:5070"
  raster_equal_area <- project(raster, target_crs,method="near")
  raster_df <- as.data.frame(raster_equal_area, xy = TRUE,)
  return(raster_df )
}


my_function_latitude_patterns=function(data)
{
  data%>%filter(variable=="all")%>%
    bind_cols(coords_present)%>%rename(lon=lon,lat=lat)->temp_data
  latitude_patterns <- temp_data%>% group_by(lat) %>%
    summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE))
  return(latitude_patterns )
}



my_function_overall_effect=function(data)
{
  data%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
    summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->data_temp
  data%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
    overall_change
  data_temp%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")->data_temp
  guild_t_test=unique(data$variable)
  T_test_result=matrix(ncol=4,nrow = 9)
  for (i in 1:9)
  {
    data%>%filter(variable==guild_t_test[i])->guild_mean
    meam_compare=t.test(guild_mean$value,mu=0)
    T_test_result[i,2:3]=meam_compare$conf.int[1:2]
    T_test_result[i,1]=meam_compare$estimate%>%as.numeric()
    T_test_result[i,4]=meam_compare$p.value
  }
  T_test_result%>%data.frame()%>%rename_all(~paste0(c("mean","low","up","pva")))%>%mutate(variable=guild_t_test)->T_test_result
  data_temp%>%left_join(T_test_result,by="variable")->data_temp
  return(data_temp)
}

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")



#load in the data stored on turbo

setwd("/Volumes/seas-zhukai/proj-soil-fungi/sdm-richness-data")
species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")
species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")

species_change_land_rcp245=readRDS("species_change_land_rcp245.rds")

species_change_land_rcp585=readRDS("species_change_land_rcp585.rds")


# when all the guilds were combined, for the impact of climate change in the scenario of rcp245

species_change_climate_rcp245%>%filter(variable=="all")%>%
  bind_cols(coords_present)%>%rename(lon=lon,lat=lat)->climate_induced_change_richness_rcp245


species_change_climate_rcp585%>%filter(variable=="all")%>%
  bind_cols(coords_present)%>%rename(lon=x,lat=y) ->climate_induced_change_richness_rcp585



species_change_land_rcp245%>%filter(variable=="all")%>%
  bind_cols(coords_present)%>%rename(lon=x,lat=y)->change_richness_land_rcp245_all

species_change_land_rcp585%>%filter(variable=="all")%>%
  bind_cols(coords_present)%>%rename(lon=x,lat=y) ->change_richness_land_rcp585_all

# to reorder the display of the guilds on the y axis

tem_df_rcp245_land$variable=factor(tem_df_rcp245_land$variable,levels=guild_type)






#create the maps

common_boundary= ggplot()+
  geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
  geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
  geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")


climate_induced_change_richness_rcp245=my_function_project(climate_induced_change_richness_rcp245)

climate_induced_change_richness_rcp585=my_function_project(climate_induced_change_richness_rcp585)

change_richness_land_rcp245_all=my_function_project(change_richness_land_rcp245_all)
change_richness_land_rcp585_all=my_function_project(change_richness_land_rcp585_all)



# create maps for the distributions



or the latitude patterns
p_land_latitude_245=ggplot()+
  geom_line(data = summary_data, aes(x = lat, y = mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
  geom_ribbon(data = summary_data, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value +sd_value), fill = "#1173EE", alpha = 0.2)+
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
  ylim(-0.08,0.08)
# for the overall effects

effect_land_overall_245=my_function_overall_effect(species_change_land_rcp245)
effect_land_overall_245$variable=factor(effect_land_overall_245$variable,
                                        levels=c("AM","EM","soilsap","littersap" ,
                                                 "woodsap","plapat" ,"para" ,"epiphy","all"))

p_land_overall_245=ggplot(data=effect_land_overall_245,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  geom_errorbar(data=effect_land_overall_245, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
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
  scale_y_discrete(breaks=guild_type,position="right",labels=c("AM (+4.1%)","EM (-4.7%)","Soil saprotroph (-4.5%)","Litter saprotroph (+3.2%)","Wood saprotroph (-0.6%)","Plant pathogen (+4.2%)","Parasite (-1.8%)","Epiphyte (-4.6%)","All (-1.2%)"))+
  geom_point(data=effect_land_overall_245,aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1",alpha=0.25)+
  geom_segment(data=effect_land_overall_245,size=0.35,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable))+
  geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
  ggtitle("RCP4.5-SSP2")+
  xlim(-0.15,0.15)+
  geom_text(data=effect_land_overall_245,size=5,color="black",
            aes(x=rep(c(0.1020640385390, -0.11087320765014266, -0.112420753516064,  0.087910320636994127, -0.03540611321197594057,  0.1027649490386, -0.0586038102421426806,
                        -0.1132754493784, -0.045710314108882),each=2),y=variable),label="***")

p_climate_245=ggplot(climate_induced_change_richness_rcp245) +
  geom_point(data = climate_induced_change_richness_rcp245, pch=15,aes(x = x, y = y, color = last*100), size = 0.175) +
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
  ylab("Climate impact")+
  ggtitle("RCP4.5-SSP2")

# for the latitude patterns

p_climate_latitude_245=ggplot()+
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
  ylim(-40,40)

# for the overal effects

climate_245_overall_effect=my_function_overall_effect(species_change_climate_rcp245)

climate_245_overall_effect$variable=factor(climate_245_overall_effect$variable,
                                      levels=c("AM","EM","soilsap","littersap" ,
                                               "woodsap","plapat" ,"para" ,"epiphy","all"))


p_climate_overall_245=ggplot(data=climate_245_overall_effect,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  geom_errorbar(data=climate_245_overall_effect, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
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
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1",alpha=0.5)+
  geom_segment(data=climate_245_overall_effect,size=0.35,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable)
  )+
  scale_y_discrete(breaks=c("AM","EM","soilsap","littersap" ,
                            "woodsap","plapat" ,"para" ,"epiphy","all"),position="right",
                   labels=c("AM (-7.8%)","EM (+1.4%)","Soil saprotroph (+1.8%)","Litter saprotroph (+3.5%)","Wood saprotroph (-1.0%)","Plant pathogen (-1.0%)","Parasite (-2.1%)","Epiphyte (-8.7%)","All (+3.7%)"))+
  ggtitle("RCP4.5-SSP2")+
  geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
  #geom_segment(data=tem_df_rcp245_climate,size=0.35,color="black",aes(x=overal_mean-overal_sd,xend=overal_mean+overal_sd,y=variable,yend=variable))+
  xlim(-0.6,0.6)+
  geom_text(data=climate_245_overall_effect,size=5,color="black",
            aes(x=rep(c(0.230632182824253641, 0.2330315022376859, -0.23403021121617,  0.231245325303022191932, 0.24201853004074219,  0.253231703030850445, -0.35423286213010492210,
                        -0.3202963021995706, -0.264321503007214357),each=2),y=variable),label="***")




          
p_land_585=ggplot(change_richness_rcp585) +
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
  ggtitle("RCP8.5-SSP5")


summary_data=my_function_latitude_patterns(species_change_land_rcp245)

summary_data_land_rcp585=my_function_latitude_patterns(species_change_land_rcp585)

p_land_latitude_585=ggplot()+
  geom_line(data = summary_data_land_rcp585, aes(x = lat, y = mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
  geom_ribbon(data = summary_data_land_rcp585, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value + sd_value), fill = "#1173EE", alpha = 0.2)+
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
  ylim(-0.08,0.08)


# for the overal effec

effect_land_overall_585=my_function_overall_effect(species_change_land_rcp585)
tem_df_rcp585_land$variable=factor(tem_df_rcp585_land$variable,levels=guild_type)

p_land_overall_585=ggplot(data=effect_land_overall_585,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  geom_errorbar(data=effect_land_overall_585, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
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
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="#1173EE",alpha=0.25)+
  xlab("")+
  geom_segment(data=effect_land_overall_585,size=0.3,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable)
  )+
  scale_y_discrete(breaks=guild_type,position="right",labels=c("AM (+0.2%)","EM (-0.3%)","Soil saprotroph (-0.3%)","Litter saprotroph (+0.3%)","Wood saprotroph (-0.05%)","Plant pathogen (+0.4%)","Parasite (-0.1%)","Epiphyte (-0.3%)","All (-0.1%)"))+
  ggtitle("RCP8.5-SSP5")+
  geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
  xlim(-0.15,0.15)+
  geom_text(data=effect_land_overall_585,size=5,color="black",
            aes(x=rep(c(0.032182824253641, -0.0315022376859, -0.03021121617,  0.03022191932, -0.03004074219,  0.03030850445, -0.03010492210,
                        -0.03021995706, -0.03007214357),each=2),y=variable),label="***")


###






climate_induced_change_richness_rcp585=my_function_project(species_change_climate_rcp585)

p_climate_585=ggplot(climate_induced_change_richness_rcp585) +
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
  ggtitle("RCP8.5-SSP5")
# latitude patterns for the rcp585

summary_data_climate_rcp585=my_function_latitude_patterns(species_change_climate_rcp585)

p_climate_latitude_585=ggplot()+
  geom_line(data = summary_data_climate_rcp585, aes(x = lat, y = mean_value), color = "#1173EE", size = 0.5) +  # Mean trend line
  geom_ribbon(data = summary_data_climate_rcp585, aes(x = lat, ymin = mean_value - sd_value, ymax = mean_value + sd_value), fill = "#1173EE", alpha = 0.2)+
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
  ylim(-0.6,0.6)

effect_climate_overall_585=my_function_overall_effect(species_change_climate_rcp585)

effect_climate_overall_585$variable=factor(effect_climate_overall_585$variable,
                                      levels=c("AM","EM","soilsap","littersap" ,
                                               "woodsap","plapat" ,"para" ,"epiphy","all"))


p_climate_overall_585=ggplot(data=effect_climate_overall_585,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  geom_errorbar(data=effect_climate_overall_585, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
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
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
  xlab("")+
  
  scale_y_discrete(breaks=guild_type,position="right",
                   labels=c("AM (-4.4%)","EM (+1.0%)","Soil saprotroph (+1.6%)","Litter saprotroph (+5.0%)","Wood saprotroph (+3.1%)","Plant pathogen (+2.6%)","Parasite (+1.0%)","Epiphyte (-7.9%)","All (+4.4%)"))+
  ggtitle("RCP 8.5-SSP5")+
  geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
  geom_segment(data=effect_climate_overall_585,size=0.35,color="black",aes(x=overal_mean-low,xend=overal_mean+up,y=variable,yend=variable))+
  geom_text(data=effect_climate_overall_585,size=5,color="black",
            aes(x=rep(c(-0.3024032182824253641, 0.235315022376859, -0.23021121617,  0.26713022191932, 0.243004074219,  0.25253030850445, -0.28543010492210,
                        0.28653021995706, 0.263007214357),each=2),y=variable),label="***")+
  xlim(-0.6,0.6)

     
       

# the overall effect of climate




p_land_245=ggplotGrob(p_land_245)
p_land_latitude_245=ggplotGrob(p_land_latitude_245)
p_land_overall_245=ggplotGrob(p_land_overall_245)
p_land_latitude_245$heights=p_land_overall_245$heights
p_land_overall_245$widths=p_climate_overall_245$widths

p1=plot_grid(p_land_245,p_land_latitude_245,p_land_overall_245,ncol=3,rel_heights = c(1,0.4,0.4),rel_widths  = c(1,0.4,0.8))

p_climate_245=ggplotGrob(p_climate_245)
p_climate_latitude_245=ggplotGrob(p_climate_latitude_245)
p_climate_overall_245=ggplotGrob(p_climate_overall_245)
p_climate_latitude_245$heights=p_climate_overall_245$heights


p2=plot_grid(p_climate_245,p_climate_latitude_245,p_climate_overall_245,ncol=3,rel_heights = c(1,0.4,0.4),rel_widths  = c(1,0.4,0.8))

p_land_585=ggplotGrob(p_land_585)
p_land_latitude_585=ggplotGrob(p_land_latitude_585)
p_land_overall_585=ggplotGrob(p_land_overall_585)
p_land_latitude_585$heights=p_land_overall_585$heights
p_climate_overall_585=ggplotGrob(p_climate_overall_585)
p_land_overall_585$widths=p_climate_overall_585$widths

p_land_overall_585$widths=p_land_overall_245$widths


p3=plot_grid(p_land_585,p_land_latitude_585,p_land_overall_585,ncol=3,rel_heights = c(1,0.4,0.4),
             rel_widths  = c(1,0.4,0.8))


p_climate_585=ggplotGrob(p_climate_585)# don't have to convert it
p_climate_latitude_585=ggplotGrob(p_climate_latitude_585)
p_climate_overall_585=ggplotGrob(p_climate_overall_585)
p_climate_latitude_585$heights=p_climate_overall_585$heights

p_climate_latitude_585$heights=p_climate_overall_245$heights

p_climate_overall_585$heights=p_land_overall_245$heights

p_climate_overall_585$widths=p_land_overall_245$widths
p_land_overall_585$widths=p_land_overall_245$widths

p4=plot_grid(p_climate_585,p_climate_latitude_585,p_climate_overall_585,ncol=3,rel_heights = c(1,0.4,0.4),
             rel_widths  = c(1,0.4,0.8))

load("~/soil-sar/plot-sar-permutation/land use effect/coords_present_new.RData")
tem_df_rcp245_climate$variable=factor(tem_df_rcp245_climate$variable,
                                      levels=c("AM","EM","soilsap","littersap" ,
                                               "woodsap","plapat" ,"para" ,"epiphy","all"))



plot_grid(p1,p2,p3,p4,ncol=1)

