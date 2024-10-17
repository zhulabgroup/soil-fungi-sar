# create the new figures when part of the regions have been removed
#for the land use effect
#species_change_land_rcp245_all
#based on the land use effect to mask the climate change effect

land_effect=bind_cols(coords_present,species_change_land_rcp245%>%filter(variable=="all"))
crs(r_climate) <- CRS("+proj=longlat +datum=WGS84")
r_land <- rasterFromXYZ(land_effect[, c("lon", "lat", "value")])


crop_data=list()
for (i in 1:9)
{
  climate_effect=bind_cols(coords_present,species_change_climate_rcp245%>%filter(variable==guild_type[i]))
  r_climate <- rasterFromXYZ(climate_effect[, c("lon", "lat", "value")])
  crs(r_climate) <- CRS("+proj=longlat +datum=WGS84")
  masked_raster <- mask(r_climate, r_land)#mask the climate change effect
  raster_df <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE)
  colnames(raster_df) <- c("x", "y", "value")
  crop_data[[i]]=raster_df%>%mutate(variable=rep(guild_type[i],dim(raster_df)[1]))
}



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

data_both_effect_rcp245=list()
 for (i in 1:9)
 {
   crop_data[[i]]%>%mutate(coordinate=paste(round(x,3),"_",round(y,3)))->climate_change_effect_crop_245
   species_change_land_rcp245%>%filter(variable==guild_type[i])%>%
     bind_cols(coords_present)%>%mutate(coordinate=paste(lon,"_",lat))->land_change_effect_245
   data_both_effect_rcp245[[i]]=land_change_effect_245%>%
     left_join(climate_change_effect_crop_245%>%dplyr::select(coordinate,value),by="coordinate")%>%
     rename(land_effect=value.x,climate_effect=value.y)
 }

data_both_effect_rcp585=list()
for (i in 1:9)
{
  crop_data_rcp585_climate[[i]]%>%mutate(coordinate=paste(round(x,3),"_",round(y,3)))->climate_change_effect_crop_585
  species_change_land_rcp585%>%filter(variable==guild_type[i])%>%
    bind_cols(coords_present)%>%mutate(coordinate=paste(lon,"_",lat))->land_change_effect_585
  data_both_effect_rcp585[[i]]=land_change_effect_585%>%
    left_join(climate_change_effect_crop_585%>%dplyr::select(coordinate,value),by="coordinate")%>%
    rename(land_effect=value.x,climate_effect=value.y)
}

common_theme=theme(legend.position = c(0.15,0.5660),
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
                   plot.margin = unit(c(0,0.3, 0, 0), "cm"),
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
      geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP4.5-SSP2",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","#E3A72F",  "#32CD92","#CD326D","gray88","white"))+
      
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
      geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP5.8-SSP5",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","#E3A72F","#32CD92","#CD326D","gray88","white"))+
      
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
      geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("RCP4.5-SSP2",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","#E3A72F","#32CD92","#CD326D","gray88","white"))+
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
      geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
      
      scale_color_manual("RCP5.8-SSP5",breaks=c("0","1","2","3","4","5"),
                         labels=c(paste0(data_percent$group ," (",round(data_percent$Freq*100,1),"%)"),""),
                         values=c("#E92316","#E3A72F","#32CD92","#CD326D","gray88","white"))+
      
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


