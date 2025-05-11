
grid_level_biomes=readRDS(file="grid_level_biomes.rds")

load(file="data_both_effect_rcp585.rds")
load(file="data_both_effect_rcp245.rds")

scenario=c("SSP2-4.5","SSP5-8.5","SSP2-4.5","SSP5-8.5")

map_data=c("data_both_effect_rcp245", "data_both_effect_rcp585")


my_function_direction_overall=function(data){
  
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
    mutate(group=case_when(joint_effect=="L-loss _ C-loss"~"Climate loss.Land loss", 
                           joint_effect=="L-gain _ C-gain"~"Climate gain.Land gain", 
                           joint_effect%in%c("L-loss _ C-NA",  "L-loss _ C-none")~"Climate NA.Land loss", 
                           joint_effect%in%c( "L-none _ C-loss", "NA _ C-loss" )~"Climate loss.Land NA",
                           joint_effect%in%c(  "L-gain _ C-NA",  "L-gain _ C-no"  )~"Climate NA.Land gain",
                           joint_effect=="L-none _ C-gain"~"Climate gain.Land NA",
                           joint_effect=="L-loss _ C-gain"~"Climate gain.Land loss", 
                           joint_effect=="L-NA _ C-NA"~"NA", 
                           joint_effect=="L-gain _ C-loss"~"Climate loss.Land gain",
                           TRUE~"other"#for example, L-none _ C-NA 
    ))->df4
  #df4=my_function(df4)
  return(df4)
}


my_function_project=function(data)
{
  if ("x"%in%colnames(data))
  {
    points <- vect(data, geom = c("x", "y"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    raster <- rasterize(points, raster_template, field = "group")
  }
  else if ("binary"%in%colnames(data))
  {
    points <- vect(data, geom = c("lon", "lat"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    raster <- rasterize(points, raster_template, field = "binary")
  }
  
  else if ("type"%in%colnames(data))
  {
    points <- vect(data, geom = c("lon", "lat"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
    raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
    raster <- rasterize(points, raster_template, field = "type")
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


full_name=c( "AM" , "EM", "Soil saprotrophs", "Litter saprotrophs","Wood saprotrophs","Plant pathogens","Parasite",  "Epiphyte" , "All") 


common_theme=theme(legend.position = c(0.02611,0.4660),
                   legend.margin = margin(t = -5, r = -5, b = -5, l = 0),
                   legend.text = element_text(size=12),
                   legend.title  = element_text(size=12),
                   text = element_text(size = 18),
                   plot.title = element_text(size = 15, hjust = 0.5), 
                   axis.text.y = element_blank(), 
                   axis.text.x = element_blank(), 
                   axis.title.y = element_text(size = 15), 
                   axis.title.x = element_text(size = 15), 
                   axis.ticks.x = element_blank(), 
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0.2, -2, 0.2), "cm"),
                   panel.background = element_rect(fill = "NA"),
                   panel.border = element_blank())

theme_grid_plot=theme(axis.text = element_text(size = 10,hjust=0),
                      plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(angle=90),
                      axis.text.y = element_text(hjust=0),
                      axis.title = element_text(size = 12),
                      axis.ticks = element_blank(),
                      panel.grid = element_blank())



map_direction_overall_guild=list()
for (m in 1:9){
  map_direction_overall=list()
  for (i in 1:2)
  {
    data=get(map_data[i])[[m]]    
    df5=my_function_direction_overall(data)  
    #to convert some categories into NA
    
    df5%>%mutate(group = if_else(grepl(c("NA|other" ),group), "NA", group))->df5
    
    df5=my_function_project(df5)
    
    df5$group=factor(df5$group,levels=c("0","1","2","3"))
    
    data_percent=table(df5$group)[c(1:4)]/sum(table(df5$group)[c(1:4)])->d
    d%>%data.frame()%>%
      mutate(group=c("CG.LG",
                     "CG.LL",
                     "CL.LG", 
                     "CL.LL"
      ))->data_percent
    
    color_palette1<- c("#1E61A5","#2D9C74","#A58FAA","#E27B60")
    
    
    #need to make sure that the names aligns with initial groups
    set.seed(42)
    data <- expand.grid(
      Variable1 = factor(c( "CL","CG")),
      Variable2 = factor(c(  "LL","LG")))
    
    data$Value <- runif(nrow(data))  #
    
    data$Combination <- interaction(data$Variable1, data$Variable2)
    
    # the order will be based on the grids
    # reorder the columns so that the cells correponsed
    data%>%left_join(data_percent%>%dplyr::rename(Combination=group),by="Combination")%>%
      arrange(Var1)%>%mutate(color=color_palette1 )%>%
      mutate(proportion=round(Freq,3)*100)->data
    
    # Plotting the grid with color mapping
    p2=ggplot(data, aes(x = Variable1, y = Variable2, fill = Combination)) +
      geom_tile(color = "black") +  # Add white borders for clarity
      scale_fill_manual("",breaks=data$Combination,values = data$color) +  # Apply the custom colors
      labs(x = "Variable 1", y = "Variable 2", fill = "Combination") +
      theme_minimal() +theme_grid_plot+
      guides(fill="none")+
      geom_text(aes(label = paste0(sprintf("%.1f", proportion) )), size = 4) +
      xlab("")+
      ylab("")+
      ggtitle("% of pixels")+
      theme(axis.text.y = element_text(size=10.5,color="black"), 
            axis.text.x = element_text(size=10.5,color="black"),
            plot.title = element_text(size = 14.5, hjust = 0.5),
            plot.margin = unit(c(0,0, 0, -1), "cm"))
    
    
    smaller_panel_grob <- ggplotGrob(p2)
    
    
    map_direction_overall[[i]]=ggplot()+
      geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
      scale_color_manual("",breaks=data$Var1,
                         labels=c(paste0(data$Combination ," (",round(data$Freq*100,1),"%)")),
                         values=data$color )+
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
      common_theme+
      guides(color="none")+
      annotation_custom(grob = smaller_panel_grob, xmin = -2740567 , xmax =  -4540567 , ymin =-707265.6, ymax = 2612262.7)+
      xlab(full_name[m])+
      #xlab("")+
      ylab("")+
      ggtitle(scenario[[i]])
  }
  map_direction_overall_guild[[m]]=map_direction_overall
}


# for the barplot showing the proportion of diversity losses and gains in each biome


select_biome=c("Temperate Grasslands, Savannas & Shrublands" , "Temperate Conifer Forests" ,  
               "Temperate Broadleaf & Mixed Forests" ,"Tropical & Subtropical Dry Broadleaf Forests")

##to get the proportional for each biome
pro_data_scenario_guild=list()
for (m in 1:9){
  
  pro_data_scenario=list()
  for ( i in 1:2)
  {
    
    data=get(map_data[i])[[m]]    
    df5=my_function_direction_overall(data)  
    
    df5%>%bind_cols(grid_level_biomes)%>%
      filter(LABEL%in%select_biome)%>%group_by(LABEL,group)%>%
      summarise(Count = n(), .groups = "drop")%>%
      filter(group%in%c("Climate gain.Land gain","Climate loss.Land loss"))->temp_data
    #add one row
    
    
    # get the total cell within each biome
    
    table(grid_level_biomes$LABEL)%>%
      data.frame()%>%filter(Var1%in%select_biome)%>%
      rename_all(~paste0(c("LABEL","number")))->number_cell_biome
    
    #to merge the data
    
    temp_data%>%left_join(number_cell_biome,by="LABEL")%>%
      mutate(pro=Count/number)->pro_data
    
    pro_data_scenario[[i]] =pro_data
  }
  pro_data_scenario_guild[[m]]=pro_data_scenario
}

#create the bar plots in the moderate-emission scenarios

biome_bar_plot_low=list()
for (m in 1:9){
  
  data_set1=pro_data_scenario_guild[[m]][[1]]
  
  data_set1$LABEL=factor(data_set1$LABEL, levels=c("Temperate Broadleaf & Mixed Forests",
                                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                                   "Tropical & Subtropical Dry Broadleaf Forests"))
  data_set1$BarWidth <- ifelse(duplicated(data_set1$LABEL) | duplicated(data_set1$LABEL, fromLast = TRUE), 0.8, 0.4)
  
  biome_bar_plot_low[[m]]=ggplot(data=data_set1,aes(fill=group,x=LABEL ,y=100*pro))+
    geom_bar(stat = "identity",position = position_dodge(width =0.8 ),width = data_set1$BarWidth) +
    scale_fill_manual("",breaks=c("Climate gain.Land gain","Climate loss.Land loss"),labels=c("Gains","Losses"),values=c( "#1E61A5","#E27B60"))+
    geom_text(aes(label = paste0("(", sprintf("%.1f", pro * 100), ")")),
              position = position_dodge(width = 0.8), 
              vjust = -1.25, hjust=0.5,size = 4)+
    theme(legend.position = c(0.05,0.9),
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(size = 12), 
          axis.text.x = element_text(size = 14,angle=45,hjust=1,vjust=1,color="black"), 
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0, 0.2, 1, 0.3), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
    scale_x_discrete(full_name[m],breaks=unique(data_set1$LABEL),
                     labels=c("Broadleaf-mixed\n forest biome","Conifer forest\n biome","Grassland\n biome","Dry forest\n biome"))+
    ylab("% of pixels")+
    
    ylim(0,55)+
    ggtitle("SSP2-4.5")
  
  
}


##create the bar plots in the high-emission scenarios
biome_bar_plot_high=list()
for (m in 1:9){
  
  data_set2=pro_data_scenario_guild[[m]][[2]]
  
  data_set2$LABEL=factor(data_set2$LABEL, levels=c("Temperate Broadleaf & Mixed Forests",
                                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                                   "Tropical & Subtropical Dry Broadleaf Forests"))
  
  data_set2$BarWidth <- ifelse(duplicated(data_set2$LABEL) | duplicated(data_set2$LABEL, fromLast = TRUE), 0.8, 0.4)
  
  biome_bar_plot_high[[m]]=ggplot(data=data_set2,aes(fill=group,x=LABEL ,y=100*pro))+
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = data_set2$BarWidth) +
    scale_fill_manual("",breaks=c("Climate gain.Land gain","Climate loss.Land loss"),labels=c("Gains","Losses"),values=c(  "#1E61A5","#E27B60"))+
    geom_text(aes(label = paste0("(", sprintf("%.1f", pro * 100), ")")),
              position = position_dodge(width = 0.8), 
              vjust = -1.25, hjust=0.5,size = 4)+
    theme(legend.position = c(0.05,0.9),
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(size = 12), 
          axis.text.x = element_text(size = 14,angle=45,hjust=1,vjust=1,color="black"), 
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0, 0.2, 1, 0.3), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
    scale_x_discrete(full_name[m],breaks=unique(data_set2$LABEL),
                     labels=c("Broadleaf-mixed\nforest biome","Conifer forest\n biome","Grassland\n biome","Dry forest\n biome"))+
    ylab("% of pixels")+
    
    ylim(0,55)+
    ggtitle("SSP5-8.5")
}

#to bind all the plots
plot_grid(
  map_direction_overall_guild[[9]][[1]],
  map_direction_overall_guild[[9]][[2]],
  biome_bar_plot_low[[9]],
  biome_bar_plot_high[[9]],
  ncol=2,
  rel_widths = c(1,1,0.5,0.5),
  label_x = 0.1,label_y = 1,label_size = 16,
  labels = paste0("(", letters[1:4], ")"))

# the final dimension for the output figure is 9 x 9 

plot_grid(map_direction_overall_guild[[1]][[1]],
          map_direction_overall_guild[[1]][[2]],
          map_direction_overall_guild[[2]][[1]],
          map_direction_overall_guild[[2]][[2]],
          map_direction_overall_guild[[3]][[1]],
          map_direction_overall_guild[[3]][[2]],
          map_direction_overall_guild[[6]][[1]],
          map_direction_overall_guild[[6]][[2]],
          label_x = 0.1,
          label_y = 0.85,
          label_size = 16,ncol=4,
          labels = paste0("(", letters[1:8], ")"))

# for the bar plots

plot_grid(biome_bar_plot_low[[1]],
          biome_bar_plot_high[[1]],
          biome_bar_plot_low[[2]],
          biome_bar_plot_high[[2]],
          biome_bar_plot_low[[3]],
          biome_bar_plot_high[[3]],
          biome_bar_plot_low[[6]],
          biome_bar_plot_high[[6]],
          label_x = 0.1,
          label_y = 1.02,
          label_size = 16,ncol=4,
          labels = paste0("(", letters[1:8], ")"))

# the final dimension for the output figure is 10 x 20 
