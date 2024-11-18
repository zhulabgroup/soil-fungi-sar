
# the color selected
#color_palette<- c("#E3A72F","#CD326D", "#F8A488",  "#4A90E2", "#32CD92", "#FFD166", "#C7F9CC", "#E6E6FA", "#B0BEC5")



color_palette=c( "#CD326D","#B0BEC5" ,"#32CD92", "#4A90E2", "#E6E6FA","#C7F9CC","#E3A72F","#F8A488" , "mediumpurple")

theme_grid_plot=theme(axis.text = element_text(size = 10,hjust=0),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle=90),
      axis.text.y = element_text(hjust=0),
      axis.title = element_text(size = 12),
      axis.ticks = element_blank(),
      panel.grid = element_blank())

# when only the direction of the impact was considered
# nine categories were classified
# when both do not have data, they were classified as NA
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




map_direction_overall=list()
  for (i in 1:2)
  {
    data=get(map_data[i])[[m]]    
    df5=my_function_direction_overall(data)  
    df5=my_function_project(df5)
    df5$group=factor(df5$group,levels=c("0","1","2","3","4","5","6","7","8"))
    
    data_percent=table(df5$group)[c(1:9)]/sum(table(df5$group)[c(1:9)])->d
    d%>%data.frame()%>%
      mutate(group=c("Climate gain.Land gain",
                     "Climate gain.Land loss",
                     "Climate gain.Land other",
                     "Climate loss.Land gain", 
                     "Climate loss.Land loss",
                     "Climate loss.Land other",
                     "Climate other.Land gain",
                     "Climate other.Land loss",
                     "Climate other.Land other"))->data_percent
   
  
    #need to make sure that the names aligns with initial groups
    set.seed(42)
    data <- expand.grid(
      Variable1 = factor(c( "Climate loss","Climate other","Climate gain")),
      Variable2 = factor(c(  "Land loss","Land other","Land gain")))
    
    data$Value <- runif(nrow(data))  #
    
    data$Combination <- interaction(data$Variable1, data$Variable2)
    
    # the order will be based on the grids
    # reorder the columns so that the cells correponsed
    data%>%left_join(data_percent%>%dplyr::rename(Combination=group),by="Combination")%>%
      arrange(Var1)%>%mutate(color=color_palette )%>%
      mutate(proportion=round(Freq,3)*100)->data
 
    
    # Plotting the grid with color mapping
    p2=ggplot(data, aes(x = Variable1, y = Variable2, fill = Combination)) +
      geom_tile(color = "black") +  # Add white borders for clarity
      scale_fill_manual("",breaks=data$Combination,values = data$color) +  # Apply the custom colors
      labs(x = "Variable 1", y = "Variable 2", fill = "Combination") +
      theme_minimal() +theme_grid_plot+
      guides(fill="none")+
      geom_text(aes(label = paste0(sprintf("%.1f", proportion) )), size = 3) +
      xlab("")+
      ylab("")+
      ggtitle("% of pixels")
    smaller_panel_grob <- ggplotGrob(p2)
    if(i<2){
      
      #group=c("Climate gain.Land gain","Climate gain.Land loss","Climate gain.Land other","Climate loss.Land gain", "Climate loss.Land loss","Climate loss.Land other","Climate other.Land gain","Climate other.Land loss")
     
      map_direction_overall[[i]]=ggplot()+
        geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
        scale_color_manual("",breaks=data$Var1,
                           labels=c(paste0(data$Combination ," (",round(data$Freq*100,1),"%)")),
                           values=data$color )+
        #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
        ggtitle("Direction of richness change \n(RCP4.5-SSP2)")+
        guides(color = guide_legend(override.aes = list(size = 2)))+
        #xlab(paste(guild_names[m]))+
        ylab("")+
        xlab("")+
        geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
        geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
        geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
        geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
        geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
        geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
        geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
        coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
        common_theme+
        guides(color="none")+
        annotation_custom(grob = smaller_panel_grob, xmin = -5640567 , xmax = -2509970, ymin =-727265.6, ymax = 3012262.7)
    }
    
    else{
      map_direction_overall[[i]]=ggplot()+
        geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
        scale_color_manual("",breaks=data$Var1,
                           labels=c(paste0(data$Combination ," (",round(data$Freq*100,1),"%)")),
                           values=data$color )+
        #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
        ggtitle("Direction of richness change \n(RCP8.5-SSP5)")+
        guides(color = guide_legend(override.aes = list(size = 2)))+
        #xlab(paste(guild_names[m]))+
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
        common_theme+
        guides(color="none")+
        annotation_custom(grob = smaller_panel_grob, xmin = -5640567 , xmax = -2509970, ymin =-727265.6, ymax = 3012262.7)
    }
  }
  

  #when considerting the magnitude of both effects


  
  my_function_guild_loss=function(data){
    data%>%filter(climate_effect<0)->data_loss_climate 
    data%>%filter(climate_effect>0)->data_gain_climate
    summary(data_loss_climate$climate_effect,na.rm=TRUE)->loss_percentile_climate 
    summary(data_gain_climate$climate_effect,na.rm=TRUE)->gain_percentile_climate# to assign a category for each cell based on the magnitude
    data%>%mutate(type_C=case_when(climate_effect<=loss_percentile_climate[2]%>%as.numeric()~"C-loss-high", 
                                   climate_effect>loss_percentile_climate[5]%>%as.numeric()&climate_effect<0~"C-loss-low", 
                                   climate_effect>=gain_percentile_climate[5]%>%as.numeric()~"C-gain-high",
                                   climate_effect<gain_percentile_climate[2]%>%as.numeric()&climate_effect>0~"C-gain-low",
                  is.na(climate_effect)~"NA",
                  TRUE~"other" ))->df2# for the impact of land use change
    
    data%>%filter(land_effect<0)->data_loss_land 
    data%>%filter(land_effect>0)->data_gain_land
    summary(data_loss_land$land_effect,na.rm=TRUE)->loss_percentile 
    summary(data_gain_land$land_effect,na.rm=TRUE)->gain_percentile
    data%>%mutate(type_L=case_when(land_effect<=loss_percentile[2]%>%as.numeric()~"L-loss-high", land_effect>loss_percentile[5]%>%as.numeric()&land_effect<0~"L-loss-low", 
                                   land_effect>=gain_percentile[5]%>%as.numeric()~"L-gain-high",
                                   land_effect<gain_percentile[2]%>%as.numeric()&land_effect>0~"L-gain-low",
                                   is.na(land_effect)~"NA",
                                   TRUE~"other" 
                                   ))->df1# combing both effects of land use conversion and climate change
   
     bind_cols(df2%>%dplyr::select(lon,lat,type_C),df1%>%dplyr::select(type_L))->df3 
    
    df3%>%mutate(joint_effect=paste(df3$type_L,"_",df3$type_C))%>%
      mutate(group=case_when(joint_effect=="L-loss-high _ C-loss-high"~"Climate loss high.Land loss high", 
                             
                            joint_effect=="L-gain-high _ C-loss-high"~"Climate loss high.Land gain high", 
                             
                            joint_effect=="L-gain-high _ C-gain-high"~"Climate gain high.Land gain high", 
                            
                            joint_effect=="L-gain-high _ other"~"Climate other.Land gain high", 
                            
                            joint_effect=="L-loss-high _ C-gain-high"~"Climate gain high.Land loss high", 
                           
                             joint_effect=="L-loss-high _ other"~"Climate other.Land loss high", 
                            
                            joint_effect=="other _ C-gain-high"~"Climate gain high.Land other", 
                            joint_effect=="other _ C-loss-high"~"Climate loss high.Land other", 
                            # joint_effect=="other _ other"~"Climate other.Land other", 
                                 joint_effect=="NA _ NA"~"NA",
                             TRUE~"other"))->df4
    #df4=my_function(df4)
    return(df4)
  }
  
  
  
  ggplot()+geom_point(data=df5,pch=15,aes(x=lon,y=lat,color=group),size=0.01)+
    guides(color = guide_legend(override.aes = list(size = 2)))+
  
    scale_color_manual(breaks=c(unique(df5$group)),
                       labels=c(unique(df5$group)),
                       values=c("white","brown","yellow","green","purple","blue","gold","pink","red","seagreen"))
  
  
  #data=get(map_data[i])[[9]]    
  #df5=my_function_guild_loss(data) 
  
  #df6=my_function_direction_overall(data) # the resulting na-na is correct
  
  #color_palette1=c("mediumpurple", "#32CD92","#C7F9CC","#4A90E2","#CD326D","#E3A72F", "#F8A488","#B0BEC5" ,"#E6E6FA")
  
  
  
  map_data=c("data_both_effect_rcp245","data_both_effect_rcp585")
  
  map_loss=list()
  for (i in 1:2)
  {
    data=get(map_data[i])[[9]]    
    df5=my_function_guild_loss(data) 
    df5=my_function_project(df5)
    df5$group=as.factor(df5$group)
    data_percent=table(df5$group)[1:9]/sum(table(df5$group)[1:9])->d
    d%>%data.frame()%>%
      mutate(group=c("Climate gain high.Land gain high",
                     "Climate gain high.Land loss high",
                     "Climate gain high.Land other",
                     "Climate loss high.Land gain high",
                     "Climate loss high.Land loss high",
                     "Climate loss high.Land other",
                     "Climate other.Land gain high",
                     "Climate other.Land loss high",
                     "Climate other.Land other"
                     ))->data_percent
    
    set.seed(42)
    data <- expand.grid(
      Variable1 = factor(c( "Climate loss high","Climate other","Climate gain high")),
      Variable2 = factor(c( "Land loss high","Land other",  "Land gain high"))
    )
    data$Value <- runif(nrow(data))  #
    
    data$Combination <- interaction(data$Variable1, data$Variable2)
  
    data_percent%>%dplyr::rename(Combination=group)%>%
      left_join(data,by="Combination")%>%
      mutate(color=color_palette1 )%>%
      mutate(proportion=round(Freq,3)*100)->data
    
    # Plotting the grid with color mapping
    p2=ggplot(data, aes(x = Variable1, y = Variable2, fill = Combination)) +
      geom_tile(color = "black") +  # Add white borders for clarity
      scale_fill_manual("",breaks=data$Combination,values = data$color) +  # Apply the custom colors
      labs(x = "Variable 1", y = "Variable 2", fill = "Combination") +
      theme_minimal() +
      theme_grid_plot+
      guides(fill="none")+
      geom_text(aes(label = paste0(sprintf("%.1f", proportion) )), size = 3) +
      xlab("")+
      ylab("")+
      ggtitle("% of pixels")
    
    smaller_panel_grob <- ggplotGrob(p2)
    
    if(i<2){
      
      map_loss[[i]]=ggplot()+
        geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
        scale_color_manual("",breaks=data$Var1,
                           labels=c(paste0(data$Combination ," (",round(data$Freq*100,1),"%)")),
                           values=data$color )+
        
        #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
        ggtitle("Magnitude of richness change \n(RCP4.5-SSP2)")+
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
      
        common_theme+
        guides(color="none")+
        annotation_custom(grob = smaller_panel_grob, xmin = -6040567 , xmax = -2509970, ymin =-907265.6, ymax = 3212262.7)
      
      
    }
    else{
      map_loss[[i]]=ggplot()+
        geom_point(data=df5,pch=15,aes(x=x,y=y,color=group),size=0.01)+
        scale_color_manual("",breaks=data$Var1,
                           labels=c(paste0(data$Combination ," (",round(data$Freq*100,1),"%)")),
                           values=data$color )+
        
        #geom_sf(data=st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0),linetype = "solid")+
        ggtitle("Magnitude of richness change \n(RCP8.5-SSP5)")+
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
        common_theme+
        guides(color="none")+
        annotation_custom(grob = smaller_panel_grob, xmin = -6040567 , xmax = -2509970, ymin =-907265.6, ymax = 3212262.7)
      
    }
  }
  
  colors <- c("#2C3E50", "#1ABC9C", "#F1C40F", "#E74C3C", "#95A5A6", 
              "#9B59B6", "#A2D9CE", "#D5B895", "#8E44AD")

 plot_grid(map_direction_overall[[1]],map_direction_overall[[2]],map_loss[[1]],map_loss[[2]],ncol=2,
           label_x = 0.08,label_y = 0.95,label_size = 20,labels = paste0("(", letters[1:4], ")"))
  
  