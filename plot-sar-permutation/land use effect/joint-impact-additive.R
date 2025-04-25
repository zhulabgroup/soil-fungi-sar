# create the map showing the additive effect of both factors


my_function_project=function(data)
{
  points <- vect(data, geom = c("lon", "lat"), crs = "EPSG:4326")  # Assuming WGS84 coordinates
  raster_template <- rast(ext(points), resolution = 0.17, crs = "EPSG:4326")  # Resolution of 1 degree
  raster <- rasterize(points, raster_template, field = "joint_impact")
  target_crs <- "EPSG:5070"
  raster_equal_area <- project(raster, target_crs,method="near")
  raster_df <- as.data.frame(raster_equal_area, xy = TRUE,)
  return(raster_df)
}
## get the data for the joint impact



select_biome=c("Temperate Grasslands, Savannas & Shrublands",
               "Temperate Conifer Forests",
               "Temperate Broadleaf & Mixed Forests",
               "Tropical & Subtropical Dry Broadleaf Forests")

guild_type=c("AM"  , "EM"  , "soilsap" ,  "littersap", "woodsap",   "plapat" ,   "para" ,     "epiphy" ,   "all" )

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

load(file="data_both_effect_rcp245.rds")

load(file="data_both_effect_rcp585.rds")


map_data=c("data_both_effect_rcp245","data_both_effect_rcp585")


### add a column to indicate the joint impact of both factors
## assuming both effects are additive

joint_impact=list()
for(i in 1:2)
{
  joint_impact_temp=list()
  for(m in 1:9)
  {
    get(map_data[i])[[m]]%>%
      mutate(joint_impact=land_effect+climate_effect)->joint_impact_temp[[m]]# for one scenario
    
  }
  joint_impact[[i]]=joint_impact_temp
}

# indicates the climate change scenarios
# to bind all the guilds
joint_data_for_pie=list()

for (i in 1:2)
{
  combined_df <- do.call(rbind, joint_impact[[i]])%>%dplyr::select(variable,joint_impact)%>%
    dplyr::rename(value=joint_impact)#in the low emission scenario
  joint_data_for_pie[[i]]=combined_df
}



# assuming them being additive
#mapping the 

joint_impact[[1]][[9]]->temp

joint_impact[[2]][[9]]->temp_585

#need to make projections
joint_all_guild=my_function_project(temp)

joint_all_guild_585=my_function_project(temp_585)


joint_all_guild%>%mutate(change=100*last)->joint_all_guild

joint_all_guild_585%>%mutate(change=100*last)->joint_all_guild_585



p_joint_rcp245=ggplot(joint_all_guild) +
  geom_tile(data =joint_all_guild ,aes(x = x, y = y, fill = change), size = 0.175)+
  scale_fill_gradientn(NULL,
                       limits = c(-70, 70),
                       colors = brewer.pal(11, "RdBu"),
                       guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(joint_all_guild, change < -70), mapping = aes(x=x,y=y,fill = last > 0.7)) +
  scale_fill_manual("Change (%)", values = "navy", labels = "> 70", 
                    guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
  ) +
  new_scale_fill() +
  geom_tile(data = filter(joint_all_guild, change < -70), mapping = aes(x=x,y=y,fill = last< -0.7)) +
  scale_fill_manual(NULL, values = "gold", labels = "< -70", 
                    guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
  theme_map+ 
  add_sf_layers()+
  coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
  xlab("")+
  ylab("Combined impact")+
  ggtitle("SSP2−4.5")


###


p_joint_rcp585=ggplot(joint_all_guild_585) +
  geom_tile(data =joint_all_guild_585 ,aes(x = x, y = y, fill = change), size = 0.175)+
  scale_fill_gradientn(NULL,
                       limits = c(-70, 70),
                       colors = brewer.pal(11, "RdBu"),
                       guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(joint_all_guild_585, change < -70), mapping = aes(x=x,y=y,fill = last > 0.7)) +
  scale_fill_manual("Change (%)", values = "navy", labels = "> 70", 
                    guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
  ) +
  new_scale_fill() +
  geom_tile(data = filter(joint_all_guild_585, change < -70), mapping = aes(x=x,y=y,fill = last< -0.7)) +
  scale_fill_manual(NULL, values = "gold", labels = "< -70", 
                    guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
  theme_map+ 
  add_sf_layers()+
  coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
  xlab("")+
  ylab("Combined impact")+
  ggtitle("SSP5−8.5")


#the latitude patterns of the change rate

my_function_latitude_patterns=function(data)
{
  data%>% group_by(lat) %>%
    summarise(mean_value = mean(joint_impact,na.rm=TRUE),sd_value = sd(joint_impact,na.rm=TRUE))->latitude_patterns
  return(latitude_patterns )
}


summary_data_joint_rcp245=my_function_latitude_patterns(temp)
summary_data_joint_rcp585=my_function_latitude_patterns(temp_585)

p_joint_latitude_245=ggplotGrob(ggplot()+
                                    geom_ribbon(data = summary_data_joint_rcp245, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "gray80")+
                                    geom_line(data = summary_data_joint_rcp245, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) +  # Mean trend line
                                    theme_latitude+
                                    xlab("Latitude")+
                                    ylab("")+
                                    ggtitle("")+
                                    coord_flip()+
                                    geom_hline(yintercept = 0,linetype="dashed")+
                                    ggtitle("SSP2−4.5")+
                                    xlim(15,65)+
                                    ylim(-40,40))


p_joint_latitude_585=ggplotGrob(ggplot()+
                                  geom_ribbon(data = summary_data_joint_rcp585, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "gray80")+
                                  geom_line(data = summary_data_joint_rcp585, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) +  # Mean trend line
                                  theme_latitude+
                                  xlab("Latitude")+
                                  ylab("")+
                                  ggtitle("")+
                                  coord_flip()+
                                  geom_hline(yintercept = 0,linetype="dashed")+
                                  ggtitle("SSP5−8.5")+
                                  xlim(15,65)+
                                  ylim(-40,40))







#(1) get the proportions of pixels that show either species gain or loss


biomes_four=c( "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests" ,                    
               "Temperate Grasslands, Savannas & Shrublands" ,"Tropical & Subtropical Dry Broadleaf Forests")

joint_change_rate_guild_scenario=list()

for (j in 1:2){#climate change scenarios
  
# each of the guild under the low 
  change_rate_guild=list()
  for (m in 1:9)#guilds
  {
    
    change_rate=matrix(ncol=3,nrow=5)
    
    for (i in 1:5)
    {
      if(i==5)#all biomes were included
      {
        joint_data_for_pie[[j]]%>%filter(variable==guild_type[m])%>%
          bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
          filter(LABEL%in%biomes_four&!is.na(value))->df1
        
        df1%>%filter(value>0)%>%dim()%>%head(1)->gain
        df1%>%filter(value==0)%>%dim()%>%head(1)->no_change
        df1%>%filter(value<0)%>%dim()%>%head(1)->loss
      }
      else{# select just one biome
        joint_data_for_pie[[j]] %>%filter(variable==guild_type[m])%>%
          bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
          filter(LABEL== biomes_four[i]&!is.na(value))->df1
        
        df1%>%filter(value>0)%>%dim()%>%head(1)->gain
        df1%>%filter(value==0)%>%dim()%>%head(1)->no_change
        df1%>%filter(value<0)%>%dim()%>%head(1)->loss
      }
      
      change_rate[i,1]=gain
      change_rate[i,2]=no_change
      change_rate[i,3]=loss
    }
    change_rate_guild[[m]]=change_rate
  }
  joint_change_rate_guild_scenario[[j]]=change_rate_guild
}
  

joint_change_rate_guild_scenario[[1]]#means joint effect in the low-emission scenario for all nine guilds

saveRDS(joint_change_rate_guild_scenario,file="joint_change_rate_guild_scenario.rds")



#for each guild, to get the data to create the pie plot
pie_data_guild_joint=list()
for (m in 1:9){
  pie_data=list()
  for (i in 1:2)# for both climate and land cover
  {
    joint_change_rate_guild_scenario[[i]][[m]]%>%
      data.frame()%>%rename_all(~paste0(c("gain","no","loss")))%>%
      mutate(biome=c("Temperate","Conifer","Grassland","Dry","All"))%>%
      melt()%>%group_by(biome)%>%mutate(percent = value / sum(value))->temp_pie_data
    temp_pie_data$biome=factor( levels=c("All" ,"Temperate", "Conifer" ,  "Grassland", "Dry" ),temp_pie_data$biome)
    
    pie_data[[i]]=temp_pie_data
  }
  pie_data_guild_joint[[m]] =pie_data
}

save(pie_data_guild_joint,file="pie_data_guild_joint.rds")



# for the climate change impact, 41394 cells were included
# estimate the mean diversity change rate among biomes 

observations=c(41140,41140,41394,41394)


guild_type=c("AM"  , "EM"  , "soilsap" ,  "littersap", "woodsap",   "plapat" ,   "para" ,     "epiphy" ,   "all" )

com_data_climate_guild_joint=list()
for (m in 1:9){
  
  com_data_climate=list()
  
  for (j in 1:2){
    
    joint_data_for_pie[[j]]%>%
      bind_cols(rep(grid_level_biomes$LABEL,9))%>%
      rename_all(~paste0(c("guild","value","LABEL")))%>%
      filter(LABEL%in%biomes_four&!is.na(value))%>%
      mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered
    
    all_biome=dim(df1)[1]
    
    # across all the biomes 
    
    df1%>%mutate(LABEL_all=rep("all",all_biome))%>%
      dplyr::select(guild,value,LABEL_all,plotid)%>%
      dplyr::rename(LABEL=LABEL_all)%>%bind_rows(df1)->df1# when it is all, it means across the biomes
    
   
    
    #to give a plot id for each cell
    # to test difference for species gains for the overall fungal diversity
    
    df1%>%filter(value>0&guild%in%c(guild_type[m]))->df_gain
    
    #need to transform the response variable
    
    df_gain$log_response <- log(df_gain$value)
    model_all <- gls(log_response ~ LABEL, data = df_gain)# performe the model and get the mean for all groups
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="all")%>%mutate(.group="")%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))->temp_mean# the predicted biomes-wide mean
    
    # perform the model just for the four biomes
    
    model_four <- gls(log_response ~ LABEL, data = df_gain%>%filter(LABEL!="all"))
    
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    # to bind the letters with the mean and add a column specifying the different letters
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))%>%
      bind_rows(temp_mean)->data_gain # 
   # the biome-wide mean was not compared with any group
    # the data used for species gains
    
    # for species loss
    
    df1%>%filter(value<0&guild%in%c(guild_type[m]))->df_loss
    
    df_loss$sqrt_response <- sqrt(abs(df_loss$value))
    
    model_all <- gls(sqrt_response~ LABEL, data = df_loss)# perform the model and get the mean for all groups
    
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    
    
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="all")%>%mutate(.group="")%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)->temp_mean
    
    model_four <- gls(sqrt_response ~ LABEL, data = df_loss%>%filter(LABEL!="all"))
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    #back transforming the rate
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)%>%
      bind_rows(temp_mean) ->data_loss
   
    com_data=rbind(data_gain,data_loss)
    com_data%>%mutate(change=ifelse(origin_mean > 0, "Positive", "Negative"))->com_data
    com_data$LABEL=factor(com_data$LABEL,
                          levels=c("all",
                                   "Temperate Broadleaf & Mixed Forests",
                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                   "Tropical & Subtropical Dry Broadleaf Forests")%>%rev())#display order of the texts
    
    com_data_climate[[j]]=com_data%>%dplyr::rename(biome=LABEL)
  }
  com_data_climate_guild_joint[[m]]=com_data_climate
}

save(com_data_climate_guild_joint,file="com_data_climate_guild_joint.rds")



net_effect_joint_guild=list()
for (m in 1:9){
  
net_effect_joint=list()
for (i in 1:2)
{
  joint_data_for_pie[[i]]%>%filter(variable==guild_type[m])%>%bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%select_biome)%>%group_by(LABEL)%>%summarise(overal_mean=mean(value,na.rm=TRUE))%>%
    dplyr::rename(biome=LABEL)->biome_specific_mean
  
  joint_data_for_pie[[i]]%>%filter(variable==guild_type[m])%>%bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%select_biome)%>%filter(!is.na(value))%>%pull(value)%>%mean()%>%
    data.frame()%>%bind_cols("all")%>%rename_all(~paste0(c("overal_mean","biome")))%>%
    dplyr::select(biome,overal_mean)->biome_wide_mean
  net_effect_joint[[i]]=bind_rows(biome_specific_mean,biome_wide_mean)

  }
net_effect_joint_guild[[m]]=net_effect_joint
}


#create the plots
#change the order of the plots
joint_effect_data=list()
for (m in 1:9)
{
  temp=list()
  for (i in 1:2)
  {
    com_data_climate_guild_joint[[m]][[i]]%>%
      left_join(net_effect_joint_guild[[m]][[i]],by="biome")->d
    
    d$biome=factor(d$biome,levels = c("all","Temperate Broadleaf & Mixed Forests",
                                      "Temperate Conifer Forests",
                   "Temperate Grasslands, Savannas & Shrublands",
                   "Tropical & Subtropical Dry Broadleaf Forests")%>%rev()) 
    d->temp[[i]]
  }
  joint_effect_data[[m]]=temp
  }

save(joint_effect_data,file="joint_effect_data.rds")#


joint_effect_data=load(file="joint_effect_data.rds")






pp_pie_guild=list()
for (m in 1:9){
  if(m==1){
    pp_pie=list()
    for (i in 1:2){
      pp_pie[[i]]=ggplotGrob(ggplot(pie_data_guild_joint[[m]][[i]], aes(x=1, y=percent, fill=variable)) +
                               ## geom_col is geom_bar(stat = "identity")(bit shorter)
                               ## use color = black for the outline
                               geom_col(width=5,position="fill", color = "black",size=0.25)+
                               coord_polar("y", start=pi) +
                               #geom_text(aes(x = 5, label = paste0(round(percent*100), "%")), size=2.5, 
                               #position = position_stack(vjust = 0.1))+
                               #labs(x = NULL, y = NULL, fill = NULL, title = "Biomes")+
                               theme(legend.position = c(0.5,0.5),
                                     legend.key.size = unit(0.3, "cm"),
                                     axis.line = element_blank(),
                                     axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     strip.text = element_blank(),
                                     plot.margin = unit(c(0, -3, 0, -3), "cm"),
                                     strip.background = element_blank(),
                                     panel.spacing.y  = unit(0.2, "lines"),
                                     panel.background = element_rect(fill = "NA"),
                                     panel.border = element_blank(),
                                     plot.title = element_text(size = 15, hjust = 0.5))+
                               facet_wrap(~biome,ncol=1)+
                               xlab("")+
                               ylab("")+
                               
                               #ylab(paste(full_name[m]))+
                               #ggtitle("")+
                               scale_fill_manual("",breaks=c("gain","no","loss"),
                                                 labels=c("Gain","No change","Loss"),values=c( "#1E61A5","gray","#E27B60")))}
  }
  else{
    pp_pie=list()
    for (i in 1:2){
      pp_pie[[i]]=ggplotGrob(ggplot(pie_data_guild_joint[[m]][[i]], aes(x=1, y=percent, fill=variable)) +
                               ## geom_col is geom_bar(stat = "identity")(bit shorter)
                               ## use color = black for the outline
                               geom_col(width=5,position="fill", color = "black",size=0.25)+
                               coord_polar("y", start=pi) +
                               #geom_text(aes(x = 5, label = paste0(round(percent*100), "%")), size=2.5, 
                               #position = position_stack(vjust = 0.1))+
                               #labs(x = NULL, y = NULL, fill = NULL, title = "Biomes")+
                               theme(legend.position ="none",
                                     legend.key.size = unit(0.3, "cm"),
                                     axis.line = element_blank(),
                                     axis.text = element_blank(),
                                     axis.ticks = element_blank(),
                                     strip.text = element_blank(),
                                     plot.margin = unit(c(0, -3, 0, -3), "cm"),
                                     strip.background = element_blank(),
                                     panel.spacing.y  = unit(0.2, "lines"),
                                     panel.background = element_rect(fill = "NA"),
                                     panel.border = element_blank(),
                                     plot.title = element_text(size = 15, hjust = 0.5))+
                               facet_wrap(~biome,ncol=1)+
                               xlab("")+
                               ylab("")+
                               
                               #ylab(paste(full_name[m]))+
                               #ggtitle("")+
                               scale_fill_manual("",breaks=c("gain","no","loss"),
                                                 labels=c("Gain","No change","Loss"),values=c("#1E61A5","gray","#E27B60")))}
  }
  pp_pie_guild[[m]]=pp_pie
}

# for the bar plot showing 

full_name=c("AM","EM","Soil saprotrophs","Litter saprotrophs","Wood saprotrophs", "Plant pathogens", "Parasite","Epiphyte","All")


pp_guild=list()
for(i in c(1,2,3,6,9))
  
{
  
  if(i>1)
  {
    pp=list()
    for (m in 1:2){
      
      pp[[m]]=ggplot(data=joint_effect_data[[i]][[m]],aes(fill=change,y=biome ,x=100*origin_mean))+
        geom_col(width = 0.3,color="black",size=0.2)+
        geom_segment(data=joint_effect_data[[i]][[m]], 
                     aes(y=biome,yend=biome,xend=  100*ori_low, x = 100*ori_up  ),
                     arrow = arrow(length = unit(0.1, "cm"),"both",angle=90))+
        
        scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#E27B60", "#1E61A5"))+
        theme(legend.position = "none",
              legend.text = element_text(size=8),
              legend.title  = element_text(size=10),
              text = element_text(size = 18),
              plot.title = element_text(size = 12, hjust = 0.5), 
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.title.x = element_text(size = 10), 
              legend.key.size = unit(0.3, "cm"),
              plot.margin = unit(c(0, 0, 0.3, 0), "cm"),
              panel.background = element_rect(fill = "NA"),
              panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
        geom_vline(xintercept =0,color="gray",linetype="dashed")+
        ylab("")+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -7),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 7),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -7),size=3.5,vjust=-2.8,
                  label=joint_effect_data[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 7),size=3.5,vjust=-2.8,
                  label=joint_effect_data[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
        geom_point(data=joint_effect_data[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="black")+
        
        scale_y_discrete(breaks=unique(joint_effect_data[[i]][[m]]$biome),position="right",
                         labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
        #xlab(full_name[i])+
        xlab("")+
        ggtitle(scenario[m])+
        xlim(xlim_list_scenario[[m]])
    }
  }
  else{
    pp=list()
    for (m in 1:2){
      
      pp[[m]]=ggplot(data=joint_effect_data[[i]][[m]],aes(fill=change,y=biome ,x=100*origin_mean))+
        geom_col(width = 0.3,color="black",size=0.2)+
        geom_segment(data=joint_effect_data[[i]][[m]], 
                     aes(y=biome,yend=biome,xend=  100*ori_low, x = 100*ori_up  ),
                     arrow = arrow(length = unit(0.1, "cm"),"both",angle=90))+
        
        scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#E27B60", "#1E61A5"))+
        theme(legend.position = "none",
              legend.text = element_text(size=8),
              legend.title  = element_text(size=10),
              text = element_text(size = 18),
              plot.title = element_text(size = 12, hjust = 0.5), 
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.title.x = element_text(size = 10), 
              legend.key.size = unit(0.3, "cm"),
              plot.margin = unit(c(0, 0, 0.3, 0), "cm"),
              panel.background = element_rect(fill = "NA"),
              panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
        geom_vline(xintercept =0,color="gray",linetype="dashed")+
        ylab("")+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -10),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 10),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -10),size=3.5,vjust=-2.8,
                  label=joint_effect_data[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
        geom_text(data=joint_effect_data[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 10),size=3.5,vjust=-2.8,
                  label=joint_effect_data[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
        geom_point(data=joint_effect_data[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=1,fill="black")+
        
        scale_y_discrete(breaks=unique(joint_effect_data[[i]][[m]]$biome),position="right",
                         labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",joint_effect_data[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
        #xlab(full_name[i])+
        xlab("")+
        ggtitle(scenario[m])+
        xlim(xlim_list[[i]])
    }
  }
  pp_guild[[i]]=pp 
}




pp_combine_effect_guild=list()
for (m in c(1,2,3,6,9)){
  if(m>1){
    pp_combine_effect=list()
    for (i in 1:2)
    {
      pp_combine_effect[[i]]=ggplotGrob(pp_guild[[m]][[i]]+
                                          annotation_custom(
                                            grob = pp_pie_guild[[m]][[i]],
                                            xmin = location_xmin[i], xmax =location_xmin[i], # Adjust x-axis position of the circle
                                            ymin = 0.21, ymax = 5.5)+
                                          ggtitle(scenario[i])+
                                          theme(
                                            legend.text = element_text(size=8),
                                            legend.title  = element_text(size=10),
                                            text = element_text(size = 18),
                                            plot.title = element_text(size = 15, hjust = 0.5), 
                                            axis.text.y = element_text(size = 12), 
                                            axis.text.x = element_text(size = 12), 
                                            axis.title.y = element_text(size = 15), 
                                            axis.title.x = element_text(size = 15), 
                                            plot.margin = unit(c(0.3, 0.1, -0.5, 0.1), "cm"),
                                            panel.background = element_rect(fill = "NA"),
                                            panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                          xlab(full_name[m]))
    }
  }
  else{
    pp_combine_effect=list()
    for (i in 1:2)
    {
      pp_combine_effect[[i]]=ggplotGrob(pp_guild[[m]][[i]]+
                                          annotation_custom(
                                            grob = pp_pie_guild[[m]][[i]],
                                            xmin = -38, xmax =-38, # Adjust x-axis position of the circle
                                            ymin = 0.21, ymax = 5.5)+
                                          ggtitle(scenario[i])+
                                          theme(
                                            legend.text = element_text(size=8),
                                            legend.title  = element_text(size=10),
                                            text = element_text(size = 18),
                                            plot.title = element_text(size = 15, hjust = 0.5), 
                                            axis.text.y = element_text(size = 12), 
                                            axis.text.x = element_text(size = 12), 
                                            axis.title.y = element_text(size = 15), 
                                            axis.title.x = element_text(size = 15), 
                                            plot.margin = unit(c(0.3, 0.1, -0.5, 0.1), "cm"),
                                            panel.background = element_rect(fill = "NA"),
                                            panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                          xlab(full_name[m]))
    }
  }
  
  pp_combine_effect_guild[[m]]=pp_combine_effect
}




p1=plot_grid(p_joint_rcp245,p_joint_latitude_245, pp_combine_effect_guild[[9]][[1]],ncol=3,rel_widths = c(1,0.35,0.5))
p2=plot_grid(p_joint_rcp585,p_joint_latitude_585, pp_combine_effect_guild[[9]][[2]],ncol=3,rel_widths = c(1,0.35,0.5))


plot_grid(p1,p2,ncol=1,
          label_x = 0.1,label_y = 1,label_size = 16)
        
## for the species loss and gains for each biome


joint_impact[[1]][[9]]->temp

joint_impact[[2]][[9]]->temp_585



temp%>%
  bind_cols(grid_level_biomes)%>%filter(!is.na(joint_impact))%>%
  mutate(type=ifelse(joint_impact>0,"gain","loss"))%>%
  group_by(LABEL)%>%count(type)->temp_data
  
#get the proportion for each biome

data=c("temp","temp_585")

bar_plot_joint_data=list()
for (k in 1:2)
{

change_type=matrix(ncol=2,nrow=4)
for (i in 1:4){
  get(data[k])%>%
    bind_cols(grid_level_biomes)%>%filter(!is.na(joint_impact))%>%
    mutate(type=ifelse(joint_impact>0,"gain","loss"))%>%
    group_by(LABEL)%>%count(type)%>%
    filter(LABEL==select_biome[i])%>%data.frame()->d
  d$n/sum(d$n)->d1
  change_type[i,]=d1
}

bar_plot_joint_data[[k]]=change_type%>%data.frame()%>%
  bind_cols(select_biome)%>%
  rename_all(~paste0(c("gain","loss","LABEL")))%>%melt()
}

## to see the proportions of species loss and gain among different biomes


biome_bar_plot_joint=list()
for (i in 1:2){
  
  data_set2= bar_plot_joint_data[[i]]
  
  data_set2$LABEL=factor(data_set2$LABEL, levels=c("Temperate Broadleaf & Mixed Forests",
                                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                                   "Tropical & Subtropical Dry Broadleaf Forests"))
  
  data_set2$BarWidth <- ifelse(duplicated(data_set2$LABEL) | duplicated(data_set2$LABEL, fromLast = TRUE), 0.8, 0.4)
  
  biome_bar_plot_joint[[i]]=ggplot(data=data_set2,aes(fill=variable,x=LABEL ,y=100*value))+
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = data_set2$BarWidth) +
    scale_fill_manual("",breaks=c("gain","loss"),labels=c("Gain","Loss"),values=c(  "#1E61A5","#E27B60"))+
    geom_text(aes(label = paste0("(", sprintf("%.1f", value * 100), ")")),
              position = position_dodge(width = 0.8), 
              vjust = -1.25, hjust=0.5,size = 4)+
    theme(legend.position = c(0.8,0.8),
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
  
    scale_x_discrete("",breaks=unique(data_set2$LABEL),
                     labels=c("Grassland\n biome", "Conifer forest\n biome","Broadleaf-mixed\n forest biome","Dry forest\n biome"))+
    ylab("% of pixels")+
    
    ylim(0,70)+
    ggtitle(scenario[i])
}

###

plot_grid(plotlist = biome_bar_plot_joint,
          ncol = 2,
          label_x = 0.1,label_y = 1,label_size = 16,
          labels = paste0("(", letters[1:2], ")") )


