
setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-climate-historical")


climate_induced_change_richness_rcp245=readRDS(file="climate_induced_change_richness_rcp245.rds")
climate_induced_change_richness_rcp585=readRDS(file="climate_induced_change_richness_rcp585.rds")
summary_data_climate_rcp245=readRDS(file="summary_data_climate_rcp245.rds")
summary_data=readRDS(file="summary_data.rds")

species_change_land_rcp245=readRDS(file="species_change_land_rcp245.rds")
species_change_land_rcp585=readRDS(file="species_change_land_rcp585.rds")

species_change_land_rcp245%>%filter(variable=="all")->species_change_land_rcp245_all

load(file="coords_present_new.RData")


species_change_land_rcp245_all%>%bind_cols(coords_present)->temp

species_change_land_rcp245_all=my_function_project(temp)

##
species_change_land_rcp585%>%filter(variable=="all")->species_change_land_rcp585_all

load(file="coords_present_new.RData")


species_change_land_rcp585_all%>%bind_cols(coords_present)->temp

species_change_land_rcp585_all=my_function_project(temp)


save(species_change_land_rcp585_all,file="species_change_land_rcp585_all.rds")

save(species_change_land_rcp245_all,file="species_change_land_rcp245_all.rds")



summary_data_climate_rcp585=readRDS(file="summary_data_climate_rcp585.rds")
summary_data_climate_rcp245=readRDS("summary_data_climate_rcp245.rds")
summary_data_land_rcp585=readRDS("summary_data_land_rcp585.rds")



#add a new column to transforme the loss and gain value
#need to update the package
climate_induced_change_richness_rcp245%>%mutate(change=100*last)->climate_induced_change_richness_rcp245

p_climate_245=ggplot(climate_induced_change_richness_rcp245) +
  geom_tile(data =climate_induced_change_richness_rcp245 ,aes(x = x, y = y, fill = change), size = 0.175)+
  scale_fill_gradientn(NULL,
                       limits = c(-70, 70),
                       colors = brewer.pal(11, "RdBu"),
                       guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(climate_induced_change_richness_rcp245, change < -70), mapping = aes(x=x,y=y,fill = last > 0.7)) +
  scale_fill_manual("Change (%)", values = "navy", labels = "> 70", 
                    guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
  ) +
  new_scale_fill() +
  geom_tile(data = filter(climate_induced_change_richness_rcp245, change < -70), mapping = aes(x=x,y=y,fill = last< -0.7)) +
  scale_fill_manual(NULL, values = "gold", labels = "< -70", 
                    guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
  theme_map+ 
  add_sf_layers()+
  coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
  xlab("")+
  ylab("Climate impact")+
  ggtitle("SSP2−4.5")


p_climate_latitude_245=ggplotGrob(ggplot()+
                                    geom_ribbon(data = summary_data_climate_rcp245, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "gray80")+
                                    geom_line(data = summary_data_climate_rcp245, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) +  # Mean trend line
                                    theme_latitude+
                                    xlab("Latitude")+
                                    ylab("")+
                                    ggtitle("")+
                                    coord_flip()+
                                    geom_hline(yintercept = 0,linetype="dashed")+
                                    ggtitle("SSP2−4.5")+
                                    xlim(15,65)+
                                    ylim(-40,40))



climate_induced_change_richness_rcp585%>%mutate(change=100*last)->climate_induced_change_richness_rcp585

p_climate_585=ggplotGrob(ggplot(climate_induced_change_richness_rcp585) +
                           geom_tile(data =climate_induced_change_richness_rcp585 ,aes(x = x, y = y, fill = change), size = 0.175)+
                           scale_fill_gradientn(NULL,
                                                limits = c(-70, 70),
                                                colors = brewer.pal(11, "RdBu"),
                                                guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
                           ggnewscale::new_scale_fill() +
                           geom_tile(data = filter(climate_induced_change_richness_rcp585, change > 70), mapping = aes(x=x,y=y,fill = last > 0.7)) +
                           scale_fill_manual("Change (%)", values = "navy", labels = "> 70", 
                                             guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
                           ) +
                           ggnewscale::new_scale_fill() +
                           geom_tile(data = filter(climate_induced_change_richness_rcp585, change < -70), mapping = aes(x=x,y=y,fill = last< -0.7)) +
                           scale_fill_manual(NULL, values = "gold", labels = "< -70", 
                                             guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
                           theme_map+
                           add_sf_layers()+
                           coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                           xlab("")+
                           ylab("Climate impact")+
                           ggtitle("SSP5−8.5"))



p_climate_latitude_585=ggplotGrob(ggplot()+
                                    geom_ribbon(data = summary_data_climate_rcp585, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "gray80")+
                                    geom_line(data = summary_data_climate_rcp585, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) +  # Mean trend line
                                    theme_latitude+
                                    xlab("Latitude")+
                                    ylab("")+
                                    ggtitle("SSP5−8.5")+
                                    coord_flip()+
                                    geom_hline(yintercept = 0,linetype="dashed")+
                                    xlim(15,65)+
                                    ylim(-60,60))

species_change_land_rcp245_all%>%mutate(change=last*100)->species_change_land_rcp245_all


p_land_245=ggplotGrob(ggplot(species_change_land_rcp245_all) +
                        geom_tile(data =species_change_land_rcp245_all ,aes(x = x, y = y, fill = change), size = 0.175)+
                        scale_fill_gradientn(NULL,
                                             limits = c(-15, 15),
                                             colors = brewer.pal(11, "RdBu"),
                                             
                                             guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
                        ggnewscale::new_scale_fill() +
                        geom_tile(data = species_change_land_rcp245_all%>%filter(change > 5), mapping = aes(x=x,y=y,fill = last > 0.7)) +
                        scale_fill_manual("Change (%)", values = "navy", labels = "> 15", 
                                          guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
                        ) +
                        
                        ggnewscale::new_scale_fill() +
                        geom_tile(data = species_change_land_rcp245_all%>%filter(change < -15), mapping = aes(x=x,y=y,fill = last< -0.7)) +
                        scale_fill_manual(NULL, values = "gold", labels = "< -15", 
                                          guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
                        
                        theme_map+
                        add_sf_layers()+
                        coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                        
                        xlab("")+
                        ylab("Land-use impact")+
                        ggtitle("SSP2−4.5"))


p_land_latitude_245=ggplotGrob(ggplot()+
                                 # Mean trend line
                                 geom_ribbon(data = summary_data, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value +100*sd_value), fill = "gray80")+
                                 geom_line(data = summary_data, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) + 
                                 theme_latitude+
                                 xlab("Latitude")+
                                 ylab("")+
                                 ggtitle("SSP2−4.5")+
                                 coord_flip()+
                                 geom_hline(yintercept = 0,linetype="dashed")+
                                 xlim(15,65)+
                                 ylim(-20,20))

species_change_land_rcp585_all%>%mutate(change=100*last)->species_change_land_rcp585_all


p_land_585=ggplotGrob(ggplot(species_change_land_rcp585_all) +
                        geom_tile(data =species_change_land_rcp585_all ,aes(x = x, y = y, fill = change), size = 0.175)+
                        scale_fill_gradientn(NULL,
                                             limits = c(-15, 15),
                                             colors = brewer.pal(11, "RdBu"),
                                             guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+
                        ggnewscale::new_scale_fill() +
                        geom_tile(data = filter(species_change_land_rcp585_all, change > 15), mapping = aes(x=x,y=y,fill = last > 0.7)) +
                        scale_fill_manual("Change (%)", values = "navy", labels = "> 15", 
                                          guide = guide_legend(order = 1,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))
                        ) +
                        
                        ggnewscale::new_scale_fill() +
                        geom_tile(data = filter(species_change_land_rcp585_all, change < -15), mapping = aes(x=x,y=y,fill = last< -0.7)) +
                        scale_fill_manual(NULL, values = "gold", labels = "< -15", 
                                          guide = guide_legend(order = 3,keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))+
                        theme_map+
                        add_sf_layers()+
                        coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
                        
                        xlab("")+
                        ylab("Land-use impact")+
                        ggtitle("SSP5−8.5"))


p_land_latitude_585=ggplotGrob(ggplot()+
                                 # Mean trend line
                                 geom_ribbon(data = summary_data_land_rcp585, aes(x = lat, ymin = 100*mean_value - 100*sd_value, ymax = 100*mean_value + 100*sd_value), fill = "gray80")+
                                 geom_line(data = summary_data_land_rcp585, aes(x = lat, y = 100*mean_value), color = "#1E61A5", size = 0.5) +  
                                 theme_latitude+
                                 xlab("Latitude")+
                                 ylab("")+
                                 ggtitle("SSP5−8.5")+
                                 coord_flip()+
                                 geom_hline(yintercept = 0,linetype="dashed")+
                                 ylim(-20,20)+
                                 xlim(15,65))

# to add the pie_plot and the mean diversity loss and gain rate



pie_data_guild=readRDS("pie_data_guild.rds")

com_data_climate_guild=readRDS("com_data_climate_guild.rds")

# for different guilds
#set the space for different plot

scenario=c("SSP2-4.5","SSP5-8.5","SSP2-4.5","SSP5-8.5")

full_name=c("AM","EM","Soil saprotrophs","Litter saprotrophs","Wood saprotrophs", "Plant pathogens", "Parasite","Epiphyte","All")


# for different guilds
#set the space for different plot

xlim_list=list(c(-40,40),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20))

# set the limit for m

xlim_list_scenario=list(c(-20,20),c(-20,20),c(-25,25),c(-25,25))


#when m=9, the plots is about the whole fungal community
# when i=1, the plot shows the effect of land use change
pp_guild=list()
for(i in c(1,2,3,6,9))
  
{
  
  if(i>1)
  {
    pp=list()
    for (m in 1:4){
      
      pp[[m]]=ggplot(data=com_data_climate_guild[[i]][[m]],aes(fill=change,y=biome ,x=100*origin_mean))+
        geom_col(width = 0.3,color="black",size=0.2)+
        geom_segment(data=com_data_climate_guild[[i]][[m]], 
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
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -7),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 7),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -7),size=3.5,vjust=-2.8,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 7),size=3.5,vjust=-2.8,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
        geom_point(data=com_data_climate_guild[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="black")+
        scale_y_discrete(breaks=unique(com_data_climate_guild[[i]][[m]]$biome),position="right",
                         labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
        #xlab(full_name[i])+
        xlab("")+
        ggtitle(scenario[m])+
        xlim(xlim_list_scenario[[m]])
    }
  }
  else{
    pp=list()
    for (m in 1:4){
      
      
      pp[[m]]=ggplot(data=com_data_climate_guild[[i]][[m]],aes(fill=change,y=biome ,x=100*origin_mean))+
        geom_col(width = 0.3,color="black",size=0.2)+
        geom_segment(data=com_data_climate_guild[[i]][[m]], 
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
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -10),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 10),size=3.5,vjust=-1.4,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -10),size=3.5,vjust=-2.8,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 10),size=3.5,vjust=-2.8,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
        geom_point(data=com_data_climate_guild[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=1,fill="black")+
        scale_y_discrete(breaks=unique(com_data_climate_guild[[i]][[m]]$biome),position="right",
                         labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
        #xlab(full_name[i])+
        xlab("")+
        ggtitle(scenario[m])+
        xlim(xlim_list[[i]])
    }
  }
  pp_guild[[i]]=pp 
}

## ad the pie plot

pp_pie_guild=list()
for (m in 1:9){
  if(m==1){
    pp_pie=list()
    for (i in 1:4){
      pp_pie[[i]]=ggplotGrob(ggplot(pie_data_guild[[m]][[i]], aes(x=1, y=percent, fill=variable)) +
                               ## geom_col is geom_bar(stat = "identity")(bit shorter)
                               ## use color = black for the outline
                               geom_col(width=5,position="fill", color = "black",size=0.25)+
                               coord_polar("y", start=0) +
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
    for (i in 1:4){
      pp_pie[[i]]=ggplotGrob(ggplot(pie_data_guild[[m]][[i]], aes(x=1, y=percent, fill=variable)) +
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

location_xmin=c(-18.7,-18.7,-23.5,-23.5)
location_xmax=c(-18.7,-18.7,-23.5,-23.5)

#for the whole fungal community and individual guilds, need to set different plot margin.

pp_combine_effect_guild=list()
for (m in c(1,2,3,6,9)){
  if(m>1){
    pp_combine_effect=list()
    for (i in 1:4)
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
    for (i in 1:4)
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

#bind all the plots

# do not set the height of the map,
# it is important to check the plot margin of the plots for better alignment

#p_land_245$heights=pp_combine_effect[[1]]$heights

p_land_latitude_245$heights=pp_combine_effect[[1]]$heights
p_land_latitude_585$heights=pp_combine_effect[[2]]$heights
p_climate_latitude_585$heights=pp_combine_effect[[4]]$heights
p_climate_latitude_245$heights=pp_combine_effect[[3]]$heights
pp_combine_effect[[3]]$widths=pp_combine_effect[[1]]$widths

#p_climate_245$widths=pp_combine_effect[[1]]$widths


p1=plot_grid(p_climate_245,p_climate_latitude_245, pp_combine_effect[[3]],ncol=3,rel_widths = c(1,0.35,0.5))
p2=plot_grid(p_land_245,p_land_latitude_245, pp_combine_effect[[1]],ncol=3,rel_widths = c(1,0.35,0.5))
p3=plot_grid(p_climate_585,p_climate_latitude_585, pp_combine_effect[[4]],ncol=3,rel_widths = c(1,0.35,0.5))

p4=plot_grid(p_land_585,p_land_latitude_585, pp_combine_effect[[2]],ncol=3,rel_widths = c(1,0.35,0.5))
plot_grid(p1,p2,p3,p4,ncol=1)
#dimension is 11x12.5



