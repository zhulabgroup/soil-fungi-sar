
setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")
rm(list = ls())

data_mean_richness_biome_model_guild=readRDS("data_mean_richness_biome_model_guild.rds")
# the diversity ratio for different fungal guilds
biome_site_level_richness_ratio_consider_nature_history=readRDS("biome_site_level_richness_ratio_consider_nature_history.rds")
df_significance=readRDS("df_significance.rds")
data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")

#if we focus on the four main fungal guilds

guild=c( "all", "AM" ,"EM","plapat","soilsap" )


biome_select=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
               "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

#set different limits for each panel
set_limits=list(c(200,2800),
                c(500,1600),
                c(20,1500),
                c(300,1000))

pp_compare=list()
for(i in 1:4)
{
  d=data_mean_richness_biome_model_guild[[1]]%>%filter(biome==biome_select[i])
  d$type=factor(d$type,levels = c("nature","crop"))
  d%>%group_by(type)%>%summarise(all_mean=mean(mean_rich),all_sd=sd(mean_rich))%>%data.frame()->all_data
  pp_compare[[i]]= ggplotGrob(ggplot()+
                                geom_point(data=d, pch=21,aes(x = type,y = mean_rich, group = plotIDM),color="gray",fill=rep(c("#32CD92","#CD326D"),each=dim(d)[1]/2),size=4,alpha=0.4) +
                                geom_line(data=d, aes(x = type, y = mean_rich, group = plotIDM),color=rep(c("gray","gray"),each=dim(d)[1]/2),alpha=0.5) +
                                guides(color="none")+
                                ylab("Diversity")+
                                xlab("Land use")+
                                scale_x_discrete(labels=c("Natural","Modified"))+
                                geom_segment(data=all_data, aes(x = type[1], y = all_mean[1],xend=type[2],yend=all_mean[2]),color="#E3A72F",linetype="dashed",size=1)+
                                geom_errorbar(data=all_data,aes(x=type,ymin = all_mean-all_sd, ymax = all_mean+all_sd), width = 0.061)+
                                geom_point(data=all_data,pch=24,aes(x=type,y=all_mean),color="black",fill=c("#32CD92","#CD326D"),size=5)+
                                
                                theme(legend.position = c(0.75,0.28),
                                      legend.text = element_text(size=8),
                                      legend.title  = element_text(size=10),
                                      #text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size=10), 
                                      axis.text.x = element_text(size = 10), 
                                      axis.title.y = element_text(size = 12), 
                                      axis.title.x = element_text(size = 12), 
                                      plot.margin = unit(c(0.3, 0.1, -.5, 0), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
                                ylim(set_limits[[i]]))
}



pp_response=list()
for (i in 1:4)
{
  data=biome_site_level_richness_ratio_consider_nature_history%>%filter( biome==biome_select[i])%>%filter(guild%in%c("all","AM","EM","soilsap","plapat"))
  data$guild=factor(data$guild,levels=rev(c("all","AM" ,"EM","plapat","soilsap")))
  
  pp_response[[i]]=ggplotGrob(ggplot(data, aes(x=guild,y = mean_ratio)) +
                                geom_bar(stat = "identity",position = "dodge",width =0.6) +
                                guides(fill="none")+
                                geom_errorbar(aes(ymin = mean_ratio-sd_ratio, ymax = mean_ratio+sd_ratio), width = 0.1)+
                                geom_hline(yintercept = 1,color="red",linetype="dashed")+
                                scale_x_discrete (breaks=guild,position = "top",labels = c("All", "AM","EM","Plant pathogens","Soil saprotrophs"))+
                                #ggtitle("Diversity ratio")+
                                xlab("")+
                                ylab("Diversity ratio")+
                                theme(legend.position = c(0.8,0.75),
                                      legend.text = element_text(size=8),
                                      legend.title  = element_text(size=10),
                                      text = element_text(size = 18),
                                      plot.title = element_text(size = 15, hjust = 0.5), 
                                      axis.text.y = element_text(size=10,margin = margin(r=-30),hjust = -3), 
                                      axis.text.x = element_text(size=10), 
                                      axis.title.y = element_text(size = 12), 
                                      axis.title.x = element_text(size = 12), 
                                      legend.key.size = unit(0.3, "cm"),
                                      plot.margin = unit(c(0.3, 0.5, 0, 0.1), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
                                geom_text(aes(x=guild,y=(mean_ratio+sd_ratio)*1.1),label =df_significance[c(1:5),][[i]],size=5)+
                                coord_flip()+
                                ylim(0,8))
  
}


pp_compare[[1]]$heights=pp_response[[1]]$heights
pp_compare[[2]]$heights=pp_response[[2]]$heights
pp_compare[[3]]$heights=pp_response[[3]]$heights
pp_compare[[4]]$heights=pp_response[[4]]$heights
#dimension 10 by 10


plot_grid(pp_compare[[1]],pp_response[[1]],
          pp_compare[[2]],pp_response[[2]],
          pp_compare[[3]],pp_response[[3]],
          pp_compare[[4]],pp_response[[4]],ncol=2,rel_widths = c(0.6,0.9))



