# the updated figure 1
# get the species composition among different land use types
# this could be done based on different biomes

sample_data(rare_all_guild_biome)%>%data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%left_join(df4,by="Site")%>%filter(!is.na(plotIDM))%>%head(45)->temp_data


biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

biome_group=list()
for (i in 1:4){
  df4%>%left_join(temp_data%>%dplyr::select(plotIDM,LABEL),by="plotIDM")%>%mutate(code=1:45)%>%
    group_by(LABEL)%>%
    filter(LABEL==biomes[i])%>%pull(code)->biome_group[[i]]
}

species_com_guild=readRDS("species_com_guild.rds")
 
# the above codes should be deleted

# based on the biomes to look at the difference in species among land use types
# all the data in the first biomes
ordination_data_guild=list()
for (m in 1:9){
  
  ordination_data=list()
  for(i in 1:4)
    {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
    species_com_guild[[m]][biome_group[[i]]]->phylo_list
    
    merged_phyloseq <- Reduce(function(x, y) merge_phyloseq(x, y), phylo_list)
    # select the samples with richness
    phylo_temp <- prune_samples(sample_sums(merged_phyloseq ) > 0, merged_phyloseq )
    phylo_temp%>%sample_data()%>%data.frame()%>%dplyr::select(type)->land_cover_type
    vege_com=otu_table(phylo_temp)%>%data.frame()%>%bind_cols(land_cover_type)
    ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
    ordination_data[[i]]=ordination$points%>%data.frame()%>%mutate(land_cover_type)
    
  }
  
  
  ordination_data_guild[[m]]=ordination_data
}


saveRDS(ordination_data_guild,file="ordination_data_guild.rds")

ordination_data_guild=readRDS("ordination_data_guild.rds")

#compare the richness among covers

compare_richness_guild=readRDS("compare_richness_guild.rds")




type=c("cultivatedCrops", "woodyWetlands",   "evergreenForest" ,"deciduousForest",
 "mixedForest","grasslandHerbaceous")
type%>%data.frame()%>%mutate(color=c("#c94e65","royalblue","forestgreen", "#037f77","chocolate1","#7c1a97"))%>%
  rename_all(~paste0(c("type","color")))->color_match

ordination_data_guild=readRDS("ordination_data_guild.rds")

pp_spcom=list()
for(i in 1:4)
  {
  ordination_data_guild[[1]][[i]]%>%distinct(type)%>%left_join(color_match,by="type")->color_select
   if(i==1)
     {
     
     pp_spcom[[i]]=ggplotGrob(ggplot(data=ordination_data_guild[[1]][[i]],aes(x=MDS1,  y= MDS2,color=type ))+
                          geom_point(data=ordination_data_guild[[1]][[i]],pch=21,color="black",aes(x=MDS1,y= MDS2,fill=type ),size=3,alpha=0.75)+
                          stat_ellipse(size=0.8,linetype="dashed")+
                          theme(legend.position = c(0.3,0.4), 
                                legend.title = element_text(size=10),
                                #text = element_text(size = 18), 
                                legend.text = element_text(size=11),
                                plot.title = element_text(size = 15, hjust = 0.5), 
                                axis.text.y = element_text(hjust = 1,size=12), 
                                axis.text.x = element_text(hjust = 0.5,size=12), 
                                axis.title.y = element_text(size = 15), 
                                axis.title.x = element_text(size = 15),
                                panel.background = element_rect(fill = "NA"), 
                                panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
                          ylim(-0.5,1)+
                          scale_fill_manual("",breaks = color_select$type,values = color_select$color)+
                          scale_color_manual("",breaks = color_select$type,values = color_select$color))
   }
  else{
    
    pp_spcom[[i]]=ggplotGrob(ggplot(data=ordination_data_guild[[1]][[i]],aes(x=MDS1,  y= MDS2,color=type ))+
    geom_point(data=ordination_data_guild[[1]][[i]],pch=21,color="black",aes(x=MDS1,y= MDS2,fill=type ),size=3,alpha=0.75)+
    stat_ellipse(size=0.8,linetype="dashed")+
      theme(legend.position = c(0.3,0.4), 
            legend.title = element_text(size=10),
            #text = element_text(size = 18), 
            legend.text = element_text(size=11),
            plot.title = element_text(size = 15, hjust = 0.5), 
            axis.text.y = element_text(hjust = 1,size=12), 
            axis.text.x = element_text(hjust = 0.5,size=12), 
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15),
            panel.background = element_rect(fill = "NA"), 
            panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
      scale_fill_manual("",breaks = color_select$type,values = color_select$color)+
    scale_color_manual("",breaks = color_select$type,values = color_select$color))
  }
}


plot_grid(pp_spcom[[1]],pp_spcom[[2]],pp_spcom[[3]],pp_spcom[[4]],ncol=2)

# add the pairwise comparison of the richness
##the function to merge the data
data_mean_richness_biome_model_guild=list()
for(m in 1:9)
{
  do.call(rbind,richness_ratio_with_rarefaction_guild[[m]])%>%data.frame()%>%
    mutate(plotid=unlist(biome_group))%>%mutate(biomes=rep(biome_select,times=plot_number))%>%
    left_join(temp,by="plotid")%>%rename_all(~paste0(c("low","mean_nature","up","mean_crop","plotid","biome","site","plotIDM")))->df1
  
  df1%>%dplyr::select(mean_nature,biome, site,plotIDM)%>%
    rename(mean_rich=mean_nature)%>%
    rbind(df1%>%dplyr::select(mean_crop,biome, site,plotIDM)%>%rename(mean_rich=mean_crop))%>%
    mutate(type=rep(c("nature","crop"),each=45))->data_mean_richness_biome_model_guild[[m]] 
}

saveRDS(data_mean_richness_biome_model_guild,file="data_mean_richness_biome_model_guild.rds")


data_mean_richness_biome_model_guild=readRDS("data_mean_richness_biome_model_guild.rds")

biome_select=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
               "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")



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
    ylab("Taxonomic diversity")+
      xlab("Land use")+
    scale_x_discrete(labels=c("Natural","Modified"))+
      geom_segment(data=all_data, aes(x = type[1], y = all_mean[1],xend=type[2],yend=all_mean[2]),color="#E3A72F",linetype="dashed",size=1)+
    geom_errorbar(data=all_data,aes(x=type,ymin = all_mean-all_sd, ymax = all_mean+all_sd), width = 0.061)+
    geom_point(data=all_data,pch=24,aes(x=type,y=all_mean),color="black",fill=c("#32CD92","#CD326D"),size=5)+
  
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
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
    ylim(0,2600))
}




plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],ncol=2,label_size = 8)


# the response ratio of different guilds

biome_site_level_richness_ratio_consider_nature_history=readRDS("biome_site_level_richness_ratio_consider_nature_history.rds")


df_significance=readRDS("df_significance.rds")

data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")

guild=c( "all", "AM" ,"EM","plapat","soilsap","littersap", "woodsap","epiphy", "para" )

pp_response=list()
for (i in 1:4)
{
  data=biome_site_level_richness_ratio_consider_nature_history%>%filter( biome==biome_select[i])
  data$guild=factor(data$guild,levels=rev(c("all","AM" ,"EM","plapat","soilsap","littersap","woodsap","epiphy","para")))
  
  pp_response[[i]]=ggplotGrob(ggplot(data, aes(x=guild,y = mean_ratio)) +
    geom_bar(stat = "identity",position = "dodge",width =0.5) +
    guides(fill="none")+
    geom_errorbar(aes(ymin = mean_ratio-sd_ratio, ymax = mean_ratio+sd_ratio), width = 0.1)+
    geom_hline(yintercept = 1,color="red",linetype="dashed")+
    scale_x_discrete (breaks=guild,position = "top",labels = c("All", "AM","EM","Pla. patho.", "Soil sapro.","Litter sapro.","Wood sapro.","Epiphyte","Parasite"))+
    #ggtitle("Response ratio")+
    xlab("")+
    ylab("Response ratio")+
    theme(legend.position = c(0.8,0.75),
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(size=12), 
          axis.text.x = element_text(size=12), 
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0.3, 0.5, 0, 0.1), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
    geom_text(aes(x=guild,y=(mean_ratio+sd_ratio)*1.1),label =df_significance[[i]],size=5)+
    coord_flip()+
    ylim(0,8))
  
}

plot_grid(pp_response[[1]],pp_response[[2]],pp_response[[3]],pp_response[[4]],ncol=1)

#if we focused on the four main fungal guilds

guild=c( "all", "AM" ,"EM","plapat","soilsap" )



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
                                      axis.text.y = element_text(size=12), 
                                      axis.text.x = element_text(size=12), 
                                      axis.title.y = element_text(size = 15), 
                                      axis.title.x = element_text(size = 15), 
                                      legend.key.size = unit(0.3, "cm"),
                                      plot.margin = unit(c(0.3, 0.5, 0, 0.1), "cm"),
                                      panel.background = element_rect(fill = "NA"),
                                      panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
                                geom_text(aes(x=guild,y=(mean_ratio+sd_ratio)*1.1),label =df_significance[c(1:5),][[i]],size=5)+
                                coord_flip()+
                                ylim(0,8))
  
}




##we focus on five land use types


 # the response ratio quantified
 
 df2$guild=factor(df2$guild,levels=rev(c("all","AM" ,"EM","epiphy","littersap","para","plapat","soilsap","woodsap")))

 biomes=c("Temperate Broadleaf & Mixed Forests",
          "Temperate Conifer Forests",
          "Temperate Grasslands, Savannas & Shrublands",
          "Tropical & Subtropical Moist Broadleaf Forests")
 
 #pp_response=list()
 for (i in 1:4)
   {
   pp_response[[i]]=ggplotGrob(ggplot(df2%>%filter(LABEL==biomes[i]), aes(x = guild, y = mean_value)) +
     geom_bar(stat = "identity", fill = "black",width=0.5)+
     geom_errorbar(aes(ymin = mean_value-sd, ymax = mean_value+sd), width = 0.1)+
     theme(legend.position = c(0.25,0.2), 
           legend.title = element_text(size=10),
           #plot.margin = unit(c(1,1, 0, 1), "cm"),
           #text = element_text(size = 18), 
           legend.text = element_text(size=11),
           plot.title = element_text(size = 12, hjust = 0.5), 
           axis.text.y = element_text(hjust = 1,size=12), 
           axis.text.x = element_text(hjust = 0.5), 
           axis.title.y = element_text(size = 18), 
           axis.title.x = element_text(size = 18),
           panel.background = element_rect(fill = "NA"), 
           panel.border = element_rect(color = "black", size = 0.7, fill = NA))+
     ylab("Response ratio")+
     geom_hline(yintercept = 1,color="red",linetype="dashed")+
     coord_flip()+
     ylim(0,6)+
     xlab("")+
     scale_x_discrete(breaks=c("all","AM" ,"EM","epiphy","littersap","para","plapat","soilsap","woodsap"), position="bottom",labels = c("All", "AM","EM","Epiphyte","Litter sapro.","Parasite","Pla. patho.","Soil sapro.","Wood sapro.")))
 }
 
 
 pp_spcom1=ggplotGrob(pp_spcom1)
 pp_compare1=ggplotGrob(pp_compare1)
 pp_response1=ggplotGrob(pp_response1)
 
 pp_spcom1$heights=pp_compare1$heights
 
 pp_spcom1$heights=pp_response1$heights
 
 pp_compare1$heights=pp_response1$heights
 
 plot_grid( pp_compare1,pp_response1)
 
 plot_grid( pp_compare1,pp_response1,ncol=3)
 
 
 
 plot_grid(pp_spcom1,pp_compare1,pp_response1,ncol=3)
 
 
 
 
 pp_spcom[[1]]$heights=pp_compare[[1]]$heights
 pp_spcom[[1]]$heights=pp_response[[1]]$heights
 pp_compare[[1]]$heights=pp_response[[1]]$heights
 
 pp_spcom[[2]]$heights=pp_compare[[2]]$heights
 pp_spcom[[2]]$heights=pp_response[[2]]$heights
 pp_compare[[2]]$heights=pp_response[[2]]$heights
 
 pp_spcom[[3]]$heights=pp_compare[[3]]$heights
 pp_spcom[[3]]$heights=pp_response[[3]]$heights
 pp_compare[[3]]$heights=pp_response[[3]]$heights
 
 pp_spcom[[4]]$heights=pp_compare[[4]]$heights
 pp_spcom[[4]]$heights=pp_response[[4]]$heights
 pp_compare[[4]]$heights=pp_response[[4]]$heights
 
 
 pp_spcom[[2]]$widths=pp_spcom[[3]]$widths
 
 pp_spcom[[3]]$widths=pp_spcom[[4]]$widths
 pp_spcom[[1]]$widths=pp_spcom[[2]]$widths
 
 
 plot_grid(pp_spcom[[1]],pp_compare[[1]],pp_response[[1]],
           pp_spcom[[2]],pp_compare[[2]],pp_response[[2]],
           pp_spcom[[3]],pp_compare[[3]],pp_response[[3]],
           pp_spcom[[4]],pp_compare[[4]],pp_response[[4]],ncol=3,rel_widths = c(0.8,0.6,1))
 #dimension 10 by 10
 