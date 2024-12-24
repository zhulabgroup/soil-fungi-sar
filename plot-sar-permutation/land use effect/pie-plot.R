##get species change rate under two scenarios
##both climate and land cover were considered with four data

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")
species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")


data=c("species_change_land_rcp245","species_change_land_rcp585","species_change_climate_rcp245","species_change_climate_rcp585")

scenario=c("SSP2-4.5","SSP5-8.5","SSP2-4.5","SSP5-8.5")

#(1) get the proportions of pixels that show either species gain or loss

change_rate_guild_scenarios=list()# the j corresponds to four scenarios for both land cover and climate effect
for (j in 1:4)
  {
change_rate_guild=list()
for (m in 1:9)
{
  change_rate=matrix(ncol=3,nrow=5)
  
  for (i in 1:5)
  {
    if(i==5)#all biomes were included
    {
      get(data[j])%>%filter(variable==guild_type[m])%>%
        bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
        filter(LABEL%in%biomes_four&!is.na(value))->df1
      
      df1%>%filter(value>0)%>%dim()%>%head(1)->gain
      df1%>%filter(value==0)%>%dim()%>%head(1)->no_change
      df1%>%filter(value<0)%>%dim()%>%head(1)->loss
    }
    else{# select just one biome
      get(data[j])%>%filter(variable==guild_type[m])%>%
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
change_rate_guild_scenarios[[j]]=change_rate_guild

}

change_rate_guild_scenarios[[1]]#means land cover effect in the low-emission scenario for all guilds

saveRDS(change_rate_guild_scenarios,file="change_rate_guild_scenarios.rds")

#for each guild, to get the data to create the pie plot
pie_data_guild=list()
for (m in 1:9){
  pie_data=list()
  for (i in 1:4)# for both climate and land cover
  {
    change_rate_guild_scenarios[[i]][[m]]%>%
      data.frame()%>%rename_all(~paste0(c("gain","no","loss")))%>%
      mutate(biome=c("Temperate","Conifer","Grassland","Dry","All"))%>%
      melt()%>%group_by(biome)%>%mutate(percent = value / sum(value))->temp_pie_data
    temp_pie_data$biome=factor( levels=c("All" ,"Temperate", "Conifer" ,  "Grassland", "Dry" ),temp_pie_data$biome)
    
    pie_data[[i]]=temp_pie_data
  }
  pie_data_guild[[m]] =pie_data
}
#pie_data_guild[[1]], means the first fungal guild
#with four list, each correspond to land and cliamte effects under two scenarios

saveRDS(pie_data_guild,file="pie_data_guild.rds")

#to get the overall effect of both factors
# it is based on the species change rate

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

biomes_four=c( "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests" ,                    
               "Temperate Grasslands, Savannas & Shrublands" ,"Tropical & Subtropical Dry Broadleaf Forests")


change_rate_mean_guild_scenario=list()
for (j in 1:4)
  {
change_rate_mean_guild=list()
for (m in 1:9)
{
change_rate=list()
for (i in 1:5)
{
  if(i ==5)
  {
    get(data[j])%>%filter(variable==guild_type[m])%>%
      bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
      filter(LABEL%in%biomes_four&!is.na(value))->df1# na cells were filtered
    
    df1%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
      summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->data_temp
    
    df1%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
      overall_change
    data_temp%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")%>%
      data.frame()->change_rate[[i]]
  }
  else
  {
    get(data[j])%>%filter(variable==guild_type[m])%>%
      bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
      filter(LABEL== biomes_four[i]&!is.na(value))->df1
    
    df1%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
      summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->data_temp
    
    df1%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
      overall_change
    data_temp%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")%>%
      data.frame()->change_rate[[i]]
    
  }
  
}
change_rate_mean_guild[[m]]=change_rate
}
change_rate_mean_guild_scenario[[j]]=change_rate_mean_guild
}

change_rate_mean_guild_scenario[[1]]#means for land cover effects in the low-emission scenarios

#to combine the list into a dataframe to create the plots
#for the two scenarios to get the replicates for the mean change
# we only select the 9 th guild when all the guilds were combined
data_biome_change_rate_guild=list()
for (m in 1:9)
  {
  data_biome_change_rate=list()
  for (i in 1:4)#i corresponds to the four data 2 scenarios x 2 effects(land and climate)
  {
    dimensions_matrix <- sapply(change_rate_mean_guild_scenario[[i]][[m]], dim)[1,]
    # the number of rows selected 
    do.call(rbind,change_rate_mean_guild_scenario[[i]][[m]])%>%
      mutate(biome=rep(c(biomes_four,"all"),time=dimensions_matrix))->data_mean_change
    data_mean_change$biome=factor(levels = rev(c("all",biomes_four)),data_mean_change$biome)
    data_biome_change_rate[[i]]=data_mean_change
  }
  data_biome_change_rate_guild[[m]]=data_biome_change_rate
}

saveRDS(data_biome_change_rate_guild,file="data_biome_change_rate_guild.rds")
#each list has 9 lists with found data sets


pie_data_guild=readRDS("pie_data_guild.rds")

## to create the pie plots
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
                                                 labels=c("Gain","No change","Loss"),values=c( "#1b9e77","gray","#d95f02")))}
  }
  else{
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
                                                 labels=c("Gain","No change","Loss"),values=c("#1b9e77","gray","#d95f02")))}
  }
pp_pie_guild[[m]]=pp_pie
}

# load in the data
data_biome_change_rate_guild=readRDS("data_biome_change_rate_guild.rds")


full_name=c("AM","EM","Soil saprotroph","Litter saprotroph","Wood saprotroph", "Plant pathogen", "Parasite","Epiphyte","All")

#set different
xlim_list=list(c(-30,30),c(-30,30),c(-60,60),c(-60,60))

xlim_list_am=list(c(-100,100),c(-100,100),c(-60,60),c(-60,60))

#set different space for the text
x_space=c(0.27,0.27,5.7,5.7)

pp_mean_effect_guild=list()
for (m in 1:9)
  {
  if(m<2){
    pp_mean_effect=list()
    for (i in 1:4)
    {
      if(i<3)
      {
        
        pp_mean_effect[[i]]=ggplotGrob(ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
          geom_col(width = 0.3)+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                       aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                       aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#ee8c11","#1173EE"))+
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
                plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
          geom_vline(xintercept =0,color="gray",linetype="dashed")+
          ylab("")+
          xlab(paste(full_name[m]))+
          xlim(xlim_list_am[[i]])+
         
          scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                           labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
          geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
          #ggtitle(paste0(scenario[i]))+
          geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
          
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                    aes(y=biome,x = -50),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                    aes(y=biome,x = 50),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                      aes(y=biome,x = -50),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                      aes(y=biome,x = 50),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        )
      }
      else{
        pp_mean_effect[[i]]=ggplotGrob(ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
          geom_col(width = 0.3)+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                       aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                       aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#ee8c11","#1173EE"))+
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
                plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
          geom_vline(xintercept =0,color="gray",linetype="dashed")+
          ylab("")+
          xlab(paste(full_name[m]))+
          xlim(xlim_list_am[[i]])+
          
          scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                           labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
          geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
          #ggtitle(paste0(scenario[i]))+
          geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
          
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                    aes(y=biome,x = -30),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                    aes(y=biome,x = 30),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                      aes(y=biome,x = -30),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                      aes(y=biome,x = 30),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        )
        
      }
      
    }
  }
  
  else{
    pp_mean_effect=list()
    for (i in 1:4)
    {
      if(i<3)
      {
        
        pp_mean_effect[[i]]=ggplotGrob(ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
          geom_col(width = 0.3)+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                       aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                       aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#ee8c11","#1173EE"))+
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
                plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
          geom_vline(xintercept =0,color="gray",linetype="dashed")+
          ylab("")+
          xlab(paste(full_name[m]))+
          xlim(xlim_list[[i]])+
         
          scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                           labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
          geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
          #ggtitle(paste0(scenario[i]))+
          geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
          
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                    aes(y=biome,x = -10),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                    aes(y=biome,x = 10),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                      aes(y=biome,x = -10),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                      aes(y=biome,x = 10),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        )
      }
      else{
        pp_mean_effect[[i]]=ggplotGrob(ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
          geom_col(width = 0.3)+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                       aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                       aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#ee8c11","#1173EE"))+
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
                plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
          geom_vline(xintercept =0,color="gray",linetype="dashed")+
          ylab("")+
          xlab(paste(full_name[m]))+
          xlim(xlim_list[[i]])+
          #ggtitle(paste0(scenario[i]))+
          scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                           labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
          geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
          #ggtitle(paste0(scenario[i]))+
          geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
          
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                    aes(y=biome,x = -30),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
          geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                    aes(y=biome,x = 30),size=3,vjust=-1.6,
                    label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                      aes(y=biome,x = -30),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
            geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                      aes(y=biome,x = 30),size=3,vjust=-3,
                      label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        )
        
      }
      
    }
  }
  
  
pp_mean_effect_guild[[m]]=pp_mean_effect
}

###for the overal fungal diversity

pp_mean_effect_guild=list()
for (m in 1:9)
{
  if(m<2){
    pp_mean_effect=list()
    for (i in 1:4)
    {
      if(i<3)
      {
        
        pp_mean_effect[[i]]=ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
                                         geom_col(width = 0.3)+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                                                      aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                                                      aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#d95f02", "#1b9e77"))+
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
                                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                                               panel.background = element_rect(fill = "NA"),
                                               panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                         geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                         ylab("")+
                                         #xlab(paste(full_name[m]))+
                                         xlim(xlim_list_am[[i]])+
                                          xlab("")+
                                         scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                                                          labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
                                         geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                         #ggtitle(paste0(scenario[i]))+
                                         geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="#7570b3")+
                                         
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -50),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 50),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -50),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 50),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        
      }
      else{
        pp_mean_effect[[i]]=ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
                                         geom_col(width = 0.3)+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                                                      aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                                                      aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#d95f02", "#1b9e77"))+
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
                                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                                               panel.background = element_rect(fill = "NA"),
                                               panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                         geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                         ylab("")+
                                         #xlab(paste(full_name[m]))+
                                         xlim(xlim_list_am[[i]])+
                                         xlab("")+
                                         scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                                                          labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
                                         geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                         #ggtitle(paste0(scenario[i]))+
                                         geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="#7570b3")+
                                         
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -30),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 30),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -30),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 30),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        
        
      }
      
    }
  }
  
  else{
    pp_mean_effect=list()
    for (i in 1:4)
    {
      if(i<3)
      {
        
        pp_mean_effect[[i]]=ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
                                         geom_col(width = 0.3)+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                                                      aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                                                      aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#d95f02", "#1b9e77"))+
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
                                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                                               panel.background = element_rect(fill = "NA"),
                                               panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                         geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                         ylab("")+
                                         #xlab(paste(full_name[m]))+
                                         xlim(xlim_list[[i]])+
                                         xlab("")+
                                         scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                                                          labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
                                         geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                         #ggtitle(paste0(scenario[i]))+
                                         geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="#7570b3")+
                                         
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -10),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 10),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -10),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 10),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        
      }
      else{
        pp_mean_effect[[i]]=ggplot(data=data_biome_change_rate_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
                                         geom_col(width = 0.3)+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                                                      aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                                                      aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                                                      arrow = arrow(length = unit(0.1, "cm"),angle=90))+
                                         scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#d95f02", "#1b9e77"))+
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
                                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                                               panel.background = element_rect(fill = "NA"),
                                               panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
                                         geom_vline(xintercept =0,color="gray",linetype="dashed")+
                                         ylab("")+
                                         #xlab(paste(full_name[m]))+
                                         xlim(xlim_list[[i]])+
                                           xlab("")+
                                         #ggtitle(paste0(scenario[i]))+
                                         scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                                                          labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
                                         geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
                                         #ggtitle(paste0(scenario[i]))+
                                         geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="#7570b3")+
                                         
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -30),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 30),size=3,vjust=-1.6,
                                                   label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
                                                   aes(y=biome,x = -30),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
                                         geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
                                                   aes(y=biome,x = 30),size=3,vjust=-3,
                                                   label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
        
        
      }
      
    }
  }
  
  
  pp_mean_effect_guild[[m]]=pp_mean_effect
}




# to combined the range and the magnitude of both effects
# each i corresponds to land and climate effect with two climate change scenarios

pp_combine_effect_guild=list()
for (m in 1:9){
  
pp_combine_effect=list()
for (i in 1:4)
{
  if(i<3){
    pp_combine_effect[[i]]=ggplotGrob(pp_mean_effect_guild[[m]][[i]]+
      annotation_custom(
        grob = pp_pie_guild[[m]][[i]],
        xmin = -32, xmax =-27, # Adjust x-axis position of the circle
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
          panel.border = element_rect(color = "black", size = 0.6, fill = NA)))
  }
  else{
    pp_combine_effect[[i]]=ggplotGrob(pp_mean_effect_guild[[m]][[i]]+
      annotation_custom(
        grob = pp_pie_guild[[m]][[i]],
        xmin = -68, xmax =-50, # Adjust x-axis position of the circle
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
          panel.border = element_rect(color = "black", size = 0.6, fill = NA)))
  }
}
pp_combine_effect_guild[[m]]=pp_combine_effect
}

# do not set the height of the map,
# it is important to check the plot margin of the plots for better alignment

#p_land_245$heights=pp_combine_effect[[1]]$heights

p_land_latitude_245$heights=pp_combine_effect[[1]]$heights

p_land_latitude_585$heights=pp_combine_effect[[2]]$heights

p_climate_latitude_585$heights=pp_combine_effect[[4]]$heights



pp_combine_effect[[3]]$widths=pp_combine_effect[[1]]$widths

#p_climate_245$widths=pp_combine_effect[[1]]$widths


p1=plot_grid(p_land_245,p_land_latitude_245, pp_combine_effect[[1]],ncol=3,rel_widths = c(1,0.3,0.5))

p2=plot_grid(p_climate_245,p_climate_latitude_245, pp_combine_effect[[3]],ncol=3,rel_widths = c(1,0.3,0.5))

p3=plot_grid(p_land_585,p_land_latitude_585, pp_combine_effect[[2]],ncol=3,rel_widths = c(1,0.3,0.5))
 
p4=plot_grid(p_climate_585,p_climate_latitude_585, pp_combine_effect[[4]],ncol=3,rel_widths = c(1,0.3,0.5))
 
 
plot_grid(p1,p2,p3,p4,ncol=1)



pp_combine_effect_guild[[m]]

plot_grid(pp_combine_effect_guild[[1]][[1]],pp_combine_effect_guild[[2]][[1]],pp_combine_effect_guild[[3]][[1]],
          pp_combine_effect_guild[[4]][[1]],ncol=2 )

pp_mean_effect_guild[[1]][[1]]$widths=pp_mean_effect_guild[[2]][[1]]$widths
pp_mean_effect_guild[[1]][[2]]$widths=pp_mean_effect_guild[[2]][[2]]$widths

#to align the panels
pp_mean_effect_guild[[1]][[1]]$widths=pp_mean_effect_guild[[1]][[3]]$widths

pp_mean_effect_guild[[2]][[1]]$widths=pp_mean_effect_guild[[2]][[3]]$widths
pp_mean_effect_guild[[2]][[2]]$widths=pp_mean_effect_guild[[2]][[3]]$widths

pp_mean_effect_guild[[3]][[1]]$widths=pp_mean_effect_guild[[3]][[3]]$widths
pp_mean_effect_guild[[4]][[1]]$widths=pp_mean_effect_guild[[4]][[3]]$widths

pp_mean_effect_guild[[2]][[4]]$widths=pp_mean_effect_guild[[2]][[3]]$widths

pp_mean_effect_guild[[3]][[2]]$widths=pp_mean_effect_guild[[3]][[1]]$widths
pp_mean_effect_guild[[3]][[4]]$widths=pp_mean_effect_guild[[3]][[1]]$widths

pp_mean_effect_guild[[5]][[3]]$widths=pp_mean_effect_guild[[5]][[2]]$widths
pp_mean_effect_guild[[5]][[1]]$widths=pp_mean_effect_guild[[5]][[2]]$widths
pp_mean_effect_guild[[5]][[4]]$widths=pp_mean_effect_guild[[5]][[2]]$widths

pp_mean_effect_guild[[6]][[2]]$widths=pp_mean_effect_guild[[6]][[1]]$widths
pp_mean_effect_guild[[7]][[2]]$widths=pp_mean_effect_guild[[7]][[1]]$widths

pp_mean_effect_guild[[2]][[2]]$widths=pp_mean_effect_guild[[2]][[1]]$widths

plot_grid(pp_mean_effect_guild[[1]][[1]],
          pp_mean_effect_guild[[2]][[1]],
          pp_mean_effect_guild[[3]][[1]],
          pp_mean_effect_guild[[4]][[1]],
          pp_mean_effect_guild[[5]][[1]],
          pp_mean_effect_guild[[6]][[1]],
          pp_mean_effect_guild[[7]][[1]],
          pp_mean_effect_guild[[8]][[1]],
          
          pp_mean_effect_guild[[1]][[3]],
          pp_mean_effect_guild[[2]][[3]],
          pp_mean_effect_guild[[3]][[3]],
          pp_mean_effect_guild[[4]][[3]],
          pp_mean_effect_guild[[5]][[3]],
          pp_mean_effect_guild[[6]][[3]],
          pp_mean_effect_guild[[7]][[3]],
          pp_mean_effect_guild[[8]][[3]],
          
          pp_mean_effect_guild[[1]][[2]],
          pp_mean_effect_guild[[2]][[2]],
          pp_mean_effect_guild[[3]][[2]],
          pp_mean_effect_guild[[4]][[2]],
          pp_mean_effect_guild[[5]][[2]],
          pp_mean_effect_guild[[6]][[2]],
          pp_mean_effect_guild[[7]][[2]],
          pp_mean_effect_guild[[8]][[2]],
          
          
          
          pp_mean_effect_guild[[1]][[4]],
          pp_mean_effect_guild[[2]][[4]],
          pp_mean_effect_guild[[3]][[4]],
          pp_mean_effect_guild[[4]][[4]],
          pp_mean_effect_guild[[5]][[4]],
          pp_mean_effect_guild[[6]][[4]],
          pp_mean_effect_guild[[7]][[4]],
          pp_mean_effect_guild[[8]][[4]],
          ncol=8,label_x = -0.1,label_y = 1.08,
          label_size = 10,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:32])

##12x16



plot_grid(pp_pie_guild[[1]][[1]],
          pp_pie_guild[[2]][[1]],
          pp_pie_guild[[3]][[1]],
          pp_pie_guild[[4]][[1]],
          pp_pie_guild[[5]][[1]],
          pp_pie_guild[[6]][[1]],
          pp_pie_guild[[7]][[1]],
          pp_pie_guild[[8]][[1]],
          
          pp_pie_guild[[1]][[3]],
          pp_pie_guild[[2]][[3]],
          pp_pie_guild[[3]][[3]],
          pp_pie_guild[[4]][[3]],
          pp_pie_guild[[5]][[3]],
          pp_pie_guild[[6]][[3]],
          pp_pie_guild[[7]][[3]],
          pp_pie_guild[[8]][[3]],
          
          pp_pie_guild[[1]][[2]],
          pp_pie_guild[[2]][[2]],
          pp_pie_guild[[3]][[2]],
          pp_pie_guild[[4]][[2]],
          pp_pie_guild[[5]][[2]],
          pp_pie_guild[[6]][[2]],
          pp_pie_guild[[7]][[2]],
          pp_pie_guild[[8]][[2]],
          
          pp_pie_guild[[1]][[4]],
          pp_pie_guild[[2]][[4]],
          pp_pie_guild[[3]][[4]],
          pp_pie_guild[[4]][[4]],
          pp_pie_guild[[5]][[4]],
          pp_pie_guild[[6]][[4]],
          pp_pie_guild[[7]][[4]],
          pp_pie_guild[[8]][[4]],
          ncol=8,label_x = -0.01251,label_y = 1.07,
          label_size = 10,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:32]
      )

