
# for different guilds
#set the space for different plot

xlim_list=list(c(-40,40),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20),c(-20,20))

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
        
        scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#E1BE6A", "#40B0A6"))+
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
                  aes(y=biome,x = -5),size=4,vjust=-1.6,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 5),size=4,vjust=-1.6,
                  label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                  aes(y=biome,x = -5),size=4,vjust=-3,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
        geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                  aes(y=biome,x = 5),size=4,vjust=-3,
                  label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
        geom_point(data=com_data_climate_guild[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="black")+
        scale_y_discrete(breaks=unique(com_data_climate_guild[[i]][[m]]$biome),position="right",
                         labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
        xlab(full_name[i])+
        ggtitle(scenario[m])+
        xlim(xlim_list[[i]])
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
      
      scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#E1BE6A", "#40B0A6"))+
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
                aes(y=biome,x = -10),size=4,vjust=-1.6,
                label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull( origin_mean )*100,")"),")")))+
      geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                aes(y=biome,x = 10),size=4,vjust=-1.6,
                label=c(paste0("(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull( origin_mean)*100,")"),")")))+
      geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative"),
                aes(y=biome,x = -10),size=4,vjust=-3,
                label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Negative")%>%pull(.group))+
      geom_text(data=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive"),
                aes(y=biome,x = 10),size=4,vjust=-3,
                label=com_data_climate_guild[[i]][[m]]%>%filter(change=="Positive")%>%pull(.group))+
      geom_point(data=com_data_climate_guild[[i]][[m]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="black")+
      scale_y_discrete(breaks=unique(com_data_climate_guild[[i]][[m]]$biome),position="right",
                       labels=paste0(rev(c("","","","","")),"(",sprintf("%.2f",com_data_climate_guild[[i]][[m]]%>%distinct( overal_mean )%>%pull()*100,")"),")"))+
      xlab(full_name[i])+
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
                                                 labels=c("Gain","No change","Loss"),values=c( "#40B0A6","gray","#E1BE6A")))}
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
                                                 labels=c("Gain","No change","Loss"),values=c("#40B0A6","gray","#E1BE6A")))}
  }
  pp_pie_guild[[m]]=pp_pie
}


pp_combine_effect_guild=list()
for (m in c(1,2,3,6,9)){
  if(m>1){
    pp_combine_effect=list()
    for (i in 1:4)
    {
      pp_combine_effect[[i]]=ggplotGrob(pp_guild[[m]][[i]]+
                                          annotation_custom(
                                            grob = pp_pie_guild[[m]][[i]],
                                            xmin = -19, xmax =-19, # Adjust x-axis position of the circle
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
                                            plot.margin = unit(c(0.3, 0.1, 0.2, 0.1), "cm"),
                                            panel.background = element_rect(fill = "NA"),
                                            panel.border = element_rect(color = "black", size = 0.6, fill = NA)))
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
                                            plot.margin = unit(c(0.3, 0.1, 0.2, 0.1), "cm"),
                                            panel.background = element_rect(fill = "NA"),
                                            panel.border = element_rect(color = "black", size = 0.6, fill = NA)))
    }
  }
  
  pp_combine_effect_guild[[m]]=pp_combine_effect
}

plot_grid(pp_combine_effect_guild[[1]][[1]],
          pp_combine_effect_guild[[2]][[1]],
          pp_combine_effect_guild[[3]][[1]],
          pp_combine_effect_guild[[6]][[1]],
          
          pp_combine_effect_guild[[1]][[3]],
          pp_combine_effect_guild[[2]][[3]],
          pp_combine_effect_guild[[3]][[3]],
          pp_combine_effect_guild[[6]][[3]],
          
          pp_combine_effect_guild[[1]][[2]],
          pp_combine_effect_guild[[2]][[2]],
          pp_combine_effect_guild[[3]][[2]],
          pp_combine_effect_guild[[6]][[2]],
          
          pp_combine_effect_guild[[1]][[4]],
          pp_combine_effect_guild[[2]][[4]],
          pp_combine_effect_guild[[3]][[4]],
          pp_combine_effect_guild[[6]][[4]],
          label_x = 0,label_y = 0.99,
          label_size = 15,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:16],ncol=4)
          

####
pp_combine_effect_guild[[1]][[1]]$widths=pp_combine_effect_guild[[1]][[2]]$widths

pp_combine_effect_guild[[2]][[1]]$widths=pp_combine_effect_guild[[2]][[2]]$widths
pp_combine_effect_guild[[2]][[3]]$widths=pp_combine_effect_guild[[2]][[2]]$widths

pp_combine_effect_guild[[2]][[4]]$widths=pp_combine_effect_guild[[2]][[2]]$widths

pp_combine_effect_guild[[3]][[1]]$widths=pp_combine_effect_guild[[3]][[2]]$widths

pp_combine_effect_guild[[3]][[3]]$widths=pp_combine_effect_guild[[3]][[2]]$widths

pp_combine_effect_guild[[6]][[3]]$widths=pp_combine_effect_guild[[6]][[2]]$widths


## get the mean value for different guilds 

species_change_land_rcp245%>%mutate(biome=rep(grid_level_biomes$LABEL,times=9))%>%
  filter(value>0)%>%
  mutate(trans_value=log(value))->df1

species_change_land_rcp245%>%mutate(biome=rep(grid_level_biomes$LABEL,times=9))%>%
  filter(value<0)%>%
  mutate(trans_value=sqrt(abs(value)))->df2

df1%>%bind_rows(df2)->df3

%>%group_by(variable,biome)%>%summarise(mean=mean(trans_value))%>%filter(variable=="all")



