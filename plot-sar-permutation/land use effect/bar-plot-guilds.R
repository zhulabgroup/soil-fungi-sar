##the below codes are disregard
pp_guild=list()

for(m in c(1,2,3,6))
{
  

pp=list()

for(i in 1:4)
{
  

pp[[i]]=ggplot(data=data_model_mean_guild[[m]][[i]],aes(fill=type,y=LABEL ,x=100*origin_mean))+
  geom_col(width = 0.3)+
  geom_segment(data=data_model_mean_guild[[m]][[i]], 
               aes(y=LABEL,yend=LABEL,xend= 100*ori_low , x = 100*ori_up),
               arrow = arrow(length = unit(0.1, "cm"),ends = "both",angle=90))+

  scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#fdc086", "#7fc97f"))+
  
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
  #xlim(xlim_list_am[[i]])+
  xlab("")+
  scale_y_discrete(breaks=unique(data_model_mean_guild[[m]][[i]]$LABEL),position="right",
                   labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_model_mean_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(5)*100,")"),")"))+
  geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
  #ggtitle(paste0(scenario[i]))+
  geom_point(data=data_model_mean_guild[[m]][[i]],aes(y=LABEL,x=100*overal_mean),pch=23,color="black",size=2,fill="#beaed4")+
  
  geom_text(data=data_model_mean_guild[[m]][[i]]%>%filter(type=="Negative"),
            aes(y=LABEL,x = -25),size=3,vjust=-1.6,
            label=c(paste0("(",sprintf("%.2f",data_model_mean_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( origin_mean )%>%round(5)*100,")"),")")))+
  geom_text(data=data_model_mean_guild[[m]][[i]]%>%filter(type=="Positive"),
            aes(y=LABEL,x = 25),size=3,vjust=-1.6,
            label=c(paste0("(",sprintf("%.2f",data_model_mean_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( origin_mean )%>%round(5)*100,")"),")")))+
  geom_text(data=data_model_mean_guild[[m]][[i]]%>%filter(type=="Negative"),
            aes(y=LABEL,x = -25),size=3,vjust=-3,
            label=data_model_mean_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
  geom_text(data=data_model_mean_guild[[m]][[i]]%>%filter(type=="Positive"),
            aes(y=LABEL,x = 25),size=3,vjust=-3,
            label=data_model_mean_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))+
  xlim(-60,60)
}

pp_guild[[m]]=pp
}

#modified the first data list

data_compare_LETTER[[1]]%>%filter(LABEL=="Tropical & Subtropical Dry Broadleaf Forests")%>%
  mutate(across(where(is.numeric), ~ 0))->df

df$type=gsub("Negative","Positive",df$type)

data_compare_LETTER[[1]]%>%bind_rows(df)->df

data_compare_LETTER[[1]]=df
##
impact=c("Land use impact","Land use impact","Climate impact","Climate impact")



pp=list()
for (j in 1:4)
  {
  data_compare_LETTER[[j]]$LABEL=factor(data_compare_LETTER[[j]]$LABEL,
                                        levels=rev(c("All",
                                                     "Temperate Broadleaf & Mixed Forests",
                                                     "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                                     "Tropical & Subtropical Dry Broadleaf Forests")))
pp[[j]]=ggplot(data=data_compare_LETTER[[j]],aes(fill=guild,y=LABEL ,x=ori_mean*100,group=guild))+
          geom_col(width = 0.8,position = "dodge",stat = "identity",color="black",size=0.1,space=0.5)+
      geom_errorbar(aes(xmin = ori_lower*100, xmax = ori_upper*100), 
                    position = position_dodge(width = 0.8), 
                    width = 0.15,size=0.3) +
      geom_text(aes(label=LETTER),position = position_dodge(width = 0.8),
                hjust = ifelse(data_compare_LETTER[[j]]$ori_mean> 0, -5.5, 8),size=2)+
      geom_text(aes(label= c(paste0("(",sprintf("%.2f",ori_mean*100,")"),")"))),
                position = position_dodge(width = 0.8),
                hjust = ifelse(data_compare_LETTER[[j]]$ori_mean > 0, -.25, 1.25),size=2)+
      xlim(-30,30)+
      scale_y_discrete(breaks=unique(data_compare_LETTER[[j]]$LABEL),position="right",
                       labels=paste0(c("",""," ","","")))+
      xlab("Richness change rate (%)")+
      ylab("")+
      ggtitle(scenario[j])+
    theme(legend.position = c(0.8,0.25),
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10), 
          axis.title.y = element_text(size = 12), 
          axis.title.x = element_text(size = 12), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0, 0, 0.1, 0.1), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
    geom_vline(xintercept =0,color="gray",linetype="dashed")+
    scale_fill_manual("",breaks=c("AM","EM","plapat","soilsap"),
                      labels=c("AM","EM","Plant pathogens","Soil saprotrophs"),
                      values=c("#c94e65","#037f77","royalblue","chocolate1"))
    
} 


#to arrange the four plots           
  plot_grid(pp[[3]],pp[[1]],pp[[4]],pp[[2]],
            ncol=2,
            label_x = 0,label_y = 1.03,
            labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:4])           
      
     
    
       
      
     
     