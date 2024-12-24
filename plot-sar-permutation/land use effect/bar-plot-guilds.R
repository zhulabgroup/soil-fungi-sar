##
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



##



data_compare_LETTER[[i]]$LABEL=factor(data_compare_LETTER[[i]]$LABEL,
                      levels=rev(c("All",
                               "Temperate Broadleaf & Mixed Forests",
                               "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                               "Tropical & Subtropical Dry Broadleaf Forests")))

  
      
    ggplot(data=data_compare_LETTER[[i]],aes(fill=guild,y=LABEL ,x=ori_mean,group=guild))+
          geom_col(width = 0.6,position = "dodge",stat = "identity")+
      geom_errorbar(aes(xmin = ori_lower, xmax = ori_upper), 
                    position = position_dodge(width = 0.6), 
                    width = 0.15,size=0.3) +
          theme(legend.position = c(0.85,0.15),
                legend.text = element_text(size=8),
                legend.title  = element_text(size=10),
                text = element_text(size = 18),
                plot.title = element_text(size = 12, hjust = 0.6), 
                axis.text.y = element_text(size = 10), 
                axis.text.x = element_text(size = 10), 
                axis.title.y = element_text(size = 10), 
                axis.title.x = element_text(size = 10), 
                legend.key.size = unit(0.3, "cm"),
                plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                panel.background = element_rect(fill = "NA"),
                panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
          geom_vline(xintercept =0,color="gray",linetype="dashed")+
          geom_text(aes(label=LETTER),position = position_dodge(width = 0.6),
                    hjust = ifelse(data_compare_LETTER[[i]]$ori_mean > 0, -5.25, 8.25),size=2)+
                   
          geom_text(aes(label= c(paste0("(",sprintf("%.2f",data_compare_LETTER[[i]]$ori_mean*100,")"),")"))),
                    position = position_dodge(width = 0.6),
                    hjust = ifelse(data_compare_LETTER[[i]]$ori_mean > 0, -.25, 1.25),size=2)
        
                    
      
       
       
      
     plot_grid(pk,pp[[2]],pp[[3]],pp[[4]],ncol=2)           
      
     plot_grid(pp[[2]],ncol=2)           
     
                
                
     pk=ggplot(data=data_compare_LETTER[[1]],aes(fill=guild,y=LABEL ,x=ori_mean))+
         geom_col(width = 0.6,position = "dodge")+
         xlim(-0.5,0.5)+
         theme(legend.position = c(0.85,0.15),
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
         
         geom_text(aes(label=LETTER),position = position_dodge(width = 0.6),
                   hjust = ifelse(data_compare_LETTER[[1]]$ori_mean > 0, -5.25, 8.25),
                   size=2)+
         geom_text(aes(label= c(paste0("(",sprintf("%.2f",data_compare_LETTER[[1]]$ori_mean*100,")"),")"))),
                   position = position_dodge(width = 0.6),
                   hjust = ifelse(data_compare_LETTER[[1]]$ori_mean > 0, -.25, 1.25),
                   size=2)
       
     
     
