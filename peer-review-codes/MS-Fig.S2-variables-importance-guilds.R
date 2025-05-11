#just select the four guilds included in the main text

d_importance%>%filter(guild%in%c("arbuscular_mycorrhizal", "ectomycorrhizal" , 
                                 "plant_pathogen" ,        "soil_saprotroph"   ))->d_importance_sub




ggplot(data=d_importance_sub, aes(x = guild, y = mean_value, fill = variable)) +
  geom_bar(stat = "identity",color="black",width = 0.8,position = "dodge")+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.5,position = position_dodge(0.8))+
  scale_fill_manual("Variables",breaks=c("temp_seasonality","soilInCaClpH",
                                         "mat_celsius","organicCPercent","mat_celsius_2","map_mm_2","map_mm"),
                    labels=c("Temp.seas.","pH", "MAT","SoilC",expression(MAT^2),expression(MAP^2),"MAP"),
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan"))+
  theme(legend.position ="right",
        legend.text = element_text(size=10),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 0.5,size=15,angle=0,color="black"), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("Importance(%)")+
  scale_x_discrete(breaks=unique(d_importance_sub$guild),
                   labels = c("AM","ECM","Plant pathogens",
                              "Soil saprotrophs" ))



