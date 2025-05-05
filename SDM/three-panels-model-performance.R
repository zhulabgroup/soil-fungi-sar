# figures for the model performance and the diversity gradient

############### for the model performance


variable_importance_present=load("variable_importance_present.RData")
#need to open this file in the repository
load("/Volumes/seas-zhukai/proj-soil-fungi/variable_importance_present.RData")

variable_importance_present%>%melt()%>%group_by(variable)%>%summarise(mean_value=mean(value),sd_value=sd(value))->
  variable_importance_present_data

variable_importance_present_data$variable=factor(variable_importance_present_data$variable,
                                                 levels=c( "temp_seasonality", "soilInCaClpH" , "organicCPercent" , "mat_celsius_2" ,   "mat_celsius"  ,   
                                                           "map_mm"  , "map_mm_2"  ))

p1=ggplot(variable_importance_present_data, aes(x = variable, y = mean_value)) +
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),size=0.35,width=0.15)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=8,angle=60), 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
  ylab("Importance (%)")+
  xlab("")+
  scale_x_discrete(breaks=c("temp_seasonality", "soilInCaClpH" ,    "organicCPercent",  "mat_celsius_2",    "mat_celsius"  ,   
                            "map_mm", "map_mm_2" ),labels=c("Temp.seas.","pH","SoilC",expression(MAT^2),"MAT","MAP",expression(MAP^2)))



model_evaluation_present=readRDS("/Volumes/seas-zhukai/proj-soil-fungi/model_evaluation_present.rds")



model_evaluation_present%>%melt()%>%group_by(variable)%>%summarise(mean_value=mean(value),sd_value=sd(value))%>%
  filter(variable%in%c("prop.correct","sensitivity","specificity","AUC"))->model_evaluation_present_data

model_evaluation_present_data$variable=factor(model_evaluation_present_data$variable,
                                              levels=c("specificity","prop.correct","sensitivity","AUC"))

p2=ggplot(model_evaluation_present_data, aes(x = variable, y = mean_value)) +
  geom_bar(stat = "identity",width=0.3)+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.2)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10,angle=60), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, -0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("Model performance")+
  scale_x_discrete(breaks=c("AUC","sensitivity","specificity","prop.correct"),
                   labels=c("AUC","Sensitivity","Specificity","Prop.correct"))+coord_flip()



###########

load(file="present_richness_all.RData")

load("/Volumes/seas-zhukai/proj-soil-fungi/present_richness_all.RData")

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

select_biome=c("Temperate Grasslands, Savannas & Shrublands",
               "Temperate Conifer Forests",
               "Temperate Broadleaf & Mixed Forests",
               "Tropical & Subtropical Dry Broadleaf Forests")

present_richness_all%>%bind_cols(grid_level_biomes%>%select(LABEL))%>%
  data.frame()%>%
  dplyr::rename(richness=...3)%>%filter(LABEL%in%select_biome)->present_richness_all


## need to make projections for the values
present_richness_all_project=my_function_project(present_richness_all)





p2=ggplot(present_richness_all_project) +
  geom_tile(data =present_richness_all_project ,aes(x = x, y = y, fill = last), size = 0.175)+
  scale_fill_gradientn("Taxonomic diversity", colors = brewer.pal(9, "Oranges"),
                       guide = guide_colorbar(order = 2,barwidth = unit(0.5, "cm"), barheight = unit(2, "cm")))+

  add_sf_layers()+
  coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
  xlab("")+
  ylab("")+
  ggtitle("")+
  theme(legend.spacing.y = unit(32, "pt"), 
        legend.position = c(0.15,0.35),
        legend.margin = margin(t = -30, r = 0, b = -1, l = 0),
        legend.text = element_text(size=12,angle=0),
        legend.box = "vertical",
        legend.justification = "center",
        legend.title = element_text(margin = margin(b = 4),size=13),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, -0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank())




p1=ggplotGrob(p1)
p2=ggplotGrob(p2)


p1$widths=p2$widths

p3=ggplotGrob(p3)


p1$heights=p3$heights/2
p2$heights=p3$heights/2

plot_grid(p1,p2,ncol=1)->p4

p4=ggplotGrob(p4)

p4$heights=p3$heights

plot_grid(p4,p3,ncol=2,rel_heights = c(1,2))

plot(p4)

p1$heights=p2$heights
