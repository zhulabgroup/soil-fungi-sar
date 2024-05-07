# changes in the land cover types
save(data_mean,file="data_mean_present.RData")
save(data_mean_future,file="data_mean_future.RData")
# compute changes in the land use, the ratio in land use cover between the two time periods

ratio=matrix(ncol=35,nrow=103567)
for (i in 4:35){
  ratio[,i]=data_mean_future[,i]/data_mean[,i]
  
}


# can we just look at part of the habitats

data_mean=mutate(data_mean,loca=paste(data_mean$lon,"_",-1*data_mean$lat))

data_mean_future=mutate(data_mean_future,loca=paste(data_mean_future$lon,"_",-1*data_mean_future$lat))

richness_zva=cbind(pred_plot,pred_zvalue[,3])%>%mutate(loca=paste(lon,"_",lat))

names(richness_zva)[4]="zvalue"
# to see some change in the land cover types for some grids
#
ratio=cbind(data_mean$loca,ratio[,10])%>%data.frame()# to see this type of the richness

ratio$X2=as.numeric(ratio$X2)
names(ratio)=c("loca","change")

species_change=merge(richness_zva,ratio,by="loca")
# the ratio between the future and present species richness
species_change=species_change%>%mutate(ratio_rich=change^zvalue)%>%mutate(future_rich=ratio_rich*rich)%>%mutate(change_rich=future_rich-rich)%>%mutate(change_rate=change_rich/rich)




# changes in the land cover of interest

ggplot(data_mean) +
  geom_point(data=data_mean,pch=15,aes(x=lon, y=-lat,color=X8), size=0.275)+
  scale_color_gradient(expression("Cover %"),low = "blue", high = "yellow")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted changes in PFT15")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

###

ggplot(data_mean_future) +
  geom_point(data=data_mean_future,pch=15,aes(x=lon, y=-lat,color=X8), size=0.275)+
  scale_color_gradient(expression("Cover %"),low = "blue", high = "yellow")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted changes in PFT15")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

# changes in the richness
ggplot(species_change) +
  geom_point(data=species_change,pch=15,aes(x=lon, y=lat,color=change_rate), size=0.275)+
  scale_color_gradient(expression("Species change %"),low = "purple", high = "seagreen1",limits = c(-1, 8))+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted species changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

# mapping the richness

ggplot(species_change) +
  geom_point(data=species_change,pch=15,aes(x=lon, y=lat,color=future_rich), size=0.275)+
  scale_color_gradient(expression("Species change %"),low = "purple", high = "seagreen1")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted species changes")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

