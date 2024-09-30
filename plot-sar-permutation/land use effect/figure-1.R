# codes for figure 1

load("~/soil-sar/plot-sar-permutation/temp_land_rich.RData")

temp_land_rich$type=factor(temp_land_rich$type,levels=c("grasslandHerbaceous", "deciduousForest"  , "mixedForest" , "evergreenForest" ,"cultivatedCrops" ,   
                                                        "shrubScrub", "pastureHay", "woodyWetlands" ))

p2=ggplot(temp_land_rich,aes(x=type,y=value,fill=type),alpha=0.5)+
  sm_raincloud()+
  geom_boxplot(width = 0.1, color = "black", outlier.size = 2) +
  theme(legend.position = c(0.5,0.15), 
        legend.title = element_text(size=10),
        text = element_text(size = 18), 
        legend.text = element_text(size=11),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=12), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA)) +
  ylim(70,260)+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) + # make four rows for the legend
  ylab("Alpha diversity") +
  xlab("Land use type")+
  scale_fill_manual("",breaks=type0,labels=c("Cultivated crops" ,    "Deciduous forest" , "EvergreenForest" , "GrasslandHerbaceous",
                                             "MixedForest", "PastureHay"  ,"ShrubScrub"  , "WoodyWetlands"),             
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+
  annotate("text",x=1,y=260,label="a",size=6)+
  annotate("text",x=2,y=230,label="b",size=6)+
  annotate("text",x=3,y=225,label="c",size=6)+
  annotate("text",x=4,y=235,label="c",size=6)+
  annotate("text",x=5,y=180,label="d",size=6)+
  annotate("text",x=6,y=200,label="e",size=6)+
  annotate("text",x=7,y=170,label="f",size=6)+
  annotate("text",x=8,y=180,label="g",size=6)


p1=ggplotGrob(p1)
p2=ggplotGrob(p2)

p1$heights=p2$heights

p1$heights=c$heights
## samples from the same season, try to see the indicative species



#to test the difference in the community composition
#for each land use type, we can selected the same number of samples and combine them as a virtual community.

set.seed(544)
sub_green=list()
for (i in 1:8)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
  dk=subset_samples(green_sample,type==type0[i])
  
  sample_names=sample_names(dk)
  sampled_names <- sample(sample_names, 50,replace = FALSE)
  sampled_physeq <- prune_samples(sampled_names, dk)
  
  sub_green[[i]]=otu_table(sampled_physeq )%>%data.frame()
}

# bind all the community data; the presented results were based on one permutation

vege_com=dplyr::bind_rows(sub_green[[1]],sub_green[[2]],sub_green[[3]],sub_green[[4]],sub_green[[5]],sub_green[[6]],sub_green[[7]],sub_green[[8]])%>%mutate(type=rep(type0,each=50))

ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
plot(ordination, type = "points")
com_score=ordination$species
com_score=ordination$points%>%data.frame()%>%mutate(type=rep(type0,each=50))

# change the 0-1 data

vege_com[vege_com>0]=1

adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)


ggplot(data=com_score,aes(x=MDS1,  y= MDS2,color=type ))+
  geom_point(data=com_score,pch=21,color="black",aes(x=MDS1,y= MDS2,fill=type ),size=3,alpha=0.75)+
  stat_ellipse(size=0.8,linetype="dashed")+
  
  theme(legend.position = c(0.449,0.213), 
        legend.title = element_text(size=10),
        text = element_text(size = 18), 
        
        legend.text = element_text(size=11),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1,margin = margin(t = -17)), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1, fill = NA)) +
  annotate("text",x=0.35,y=0.705,label=expression(italic(F)*" = 5.12"),size=7)+
  annotate("text",x=0.84,y=0.705,label=expression(italic(P)*" <0.001"),size=7)+
  ylim(-1.7,0.75)+
  xlab("MDS1")+
  
  scale_fill_manual("",breaks=type0,labels=c("CultivatedCrops" ,    "DeciduousForest" , "EvergreenForest" , "GrasslandHerbaceous",
                                             "MixedForest", "PastureHay"  ,"ShrubScrub"  , "WoodyWetlands"),             
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+
  
  scale_color_manual("",breaks=type0,labels=c("CultivatedCrops" ,    "DeciduousForest" , "EvergreenForest" , "GrasslandHerbaceous",
                                              "MixedForest"   , "PastureHay"  ,"ShrubScrub"  , "WoodyWetlands"),
                     values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # make four rows for the legend



p1=ggplotGrob(p1)
p2=ggplotGrob(p2)

p1$heights=p2$heights

p1$heights=c$heights

