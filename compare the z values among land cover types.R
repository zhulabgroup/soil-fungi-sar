# comparing the z values among different land cover types
#the land cover is plot-based

library(neonUtilities)
library(ggplot2)
library(car)
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding="euc-jp")

plant.data <- loadByProduct(dpID="DP1.10058.001")
# get the plant diversity and cover of the plot
plant.coar=plant.data[[4]]
plant.fine=plant.data[[5]]
plotID=substr(plant.fine$namedLocation,1,8)
vegetation_type=data.frame(cbind(plotID,plant.fine$nlcdClass))# 11 types of vegetation
vegetation_type=unique(vegetation_type)
comp_vege=merge(all_z_30,vegetation_type,by="plotID")# the former data set is all the estimated z for each plot(30 z values per plot)
names(comp_vege)=c("plotID","z","type")
comp_vege=subset(comp_vege,z<10)

write.csv(comp_vege,"comp_vege.csv")# for some reasons, there are NA for two plots as to their land cover

comp_vege=subset(comp_vege,plotID!="UNDE_044")

comp_vege=subset(comp_vege,plotID!="CPER_048")

#test the hormogenety of the variance


leveneTest(z~type,data=comp_vege)

oneway.test(z~type,data=comp_vege,na.action=na.omit, var.equal=FALSE)

tukey(comp_vege$z, comp_vege$type, method="G")

source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding="euc-jp")

# create the plots
k=aggregate(z~type,data=comp_vege,FUN=mean)
od=k[order(k$z),]# with the increase trend to display the box plots
comp_vege$type=factor(comp_vege$type,levels=od$type)
ggboxplot(comp_vege,x="type",y="z",fill = "type")+
  xlab("")+
  theme(legend.position = "bottom", legend.text = element_text(size=14), text = element_text(size=15), axis.text.x = element_blank(),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank())+
  ylab("z")+
  guides(fill=guide_legend(nrow=4, byrow=TRUE))+# make four rows for the legend
  geom_hline(yintercept =0.7559199,linetype="dashed",color="red" )+
  scale_fill_discrete("",breaks=od$type,labels=c("dwarfScrub(N=3)" , "emergentHerbaceousWetlands(N=4)", "grasslandHerbaceous(N=64)" , "pastureHay(N=17)", "shrubScrub(N=51)" , "evergreenForest(N=93)","woodyWetlands(N=29)","deciduousForest(N=96)",  "cultivatedCrops(N=22)"  ,  "mixedForest(N=20)" , "sedgeHerbaceous(N=4)"   ))
