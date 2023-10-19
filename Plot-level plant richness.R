#get the plot level plant diversity
library(neonUtilities)

plant.data <- loadByProduct(dpID="DP1.10058.001")
# get the plant diversity and cover of the plot

plant.coar=plant.data[[4]]
plant.fine=plant.data[[5]]
plotID=substr(plant.fine$namedLocation,1,8)
vegetation_type=data.frame(cbind(plotID,plant.fine$nlcdClass))# 11 types of vegetation
# at the 20*20 m plot, species were recorded with abundance?
# why thoses recorded in the 10 m2 plot still be recored in the 100 m2
# just focused on the plotID and all species recorded in the plot are included for richness
plant.rich.fine=plant.fine[,c("plotID","scientificName","subplotID")]# 1m2 scale data
plant.rich.coar=plant.coar[,c("plotID","scientificName","subplotID")]# 10m2 and 100m2 scale data
plant.rich=rbind(plant.rich.fine,plant.rich.coar)
plot_id=unique(plant.rich$plotID)

rich=numeric()
for (i in 1:length(plot_id)){
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  a=subset(plant.rich,plotID==plot_id[i])
  rich[i]=length(unique(a$scientificName))# 'sp'with unkonwn identity was treated as a species
}

plot_plant_rich=data.frame(cbind(plot_id,rich))

write.csv(plot_plant_rich,"plot_plant_rich.csv")

