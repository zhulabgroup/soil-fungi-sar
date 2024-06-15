# the land cover type is plot-based
library(neonUtilities)
library(ggplot2)
library(ggpubr)
library(car)
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

# this analysis requires the data generated in step 1

1.# create a data frame with each plot with 30 estimated z values
all_z_ranall=read.csv("all_z.ranall.csv")
a=list()
for (i in 1:dim(all_z_ranall)[1])
  {
 a[[i]] =t(all_z_ranall[i,3:32])
}

b=a[[1]]
for (i in 2:dim(all_z_ranall)[1])
  {
 b=rbind(b,a[[i]]) 
}
plotID=rep(all_z_ranall$a1,each=30)
all_z_30=cbind(plotID,b)
all_z_30 <- data.frame(all_z_30)
all_z_30$V2 <- as.numeric(all_z_30$V2)

2. # extract the land cover type data
plant.data <- loadByProduct(dpID = "DP1.10058.001")
# get the plant diversity and cover of the plot
plant.coar <- plant.data[[4]]
plant.fine <- plant.data[[5]]
plotID <- substr(plant.fine$namedLocation, 1, 8)
vegetation_type <- data.frame(cbind(plotID, plant.fine$nlcdClass)) # 11 land cover types
vegetation_type <- unique(vegetation_type)

3. #cbind the z values and land cover type
comp_vege <- merge(all_z_30, vegetation_type, by = "plotID") # the former data set is all the estimated z for each plot(30 z values per plot)
names(comp_vege) <- c("plotID", "z", "type")
comp_vege <- subset(comp_vege, z < 10)
comp_vege <- subset(comp_vege, plotID != "UNDE_044")
comp_vege <- subset(comp_vege, plotID != "CPER_048")
write.csv(comp_vege, "comp_vege.csv") # for some reasons, there are NA for two plots with their land cover


4. #test the homogeneity of the variance
leveneTest(z ~ type, data = comp_vege)
oneway.test(z ~ type, data = comp_vege, na.action = na.omit, var.equal = FALSE)
tukey(comp_vege$z, comp_vege$type, method = "G")
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

5.#create the plot
k <- aggregate(z ~ type, data = comp_vege, FUN = mean)
od <- k[order(k$z), ] # with the increase trend to display the box plots
comp_vege$type <- factor(comp_vege$type, levels = od$type)

ggboxplot(comp_vege, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.8058), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  scale_fill_manual("", breaks = od$type, values = c("chocolate1", "gray", "cadetblue1", "lavender", "greenyellow", "mediumseagreen", "cornsilk", "tan", "wheat", "gold", "mediumpurple"), labels = c("dwarfScrub(N=3)", "emergentHerbaceous\nWetlands(N=4)", "grassland\nHerbaceous(N=64)", "pastureHay(N=17)", "shrubScrub(N=51)", "evergreenForest(N=93)", "woodyWetlands(N=29)", "deciduousForest(N=96)", "cultivatedCrops(N=22)", "mixedForest(N=20)", "sedgeHerbaceous(N=4)")) +
  annotate("text", x = 1, y = 0.96, label = "b", size = 6) +
  annotate("text", x = 2, y = 0.94, label = "b", size = 6) +
  annotate("text", x = 3, y = 1.05, label = "b", size = 6) +
  annotate("text", x = 4, y = 1.06, label = "b", size = 6) +
  annotate("text", x = 5, y = 1.07, label = "b", size = 6) +
  annotate("text", x = 6, y = 1.15, label = "a", size = 6) +
  annotate("text", x = 7, y = 1.13, label = "a", size = 6) +
  annotate("text", x = 8, y = 1.08, label = "a", size = 6) +
  annotate("text", x = 9, y = 1.21, label = "a", size = 6) +
  annotate("text", x = 10, y = 1.11, label = "a", size = 6) +
  annotate("text", x = 11, y = 1.5, label = "ab", size = 6)
