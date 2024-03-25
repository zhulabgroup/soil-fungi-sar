library(ggplot2)
library(ggpubr)
# date required
# comp_vege#each plot has 30 estimated z values
# com_guild#each fungal guil has 30 estimated z values

1. # test the homogeneity of the variance
# excluded three land cover types that have less than 5 plots 
comp_vege=subset(comp_vege,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
leveneTest(z ~ type, data = comp_vege)
oneway.test(z ~ type, data = comp_vege, na.action = na.omit, var.equal = FALSE)
tukey(comp_vege$z, comp_vege$type, method = "G")

source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

2. # create the plot
k <- aggregate(z ~ type, data = comp_vege, FUN = mean)
od <- k[order(k$z), ] # with the increase trend to display the box plots
comp_vege$type <- factor(comp_vege$type, levels = od$type)

aa=ggboxplot(comp_vege, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  ylab("z")+
  theme(legend.position = c(0.48, 0.8058), 
        legend.text = element_text(size = 14),
        text = element_text(size = 15), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "italic", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199,
             linetype = "dashed", 
             color = "red")+ 
  scale_fill_manual("", breaks = od$type, 
  values = c(  "cadetblue1", "tan", "greenyellow", "mediumseagreen", "wheat", "cornsilk","lavender","mediumpurple"), 
  labels = c( "grassland\nHerbaceous(N=64)", "pastureHay(N=17)", "shrubScrub(N=51)", "evergreenForest(N=93)", "woodyWetlands(N=29)", "deciduousForest(N=96)", "cultivatedCrops(N=22)", "mixedForest(N=20)"))+
  annotate("text", x = 1, y = 1.03, label = "c", size = 6) +
  annotate("text", x = 2, y = 1.05, label = "c", size = 6) +
  annotate("text", x = 3, y = 1.07, label = "c", size = 6) +
  annotate("text", x = 4, y = 1.15, label = "b", size = 6) +
  annotate("text", x = 5, y = 1.13, label = "ab", size = 6) +
  annotate("text", x = 6, y = 1.125, label = "ab", size = 6) +
  annotate("text", x = 7, y = 1.218, label = "ab", size = 6) +
  annotate("text", x = 8, y = 1.1108, label = "a", size = 6)

# load()

3. # test the homogeneity of the variance
leveneTest(z ~ guild, data = com_guild) # testing variance homogenety, unblanced
oneway.test(z ~ guild, data = com_guild, na.action = na.omit, var.equal = FALSE)
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

tukey(com_guild$z, com_guild$guild, method = "G") # multiple comparison with 'Games.Howell' Post-hoc Test

4. # create a plot
k <- aggregate(z ~ guild, data = com_guild, FUN = mean)
od1 <- k[order(k$z), ] # with the increase trend to display the box plots

com_guild$guild <- factor(com_guild$guild, levels = od1$guild)

 bb=ggboxplot(com_guild, x = "guild", y = "z", fill = "guild", outlier.colour = "gray", outlier.shape = NA) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  geom_hline(yintercept = 0.787, linetype = "dashed", color = "red") +
  xlab("") +
  theme(legend.position = c(0.51, 0.85), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  scale_fill_manual("", breaks = od1$guild, values = c("chocolate1", "gray", "cadetblue1", "greenyellow", "forestgreen", "purple","lavender", "tan"), labels = c("soil saprotroph(N=493)", "parasite(N=459)", "wood saprotroph(N=482)", "plant pathogen(N=482)", "litter saprotroph(N=493)", "ECM(N=493)","epiphyte(N=482)",  "ACM(N=482)")) +
  annotate(x = 1, y = 1.2, "text", label = "g", size = 6) +
  annotate(x = 2, y = 1.5, "text", label = "f", size = 6) +
  annotate(x = 3, y = 1.56, "text", label = "e", size = 6) +
  annotate(x = 4, y = 1.53, "text", label = "e", size = 6) +
  annotate(x = 5, y = 1.45, "text", label = "d", size = 6) +
  annotate(x = 6, y = 1.4, "text", label = "c", size = 6) +
  annotate(x = 7, y = 1.705, "text", label = "b", size = 6) +
  annotate(x = 8, y = 1.8, "text", label = "a", size = 6) +
  ylim(0, 3)
 
5. # arrange the two plots

plot_grid(aa, bb, ncol = 1, labels = c("(a)", "(b)"))
