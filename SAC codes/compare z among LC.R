# compare the turnover rate among land cover types
library(dplyr)
library(neonUtilities)
library(ggplot2)
library(ggpubr)
library(car)
library(tidyverse)
library(ggthemes)
library(multcompView)

load("~/soil-sar/data/comp_guild.RData")# this should be the updated data
load("~/soil-sar/data/comp_vege.RData")
d=comp_vege[,c(2,4)]%>%unique()
d1=merge(com_guild,d,by="plotID")
gu=unique(d1$guild)

# for each guild, i examined if turnover rate differ among land cover types
# land cover types with less than 5 plots were excluded

source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

comp=list()
od=list()
ann=list()# to test if species turn over rate differ among land cover types
for (i in 1:length(gu))
     {
       guild_spe=subset(d1,guild==gu[i])
       guild_spe=subset(guild_spe,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
       leveneTest(z ~ type, data = guild_spe)
       ann[[i]]=oneway.test(z ~ type, data = guild_spe, na.action = na.omit, var.equal = FALSE)
       k <- aggregate(z ~ type, data = guild_spe, FUN = mean)
       od[[i]]<- k[order(k$z), ] # with the increase trend to display the box plots
       guild_spe$type <- factor(guild_spe$type, levels = od[[i]]$type)
       tk=tukey(guild_spe$z, guild_spe$type, method = "G")
       comp[[i]]=tk[[2]]
}


data.frame(round(p.adjust(comp[[1]][,"p"]),2))



# ACM fungi

round(data.frame(p.adjust(comp[[5]][,3])),digits = 3)
data_acm=subset(d1,guild==gu[5])
data_acm=subset(data_acm,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
data_acm$type <- factor(data_acm$type, levels = od[[5]]$type)

p1=ggboxplot(data_acm, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 1.875, label = "c", size = 6) +
  annotate("text", x = 2, y = 1.648485, label = "bc", size = 6) +
  annotate("text", x = 3, y = 1.89325403285, label = "abc", size = 6) +
  annotate("text", x = 4, y = 1.7038526, label = "abc", size = 6) +
  annotate("text", x = 5, y = 1.86427, label = "a", size = 6) +
  annotate("text", x = 6, y = 1.70630584860525, label = "abc", size = 6) +
  annotate("text", x = 7, y = 1.8576023, label = "a", size = 6) +
  annotate("text", x = 8, y = 1.850652, label = "ab", size = 6) +
  xlab("ACM (N=5 pairs)")+
  ylim(0,4)+
  scale_fill_manual("", breaks = od[[5]]$type, values = c( "greenyellow", "cadetblue1", "tan" ,"cornsilk",  "mediumseagreen","wheat","lavender","purple"), 
                    labels = c( "shrubScrub(N=42)","grassland\nHerbaceous(N=63)","pastureHay(N=15)" ,"deciduousForest(N=87)", "evergreenForest(N=56)", "woodyWetlands(N=20)","cultivatedCrops(N=21)", "mixedForest(N=12)")) 


## for ecm fungi
round(data.frame(p.adjust(comp[[4]][,3])),digits = 3)
data_ecm=subset(d1,guild==gu[4])
data_ecm=subset(data_ecm,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
data_ecm$type <- factor(data_ecm$type, levels = od[[4]]$type)

# create the plots

p2=ggboxplot(data_ecm, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 1.35, label = "d", size = 6) +
  annotate("text", x = 2, y = 1.3309285, label = "c", size = 6) +
  annotate("text", x = 3, y = 1.43285, label = "c", size = 6) +
  annotate("text", x = 4, y = 1.4038526, label = "c", size = 6) +
  annotate("text", x = 5, y = 1.427, label = "c", size = 6) +
  annotate("text", x = 6, y = 1.54860525, label = "c", size = 6) +
  annotate("text", x = 7, y = 1.576023, label = "b", size = 6) +
  annotate("text", x = 8, y = 1.652, label = "a", size = 6) +
  xlab("ECM (N=18 pairs)")+
  ylim(0,4)+
  scale_fill_manual("N=18 pairs", breaks = od[[4]]$type, values = c("tan", "mediumseagreen","wheat", "purple","cornsilk", "cadetblue1", "greenyellow", "lavender"), 
                    labels = c( "pastureHay(N=17)", "evergreenForest(N=93)","woodyWetlands(N=29)" ,"mixedForest(N=20)","deciduousForest(N=96)","grassland\nHerbaceous(N=64)","shrubScrub(N=51)","cultivatedCrops(N=22)" )) 

# for wood sap
round(data.frame(p.adjust(comp[[3]][,3])),digits = 3)

data_woosap=subset(d1,guild==gu[3])
data_woosap=subset(data_woosap,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
data_woosap$type <- factor(data_woosap$type, levels = od[[3]]$type)

p3=ggboxplot(data_woosap, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 1.105, label = "e", size = 6) +
  annotate("text", x = 2, y = 1.325, label = "dc", size = 6) +
  annotate("text", x = 3, y = 1.25243285, label = "d", size = 6) +
  annotate("text", x = 4, y = 1.1238526, label = "d", size = 6) +
  annotate("text", x = 5, y = 1.1427, label = "bc", size = 6) +
  annotate("text", x = 6, y = 1.16525, label = "b", size = 6) +
  annotate("text", x = 7, y = 1.1623, label = "ab", size = 6) +
  annotate("text", x = 8, y = 1.12562, label = "a", size = 6) +
  xlab("Wood saprotroph (N=20 pairs)")+
  ylim(0,4)+
  scale_fill_manual("", breaks = od[[3]]$type, values = c( "cadetblue1", "lavender", "greenyellow", "tan", "cornsilk", "mediumseagreen", "wheat", "purple"), 
                    labels = c(  "grassland\nHerbaceous(N=64)","cultivatedCrops(N=22)",  "shrubScrub(N=51)", "pastureHay(N=17)","deciduousForest(N=95)","evergreenForest(N=92)", "woodyWetlands(N=29)", "mixedForest(N=20)")) 





# for the litter saprotroph

od[[6]]# order of the z among types

round(data.frame(p.adjust(comp[[6]][,3])),digits = 3)
data_litsap=subset(d1,guild==gu[6])
data_litsap=subset(data_litsap,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
data_litsap$type <- factor(data_litsap$type, levels = od[[6]]$type)

# create the plot

p4=ggboxplot(data_litsap, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 1.25, label = "e", size = 6) +
  annotate("text", x = 2, y = 1.25, label = "de", size = 6) +
  annotate("text", x = 3, y = 1.43285, label = "d", size = 6) +
  annotate("text", x = 4, y = 1.4038526, label = "cd", size = 6) +
  annotate("text", x = 5, y = 1.427, label = "c", size = 6) +
  annotate("text", x = 6, y = 1.54860525, label = "bc", size = 6) +
  annotate("text", x = 7, y = 1.576023, label = "b", size = 6) +
  annotate("text", x = 8, y = 1.52, label = "a", size = 6) +
  xlab("Litter saprotroph (N=20 pairs)")+
  ylim(0,4)+
  scale_fill_manual("", breaks = od[[6]]$type, values = c("tan", "cadetblue1",  "greenyellow", "lavender", "cornsilk", "wheat","mediumseagreen", "purple"), 
                    labels = c("pastureHay(N=17)",  "grassland\nHerbaceous(N=64)","shrubScrub(N=51)","cultivatedCrops(N=22)"   ,"deciduousForest(N=96)","woodyWetlands(N=29)","evergreenForest(N=93)","mixedForest(N=20)")) 



# for soil saprotroph
round(data.frame(p.adjust(comp[[1]][,3])),digits = 3)

data_soilsap=subset(d1,guild==gu[1])
data_soilsap=subset(data_soilsap,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
data_soilsap$type <- factor(data_soilsap$type, levels = od[[1]]$type)


p5=ggboxplot(data_soilsap, x = "type", y = "z", fill = "type", outlier.shape = NA) +
  xlab("") +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  ylab("z") +
  theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 1.1325, label = "c", size = 6) +
  annotate("text", x = 2, y = 1.1285, label = "bc", size = 6) +
  annotate("text", x = 3, y = 1.143285, label = "bc", size = 6) +
  annotate("text", x = 4, y = 1.14038526, label = "bc", size = 6) +
  annotate("text", x = 5, y = 1.15427, label = "bc", size = 6) +
  annotate("text", x = 6, y = 1.25160584860525, label = "b", size = 6) +
  annotate("text", x = 7, y = 1.16576023, label = "ab", size = 6) +
  annotate("text", x = 8, y = 1.321652, label = "a", size = 6) +
  xlab("Soil saprotroph (N=8 pairs)")+
  ylim(0,4)+
  scale_fill_manual("", breaks = od[[1]]$type, values = c("cornsilk","purple",  "cadetblue1", "mediumseagreen", "wheat", "greenyellow",  "tan" , "lavender"), 
                    labels = c("deciduousForest(N=96)", "mixedForest(N=20)", "grassland\nHerbaceous(N=64)", "evergreenForest(N=93)", "woodyWetlands(N=29)", "shrubScrub(N=51)","pastureHay(N=17)", "cultivatedCrops(N=22)" )) 





  # for epiphyte fungi
  round(data.frame(p.adjust(comp[[7]][,3])),digits = 3)
  data_epiphy=subset(d1,guild==gu[7])
  data_epiphy=subset(data_epiphy,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
  data_epiphy$type <- factor(data_epiphy$type, levels = od[[7]]$type)
  
  p6=ggboxplot(data_epiphy, x = "type", y = "z", fill = "type", outlier.shape = NA) +
    xlab("") +
    theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    ylab("z") +
    theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
    geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
    annotate("text", x = 1, y = 1.7625, label = "d", size = 6) +
    annotate("text", x = 2, y = 1.53985285, label = "c", size = 6) +
    annotate("text", x = 3, y = 1.62403285, label = "cd", size = 6) +
    annotate("text", x = 4, y = 1.8038526, label = "c", size = 6) +
    annotate("text", x = 5, y = 2.000427, label = "bc", size = 6) +
    annotate("text", x = 6, y = 1.654860525, label = "b", size = 6) +
    annotate("text", x = 7, y = 1.6576023, label = "ab", size = 6) +
    annotate("text", x = 8, y = 2.15, label = "a", size = 6) +
    xlab("Epiphyte (N=16 pairs)")+
    ylim(0,4)+
    scale_fill_manual("", breaks = od[[7]]$type, values = c("greenyellow", "cadetblue1", "tan","mediumseagreen","wheat","cornsilk",  "purple","lavender"), 
                      labels = c( "shrubScrub(N=51)","grassland\nHerbaceous(N=64)","pastureHay(N=17)" , "evergreenForest(N=92)","woodyWetlands(N=29)"  ,"deciduousForest(N=95)" , "mixedForest(N=20)","cultivatedCrops(N=22)")) 
  
  
  # for  parasitic fungi
  
  round(data.frame(p.adjust(comp[[8]][,3])),digits = 3)
  
  data_para=subset(d1,guild==gu[8])
  data_para=subset(data_para,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
  data_para$type <- factor(data_para$type, levels = od[[8]]$type)
  
  p7=ggboxplot(data_para, x = "type", y = "z", fill = "type", outlier.shape = NA) +
    xlab("") +
    theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    ylab("z") +
    theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
    geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
    annotate("text", x = 1, y = 1.325, label = "d", size = 6) +
    annotate("text", x = 2, y = 1.285, label = "cd", size = 6) +
    annotate("text", x = 3, y = 1.43285, label = "cd", size = 6) +
    annotate("text", x = 4, y = 1.4038526, label = "c", size = 6) +
    annotate("text", x = 5, y = 1.5427, label = "bc", size = 6) +
    annotate("text", x = 6, y = 1.60584860525, label = "b", size = 6) +
    annotate("text", x = 7, y = 1.6576023, label = "ab", size = 6) +
    annotate("text", x = 8, y = 1.652, label = "a", size = 6) +
    xlab("Parasite (N=14 pairs)")+
    ylim(0,4)+
    scale_fill_manual("", breaks = od[[8]]$type, values = c("cornsilk","tan", "purple", "cadetblue1","wheat" , "greenyellow", "mediumseagreen"  , "lavender"), 
                      labels = c("deciduousForest(N=96)","pastureHay(N=17)","mixedForest(N=18)",  "grassland\nHerbaceous(N=64)", "woodyWetlands(N=27)", "shrubScrub(N=46)","evergreenForest(N=91)", "cultivatedCrops(N=21)" )) 
  

 
  # for plant pathogens 
  
  round(data.frame(p.adjust(comp[[2]][,3])),digits = 3)
  data_plapat=subset(d1,guild==gu[2])
  data_plapat=subset(data_plapat,type!="dwarfScrub"&type!="sedgeHerbaceous"&type!="emergentHerbaceousWetlands")
  data_plapat$type <- factor(data_plapat$type, levels = od[[2]]$type)
  
  p8=ggboxplot(data_plapat, x = "type", y = "z", fill = "type", outlier.shape = NA) +
    xlab("") +
    theme(legend.position = "bottom", legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    ylab("z") +
    theme(legend.position = c(0.45, 0.7558), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
    geom_hline(yintercept = 0.7559199, linetype = "dashed", color = "red") +
    annotate("text", x = 1, y = 1.25, label = "e", size = 6) +
    annotate("text", x = 2, y = 1.43285, label = "d", size = 6) +
    annotate("text", x = 3, y = 1.2403285, label = "cd", size = 6) +
    annotate("text", x = 4, y = 1.4038526, label = "cd", size = 6) +
    annotate("text", x = 5, y = 1.43627, label = "bc", size = 6) +
    annotate("text", x = 6, y = 1.654860525, label = "b", size = 6) +
    annotate("text", x = 7, y = 1.576023, label = "ab", size = 6) +
    annotate("text", x = 8, y = 1.4552, label = "a", size = 6) +
    xlab("Plant pathogens (N=18 pairs)")+
    ylim(0,4)+
    scale_fill_manual("", breaks = od[[2]]$type, values = c( "cadetblue1","greenyellow", "lavender", "tan","cornsilk", "mediumseagreen", "wheat","purple"), 
                      labels = c(  "grassland\nHerbaceous(N=64)","shrubScrub(N=51)","cultivatedCrops(N=22)" ,"pastureHay(N=17)", "deciduousForest(N=95)", "evergreenForest(N=92)","woodyWetlands(N=29)" , "mixedForest(N=20)")) 
  
  
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,ncol=2,labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),label_size = 14)

   # 14x25 pdf

round(data.frame(p.adjust(comp[[6]][,3])),digits = 3)

nbb=numeric()
for (i in 1:8)
{
  
  signb=round(data.frame(p.adjust(comp[[i]][,3])),digits = 3) 
  
  nbb[i]=dim(subset(signb,p.adjust.comp..i.....3..<0.05))[1]
  
}



