##

library(plyr)
library(dplyr)
library(gstat)
library(raster)
library(ggplot2)
library(car)
library(classInt)
library(RStoolbox)
library(caret)
library(caretEnsemble)
library(doParallel)
library(gridExtra)

###

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
grid<-read.csv( "GP_prediction_grid_data.csv", header= TRUE)

train=read.csv("train_data.csv",sep=",")

powerTransform(train$SOC)

train$SOC.bc<-bcPower(train$SOC, 0.2523339)

train.xy<-train[,c(1,24,8:9)]
train.df<-train[,c(1,24,11:21)]

grid.xy<-grid[,c(1,2:3)]
grid.df<-grid[,c(4:14)]

RESPONSE<-train.df$SOC.bc
train.x<-train.df[3:13]


coordinates(train.xy) = ~x+y
coordinates(grid.xy) = ~x+y


mc <- makeCluster(detectCores())
registerDoParallel(mc)

myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1856)
GLM<-train(train.x,
           RESPONSE,
           method = "glm",
           trControl=myControl,
           preProc=c('center', 'scale'))
print(GLM)




# Extract residuals
train.xy$residuals.glm<-resid(GLM)
# Variogram
v.glm<-variogram(residuals.glm~ 1, data = train.xy,cutoff=300000,width=300000/15)

v.glm<-variogram(residuals.glm~ 1, data = train.xy)

# Intial parameter set by eye esitmation
m.glm<-vgm(0.15,"Exp",40000,0.05)
# least square fit
m.f.glm<-fit.variogram(v.glm, m.glm)
m.f.glm


###

leveneTest(core_rich ~ type, data = compare_richness)
oneway.test(core_rich ~ type, data = compare_richness, na.action = na.omit, var.equal = FALSE)
tukey(compare_richness$core_rich, compare_richness$type, method = "G")
###create the plots
k <- aggregate(core_rich ~ type, data = compare_richness, FUN = mean)
od <- k[order(k$core_rich), ] # with the increase trend to display the box plots
compare_richness$type <- factor(compare_richness$type, levels = od$type)



p1=ggplot(compare_richness,aes(x=type,y=core_rich,fill=type),alpha=0.5)+
  sm_raincloud(size=0.1)+
  geom_boxplot(width = 0.05, color = "black", size=0.1,outlier.size = 2) +
  theme(legend.position = c(0.5,0.2), 

        plot.title = element_text(hjust = 0.5, vjust = 2),
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(angle=90), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  
  ylab("Alpha diversity") +
  xlab("")+
  annotate("text", x = 1, y = 340, label = "c", size = 6) +
  annotate("text", x = 2, y = 260, label = "c", size = 6) +
  annotate("text", x = 3, y = 240, label = "bc", size = 6) +
  annotate("text", x = 4, y = 240, label = "bc", size = 6) +
  annotate("text", x = 5, y = 280, label = "bc", size = 6) +
  annotate("text", x = 6, y = 320, label = "bc", size = 6) +
  annotate("text", x = 7, y = 300, label = "b", size = 6) +
  annotate("text", x = 8, y = 340, label = "a", size = 6) +
  ylim(-160,360)
  
  



### for the plot-level diversity

library(dplyr)
library(neonUtilities)
library(ggplot2)
library(ggpubr)
library(car)
library(tidyverse)
library(ggthemes)
library(multcompView)

source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")


leveneTest(log(beta_mean) ~ type, data = compare_richness)
oneway.test(log(beta_mean) ~ type, data = compare_richness, na.action = na.omit, var.equal = FALSE)
tukey(log(compare_richness$beta_mean), compare_richness$type, method = "G")

###create the plots
k <- aggregate(log(beta_mean) ~ type, data = compare_richness, FUN = mean)
od <- k[order(k$`log(beta_mean)`), ] # with the increase trend to display the box plots
compare_richness$type <- factor(compare_richness$type, levels = od$type)
tukey(log(compare_richness$beta_mean), compare_richness$type, method = "G")
tukey(log(compare_richness$beta_mean), compare_richness$type, method = "G")$Games.Howell%>%data.frame()%>%dplyr::select(p)%>%round(digits = 3)

p2=ggplot(compare_richness,aes(x=type,y=log(beta_mean),fill=type),alpha=0.5)+
  sm_raincloud(size=0.1)+
  geom_boxplot(width = 0.05, color = "black", size=0.1,outlier.size = 1) +
  theme(legend.position = c(0.5,0.2), 
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(angle=90), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Log(Beta diversity)") +
  xlab("")+
  annotate("text", x = 1, y = 0.03, label = "ab", size = 6) +
  annotate("text", x = 2, y = 0.03, label = "ab", size = 6) +
  annotate("text", x = 3, y = 0.03, label = "b", size = 6) +
  annotate("text", x = 4, y = 0.03, label = "ab", size = 6) +
  annotate("text", x = 5, y = 0.03, label = "ab", size = 6) +
  annotate("text", x = 6, y = 0.03, label = "ab", size = 6) +
  annotate("text", x = 7, y = 0.03, label = "a", size = 6) +
  annotate("text", x = 8, y = 0.03, label = "ab", size = 6) +
  ylim(-0.43,0.051)


## for the plot level richness

leveneTest(log(plot_rich) ~ type, data = compare_richness)

oneway.test(log(plot_rich) ~ type, data = compare_richness, na.action = na.omit, var.equal = FALSE)

tukey(log(compare_richness$plot_rich), compare_richness$type, method = "G")

###create the plots
k <- aggregate(log(plot_rich) ~ type, data = compare_richness, FUN = mean)

od <- k[order(k$`log(plot_rich)`), ] # with the increase trend to display the box plots

compare_richness$type <- factor(compare_richness$type, levels = od$type)


tukey(log(compare_richness$plot_rich), compare_richness$type, method = "G")

tukey(log(compare_richness$plot_rich), compare_richness$type, method = "G")$Games.Howell%>%data.frame()%>%dplyr::select(p)%>%round(digits = 3)

# need to adjust the plot

tukey(log(compare_richness$plot_rich), compare_richness$type, method = "G")$Games.Howell%>%data.frame()%>%dplyr::select(p)->op

round(p.adjust(op$p),digits = 3)%>%data.frame()%>%mutate(pari=rownames(op))


##

p3=ggplot(compare_richness,aes(x=type,y=log(plot_rich),fill=type),alpha=0.5)+
  sm_raincloud(size=0.1)+
  geom_boxplot(width = 0.05, color = "black", size=0.1,outlier.size = 1) +
  theme(legend.position = c(0.5,0.2), 
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Log(Gamma diversity)") +
  xlab("")+
  annotate("text", x = 1, y = 9, label = "b", size = 6) +
  annotate("text", x = 2, y = 9, label = "b", size = 6) +
  annotate("text", x = 3, y = 9, label = "b", size = 6) +
  annotate("text", x = 4, y = 9, label = "ab", size = 6) +
  annotate("text", x = 5, y = 9, label = "ab", size = 6) +
  annotate("text", x = 6, y = 9, label = "ab", size = 6) +
  annotate("text", x = 7, y = 9, label = "a", size = 6) +
  annotate("text", x = 8, y = 9, label = "a", size = 6) +
  ylim(4,9)

cowplot::plot_grid(p1,p2,p3,ncol=2)

p1=ggplotGrob(p1)

p2=ggplotGrob(p2)
p3=ggplotGrob(p3)

p1$widths=p2$widths
p1$widths=p3$widths

cowplot::plot_grid(p1,p2,p3,ncol=3)

temp_land_rich

type0=c( "cultivatedCrops",  "deciduousForest","evergreenForest", "grasslandHerbaceous","mixedForest"  ,"pastureHay",  "shrubScrub","woodyWetlands")

leveneTest(core_rich ~ type, data = compare_richness)

oneway.test(core_rich ~ type, data = compare_richness, na.action = na.omit, var.equal = FALSE)

tukey(compare_richness$core_rich, compare_richness$type, method = "G")

###create the plots
k <- aggregate(value ~ type, data = temp_land_rich, FUN = mean)
od <- k[order(k$value), ] # with the increase trend to display the box plots
temp_land_rich$type <- factor(temp_land_rich$type, levels = od$type)




ggplot(temp_land_rich,aes(x=type,y=value,fill=type),alpha=0.5)+
  sm_raincloud(size=0.51)+
  geom_boxplot(width = 0.05, color = "black", outlier.size = 2) +
  theme(legend.position =c(0.5,0.153), legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        panel.grid.major = element_line(color="white"),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Alpha diversity") +
  annotate("text", x = 1, y = 175, label = "g", size = 6) +

  annotate("text", x = 2, y = 170, label = "f", size = 6) +

  annotate("text", x = 3, y = 200, label = "e", size = 6) +

  annotate("text", x = 4, y = 178, label = "d", size = 6) +
  annotate("text", x = 5, y = 230, label = "c", size = 6) +
  annotate("text", x = 6, y = 220, label = "c", size = 6) +
  annotate("text", x = 7, y = 225, label = "b", size = 6) +
  annotate("text", x = 8, y = 255, label = "a", size = 6) +
  ylim(80,255)
  
  ### when the season was not consided 

ggplot(land_richness,aes(x=variable,y=value,fill=variable),alpha=0.5)+
  sm_raincloud(size=0.1)+
  geom_boxplot(width = 0.05, color = "black", size=0.1,outlier.size = 1) +
  theme(legend.position = c(0.5,0.2), 
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Alpha diversity") +
  xlab("")+
  annotate("text", x = 1, y = 175, label = "h", size = 6) +
  annotate("text", x = 2, y = 170, label = "g", size = 6) +
  annotate("text", x = 3, y = 180, label = "f", size = 6) +
  annotate("text", x = 4, y = 180, label = "e", size = 6) +
  annotate("text", x = 5, y = 195, label = "d", size = 6) +
  annotate("text", x = 6, y = 210, label = "c", size = 6) +
  annotate("text", x = 7, y = 218, label = "b", size = 6) +
  annotate("text", x = 8, y = 255, label = "a", size = 6) +
  ylim(40,260)

mod1=aov(log(value)~variable,data=land_richness)
summary(mod1)
  
TukeyHSD(mod1)

