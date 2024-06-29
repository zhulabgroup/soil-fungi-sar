## to include the land use data with original data

library(ggplot2)
library(dplyr)
library(multcomp)

load("~/soil-sar/data/comp_vege.RData")
load("~/soil-sar/plot-sar-permutation/rare_all.Rdata")

d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site correspondes to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(rare_all, plotIDM)# merge the new plotid with the initial data 

sample_data(d)%>%data.frame()%>%dplyr::select(plotIDM)->temp
names(temp)="plotID"

comp_vege%>%dplyr::select(plotID,type)%>%distinct()->land

left_join(temp,land,by="plotID")%>%dplyr::select(plotID,type)->temp

temp%>%mutate(type = replace_na(type, "evergreenForest"))%>%data.frame()%>%dplyr::select(type)->land

# for na rows assign that with "evergreen forest"
row.names(land)=row.names(sample_data(d))
land=sample_data(land)
d<- merge_phyloseq(rare_all, land)#
# each core is assigned with a land cover type
# to see the richness for each

#type Freq
#             cultivatedCrops  175
#            deciduousForest 1344
#                  dwarfScrub   12
#  emergentHerbaceousWetlands   46
#            evergreenForest 2602
#         grasslandHerbaceous  908
#                 mixedForest  260
#                  pastureHay  222## this is artificial land use type
#             sedgeHerbaceous   14
#                shrubScrub  515
#              woodyWetlands  280

# get the same number of cores and determine the mean value of core-level richness

type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","pastureHay","shrubScrub","woodyWetlands")
# does not consider differences in the sampling season
set.seed(123)
kk=list()
for(j in 1:500)
  {
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  
  richness=numeric()
  for (i in 1:8)
  {
  
    dk=subset_samples(d,type==type0[i])
    
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 175)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich=matrix(ncol=8,nrow=500)
for (i in 1:500)
{
  land_rich[i,]=kk[[i]]
    
}

land_rich%>%data.frame()%>%rename_all(~paste0(c(type0)))->d

land_rich=melt(land_rich)

# create a plot to show the core level richness

d$variable<- factor(d$variable, levels = c("evergreenForest","woodyWetlands","pastureHay","cultivatedCrops","shrubScrub","mixedForest","deciduousForest","grasslandHerbaceous"))

load("~/soil-sar/land_richness.RData")
land_richness=d
save(land_richness,file="land_richness.RData")

ggboxplot(d%>%filter(variable!="pastureHay"), x = "variable", y = "value", fill = "variable", outlier.shape = NA) +
  xlab("") +
  ylab("Richness") +
  theme(legend.position = c(0.25, 0.7558), 
        legend.text = element_text(size = 14), 
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15),
        axis.title.y = element_text( size = 20), 
        axis.title.x = element_text(size = 20)) +
  guides(fill = guide_legend(nrow = 7, byrow = TRUE)) + # make four rows for the legend
  xlab("Land cover types")

# consider the difference in the sampling season 
## all the samples were taken from the season of peakGreeness

sample_data(rare_all)%>%data.frame()%>% dplyr::select(sampleTiming)%>%count(sampleTiming,name="count")

green_sample=subset_samples(d,sampleTiming=="peakGreenness")

sample_data(green_sample)%>%data.frame()%>%dplyr::select(type)%>%count(type,name="count")->green_sample

##
#1             cultivatedCrops    51
#2             deciduousForest   461
#3                  dwarfScrub    12
#4  emergentHerbaceousWetlands    18
#5             evergreenForest   668
#6         grasslandHerbaceous   306
#7                 mixedForest    90
#8                  pastureHay    50
#9             sedgeHerbaceous    14
#10                 shrubScrub   213
#11              woodyWetlands   104


sub_green=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:8)
  {
    dk=subset_samples(green_sample,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 50)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  sub_green[[j]]=richness
}

# when the sampling season is the same

land_rich_green=matrix(ncol=8,nrow=500)
for (i in 1:500)
{
  land_rich_green[i,]=sub_green[[i]]
  
}

land_rich_green%>%data.frame()%>%rename_all(~paste0(c(type0)))->temp_land_rich
melt(land_rich_green)%>%mutate(type=rep(type0,each=500))->temp_land_rich
mod=aov(log(value)~type,data=temp_land_rich)

save(temp_land_rich,file="temp_land_rich.RData")

ggplot(temp_land_rich,aes(x=type,y=value,fill=type),alpha=0.5)+
  sm_raincloud()+
  geom_boxplot(width = 0.1, color = "black", outlier.size = 2) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank()) +
  ylab("Alpha diversity") 


## samples from the same season

sample_data(rare_all)%>%data.frame()%>% select(sampleTiming)%>%count(sampleTiming,name="count")

green_sample=subset_samples(d,sampleTiming=="peakGreenness")

sample_data(green_sample)%>%data.frame()%>%dplyr::select(type)%>%count(type,name="count")->green_sample

##
#1             cultivatedCrops    51
#2             deciduousForest   461
#3                  dwarfScrub    12
#4  emergentHerbaceousWetlands    18
#5             evergreenForest   668
#6         grasslandHerbaceous   306
#7                 mixedForest    90
#8                  pastureHay    50
#9             sedgeHerbaceous    14
#10                 shrubScrub   213
#11              woodyWetlands   104

#to test the difference in the community composition
#for each land use type, we can selected the same number of samples and combine them as a virtual community.

sub_green=list()
 
  for (i in 1:8)
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
    dk=subset_samples(green_sample,type==type0[i])
    
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 50)
    sampled_physeq <- prune_samples(sampled_names, dk)
    
  sub_green[[i]]=otu_table(sampled_physeq )%>%data.frame()
}
# bind all the community data

vege_com=dplyr::bind_rows(sub_green[[1]],sub_green[[2]],sub_green[[3]],sub_green[[4]],sub_green[[5]],sub_green[[6]],sub_green[[7]],sub_green[[8]])%>%mutate(type=rep(type0,each=50))

ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
plot(ordination, type = "points", col =1:8)

com_score=ordination$species
com_score=ordination$points%>%data.frame()%>%mutate(type=rep(type0,each=50))

# change the 0-1 data

vege_com[vege_com>0]=1

adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)


ggplot(data=com_score,aes(x=MDS1,  y= MDS2,color=type ))+
  geom_point(data=com_score,pch=21,color="black",aes(x=MDS1,y= MDS2,fill=type ),size=4,alpha=0.75)+
  stat_ellipse(size=0.8,linetype="dashed")+
  theme(legend.position = c(0.849,0.213), 
        legend.title = element_text(size=10),
        text = element_text(size = 18), 
        legend.text = element_text(size=11),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  annotate("text",x=0.35,y=0.705,label=expression(italic(F)*" = 5.12"),size=7)+
annotate("text",x=0.75,y=0.705,label=expression(italic(P)*" <0.001"),size=7)+
  ylim(-1.2,0.75)+
  scale_fill_manual("NLCD",breaks=type0,labels=type0,
    values=c("magenta","tomato","seagreen1","purple","seagreen","royalblue","orange","pink"))+

scale_color_manual("NLCD",breaks=type0,labels=type0,
                  values=c("magenta","tomato","seagreen1","purple","seagreen","royalblue","orange","pink"))+
  geom_vline(xintercept = 0,color="gray",linetype="dashed",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dashed",size=1)
  

  

## difference in the mean richness among different land cover types, based on plot-level measurements 

compare_richness=plot_diversity_env[,c(1,2,3,5)]%>%left_join(vege_tep)%>%mutate(type = ifelse(nchar(plotID)<4, "evergreenForest", type))%>%mutate(site=substr(plotID,1,4))%>%filter(!is.na(type))%>%filter(type%in%type0)

## compare the core level richness for the forest types
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

leveneTest(core_rich ~ type, data = compare_richness)
oneway.test(core_rich ~ type, data = compare_richness, na.action = na.omit, var.equal = FALSE)
tukey(compare_richness$core_rich, compare_richness$type, method = "G")

###create the plots
k <- aggregate(core_rich ~ type, data = compare_richness, FUN = mean)
od <- k[order(k$core_rich), ] # with the increase trend to display the box plots
compare_richness$type <- factor(compare_richness$type, levels = od$type)
ggplot(compare_richness,aes(x=type,y=core_rich,fill=type),alpha=0.5)+
  sm_raincloud()+
  geom_boxplot(width = 0.1, color = "black", outlier.size = 2) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank()) +
  ylab("Core level richness") 


save(compare_richness,file="compare_richness.RData")
tukey(compare_richness$core_rich, compare_richness$type, method = "G")$result1

tukey(compare_richness$core_rich, compare_richness$type, method = "G")$Games.Howell%>%data.frame()%>%dplyr::select(p)%>%round(digits = 3)%>%filter(p<0.05)


round(data.frame(p.adjust(comp[[3]][,3])),digits = 3)


### compare the beta diversity

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
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  ylab("Core level richness") +
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







