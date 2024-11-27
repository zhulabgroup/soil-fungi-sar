## to include the land use data with original data

library(ggplot2)
library(dplyr)
library(multcomp)
library(phyloseq)
library(tidyr)

load("~/soil-sar/data/comp_vege.RData")
load("~/soil-sar/plot-sar-permutation/rare_all.Rdata")

d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site corresponds to a plot
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

# assign the NA rows with "evergreen forest"
row.names(land)=row.names(sample_data(d))
land=sample_data(land)
d<- merge_phyloseq(rare_all, land)#
# to see the number of different soil cores among different land use types
sample_data(d)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%dplyr::summarize(count=n())

sample_data(d)%>%data.frame()%>%dplyr::select(plotID,type)

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

d=land_richness

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





## all the samples were taken from peakGreeness

sample_data(rare_all)%>%data.frame()%>% dplyr::select(sampleTiming)%>%count(sampleTiming,name="count")#

green_sample=subset_samples(d,sampleTiming=="peakGreenness")

sample_data(green_sample)%>%data.frame()%>%dplyr::select(type)%>%count(type,name="count")

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

# when the sampling season is the same
set.seed(223)
sub_green=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:8)
  {
    dk=subset_samples(green_sample,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 40)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  sub_green[[j]]=richness
}


land_rich_green=matrix(ncol=8,nrow=500)
for (i in 1:500)
{
  land_rich_green[i,]=sub_green[[i]]
  
}

land_rich_green%>%data.frame()%>%rename_all(~paste0(c(type0)))->temp_land_rich
melt(land_rich_green)%>%mutate(type=rep(type0,each=500))->temp_land_rich
mod=aov(log(value)~type,data=temp_land_rich)

TukeyHSD(mod)
# to increase the range of the two human-modified plots, i resampled 40 samples for each plot
# and save the data on turbo

setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")

save(temp_land_rich,file="temp_land_rich_40_samples.RData")

save(temp_land_rich,file="temp_land_rich.RData")

load("~/soil-sar/plot-sar-permutation/temp_land_rich.RData")

temp_land_rich$type=factor(temp_land_rich$type,levels=c("grasslandHerbaceous", "deciduousForest"  , "mixedForest" , "evergreenForest" ,"cultivatedCrops" ,   
                                                        "shrubScrub", "pastureHay", "woodyWetlands" ))
p2=ggplot(temp_land_rich,aes(x=type,y=value,fill=type),alpha=0.5)+
  sm_raincloud()+
  geom_boxplot(width = 0.1, color = "black", outlier.size = 2) +
  theme(legend.position = c(0.45,0.2), 
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
  scale_fill_manual("",breaks=type0,labels=c("Cultivated crops" ,    
                                             "Deciduous forest" , 
                                             "Evergreen forest" ,
                                             "Grassland herbaceous",
                                                 "Mixed forest", 
                                             "Pasture hay"  ,
                                             "Shrub scrub"  , "Woody wetlands"),             
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+
  annotate("text",x=1,y=260,label="a",size=6)+
annotate("text",x=2,y=240,label="b",size=6)+
annotate("text",x=3,y=225,label="c",size=6)+
annotate("text",x=4,y=245,label="c",size=6)+
annotate("text",x=5,y=195,label="d",size=6)+
annotate("text",x=6,y=210,label="e",size=6)+
annotate("text",x=7,y=190,label="f",size=6)+
annotate("text",x=8,y=185,label="g",size=6)

  
p1=ggplotGrob(p1)
p2=ggplotGrob(p2)

p1$heights=p2$heights

p1$widths=p2$widths
## samples from the same season, try to see the indicative species

plot_grid(p1,p2,ncol=2,labels = c("(a)","(b)"))

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


p1=ggplot(data=com_score,aes(x=MDS1,  y= MDS2,color=type ))+
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
  annotate("text",x=0.15,y=0.705,label=expression(italic(F)*" = 5.12"),size=7)+
annotate("text",x=0.64,y=0.705,label=expression(italic(P)*" < 0.001"),size=7)+
  ylim(-1.7,0.75)+
  xlab("MDS1")+

  scale_fill_manual("",breaks=type0,labels=c("CultivatedCrops" ,    "DeciduousForest" , "EvergreenForest" , "GrasslandHerbaceous",
                                                 "MixedForest", "PastureHay"  ,"ShrubScrub"  , "WoodyWetlands"),             
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+

scale_color_manual("",breaks=type0,labels=c("CultivatedCrops" ,    "DeciduousForest" , "EvergreenForest" , "GrasslandHerbaceous",
                                                "MixedForest"   , "PastureHay"  ,"ShrubScrub"  , "WoodyWetlands"),
                  values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # make four rows for the legend
  






### to see if the result changes with permutation


num_cores <- detectCores() - 5  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","pastureHay","shrubScrub","woodyWetlands")

result=foreach(1:8,.packages=c("phyloseq","dplyr"))%dopar% {
  type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","pastureHay","shrubScrub","woodyWetlands")
  
  sub_green=list()
for (i in 1:8)
{
  dk=subset_samples(green_sample,type==type0[i])
  sample_names=sample_names(dk)
  sampled_names <- sample(sample_names, 50,replace = FALSE)
  sampled_physeq <- prune_samples(sampled_names, dk)
  
  sub_green[[i]]=otu_table(sampled_physeq )%>%data.frame()
}

# bind all the community data; the initial results were based on one time of simulation
vege_com=dplyr::bind_rows(sub_green[[1]],sub_green[[2]],sub_green[[3]],sub_green[[4]],sub_green[[5]],sub_green[[6]],sub_green[[7]],sub_green[[8]])%>%mutate(type=rep(type0,each=50))
vege_com[vege_com>0]=1
adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)
adonis_result$aov.tab[6]$`Pr(>F)`[1]# get the pvalues for the test

}

stopCluster(cl)



######## use the pararell function

set.seed(544)

sub_green=list()
for (i in 1:8)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
  dk=subset_samples(green_sample,type==type0[i])
  sample_names=sample_names(dk)
  num_cores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  result=foreach(1:20,.packages=c("phyloseq","dplyr"))%dopar% {
  sampled_names <- sample(sample_names, 10,replace = FALSE)
  #otu_table(sampled_physeq )%>%data.frame()
  }
  sub_green[[i]]=result
}




# bind all the community data; the initial results were based on one time of simulation

vege_com=dplyr::bind_rows(sub_green[[1]],sub_green[[2]],sub_green[[3]],sub_green[[4]],sub_green[[5]],sub_green[[6]],sub_green[[7]],sub_green[[8]])%>%mutate(type=rep(type0,each=50))

ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
plot(ordination, type = "points")
com_score=ordination$species
com_score=ordination$points%>%data.frame()%>%mutate(type=rep(type0,each=50))

# change the 0-1 data

vege_com[vege_com>0]=1

adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)







## if we repeat the analysis 100 times,

pva=numeric()
fva[j]=numeric()
for (j in 1:100)
  {
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
sub_green=list()
for (i in 1:8)
{
  
  dk=subset_samples(green_sample,type==type0[i])
  
  sample_names=sample_names(dk)
  
  sampled_names <- sample(sample_names, 50,replace = FALSE)
  sampled_physeq <- prune_samples(sampled_names, dk)
  
  sub_green[[i]]=otu_table(sampled_physeq )%>%data.frame()
}

# bind all the community data; the initial results were based on one time of simulation

vege_com=dplyr::bind_rows(sub_green[[1]],sub_green[[2]],sub_green[[3]],sub_green[[4]],sub_green[[5]],sub_green[[6]],sub_green[[7]],sub_green[[8]])%>%mutate(type=rep(type0,each=50))


vege_com[vege_com>0]=1

adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)

pva[j]=adonis_result$aov.tab[6]$`Pr(>F)`[1]# 
fva[j]=adonis_result$aov.tab[4]$F.Model
}

save(green_sample,file="green_sample.RData")

# try to examine the indication analysis


vege_com[vege_com>0]=1

vege_com_trans=t(vege_com[1:67769,])%>%data.frame()




#(2)# compare the difference in the mean richness among land cover types based on the true plot-level measurements rather than virtual communities

load("~/soil-sar/plot-level-diversity/plot_diversity_env.RData")
compare_richness=plot_diversity_env[,c(1,2,3,5)]%>%left_join(vege_tep)%>%mutate(type = ifelse(nchar(plotID)<4, "evergreenForest", type))%>%
  mutate(site=substr(plotID,1,4))%>%filter(!is.na(type))%>%filter(type%in%type0)

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


### compare the beta diversity based on plot-level measurements

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


## for the plot level gamma diversity 

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


### compare the species composition with pcoa analysis

compare_richness%>%group_by(type)%>%summarize(count=n())

# if we just select the same number of plots for comparison
#
data_matrix=compare_richness%>%dplyr::select(core_rich, plot_rich, beta_mean )
# if we just look at the neon sites

data_matrix=compare_richness%>%filter(nchar(plotID)>2)%>%dplyr::select(core_rich, plot_rich, beta_mean )

compare_richness%>%filter(nchar(plotID)>2)%>%dplyr::select(core_rich, plot_rich, beta_mean,type )->tem_data


data_matrix=apply(data_matrix,2,range01)

dist_matrix <- vegdist(data_matrix,method = "bray")

pcoa_result <- cmdscale(dist_matrix, eig=TRUE, k=2)  # k=2 for 2D PCoA

grouping_factor <- tem_data$type

perm_result <- adonis2(dist_matrix ~ grouping_factor, data = tem_data, permutations = 999)


# Extracting the coordinates
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

pcoa_coords%>%mutate(type=tem_data$type)->pcoa_coords

ggplot(data=pcoa_coords,aes(x=PCoA1,y=PCoA2,color=type))+
  geom_point()+
  stat_ellipse(linewidth=0.6,linetype="solid")


##

p3=ggplot(data=sensitivity_guild,aes(x=guild,y=sensitivity))+
   geom_bar(stat = "identity",width=0.5)+
 scale_x_discrete(breaks=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"),
                      labels=c("All","AM","EM","Epiphyte","Litter sprotroph","Parasite","Plant pathogen","Soil saprotroph","Wood saprotroph"))+
  theme(legend.position = c(0.8,0.87),
           legend.text = element_text(size=8),
           legend.title  = element_text(size=10),
           text = element_text(size = 18),
           plot.title = element_text(size = 15, hjust = 0.5),
           axis.text.y = element_text(hjust = 0),
           axis.text.x = element_text(hjust = 1,angle=90),
           axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
           axis.ticks.x = element_blank(),
           legend.key.size = unit(0.3, "cm"),
           panel.background = element_rect(fill = "NA"),
           panel.border = element_rect(color = "black", size = 1, fill = NA))+
   geom_hline(yintercept = 1,color="red",linetype="dashed",size=1.1)+
   ylab("Response ratio")+
   xlab("")

p1=ggplotGrob(p1)
p2=ggplotGrob(p2)
p3=ggplotGrob(p3)

p1$heights=p2$heights
p1$heights=p3$heights
