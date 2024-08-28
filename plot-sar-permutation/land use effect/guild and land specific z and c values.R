## determining the habitat affinity for each guild and land use type


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

plot_diversity_env_land%>%dplyr::select(plotID,type)%>%distinct()->land

land%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))->land

land%>%filter(!is.na(type))->land#(if we do not filter out the NAs, the join data will be longer)

temp%>%left_join(land,by="plotID")%>%dplyr::select(plotID,type)%>%dplyr::select(type)->temp


# for na rows assign that with "evergreen forest"
row.names(temp)=row.names(sample_data(d))

temp=sample_data(temp)

d<- merge_phyloseq(rare_all, temp)#

# to see the number of different soil cores among different land use types
sample_data(d)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%summarize(count=n())

sample_data(d)%>%data.frame()%>%dplyr::select(plotID,type)


# cultivatedCrops              397
# deciduousForest             1445
# dwarfScrub                    12
# emergentHerbaceousWetlands    46
# evergreenForest             2066
# grasslandHerbaceous         1018
# mixedForest                  342
# sedgeHerbaceous               14
# shrubScrub                   556
# woodyWetlands                312
# NA                           170#some plots still do not have land use data

# get the same number of cores and determine the mean value of core-level richness

type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","shrubScrub","woodyWetlands")

# does not consider differences in the sampling season
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  
  richness=numeric()
  for (i in 1:7)
  {
    
    dk=subset_samples(d,type==type0[i])
    
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300,replace = FALSE)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich[i,]=kk[[i]]
  
}

land_rich%>%data.frame()%>%rename_all(~paste0(c(type0)))->d

land_rich=melt(d)

ggboxplot(land_rich, x = "variable", y = "value", fill = "variable", outlier.shape = NA) +
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









# compare the z and c for specific guilds

plot_diversity_env_land%>%dplyr::select(plotID,type)->temp

# combined the pasture hay and cultivated croplands

temp%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))->temp


            
  model_data_SAR%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%filter(!is.na(type))%>%
    group_by(guild,type)%>%summarise(mean_zvalue=mean(zvalue,na.rm=TRUE))
  
  
  model_data_SAR%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%filter(!is.na(type))%>%
    group_by(guild,type)%>%summarise(mean_cvalue=mean(logc,na.rm=TRUE))%>%mutate(c=2.71828^mean_cvalue)
  