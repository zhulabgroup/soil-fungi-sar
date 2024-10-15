#
library(phyloseq)
library(reshape2)
library(stringr)
library(ggplot2)
library(terra)


## the function to determine the guild-specific affinity
# this function does not work, it returns with the same
my_function_guild_level_richness=function(data,i,j)
  {
  set.seed(123)
  kk=list()
  for(j in 1:50)
  {
    cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
    type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","shrubScrub","woodyWetlands")
    
    richness=numeric()
    for (i in 1:7)
    {
      dk=subset_samples(data,type==type0[i])
      sample_names=sample_names(dk)
      sampled_names <- sample(sample_names, 30)
      sampled_physeq <- prune_samples(sampled_names, dk)
      richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
    } 
    kk[[j]]=richness
  }
  
  df=matrix(ncol=7,nrow=500)
  for (i in 1:500)
  {
    df[i,]=kk[[i]]
  }
  df%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->df
  
return(df)
}


save(land_rich_all_updated,file="land_rich_all_updated.RData")



#(1)# determining the habitat affinity for each fungal guild and land use type

#add the fungal guild data to the full sample data
ft<- read_csv("~/soil-sar/plot-sar-permutation/FungalTraits_1.2_ver_16Dec_2020.csv")
ft%>%data.frame()%>%rename(genus=GENUS)%>%select(genus,primary_lifestyle)->ft_temp
tax_table(rare_all)%>%data.frame()%>%left_join(ft_temp,by="genus")->guild_temp
taxa_matrix <- as.matrix(guild_temp)
new_taxa_table <- tax_table(taxa_matrix)
row.names(new_taxa_table )=row.names(tax_table(rare_all))
tax_table(rare_all) <- new_taxa_table
rare_all_guild=rare_all

save(rare_all_guild,file="rare_all_guild.RData")

#(2) add the land use type data to the whole dataset

load("~/soil-sar/plot-sar-permutation/plot_diversity_env_land.RData")

d=sample_data(rare_all_guild)
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
d<- merge_phyloseq(rare_all_guild, plotIDM)# merge the new plotid with the initial data 

sample_data(d)%>%data.frame()%>%dplyr::select(plotIDM)%>%rename(plotID=plotIDM)->temp
plot_diversity_env_land%>%dplyr::select(plotID,type)%>%distinct()%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))%>%
  filter(!is.na(type))->land#(if we do not filter out the NAs, the join data will be longer)

temp%>%left_join(land,by="plotID")%>%dplyr::select(plotID,type)%>%dplyr::select(type)->temp


#assign NA rows with "evergreen forest"
row.names(temp)=row.names(sample_data(d))
temp=sample_data(temp)
d<- merge_phyloseq(rare_all_guild, temp)#
# to see the number of different soil cores among different land use types
sample_data(d)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%summarize(count=n())

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

rare_all_guild=d
save(rare_all_guild,file="rare_all_guild.RData")

plot_diversity_env_land%>%dplyr::select(plotID,type)->temp

# combined the pasture hay and cultivated croplands

temp%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))%>%filter(!is.na(type))->temp


#(3) get the guild-specific c and z values 
#note, some plots have large c values

model_data_SAR%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%
  group_by(guild)%>%summarise(mean_cvalue=mean(logc,na.rm=TRUE))%>%mutate(c=2.71828^mean_cvalue)->c_temp

#guild     mean_cvalue      c
#1 AM             -3.27  0.0382
#2 EM             -0.696 0.499 
#3 all             1.88  6.52  
#4 epiphy         -3.02  0.0490
#5 littersap      -1.09  0.336 
#6 para           -1.40  0.246 
#7 plapat         -1.60  0.201 
#8 soilsap         0.324 1.38  
#9 woodsap        -1.89  0.152 

model_data_SAR%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%filter(!is.na(type))%>%
  group_by(guild)%>%summarise(mean_zvalue=mean(zvalue,na.rm=TRUE))->z_temp

#guild     mean_zvalue
#1 AM              0.763
#2 EM              0.754
#3 all             0.706
#4 epiphy          0.730
#5 littersap       0.717
#6 para            0.654
#7 plapat          0.693
#8 soilsap         0.659
#9 woodsap         0.741


#(4)# get the guild-specific habitat affinity with the formula of affinity=(S1/S2)^1/z
# where S1 and S2 indicate the core-level of fungal richness in the human-modified and the natural habitats, respectively.
# get the same number of cores from each land use type and determined the mean value of core-level richness
#####
# split the data into different guilds

data_EM <- subset_taxa(rare_all_guild, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild, primary_lifestyle == "epiphyte")

### write a function to manipulate all the data sets


###
type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","shrubScrub","woodyWetlands")

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(rare_all_guild,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_all_updated=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_all_updated[i,]=kk[[i]]
}

land_rich_all_updated%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_all_updated

save(land_rich_all_updated,file="land_rich_all_updated.RData")
####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_plapat,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_plapat=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_plapat[i,]=kk[[i]]
}

land_rich_plapat%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_plapat_updated

save(land_rich_plapat_updated,file="land_rich_plapat_updated.RData")

###
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_soilsap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_soilsap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_soilsap[i,]=kk[[i]]
}

land_rich_soilsap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_soilsap_updated

save(land_rich_soilsap_updated,file="land_rich_soilsap_updated.RData")
#####
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_littersap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_littersap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_littersap[i,]=kk[[i]]
}

land_rich_littersap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_littersap_updated

save(land_rich_littersap_updated,file="land_rich_littersap_updated.RData")

##

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_woodsap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_woodsap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_woodsap[i,]=kk[[i]]
}

land_rich_woodsap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_woodsap_updated

save(land_rich_woodsap_updated,file="land_rich_woodsap_updated.RData")
###

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_para,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_para=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_para[i,]=kk[[i]]
}

land_rich_para%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_para_updated

save(land_rich_para_updated,file="land_rich_para_updated.RData")
####
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_epiphy,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_epiphy=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_epiphy[i,]=kk[[i]]
}

land_rich_epiphy%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_epiphy_updated

save(land_rich_epiphy_updated,file="land_rich_epiphy_updated.RData")




####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_EM,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_EM=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_EM[i,]=kk[[i]]
}

land_rich_EM%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_EM_updated

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

save(land_rich_EM_updated,file="land_rich_EM_updated.RData")


####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_AM,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_AM=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_AM[i,]=kk[[i]]
}

land_rich_AM%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_AM_updated

save(land_rich_AM_updated,file="land_rich_AM_updated.RData")


####

land_rich_AM_updated=richness_among_land(data_AM)

land_rich_soilsap_updated=richness_among_land(data_soilsap)
land_rich_woodsap_updated=richness_among_land(data_woodsap)
land_rich_littersap_updated=richness_among_land(data_littersap)
land_rich_plapat_updated=richness_among_land(data_plapat)
land_rich_para_updated=richness_among_land(data_para)
land_rich_epiphy_updated=richness_among_land(data_epiphy)
##save all the results
save(land_rich_AM_updated,file="land_rich_AM_updated.RData")
save(land_rich_soilsap_updated,file="land_rich_soilsap_updated.RData")
save(land_rich_woodsap_updated,file="land_rich_woodsap_updated.RData")
save(land_rich_littersap_updated,file="land_rich_littersap_updated.RData")
save(land_rich_plapat_updated,file="land_rich_plapat_updated.RData")
save(land_rich_para_updated,file="land_rich_para_updated.RData")
save(land_rich_epiphy_updated,file="land_rich_epiphy_updated.RData")



mean_richness_guild=bind_rows(land_rich_AM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_EM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_soilsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_littersap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_woodsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_plapat_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_para_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_epiphy_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_all_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)))%>%
  mutate(guild=rep(c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all"),each=7))%>%data.frame()

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

affinity=numeric()
for (i in 1:9)
{
  df=mean_richness_guild%>%filter(guild==guild_type[i])
  affinity[i]=df[1,2]/df[2:7,2] %>%mean()
}

# if we test the significance of the response

reponse_significance=numeric()
for (i in 1:9)
{
  df=mean_richness_guild%>%filter(guild==guild_type[i])
  response=df[1,2]/df[2:7,2]
  
  reponse_significance[i]=t.test(response,mu=1,alternative = "two.sided")[[3]]
  
}

reponse_significance%>%data.frame()%>%bind_cols(guild_type)%>%
  rename_all(~paste0(c("pvalue","guild")))
# if we based on the overl richness among different guilds to get the ratio

my_function_response_ratio=function(data)
{
  data%>%mutate(type=if_else(!variable%in% c("cultivatedCrops"),"modified","Natural"))->d
  # to compare the means among the modified and the natural natural communities
  #mod=aov(value~type,data=d)
  d%>%group_by(type)%>%summarise(mean_value=mean(value,na.rm=TRUE))->d
  
  ratio=d[2,2]/d[1,2]
  return(ratio)
}

#difference in the ratio among species
ratio=numeric()

for(i in 1:9)
  {
  df=mean_richness_guild%>%filter(guild==guild_type[i])
  ratio[i]=my_function_response_ratio(df)%>%pull()
  
}

# test the significance of the richness values based on the broad classification of the difference

richness_guild=bind_rows(land_rich_AM_updated,
                              land_rich_EM_updated,
                              land_rich_soilsap_updated,
                              land_rich_littersap_updated,
                              land_rich_woodsap_updated,
                              land_rich_plapat_updated,
                              land_rich_para_updated,
                              land_rich_epiphy_updated,
                              land_rich_all_updated)%>%
mutate(guild=rep(c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all"),each=3500))%>%data.frame()


guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")




ratio=numeric()
for(i in 1:9)
{
  df=richness_guild%>%filter(guild==guild_type[i])
  df%>%mutate(type=if_else(!variable%in% c("cultivatedCrops"),"modified","Natural"))->d
  #leveneTest(zvalue ~ type, data = type_mean)
  #k=oneway.test(value ~ type, data = d, na.action = na.omit, var.equal = FALSE)
  k=tukey(d$value, d$type, method = "G")
  ratio[i]=k$result1[1,2]/k$result1[2,2]# the richness ratio between the modified and native ones
  #ratio[i]=k$p.value
}

#create the box plot for the mean richness

pp_affinity=list()
for(i in 1:9)
{
  df=richness_guild%>%filter(guild==guild_type[i])
  df%>%mutate(type=if_else(!variable%in% c("cultivatedCrops"),"Modified","Natural"))->d
  
  pp_affinity[[i]]=ggboxplot(d, x = "type", y = "value", fill = "type", outlier.shape = NA) +
    xlab("") +
    ylab("Richness") +
    theme(legend.position = c(0.25, 0.7558), 
          legend.text = element_text(size = 14), 
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          text = element_text(size = 15),
          axis.title.y = element_text( size = 20), 
          axis.title.x = element_text(size = 20)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # make four rows for the legend
    xlab("Land cover types")+
    ggtitle(paste0(d$guild)%>%unique())
  pp_affinity[[i]]=ggplotGrob(pp_affinity[[i]])
}

plot_grid(pp_affinity[[1]],
          pp_affinity[[2]],
          pp_affinity[[3]],
          pp_affinity[[4]],
          pp_affinity[[5]],
          pp_affinity[[6]],
          pp_affinity[[7]],
          pp_affinity[[8]],
          pp_affinity[[9]],
          ncol=3)

pp_affinity[[1]]$widths=pp_affinity[[6]]$widths

pp_affinity[[1]]$widths=pp_affinity[[9]]$widths
pp_affinity[[1]]$widths=pp_affinity[[2]]$widths

pp_affinity[[9]]$widths=pp_affinity[[2]]$widths


plot_grid(pp_affinity[[1]],
          pp_affinity[[2]],
          
          pp_affinity[[6]],
          
          pp_affinity[[9]],
          ncol=2)




affinity%>%data.frame()%>%bind_cols(guild_type)%>%data.frame()%>%rename_all(~paste0(c("sensitivity","guild")))%>%
  left_join(z_temp,by="guild")%>%mutate(affinity=sensitivity^(1/mean_zvalue))%>%left_join(c_temp,by="guild")->parameters

#sensitivity     guild mean_zvalue  affinity mean_cvalue          c
#1   1.1558572        AM   0.7627419 1.2091247  -3.2654048 0.03818156
#2   0.6982490        EM   0.7541218 0.6210843  -0.6956156 0.49876755
#3   0.6899942   soilsap   0.6588153 0.5693605   0.3238652 1.38246063
#4   1.1260563 littersap   0.7171882 1.1800272  -1.0910826 0.33585294
#5   0.9118643   woodsap   0.7414621 0.8829954  -1.8857109 0.15172135
#6   1.1898728    plapat   0.6928434 1.2852038  -1.6036925 0.20115262
#7   0.8387798      para   0.6540632 0.7643021  -1.4019191 0.24612441
#8   0.6773579    epiphy   0.7304279 0.5866511  -3.0156364 0.04901473
#9   0.8797784       all   0.7064055 0.8341686   1.8750299 6.52100604

###

#parameters
#sensitivity     guild mean_zvalue  affinity mean_cvalue           c
#1   1.2114232        AM   0.7627419 1.2858961  0.07015474    1.072674
#2   0.7160444        EM   0.7541218 0.6421608  0.60695776    1.834840
#3   0.7202833   soilsap   0.6588153 0.6077259  1.72830766    5.631110
#4   1.1929157 littersap   0.7171882 1.2788499  0.46743711    1.595898
#5   0.9605541   woodsap   0.7414621 0.9471689  0.19619948    1.216769
#6   1.2578671    plapat   0.6928434 1.3925333  0.37924158    1.461176
#7   0.8884922      para   0.6540632 0.8346345  0.37564831    1.455935
#8   0.7177464    epiphy   0.7304279 0.6350613  0.07984360    1.083118
#9   0.9266827       all   0.7064055 0.8978153  7.65814512 2117.814570

save(parameters,file="parameters.RData")

