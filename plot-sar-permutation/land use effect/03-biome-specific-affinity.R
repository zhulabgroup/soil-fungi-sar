
# get the habitat affinity within the same biome
# to extract the biome data fro each plot
library(dplyr)
library(raster)
library(phyloseq)
library(sp)
library(sf)
library(terra)
library(tidyr)

setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")


biomes <- st_read("wwf_terr_ecos.shp")

biomes <- group_by(biomes, BIOME) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup() %>%
  mutate(BIOME = as.character(BIOME))

biome_labels <- data.frame(
  BIOME = as.character(c(seq(1, 14), 98, 99)),
  label = c(
    "Tropical & Subtropical Moist Broadleaf Forests",
    "Tropical & Subtropical Dry Broadleaf Forests",
    "Tropical & Subtropical Coniferous Forests",
    "Temperate Broadleaf & Mixed Forests",
    "Temperate Conifer Forests",
    "Boreal Forests/Taiga",
    "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    "Temperate Grasslands, Savannas & Shrublands",
    "Flooded Grasslands & Savannas",
    "Montane Grasslands & Shrublands",
    "Tundra",
    "Mediterranean Forests, Woodlands & Scrub",
    "Deserts & Xeric Shrublands",
    "Mangroves",
    "Undefined",
    "Undefined2"
  ),
  stringsAsFactors = FALSE
)

biomes$LABEL <- biome_labels$label[match(biomes$BIOME, biome_labels$BIOME)]
cbind(biomes$BIOME, biomes$LABEL)

biomes <- st_crop(biomes, c(xmin=-170,xmax=-55,ymin=17,ymax=72))



# just create an empty raster to save the values.
r <- rast(ext(biomes),resolution = 0.17,   # the resolution of your targeted raster
          crs = "EPSG:4326")
#r <- rasterize(biomes, r, field = "BIOME")  # 'field' can be a column name or a constant value

r <- rasterize(biomes, r, field = "LABEL")  # 'field' can be a column name or a constant value



load("rare_all_guild.RData")

full_parameter_data <- readRDS("full_parameter_data.rds")

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
rare_all_guild<- merge_phyloseq(rare_all_guild, plotIDM)# merge the new plotid with the initial data 


sample_data(rare_all_guild)%>%data.frame()%>%dplyr::select(plotIDM,lon,lat)%>%group_by(plotIDM)%>%
  summarise(lon=mean(lon),lat=mean(lat))%>%rename(plotID=plotIDM)->plot_coordinates




#full_parameter_data%>%left_join(plot_coordinates,by="plotID")->full_parameter_data

#plot_biomes=full_parameter_data%>%dplyr::select(lon,lat,plotID)%>%distinct()

extract_biomes=terra::extract(r,plot_coordinates[,c("lon","lat")])

plot_coordinates%>%bind_cols(extract_biomes%>%dplyr::select(LABEL))%>%dplyr::select(plotID,LABEL)->plot_biomes

# some plots do not have biome information 

# to see the plot level land cover types


load("~/soil-sar/data/comp_vege.RData")


sample_data(rare_all_guild)%>%data.frame()%>%dplyr::select(plotIDM)%>%
  rename(plotID=plotIDM)->temp

#add the rownames for the 

#temp%>%mutate(name=rownames(temp))->temp

#load in the updated land use data

land=read.csv("All_NEON_TOS_Plot_Centroids_V11.csv")

land%>%dplyr::select(plotID,nlcdClass)->land
land=distinct(land)

# the plot TREE_015 was defined with two land cover types

land=land[-1882,]#remove the dumplicated 

left_join(temp,land,by="plotID")%>%dplyr::select(plotID,nlcdClass)->temp

temp%>%mutate(type = replace_na(nlcdClass, "evergreenForest"))%>%data.frame()->land

#get the biome of the plots

land%>%left_join(plot_biomes,by="plotID")%>%
  mutate(joint_land=paste(type,"_",LABEL))%>%
  dplyr::select(joint_land, LABEL)->land_biome

# assign the NA rows with "evergreen forest" and then combine that with the phyloseq object

row.names(land_biome)=row.names(sample_data(d))
land_biome=sample_data(land_biome)
rare_all_guild_biome<- merge_phyloseq(rare_all_guild, land_biome)#


saveRDS(rare_all_guild_biome,"rare_all_guild_biome.rds")

# to see the number of different soil cores among different land use types
sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%dplyr::summarize(count=n())

sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(joint_land)%>%group_by(joint_land)%>%dplyr::summarize(count=n())


# to compare the richness between the natural and croplands

# where the we have different biomes

#cultivatedCrops _ Temperate Broadleaf & Mixed Forests                       33
#cultivatedCrops _ Temperate Conifer Forests                                 10
#cultivatedCrops _ Temperate Grasslands, Savannas & Shrublands              147
#cultivatedCrops _ Tropical & Subtropical Moist Broadleaf Forests            11

dk=subset_samples(d,type=="cultivatedCrops")
# need to compare the richness

#the distance between the plots 

distance_matrix <- distm(plot_coordinates[, c("lon", "lat")], fun = distHaversine)

rownames(distance_matrix) <- plot_coordinates$plotID
colnames(distance_matrix) <- plot_coordinates$plotID
diag(distance_matrix) <- NA

nearest_plot <- apply(distance_matrix, 1, function(x) {
  nearest_index <- which.min(x)  # Find the index of the minimum distance
  return(colnames(distance_matrix)[nearest_index])  # Return the corresponding plot name
})


nearest_df <- data.frame(
  plot = rownames(distance_matrix),
  nearest_plot = nearest_plot
)

# determine the sites distance

sample_data(d)%>%data.frame()%>%dplyr::select(Site,lon,lat)%>%
  group_by(Site)%>%summarise(lon=mean(lon,rm.na=TRUE),lat=mean(lat,rm.na=TRUE))->site_coordinates

distance_matrix <- distm(site_coordinates[, c("lon", "lat")], fun = distHaversine)

rownames(distance_matrix) <- site_coordinates$Site
colnames(distance_matrix) <- site_coordinates$Site
diag(distance_matrix) <- NA

nearest_site <- apply(distance_matrix, 1, function(x) {
  nearest_index <- which.min(x)  # Find the index of the minimum distance
  return(colnames(distance_matrix)[nearest_index])  # Return the corresponding plot name
})


nearest_df_site <- data.frame(
  site = rownames(distance_matrix),
  nearest_site = nearest_site
)


data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")

data_EM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild_biome, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "epiphyte")

biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
         "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

#to see how many soil samples are there in the natural and human-dominated communities

core_number=matrix(ncol=2,nrow=4)
site=list()
for (j in 1:4)
  
  if(j%in%c(1,2,4))
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) 
  dk=subset_samples(rare_all_guild_biome,LABEL==biomes[j])# select the biome of interest
  # to see which are croplands
  subset_samples(dk,type=="cultivatedCrops")->modified_data
  #soil samples
  
  dim(sample_data(modified_data)%>%data.frame())[1]->n_crop_sample
  subset_samples(dk,type!="cultivatedCrops")->natural_data
  dim(sample_data(natural_data)%>%data.frame())[1]->n_nature_sample
  sample_data(modified_data)%>%data.frame()%>%pull(plotIDM)%>%
    substr(start = 1, stop = 4)%>%unique()->modified_site
  
  sample_data(natural_data)%>%data.frame()%>%pull(plotIDM)%>%
    substr(start = 1, stop = 4)%>%unique()->natural_site
  overlap <- intersect(modified_site, natural_site)
  core_number[j,1]=n_crop_sample
  core_number[j,2]=n_nature_sample
  site[[j]]=list(modified_site,natural_site)
}

# if we just based on the site-specific comparison for the richness 

hehe=list()

for(j in 1:4)
  
{

dk=subset_samples(rare_all_guild_biome,LABEL==biomes[j])# select the biome of interest
# to see which are croplands
subset_samples(dk,type=="cultivatedCrops")->modified_data
#soil samples
dim(sample_data(modified_data)%>%data.frame())[1]->n_crop_sample
subset_samples(dk,type!="cultivatedCrops")->natural_data
dim(sample_data(natural_data)%>%data.frame())[1]->n_nature_sample
sample_data(modified_data)%>%data.frame()%>%pull(plotIDM)%>%
  substr(start = 1, stop = 4)%>%unique()->modified_site
sample_data(natural_data)%>%data.frame()%>%pull(plotIDM)%>%
  substr(start = 1, stop = 4)%>%unique()->natural_site

overlap <- intersect(modified_site, natural_site)

site_number=length(overlap)


  times=5
  
  if(length(overlap)>=1)
  {

site_response=matrix(ncol=site_number,nrow=times)#calculate the response ratio for
 
for(m in 1:site_number)
 {
  # when there is overlap between these two data sets
 
  set.seed(234)
  # based on the selected site names to select the samples
  subset_samples(natural_data,Site%in%overlap[m])->sub_natural_data
  subset_samples(modified_data,Site%in%overlap[m])->sub_modified_data
  n_sub_natural_samples=nsamples(sub_natural_data)
  n_sub_modified_samples=nsamples(sub_modified_data)

  if(n_sub_natural_samples>=n_sub_modified_samples)
  {
    set.seed(123)
    sample_names=sample_names(sub_natural_data)
    sampled_names <- replicate(times,sample(sample_names, n_sub_modified_samples,replace = FALSE))
    #richness_natural=numeric()
    
    richness_natural=numeric()
    for (i in 1:times)
    {
      # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) 
      sampled_physeq <- prune_samples(sampled_names[,i], sub_natural_data)
      richness_natural[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%
        summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
        as.numeric()
    }
    
    richness_modified=estimate_richness(sub_modified_data, measures = "Observed")%>%
      summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
      as.numeric()
    site_response[,m]=richness_modified/richness_natural
  }
  
else
  {
  set.seed(123)
  sample_names=sample_names(sub_modified_data)
  sampled_names <- replicate(times,sample(sample_names, n_sub_natural_samples,replace = FALSE))
  richness_modified=numeric()
  for (i in 1:times)
  {
    # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) 
    sampled_physeq <- prune_samples(sampled_names[,i], sub_modified_data)
    richness_modified[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%
      summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
      as.numeric()
  }

richness_natural=estimate_richness(sub_natural_data, measures = "Observed")%>%
  summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
  as.numeric()
# the ratio between the natural and modified communities
site_response[,m]=richness_modified/richness_natural
  }
  
}#wit the for

hehe[[j]]=site_response
}

  else{
    hehe[[j]]="hehe"
  }

}












mm=list()
for (k in 1:9)
{
  
response_ratio=list()
for (j in 1:4)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) 
dk=subset_samples(get(data[k]),LABEL==biomes[j])# select the biome of interest
# to see which are croplands
subset_samples(dk,type=="cultivatedCrops")->modified_data
#soil samples
dim(sample_data(modified_data)%>%data.frame())[1]->n_crop_sample

subset_samples(dk,type!="cultivatedCrops")->natural_data
dim(sample_data(natural_data)%>%data.frame())[1]->n_nature_sample


sample_data(modified_data)%>%data.frame()%>%pull(plotIDM)%>%
substr(start = 1, stop = 4)%>%unique()->modified_site
sample_data(natural_data)%>%data.frame()%>%pull(plotIDM)%>%
  substr(start = 1, stop = 4)%>%unique()->natural_site
overlap <- intersect(modified_site, natural_site)
if(length(overlap)>1)
  {
# when there is overlap between these two data sets
    set.seed(234)
  # based on the selected site names to select the samples
  subset_samples(natural_data,Site%in%overlap)->sub_natural_data
    sample_names=sample_names(sub_natural_data)
    times=100
    sampled_names <- replicate(times,sample(sample_names, n_crop_sample,replace = FALSE))
    richness_natural=numeric()
    for (i in 1:times)
    {
     # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) 
    sampled_physeq <- prune_samples(sampled_names[,i], sub_natural_data)
    richness_natural[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%
      summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
      as.numeric()
    }
    
richness_modified=estimate_richness(modified_data, measures = "Observed")%>%
  summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
  as.numeric()
# the ratio between the natural and modified communities
response_ratio[[j]]=richness_modified/richness_natural
}

else if(length(overlap)==1)
{
  # when there is overlap between these two data sets
  set.seed(234)
  
  # based on the selected site names to select the samples
  
  subset_samples(modified_data,Site%in%overlap)->sub_modified_data
  
  #more samples in the human modified sites
  
  sample_names=sample_names(sub_modified_data)
  
  times=100
  
  sampled_names <- replicate(times,sample(sample_names, n_nature_sample,replace = FALSE))
  
  richness_modified=numeric()
  for (i in 1:times)
  {
    # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) 
    
    sampled_physeq <- prune_samples(sampled_names[,i], sub_modified_data)
    
    richness_modified[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%
      summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
      as.numeric()
  }
  
  richness_natural=estimate_richness(natural_data, measures = "Observed")%>%
    summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
    as.numeric()
  # the ratio between the natural and modified communities
  response_ratio[[j]]=richness_modified/richness_natural
}

else if (length(overlap)<1)
  
  {
  nearest_df_site%>%filter(site%in%modified_site)%>%pull(nearest_site)->near_sites
  
  subset_samples(natural_data,Site%in%near_sites)->sub_natural_data
  set.seed(234)
  sample_names=sample_names(sub_natural_data)
  times=100
  sampled_names <- replicate(times,sample(sample_names, n_crop_sample,replace = FALSE))
  
  richness_natural=numeric()
  for (i in 1:times)
  {
    # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) 
    sampled_physeq <- prune_samples(sampled_names[,i], sub_natural_data)
    richness_natural[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%
      summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
      as.numeric()
  }
  
  richness_modified=estimate_richness(modified_data, measures = "Observed")%>%
    summarize(mean_A = mean(Observed, na.rm = TRUE))%>%
    as.numeric()
  # the ratio between the natural and modified communities
  response_ratio[[j]]=richness_modified/richness_natural
}
  
}

mm[[k]]=response_ratio
}

# to get the response ratio of different guilds
# get the mean of the response ration
kkm=matrix(ncol=4,nrow=9)
for(j in 1:9){
  
kk=numeric()
for (i in 1:4)
  {
  kk[i]=mm[[j]][[i]]%>%mean()
  
}
kkm[j,]=kk
}

kkm%>%data.frame()%>%
  mutate(guild=c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para"))%>%
  rename_all(~paste0(c(biomes,"guild")))%>%melt()%>%
  rename(biome=variable,response=value)%>%
  mutate(guild_biome=paste(guild,"_",biome))->guild_specific_response_ratio

# to combine the response ratio with the c and z values

full_parameter_data%>%left_join(land%>%distinct(),by="plotID")%>%
  left_join(plot_biomes)%>%group_by(guild,LABEL)%>%
  summarise(mean_zvalue=mean(zvalue,na.rm=TRUE),mean_cvalue=mean(2.71828^logc,na.rm=TRUE))%>%
  rename(biome=LABEL)%>%mutate(guild_biome=paste(guild,"_",biome))%>%
  dplyr::select( mean_zvalue, mean_cvalue,guild_biome)%>%
left_join(guild_specific_response_ratio%>%dplyr::select(response,guild_biome),by="guild_biome")%>%
  mutate(affinity=response^(1/mean_zvalue))->full_guild_affinity_data

saveRDS(full_guild_affinity_data,file="full_guild_affinity_data.rds")

full_guild_affinity_data%>%filter(affinity<1)

# to test the difference in the affinity among different land use types

ratio=data.frame(ncol=3,nrow=4)
for (i in 1:4)
  {
  ratio[i,1]=t.test(response_ratio[[i]],mu=1,alternative = "two.side")$estimate
  ratio[i,2:3]=t.test(response_ratio[[i]],mu=1,alternative = "two.side")[4]$conf.int[1:2]
}

ratio%>%rename_all(~paste0(c("mean","low","up")))->ratio

ggplot()+
  geom_point(data=ratio,aes(x=1:4,y=mean),size=3)+
  
  geom_segment(data=ratio,aes(x=1:4,y=low,xend=1:4,yend=up),color="red")+
  geom_hline(yintercept = 1,linetype="dashed",color="red")+
  ylab("response ratio")+
  scale_x_discrete(limits=1:4,labels=biomes)+
  theme(legend.position = c(0.25,0.38),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size = 12), 
        axis.text.x = element_text(hjust = 1,size = 12,angle=80), 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))
  

# non sites available in the natural data, we need to find the neighboring with the same
# determine the distance between them




