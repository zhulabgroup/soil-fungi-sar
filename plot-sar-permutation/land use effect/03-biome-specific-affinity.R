
# get the habitat affinity within the same biome
# to extract the biome data fro each plot
library(dplyr)
library(raster)
library(phyloseq)
library(sp)
library(sf)
library(terra)
library(tidyr)
library(geosphere)

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


#load("~/soil-sar/data/comp_vege.RData")


sample_data(rare_all_guild)%>%data.frame()%>%dplyr::select(plotIDM)%>%
  rename(plotID=plotIDM)->temp

#add the rownames for the 

#temp%>%mutate(name=rownames(temp))->temp

#load in the updated land use data

land=read.csv("All_NEON_TOS_Plot_Centroids_V11.csv")

land%>%dplyr::select(plotID,nlcdClass)->land
land=distinct(land)

#not sure why some plot still do not have landcover data need to consult the neon falculty
# the plot TREE_015 was defined with two land cover types

land=land[-1882,]#remove the dumplicated 

left_join(temp,land,by="plotID")%>%dplyr::select(plotID,nlcdClass)->temp

temp%>%mutate(type = replace_na(nlcdClass, "evergreenForest"))%>%data.frame()->land

#get the biome of the plots

land%>%left_join(plot_biomes,by="plotID")%>%
  mutate(joint_land=paste(type,"_",LABEL))%>%
  dplyr::select(type,joint_land, LABEL)->land_biome

# replace the 
land_biome$type=gsub("pastureHay","cultivatedCrops",land_biome$type)

# assign the NA rows with "evergreen forest" and then combine that with the phyloseq object


sample_data_df <- as.data.frame(sample_data(rare_all_guild))

# Delete the column you want (replace "column_name" with the actual column name)
sample_data_df$type <- NULL

# Assign the modified sample data back to the phyloseq object
sample_data(rare_all_guild) <- sample_data(sample_data_df)



row.names(land_biome)=row.names(sample_data(d))
land_biome=sample_data(land_biome)
rare_all_guild_biome<- merge_phyloseq(rare_all_guild, land_biome)#


saveRDS(rare_all_guild_biome,"rare_all_guild_biome.rds")

rare_all_guild_biome=readRDS("rare_all_guild_biome.rds")

# to see the number of different soil cores among different land use types
sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%dplyr::summarize(count=n())

sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(joint_land)%>%group_by(joint_land)%>%dplyr::summarize(count=n())


# to compare the richness between the natural and croplands

# where the we have different biomes

#cultivatedCrops _ Temperate Broadleaf & Mixed Forests                       33
#cultivatedCrops _ Temperate Conifer Forests                                 10
#cultivatedCrops _ Temperate Grasslands, Savannas & Shrublands              147
#cultivatedCrops _ Tropical & Subtropical Moist Broadleaf Forests            11


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
# and then get the mean for the site level response ratio

all_data=list()
for (n in 1:9)
{
  
data_one_guild=list()
for(j in 1:4)
{
dk=subset_samples(get(data[n]),LABEL==biomes[j])# select the biome of interest
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
  times=500
  
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
  
}

data_one_guild[[j]]=site_response
}

  else{
    nearest_df_site%>%filter(site%in%modified_site)%>%pull(nearest_site)->near_sites
    site_response=matrix(ncol=2,nrow=times)
    for (m in 1:2)
    {
      subset_samples(natural_data,Site%in%near_sites[m])->sub_natural_data
      subset_samples(modified_data,Site%in%modified_site[m])->sub_modified_data
      n_sub_natural_samples=nsamples(sub_natural_data)
      n_sub_modified_samples=nsamples(sub_modified_data)
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
     data_one_guild[[j]]=site_response
  }
}

all_data[[n]]=data_one_guild
}



# for each biome to get the mean value

  
  biome_site_response=matrix(ncol=4,nrow=9)
  for (j in 1:4)
    for (i in 1:9)
  {
    biome_site_response[i,j]=apply(all_data[[i]][[j]]%>%data.frame%>%filter_all(all_vars(!is.infinite(.))),2,FUN=mean)%>%mean()
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
    sampled_names <- replicate(times,sample(sample_names, n_crop_sample,replace = FALSE))# the same number of the soil cores
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

# to test the difference between the ratio and the 
test_result_guild=list()
for (m in 1:9){

set.seed(123)
test_result=matrix(nrow=4,ncol=3)
for (j in 1:4)
{
  # get a subset data for that
  sub_data=sample(mm[[m]][[j]],100)
  
test_result[j,1]=t.test(sub_data,mu=1,alternative = "two.sided")$estimate
test_result[j,2:3]=t.test(sub_data,mu=1,alternative = "two.sided")[4]$conf.int[1:2]
}

test_result_guild[[m]]=test_result
}

do.call(rbind,test_result_guild)%>%data.frame()%>%
  rename_all(~paste0(c("mean","low","up")))%>%
  mutate(guild=rep(c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para"),each=4))%>%
  mutate(biome=rep(biomes,9))->t_test_result_guild

ggplot()+
  geom_point(data=t_test_result_guild%>%filter(guild=="all"),size=2,aes(x=1:4,y=mean))+
  geom_segment(data=t_test_result_guild%>%filter(guild=="all"),size=0.35,color="black",aes(x=1:4,xend=1:4,y=low,yend=up))+
  geom_hline(yintercept = 1,color="red",lty="dashed")+
  ylab("Response ratio")+
  xlab("Biomes")


# if we just look at the richness and compared the mean rather than the ratio






# to get the response ratio of different guilds
# get the mean of the response ration
kkm=matrix(ncol=4,nrow=9)
for(j in 1:9){
  
kk=numeric()
for (i in 1:4)
  {
  kk[i]=mm[[j]][[i]]%>%mean()#need to test the significant among
  
}
kkm[j,]=kk
}

kkm%>%data.frame()%>%
  mutate(guild=c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para"))%>%
  rename_all(~paste0(c(biomes,"guild")))%>%melt()%>%
  rename(biome=variable,response=value)%>%
  mutate(guild_biome=paste(guild,"_",biome))->guild_specific_response_ratio

op=readRDS("full_guild_affinity_data.rds")

# to combine the response ratio with the c and z values

full_parameter_data%>%left_join(land%>%distinct(),by="plotID")%>%
  left_join(plot_biomes)%>%group_by(guild,LABEL)%>%
  summarise(mean_zvalue=mean(zvalue,na.rm=TRUE),mean_cvalue=mean(2.71828^logc,na.rm=TRUE))%>%
  rename(biome=LABEL)%>%mutate(guild_biome=paste(guild,"_",biome))%>%
  dplyr::select( biome,mean_zvalue, mean_cvalue,guild_biome)%>%
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
  

# none sites available in the natural data, we need to find the neighboring with the same
# determine the distance between them


melt(biome_site_response)%>%bind_cols(melt(kkm))%>%data.frame()%>%
  dplyr::select(value...3,value...6)%>%
  rename_all(~paste0(c("site_specific","non_site_specific")))->temp

ggplot()+
  geom_point(data=temp,aes(y=site_specific,x= non_site_specific),size=3)

# determine the difference in the species composition among land use types and different biomes

all_richness_data_guild=list()
for( m in 1:9)#select three guilds
{
cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # 
all_richness_data=list()
  for(k in 1:4)
  {
    #cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # 
    biome_data=subset_samples(get(data[m]),LABEL==biomes[k])# select the biome of interest
    # to see which are croplands
    subset_samples(biome_data,type=="cultivatedCrops")->modified_data
    #soil samples
    dim(sample_data(modified_data)%>%data.frame())[1]->n_crop_sample
    subset_samples(biome_data,type!="cultivatedCrops")->natural_data
    dim(sample_data(natural_data)%>%data.frame())[1]->n_nature_sample
    sample_data(modified_data)%>%data.frame()%>%pull(plotIDM)%>%
      substr(start = 1, stop = 4)%>%unique()->modified_site
    sample_data(natural_data)%>%data.frame()%>%pull(plotIDM)%>%
      substr(start = 1, stop = 4)%>%unique()->natural_site
    overlap <- intersect(modified_site, natural_site)
    # with the same site to compare the species composition
 
    natural_data%>%sample_data%>%group_by(type)%>%count()
    modified_data%>%sample_data%>%group_by(type)%>%count()
    
   # deciduousForest  1129
    # evergreenForest   573
    #mixedForest       313
    #shrubScrub         10
    #woodyWetlands     207
    #cultivatedCrops   117
  # if we focused on the 5 forest types, we can just compare the composition
  # select 100 samples for each of the type
   
    type01=c("cultivatedCrops","deciduousForest","evergreenForest","mixedForest","woodyWetlands")
    type02=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous", "mixedForest","shrubScrub","woodyWetlands")
    type03=c("cultivatedCrops","deciduousForest","grasslandHerbaceous", "emergentHerbaceousWetlands","shrubScrub")
    type04=c("cultivatedCrops","grasslandHerbaceous", "evergreenForest")
    type=list(type01,type02,type03,type04)
    times_sample=c(100,20,25,8)
    # does not consider differences in the sampling season
      times=500
     # set.seed(123)
      number_types=type[[k]]%>%length()
      richness_guild=matrix(ncol=number_types,nrow=times)
      for (i in 1:number_types)
      {
        type_select=type[[k]][i]
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
        richness=numeric()
        for (j in 1:times){
        #cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
        dk=subset_samples(biome_data,type==type_select)
        sample_names=sample_names(dk)
        sampled_names <- replicate(times,sample(sample_names, times_sample[k]))
        sampled_physeq <- prune_samples(sampled_names[,j], dk)
        richness[j]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
      } 
        richness_guild[,i]=richness
    }
      all_richness_data[[k]]=richness_guild
  }
      
all_richness_data_guild[[m]]=all_richness_data
}
      
#save the results

saveRDS(all_richness_data,"all_richness_data.rds")

saveRDS(all_richness_data_guild,"all_richness_data_guild.rds")#data saved


#compared the richness among different land use types

# construct a color matrix to match the color

unique(unlist(type))%>%data.frame()%>%
  mutate(values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan","#f0a73a"))%>%
  rename_all(~paste0(c("landcover","color")))->color_match



pp=list()
for(i in 1:4)
{

all_richness_data_guild[[1]][[i]]%>%data.frame()%>%
  rename_all(~paste0(type[[i]]))%>%melt()->tem_data
k <- aggregate(value ~ variable, data = tem_data, FUN = mean)
od <- k[order(-k$value), ] # with the increase trend to display the box plots
tem_data$variable <- factor(tem_data$variable, levels = od$variable)

#test the significant of the difference
color_match%>%
  mutate(land_cover=c("Cultivated crops", 
                      "Deciduous forest","Evergreen forest",
                      "Mixed forest","Woody wetlands",
                      "Grassland herbaceous","Shrub scrub",
                      "Emergent herbaceous wetlands"))->color_match

#comp_result[[2]]$Letters%>%data.frame()%>%rownames_to_column()%>%
 # rename_all(~paste0(c("landcover","letter")))->significant

type[[i]]%>%data.frame()%>%rename_all(~paste0(c("landcover")))%>%left_join(color_match,by="landcover")->
  color_select
pp[[i]]=ggplot(tem_data,aes(x=variable,y=value,fill=variable),alpha=0.5)+
  sm_raincloud()+
  geom_boxplot(width = 0.1, color = "black", outlier.size = 2) +
  theme(legend.position =c(0.32,0.20), legend.text = element_text(size = 12), 
        text = element_text(size = 12), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),
        panel.background = element_rect(fill = "NA")) +
  ylab("Core-level ichness") +
  xlab(paste0(biomes[i]))+
  scale_fill_manual("",breaks=type[[i]],values = color_select$color,labels=color_select$land_cover)+
  #ylim(-10,300)+
  guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.7))+
  guides(fill="none")+
  guides(color="none")
}

plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],ncol=2)

## test the significant of the difference
comp_result=list()
for(i in 1:4)
{
  
  all_richness_data[[i]]%>%data.frame()%>%
    rename_all(~paste0(type[[i]]))%>%melt()->tem_data
  k <- aggregate(value ~ variable, data = tem_data, FUN = mean)
  od <- k[order(-k$value), ] # with the increase trend to display the box plots
  tem_data$variable <- factor(tem_data$variable, levels = od$variable)
anova_result <- aov(value ~ factor(variable), data = tem_data)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_pvalues <- tukey_result$`factor(variable)`[, 4]  # Extract p-values
comp_result[[i]] <- multcompLetters(tukey_pvalues, compare = "<", threshold = 0.05)
}


      
    
#compare the species composition among land use types
# for each land use type select the same number
# for the first biome

community_data_guild=list()
statistic_guild=list()
for(m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # 
  
      set.seed(556)
      # the number of the land use type
      
      community_data=list()
      statistic=matrix(ncol=3,nrow=4)
      for (k in 1:4)
      {
     # cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # 
        
      type01=c("cultivatedCrops","deciduousForest","evergreenForest","mixedForest","woodyWetlands")
      type02=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous", "mixedForest","shrubScrub","woodyWetlands")
      type03=c("cultivatedCrops","deciduousForest","grasslandHerbaceous", "emergentHerbaceousWetlands","shrubScrub")
      type04=c("cultivatedCrops","grasslandHerbaceous", "evergreenForest")
      type=list(type01,type02,type03,type04)
      #times_sample=c(100,20,25,8)
      times_sample=shared_sample_number_guild[[m]]
      
      number_types=type[[k]]%>%length()
      
      sub_com=list()
      for (i in 1:number_types)
      {
        #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
        type_select=type[[k]][i]
        biome_data=subset_samples(get(data[m]),LABEL==biomes[k])#
        biome_data= prune_samples(sample_sums(biome_data)> 0, biome_data)#select non-zero rows
       # number of samples need to be sampled
         dk=subset_samples(biome_data,type==type_select)
        sample_names=sample_names(dk)
        sampled_names <- sample(sample_names, times_sample[k],replace = FALSE)
        sampled_physeq <- prune_samples(sampled_names, dk)
        sub_com[[i]]=otu_table(sampled_physeq )%>%data.frame()
      }
    
      # bind all the community data; the presented results were based on one permutation
      vege_com=combined_df <- do.call(rbind, sub_com)%>%
        mutate(type=rep(type[[k]],each=times_sample[k]))#
      ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
      com_score=ordination$species
      com_score=ordination$points%>%data.frame()%>%mutate(type=rep(type[[k]],each=times_sample[k]))
      community_data[[k]]=com_score
      vege_com[vege_com>0]=1
      adonis_result <- adonis(vege_com[, -ncol(vege_com)] ~ type, data =com_score, permutations = 999)
      statistic[k,1]=ordination$stress
      statistic[k,2]=adonis_result $aov.tab[1,c(4,6)][1]%>%pull()#f-value
      statistic[k,3]=adonis_result $aov.tab[1,c(4,6)][2]%>%pull()#p-vale
      }
     
      community_data_guild[[m]]=community_data
      
      statistic_guild[[m]]=statistic
}



saveRDS(community_data_guild,file="community_data_guild.rds")
saveRDS( statistic_guild,file=" statistic_guild.rds")



   
# need to consider different fungal guilds
  # at the guild, some sample do not have any species
  
  
  
      
      
statistic_guild=round(statistic,2)
      
pp1=list()
      for (i in 1:4){
        
data_range=community_data_guild[[1]][[i]]%>%data.frame()
        y_range=range(data_range$MDS2)[2]
        x_range=range(data_range$MDS1)[1]
       
type[[i]]%>%data.frame()%>%rename_all(~paste0(c("landcover")))%>%left_join(color_match,by="landcover")->
          color_select
        
 pp1[[i]]=ggplotGrob(ggplot(data=community_data_guild[[1]][[i]],aes(x=MDS1,  y= MDS2,color=type ))+
        geom_point(data=community_data_guild[[1]][[i]],pch=21,color="black",aes(x=MDS1,y= MDS2,fill=type ),size=3,alpha=0.75)+
        stat_ellipse(size=0.8,linetype="dashed")+
        theme(legend.position = c(0.3,0.3), 
              legend.title = element_text(size=10),
              text = element_text(size = 18), 
              legend.text = element_text(size=11),
              plot.title = element_text(size = 15, hjust = 0.5), 
              axis.text.y = element_text(hjust = 1), 
              axis.text.x = element_text(hjust = 0), 
              axis.title.y = element_text(size = 18), 
              axis.title.x = element_text(size = 18),
              panel.background = element_rect(fill = "NA"), 
              panel.border = element_rect(color = "black", size = 1, fill = NA))+
   #guides(color="none",fill="none")+
        
        annotate("text",x=max(community_data_guild[[1]][[i]]$MDS1),y=max(community_data_guild[[1]][[i]]$MDS2),label=bquote(italic(F)~ "=" ~ .(statistic_guild[[1]][i, 1])),size=5)+
        #annotate("text",x=0.64,y=0.705,label=expression(italic(P)*" < 0.001"),size=5)+
        xlab("MDS1")+
        scale_fill_manual("",breaks=type[[i]],values = color_select$color,labels=color_select$land_cover)+
      scale_color_manual("",breaks=type[[i]],values = color_select$color,labels=color_select$land_cover))
 
      
        }
      



pp[[1]]$widths=pp[[2]]$widths
pp[[2]]$widths=pp[[3]]$widths
pp[[4]]$widths=pp[[3]]$widths
pp[[1]]$widths=pp[[3]]$widths

plot_grid(pp1[[1]],pp[[1]],pp1[[2]],pp[[2]], pp1[[3]], pp[[3]],pp1[[4]],pp[[4]],ncol=2)

#to creat the legend
combined_legend=cowplot::get_legend(pp[[1]])

combined_legend1=cowplot::get_legend(pp[[2]])

legend_combined <- plot_grid(combined_legend, combined_legend1, ncol = 2)


# Get legend from one of the plots


rare_all_guild%>%sample_data()%>%data.frame()%>%
  distinct(plotIDM)%>%pull(plotIDM)->plot_id
rare_all_guild%>%sample_data()%>%data.frame()%>%filter(Project=="NEON")->temp_data

dd=list()
for(i in 1:515){

dd[[i]]=temp_data%>%filter(plotIDM==plot_id[i])%>%distinct(lon)
}


#
mm=list()
mm_size=list()
for (k in 1:9)
{
  
  response_ratio=list()
  size=matrix(nrow=4,ncol=2)
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
    
    richness_natural=estimate_richness(natural_data, measures = "Observed")
    
    richness_modified=estimate_richness(modified_data, measures = "Observed")
    
    richness_natural%>%data.frame()%>%rename_all(~paste0(c("Observed")))->richness_natural
    
    richness_modified=estimate_richness(modified_data, measures = "Observed")
    
    N1=dim(richness_modified)[1]
    N2=dim(richness_natural)[1]
    
    richness_modified%>%bind_rows(richness_natural)%>%
      mutate(type=rep(c("modified","nature"),times=c(N1,N2)))->nature_modified
    
    response_ratio[[j]]=nature_modified%>%group_by(type)%>%summarise(mean_value=mean(Observed))%>%data.frame()
    size[j,1]=N1
    size[j,2]=N2
  }
  
  mm[[k]]=response_ratio
  mm_size[[k]]=size
}
