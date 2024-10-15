## get the history map for the continental US

setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")

land_use_data <- rast("nlcd_2001_land_cover_l48_20210604.img")

plot_coordinates_convert=plot_coordinates[,2:3]

plot_coordinates_convert%>%
  st_as_sf(coords = c('lon', 'lat'), crs ="EPSG:4326" )%>%
  vect()%>%terra::project(land_use_data )->temp_coords

values=terra::extract(result_raster,temp_coords)

values=terra::extract(result_raster,plot_coordinates_convert)

## to select the human modified plots

rare_all_guild_biome=readRDS("rare_all_guild_biome.rds")

sample_data(rare_all_guild_biome)%>%data.frame()%>%
filter(type=="cultivatedCrops")%>%dplyr::select(Site,lon,lat,plotIDM)%>%distinct(plotIDM)


#BLAN_033-M-20.5-2.5-20151110-GEN-DNA1 BLAN
#SCBI_006-M-9-3.5-20151113-GEN-DNA2    SCBI
#STER_031-M-25-7.5-20160420-GEN-DNA1   STER
#SERC_007-M-6.5-28-20170412-GEN-DNA1   SERC
#DSNY_041-M-36-36-20170829-GEN-DNA1    DSNY
#JERC_004-M-0.5-25-20170724-GEN-DNA1   JERC
#KONA_003-M-12-2.5-20170807-GEN-DNA1   KONA
#ORNL_021-O-28.5-33-20161102-GEN-DNA1  ORNL
#LAJA_044-M-39-2.5-20170307-GEN-DNA1   LAJA

# 45 plots from 9 sites are human-modified plots

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

plot_coordinates=saveRDS(plot_coordinates,file="plot_coordinates.rds")

plot_coordinates= readRDS("plot_coordinates.rds")

#for each site, determine its nearest site

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

# determine the spatial distance among sites

sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(Site,lon,lat)%>%
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


#get the human-modified plots
# croplands and pastures were pooled
#374 samples 

sample_data(rare_all_guild_biome)%>%data.frame()%>%
  filter(type=="cultivatedCrops")%>%dplyr::select(Site,lon,lat,plotIDM)%>%distinct()->crop_plot


# extract the pre-modified land use type for each plot based on their coordinates
# based on a publicly available data(https://doi.org/10.5281/zenodo.7055086)

#LULC type:
#0 nodata value
#1 urban
#2 crop
#3 pasture
#4 forest
#5 shrub
#6 grassland
#7 wetland
#8 water
#9 barren

historical_land_use=list()
for (i in 1630:2020)
{
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
file_name=paste("conus_lulc_",i,".tif")
file_name=gsub(" ","",file_name)
land_history_data <- rast(file_name)
result_raster_history <- project(land_history_data, 'EPSG:4326',method="near")# select a CRS
historical_land_use[[i]]=terra::extract(result_raster_history,crop_plot[,c("lon","lat")],method="simple")%>%
  data.frame()%>%dplyr::select(-ID)%>%rename_all(~paste0(c("year")))
}



saveRDS(historical_land_use,file="historical_land_use.rds")

historical_land_use=readRDS("historical_land_use.rds")


subset_columns <- historical_land_use[1630:2020]

do.call(cbind,subset_columns)%>%data.frame()->land_use_among_year

# to see when the cropland and pasture were firstly appeared


land_use_among_year%>%melt()%>%mutate(year=rep(c(1630:2020),each=45))%>%
  mutate(plotid=rep(crop_plot$plotIDM,times=391))->values_temp

values_temp%>%filter(value==2)%>%group_by(plotid)%>%
  summarise(appear_year=min(year))->df1

values_temp%>%filter(value==3)%>%group_by(plotid)%>%
  summarise(appear_year=min(year))->df2

df1%>%bind_rows(df2)

first_year=numeric()
for (i in 1:45)
{
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  values_temp%>%filter(plotid==plotid[i])->temp_data
  if(any(unique(temp_data$value)%in%c(2,3)))
  {
    first_year[i]=temp_data%>% filter(value == "2" | value == "3")%>%
      filter(year==min(year,na.rm=TRUE))%>%pull(year) 
  }
  else
    {
      first_year[i]=NA
  }
}


# to see the variation of the land use type for each plot
# the most often seen type was considered as the pre-modified land use type for each plot

type=numeric()
for (i in 1:45)
  {
 type[i]=land_use_among_year[i,]%>%t()%>%data.frame()%>%
   rename_all(~paste0(c("name")))%>%
 group_by(name)%>%count()%>%data.frame()%>%
 slice_max(n,n=1)%>%pull(name)
  
}

# the historical land use type for each plot
# all the plots in the LAJA site were assigend to forest because of none data
crop_plot%>%bind_cols(type)%>%
  rename(history=...5)%>%
  mutate(historical_type=case_when(history%in%c("4","NA")~"forest", 
                                   history=="7"~"wetland", 
                                   history=="6"~"grassland", 
           TRUE~"forest"))->df4


#the current land use type for each of the nine sites where human-modified plots occur
#there are cases that all the plots within a site are human-modified plots

site=df4$Site
d=list()
for(i in 1:9)
{
  d[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(type)
  
}
# for the STER site, all the plots are the croplands


# test if one site has only one biome？(yes!)
biome_type=list()
for(i in 1:9)
{
  biome_type[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(LABEL)
}

  site=unique(df4$Site)
  
# the number land cover types historically per site
# in most cases, there are only one type for their historical landuse type
# for KONA and DSNY, there are two land use types per site

crop_type=list()
for (i in 1:9)
  {
  crop_type[[i]]=df4%>%filter(Site==site[i])#need to make sure that the site has nine elements
}


#create a DF so that the historical and present land use type can be matched
all_land_cover=sample_data(rare_all_guild_biome)%>%data.frame()%>%distinct(type)%>%
  mutate(broad_type=case_when(type%in%c("evergreenForest", "deciduousForest","mixedForest")~"forest", 
                                   type%in%c("emergentHerbaceousWetlands","woodyWetlands")~"wetland", 
                                   type%in%c("grasslandHerbaceous","sedgeHerbaceous")~"grassland", 
                                  type==c("shrubScrub","dwarfScrub")~"shrub", 
                                   TRUE~"other"))

# to select the same land use type to get the mean richness
# determine the richness based on each plot
# the pairwise plots were pooled rather than being plot-specific

#for each human-modified plot, calculate the richness in the plot
# and find the co-site plots having the analogous land use type as the human-modified plot's historical land use type



data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")
data_EM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild_biome, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "epiphyte")

###################used the below codes to get the species ratio##############

plotid=df4$plotIDM
site=df4$Site

# for other tropic guilds
#compare_richness_guild=list()
response_ratio_guild=list()
for (m in 1:9)
{
  
cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = ""))

response_ratio=numeric()
sample_size=matrix(nrow=45,ncol=2)
#t_test_result=numeric()
#compare_richness=list()
#species_com=list()
for (i in c(1:45))
{
 cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  subset_samples(get(data[m]),plotIDM%in%df4$plotIDM)%>%
    subset_samples(plotIDM==plotid[i])->modified_data
  
  df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
  # the land cover to be selected from the adjacent natural plots
  all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
  unique(df4$Site)->modified_site
  subset_samples(get(data[m]),type!="cultivatedCrops")->natural_data
  sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
  
  #to see if the two data sets have shared sites
  #setdiff(modified_site,natural_site)->sites_no_in_nature
  
  # find the plot within the same site
  # all sites in the 
  Site=df4$Site
  #when no sites are available, we can find the plots based on the most close sites
  nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site
  
  matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)# for these plots, they do not have shared sites
  
  if(i%in%c(matching_elements,24))
  {
    subset_samples(natural_data,Site==near_site&type%in%select_type)->natural_data_sub
    
  }
  else if(i==12)#can not find analogous plots
  {
    subset_samples(natural_data,Site==site[i])->natural_data_sub
  }
  
  else{
    subset_samples(natural_data,Site==site[i]&type%in%select_type)->natural_data_sub
    
  }
  #when i=12, the history land use type was grassland but all current plots are deciduous forest
  # the nearest plots are also forest, so we used forest as the analogous land use type
  #when i=24 the history land use type was forest but currently all are woodyWetlands
  #and we get the adjacent sites
  # for the LAJA site we do not have historical data and we assumed the historical land use type for all the plots was forest
  n_natural_sample=nsamples(natural_data_sub)
  n_modified_sample=nsamples(modified_data)
  
  # need to get the land use type for the natural data
  natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type)->natural_land_use
  modified_data%>%sample_data()%>%data.frame()%>%dplyr::select(type)->modified_land_use
  
  combined_land_use=rbind(modified_land_use,natural_land_use)
  richness_modified=estimate_richness(modified_data, measures = "Observed")
  richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
  
  #otu_table(modified_data)%>%data.frame()->otu_modified
  #otu_table(natural_data)%>%data.frame()->otu_natural
  
  sample_size[i,1]=n_modified_sample
  sample_size[i,2]=n_natural_sample
  #compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)%>%
    #mutate(type = if_else(type %in% c("cultivatedCrops"), type, "natural"))
   response_ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
  #species_com[[i]]=bind_rows(otu_modified,otu_natural)

}

#compare_richness_guild[[m]]=compare_richness
response_ratio_guild[[m]]=response_ratio
}

saveRDS(response_ratio_guild,file="response_ratio_guild.rds")

# to get the richness data among land use types, which was computaed on great lakes

compare_richness_guild=readRDS("compare_richness_guild.rds")

response_ratio_guild=readRDS("response_ratio_guild.rds")



# get the mean of the response ratio for each guild
# and biome without considering the significance of the difference

biome_response_ratio=list()
for (i in 1:9)
{
  
bind_cols(df4,response_ratio_guild[[i]])%>%
  rename(ratio=...7)%>%
    bind_cols(sample_size)%>%
    filter(...8>=3)->plot_response_ratio

sample_data(rare_all_guild_biome)%>%
  data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%
  left_join(plot_response_ratio,by="Site")%>%
  filter(!is.na(ratio)) %>%
  mutate(LABEL= ifelse(is.na(as.character(LABEL)), "Tropical & Subtropical Moist Broadleaf Forests", as.character(LABEL)))%>%
  group_by(LABEL)%>%distinct()->ratio_temp

biome_response_ratio[[i]]=ratio_temp%>%summarise(mean_value=mean(ratio))
}

# in this case, we will have positive response
# Temperate Broadleaf & Mixed Forests                 0.936
# Temperate Conifer Forests                           1.06 
# Temperate Grasslands, Savannas & Shrublands         0.829
# Tropical & Subtropical Moist Broadleaf Forests      1.08 

#to test if the difference is significant
# need to select the plots with more than 3 samples
# significance test are the same if we use the aov test with unequal variance

# select the plots with more than three samples to test the significance of differene in richness
sample_size%>%data.frame()%>%
  bind_cols(df4$plotIDM)%>%
  mutate(no=1:45)%>%
  filter(X1>=3)->plot_with_three_cores

# for each guild to test the difference in richness among the natural and the human modified communities
#get the p-values

significant_test=list()
for (m in 1:9 ){
result=numeric()
for (i in plot_with_three_cores$no )
  {
test_result=t.test(compare_richness_guild[[m]][[i]]%>%filter(type=="cultivatedCrops")%>%
         dplyr::select(Observed),
         compare_richness_guild[[m]][[i]]%>%
         filter(type=="natural")%>%dplyr::select(Observed), var.equal = FALSE)
result[i]=test_result$p.value
}
significant_test[[m]]=result
}

#based on the t-test to select the plots showing significant difference

final_biome_ratio=list()
for (m in 1:9)
{
  final_biome_ratio[[m]]=df4%>%bind_cols(response_ratio_guild[[m]])%>%
    bind_cols(significant_test[[m]])%>%
    left_join(sample_data(rare_all_guild_biome)%>%
    data.frame()%>%dplyr::select(Site,LABEL),by="Site")%>%
    distinct()%>%mutate(LABEL= ifelse(is.na(as.character(LABEL)), "Tropical & Subtropical Moist Broadleaf Forests", as.character(LABEL)))%>%
    group_by(LABEL)%>%distinct()%>%
    rename(ratio=...7,pval=...8)%>%filter(pval<0.1)%>%
    group_by(LABEL)%>%summarise(mean_value=mean(ratio))
}



#if we do not do significant test and just look at the overall trend of difference between the richness

trend_biome_ratio=list()
for (m in 1:9)
{
  trend_biome_ratio[[m]]=df4%>%bind_cols(response_ratio_guild[[m]])%>%
    left_join(sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(Site,LABEL),by="Site")%>%
    distinct()%>%mutate(LABEL= ifelse(is.na(as.character(LABEL)), "Tropical & Subtropical Moist Broadleaf Forests", as.character(LABEL)))%>%
    group_by(LABEL)%>%distinct()%>%
    rename(ratio=...7)%>%
    group_by(LABEL)%>%summarise(mean_value=mean(ratio))
}

row_counts <- sapply(final_biome_ratio, nrow)
do.call(rbind,final_biome_ratio)%>%
  mutate(guild=rep(c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para"),times=row_counts))->df1


do.call(rbind,trend_biome_ratio)%>%
  mutate(guild=rep(c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para"),each=4))->df2

df3=rbind(df2,df1)%>%mutate(type=rep(c("trend","sig"),times=c(36,34)))


a=ggplot(df3%>%filter(type=="trend"), aes(x = LABEL, y = mean_value, fill = guild)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle=90))+
  geom_hline(yintercept = 1,color="red",linetype="dashed")


b=ggplot(df3%>%filter(type=="sig"), aes(x = LABEL, y = mean_value, fill = guild)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle=90))+
  geom_hline(yintercept = 1,color="red",linetype="dashed")

plot_grid(a,b,ncol=2,labels = c("trend","sig"))

# based on the significance test to get the mean value of the response ratio 



# if we make the comparison based on the pooled data rather than the individual plots
# cores from the cropland and cores from natural plots were pooled 
# soil samples taken from the same plots of the same historical land use type and same site were pooled
site_id=unique(df4$Site)
site_id=unique(df4$Site)

dd=list()
for (i in c(1:9))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  subset_samples(rare_all_guild_biome,type=="cultivatedCrops")->modified_data

  # when no sites overlap between the two, we need to find the adjacent sites
  
  subset_samples(modified_data,Site==site_id[i])->modified_temp
    
    # to see the historical land use type for these samples
    sample_data(modified_temp)%>%data.frame()%>%distinct(plotIDM)->select_plot
    
    df4%>%filter(plotIDM%in%select_plot$plotIDM)%>%pull(historical_type)->historical_type
    dd[[i]]=historical_type
    
    #based on the historical land use type to select the current plots
   # all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type#to be select in the natural data
}


# when i =5 and i=7, there are two types of types need to 



ratio=numeric()
for (i in c(1:9))
{
  if(i%in%c(1:3,6,8,9))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
    subset_samples(rare_all_guild_biome,type=="cultivatedCrops")->modified_data
    subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
    
  # when no sites overlap between the two, we need to find the adjacent sites
    
    nearest_df_site%>%filter(site==site_id[i])%>%pull(nearest_site)->near_site
    subset_samples(modified_data,Site==site_id[i])->modified_temp
    natural_data%>%sample_data()%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
    
     if(site_id[i]%in%natural_site)#they have shared sites
      {
    
    # to see the historical land use type for these samples
    sample_data(modified_temp)%>%data.frame()%>%distinct(plotIDM)->select_plot
    df4%>%filter(plotIDM%in%select_plot$plotIDM)%>%pull(historical_type)->historical_type
    #to be select in the natural data
    all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
    subset_samples(natural_data,Site==site_id[i])%>%sample_data()%>%data.frame()%>%distinct(type)%>%pull(type)->natural_type
    
    shared_type=intersect(select_type,natural_type)
    
    if(length(shared_type)<1)# when no analogous type was found for the current plots
      {
      subset_samples(natural_data,Site==near_site&type%in%select_type)->natural_data_sub
    }
    
    else{
    # since all historical land cover type is the same, we pooled the site-level data
    subset_samples(natural_data,Site==site_id[i]&type%in%select_type)->natural_data_sub
      }
    richness_modified=estimate_richness(modified_temp, measures = "Observed")
    richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
    ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed) 
    
     }
  
    # when the modified sites do not belong to the natural sites
    #based on the adjacent sites
    else{
      
      # to see the historical land use type for these samples
      
      nearest_df_site%>%filter(site==site_id[i])%>%pull(nearest_site)->near_site
      sample_data(modified_temp)%>%data.frame()%>%distinct(plotIDM)->select_plot
      df4%>%filter(plotIDM%in%select_plot$plotIDM)%>%pull(historical_type)->historical_type
      #based on the historical land use type to select the current plots
      all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type#to be select in the natural data
  
    subset_samples(natural_data,Site==near_site&type%in%select_type)->natural_data_sub
    
      richness_modified=estimate_richness(modified_temp, measures = "Observed")
      richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
      ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed) 
    }
  }
    
else if(i==7)
  {
    subset_samples(rare_all_guild_biome,type=="cultivatedCrops")->modified_data
    subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
    
    nearest_df_site%>%filter(site==site_id[i])%>%pull(nearest_site)->near_site
    
    # when no sites overlap between the two, we need to find the adjacent sites
    subset_samples(modified_data,Site==site_id[i])->modified_temp
    #based on the historical land use type to select the current plot
    
    subset_samples(natural_data,Site==near_site)->natural_data_sub#can not find related
    
    richness_modified=estimate_richness(modified_temp, measures = "Observed")
    richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
    ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed) 
    }
    else
      {
        subset_samples(rare_all_guild_biome,type=="cultivatedCrops")->modified_data
        subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
        nearest_df_site%>%filter(site==site_id[i])%>%pull(nearest_site)->near_site
        # when no sites overlap between the two, we need to find the adjacent sites
        subset_samples(modified_data,Site==site_id[i])->modified_temp
        #based on the historical land use type to select the current plot
        subset_samples(natural_data,Site==site_id[i])->natural_data_sub#can not find related
       
         richness_modified=estimate_richness(modified_temp, measures = "Observed")
        richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
        ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed) 
     }
}

# how to test the significance

sample_data(rare_all_guild_biome)%>%data.frame()%>%dplyr::select(Site,LABEL)%>%distinct()%>%
  filter(Site%in%site_id)%>%head(9)%>%bind_cols(ratio)%>%group_by(LABEL)%>%summarise(m=mean(...3))


# get the mean of the response ratio
#1 Temperate Broadleaf & Mixed Forests            0.947
#2 Temperate Conifer Forests                      0.974
#3 Temperate Grasslands, Savannas & Shrublands    0.847
#4 Tropical & Subtropical Moist Broadleaf Forests 1.07 




  












 
    
    
    
   
    
    
    
  df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
  # the land cover to be selected from the adjacent natural plots
  
  
  unique(df4$Site)->modified_site
  subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
  sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
  
  # find the sites not included in the natural data
  #to see if the two data sets have overlap
  setdiff(modified_site,natural_site)->sites_no_in_nature
  
  # find the plot within the same site
  Site=df4$Site
  
  #when no sites are available, we can find the plots based on the most close sites
  
  nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site
  
  matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)
  
  if(i%in%c(matching_elements,24))
  {
    subset_samples(natural_data,Site==near_site&type%in%select_type)->natural_data_sub
    
  }
  else if(i==12)
  {
    subset_samples(natural_data,Site==site[i])->natural_data_sub
  }
  
  else{
    subset_samples(natural_data,Site==site[i]&type%in%select_type)->natural_data_sub
    
  }
  #when i =12, the history land use type was grassland but currently all are deciduous forest
  # the nearest is also forest, so we used forest as the analogous land use type
  
  #when i=24 the history land use type was forest but currently all are woodyWetlands
  #and we get the adjacent sites
  
  # for the LAJA site we do not have historical data and we assumed the historical land use type was forest
  n_natural_sample=nsamples(natural_data_sub)
  n_modified_sample=nsamples(modified_data)
  
  # need to get the land use type for the natural data
  natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type)->natural_land_use
  modified_data%>%sample_data()%>%data.frame()%>%dplyr::select(type)->modified_land_use
  
  combined_land_use=rbind(modified_land_use,natural_land_use)
  
  richness_modified=estimate_richness(modified_data, measures = "Observed")
  richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
  
  #otu_table(modified_data)%>%data.frame()->otu_modified
  #otu_table(natural_data)%>%data.frame()->otu_natural
  
  #compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)
  response_ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
  
  sample_size[i,1]=n_modified_sample
  sample_size[i,2]=n_natural_sample
  #species_com[[i]]=bind_rows(otu_modified,otu_natural)
  
  #t_test_result <- t.test(richness_modified$Observed, richness_natural$Observed, var.equal = FALSE)
  #t_test_result[i]=t_test_result$p.value
}
