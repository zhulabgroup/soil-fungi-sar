## get the history map for the continental US

setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-effect")

land_use_data <- rast("nlcd_2001_land_cover_l48_20210604.img")

plot_coordinates_convert=plot_coordinates[,2:3]

plot_coordinates_convert%>%
  st_as_sf(coords = c('lon', 'lat'), crs ="EPSG:4326" )%>%
  vect()%>%terra::project(land_use_data )->temp_coords

values=terra::extract(result_raster,temp_coords)

values=terra::extract(result_raster,plot_coordinates_convert)

## to select the plots belong to the croplands

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

# 45 plots from 9 sites belong to human-modified land

sample_data(rare_all_guild_biome)%>%data.frame()%>%
  filter(type=="cultivatedCrops")%>%dplyr::select(Site,lon,lat,plotIDM)%>%distinct()->crop_plot


# extract the pre-modified land use type for each crop based on their coordinates for the nlcd

values=terra::extract(result_raster,crop_plot[,c("lon","lat")])
bind_cols(values,crop_plot)

# most of the data are the same as their current forms

# to use another data set to extract the land use data

values=list()
for (i in 1630:2020)
{
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
file_name=paste("conus_lulc_",i,".tif")
file_name=gsub(" ","",file_name)
land_history_data <- rast(file_name)
result_raster_history <- project(land_history_data, 'EPSG:4326',method="near")# select a CRS
values[[i]]=terra::extract(result_raster_history,crop_plot[,c("lon","lat")],method="simple")%>%
  data.frame()%>%dplyr::select(-ID)%>%rename_all(~paste0(c("year")))
}

historical_land_use=values

saveRDS(historical_land_use,file="historical_land_use.rds")

values=readRDS("historical_land_use.rds")


subset_columns <- values[1630:2020]

do.call(cbind,subset_columns)%>%data.frame()->land_use_among_year

# to see when the crop land stars to appear in the 400 years

land_use_among_year[1, land_use_among_year[1, ] == 2]
land_use_among_year[land_use_among_year == 2 | land_use_among_year == 3]


# to see the variation in the land use type for each plot
# for each plot to select the most often or year land use type

type=numeric()
for (i in 1:45)
  {
 type[i]=land_use_among_year[i,]%>%t()%>%data.frame()%>%
   rename_all(~paste0(c("name")))%>%
 group_by(name)%>%count()%>%data.frame()%>%
 slice_max(n,n=1)%>%pull(name)
  
}

# to find the latest version or based on the nearest plot
crop_plot%>%bind_cols(type)%>%
  rename(history=...5)%>%
  mutate(historical_type=case_when(history%in%c("4","NA")~"forest", 
                                   history=="7"~"wetland", 
                                   history=="6"~"grassland", 
           TRUE~"forest"))->df4



# to get the pairwise plot to see how the plots are mathces
#there are cases that the initial land use type has disappeared so we need to find the neariest
site=df4$Site

d=list()
for(i in 1:9)
{
  d[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(type)
  
}
# for the STER site, all the plots are the croplands
  #to check if one site has several biomes?

# test if one site has only one biome？
biome_type=list()
for(i in 1:9)
{
  biome_type[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(LABEL)
  
}

  

crop_type=list()
for (i in 1:9)
  {
  crop_type[[i]]=df4%>%filter(Site==site[i])
}

# in most cases, there are only one type for their historical landuse type


rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i]&type=="deciduousForest")
# for the cropland data



# need to find the plots with forest
#create a dataframe so that the historical and present land use type can be matched

all_land_cover=sample_data(rare_all_guild_biome)%>%data.frame()%>%distinct(type)%>%
  mutate(broad_type=case_when(type%in%c("evergreenForest", "deciduousForest","mixedForest")~"forest", 
                                   type%in%c("emergentHerbaceousWetlands","woodyWetlands")~"wetland", 
                                   type%in%c("grasslandHerbaceous","sedgeHerbaceous")~"grassland", 
                                  type==c("shrubScrub","dwarfScrub")~"shrub", 
                              
                                   TRUE~"other"))

# to select the same land use type to get the mean richness
# determine the richness based on each plot
# the pairwise plots were pooled


#crop_data=subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)

#crop_data%>%sample_data()%>%data.frame()%>%filter(Site==site[i])->crop_sub_data

plotid=unique(df4$plotIDM)

response_ratio=numeric()
for (i in c(1:45))
  {
 cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))

  subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)%>%
  subset_samples(plotIDM==plotid[i])->modified_data
  
  df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
  all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
  
unique(df4$Site)->modified_site
subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
 
 # find the sites not included in the natural 
setdiff(modified_site,natural_site)->sites_no_in_nature
   
# find the plot within the same site
Site=df4$Site

#to see if have overlap for the two data sets
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
  #subset_samples(natural_data,Site==site[i])%>%sample_data()%>%data.frame()%>%distinct(type)%>%pull(type)->natural_type
  
  subset_samples(natural_data,Site==site[i]&type%in%select_type)->natural_data_sub
  
}
#when i =12, the history land use type was grassland but currently all are deciduous forest
# the neariest is also forest
#when i=24 the history landtype was forest but currently all are woodyWetlands
# for the LAJA site we do not have historical data and we assumed it was forest
#n_natural_sample=nsamples(natural_data_sub)
#n_modified_sample=nsamples(modified_data)

richness_modified=estimate_richness(modified_data, measures = "Observed")
richness_natural=estimate_richness(natural_data_sub, measures = "Observed")

response_ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
}


# only select the forest type and should not select the grassland
# it might be better to do that for each plot
# to get the richness of each plot

site=df4$Site

compare_richness=list()
#species_com=list()
for (i in c(1:45))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)%>%
    subset_samples(plotIDM==plotid[i])->modified_data
  
  df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
  
  all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
  
  
  unique(df4$Site)->modified_site
  subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
  sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
  
  # find the sites not included in the natural 
  setdiff(modified_site,natural_site)->sites_no_in_nature

  # find the plot within the same site
  Site=df4$Site
  
  #to see if have overlap for the two data sets
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
    #subset_samples(natural_data,Site==site[i])%>%sample_data()%>%data.frame()%>%distinct(type)%>%pull(type)->natural_type
    
    subset_samples(natural_data,Site==site[i]&type%in%select_type)->natural_data_sub
    
  }
  #when i =12, the history landtype was grassland but currently all are deciduous forest
  # the neariest is also forest
  #when i=24 the history landtype was forest but currently all are woodyWetlands
  # for the LAJA site we do not have historical data and we assumed it was forest
  #n_natural_sample=nsamples(natural_data_sub)
  #n_modified_sample=nsamples(modified_data)
  
  # need to get the land use type for the natural data
  natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type)->natural_land_use
  modified_data%>%sample_data()%>%data.frame()%>%dplyr::select(type)->modified_land_use
  
  combined_land_use=rbind(modified_land_use,natural_land_use)
  
  richness_modified=estimate_richness(modified_data, measures = "Observed")
  richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
  
  #otu_table(modified_data)%>%data.frame()->otu_modified
  #otu_table(natural_data)%>%data.frame()->otu_natural
  
  compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)
#species_com[[i]]=bind_rows(otu_modified,otu_natural)
}

#save the work space and run on great lakes
plotid=df4$plotIDM

site=df4$Site

# for other tropic guilds
compare_richness_guild=list()
for (m in 1:9)
{
  
  #cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), i, collapse = ""))
  
response_ratio=numeric()
sample_size=matrix(nrow=45,ncol=2)

#t_test_result=numeric()
compare_richness=list()
#species_com=list()
for (i in c(1:45))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)%>%
    subset_samples(plotIDM==plotid[i])->modified_data
  
  df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
  # the land cover to be selected from the adjacent natural plots
  all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
  
  
  unique(df4$Site)->modified_site
  
  subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
  sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
  
  
  #to see if the two data sets have shared sites
  
  setdiff(modified_site,natural_site)->sites_no_in_nature
  
  # find the plot within the same site
  # all sites in the 
  Site=df4$Site
  
  #when no sites are available, we can find the plots based on the most close sites
  
  nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site
  
  matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)
  
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
  
  sample_size[i,1]=n_modified_sample
  sample_size[i,2]=n_natural_sample
  
  compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)%>%
    mutate(type = if_else(type %in% c("cultivatedCrops"), type, "natural"))
  response_ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
  
  
  #species_com[[i]]=bind_rows(otu_modified,otu_natural)
  
  #t_test_result <- t.test(richness_modified$Observed, richness_natural$Observed, var.equal = FALSE)
  #t_test_result[i]=t_test_result$p.value
}

compare_richness_guild[[m]]=response_ratio
}

# get the mean of the response ratio

bind_cols(df4,response_ratio)%>%
  rename(ratio=...7)->plot_response_ratio

sample_data(rare_all_guild_biome)%>%
  data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%
  left_join(plot_response_ratio,by="Site")%>%
  filter(!is.na(ratio)) %>%
  mutate(LABEL= ifelse(is.na(as.character(LABEL)), "Tropical & Subtropical Moist Broadleaf Forests", as.character(LABEL)))%>%
  group_by(LABEL)->ratio_temp

ratio_temp%>%summarise(mean_value=mean(ratio))

# in this case, we will have positive response


# Temperate Broadleaf & Mixed Forests                 0.936
# Temperate Conifer Forests                           1.06 
# Temperate Grasslands, Savannas & Shrublands         0.829
# Tropical & Subtropical Moist Broadleaf Forests      1.08 

#to test if the difference is significant
# need to select the plots with more than 3 samples
sample_size%>%data.frame()%>%
  bind_cols(df4$plotIDM)%>%
  mutate(no=1:45)%>%
  filter(X1>=3)->plot_with_three_cores

result=numeric()
for (i in plot_with_three_cores$no )
  {
test_result=t.test(compare_richness[[i]]%>%filter(type=="cultivatedCrops")%>%
         dplyr::select(Observed),
       compare_richness[[i]]%>%
         filter(type=="natural")%>%dplyr::select(Observed), var.equal = FALSE)
result[i]=test_result$p.value
}

# if we test the results and get the values of each plot

plot_response_ratio%>%bind_cols(result)%>%
  filter(...8<=0.1)%>%
  left_join(sample_data(rare_all_guild_biome)%>%
 data.frame()%>%dplyr::select(Site,LABEL)%>%
   distinct(),by="Site")%>%
  mutate(LABEL= ifelse(is.na(as.character(LABEL)), "Tropical & Subtropical Moist Broadleaf Forests", as.character(LABEL)))%>%
  group_by(LABEL)%>%summarise(mm=mean(ratio))

#Temperate Broadleaf & Mixed Forests            0.754
# Temperate Conifer Forests                      1.23 
# Temperate Grasslands, Savannas & Shrublands    0.754
#Tropical & Subtropical Moist Broadleaf Forests 1.11 


result1=numeric()
for (i in plot_with_three_cores$no )
{
  test_result=oneway.test(Observed~type, data=compare_richness[[i]], var.equal = FALSE)
  result1[i]=test_result$p.value
}




plot_with_three_cores


biomes=unique(ratio_temp$LABEL)

t.test(ratio_temp%>%filter(LABEL==biomes[5])%>%data.frame()%>%dplyr::select(ratio),mu=1,alternative = "two.sided")

# it looks only for the grassland the response ratio is significant
# if we do some re sampling









# if we make the comparison based on the pooled data rather than the individual plots
# cores from the cropland and cores from the plots  the land use type analogous land cover type
# if the plots have the same historical land use type and are located within the same site we can pooled the data


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

site_id=unique(df4$Site)

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




  












  else{
    
  
    # when the historical have several types
    # when i=5, we have forest and wetlands for the plots selected
    # select the samples based on the different 
    
    ratio=numeric()
    
    historical_type%>%unique()->unique_historical_land
    
    for (j in 1:2)
    {
      
    df4%>%filter(plotIDM%in%select_plot$plotIDM)%>%filter(historical_type==unique_historical_land[j])%>%pull(plotIDM)->sub_type_plot
    
    subset_samples(modified_data,plotIDM%in%sub_type_plot)->sub_modified_data
    
    all_land_cover%>%filter(broad_type%in%unique_historical_land[j])%>%pull(type)->select_type_for_nature
    
    subset_samples(natural_data,Site==site[i]&type%in%select_type_for_nature)->natural_data_sub# all the plots are wetland
    
    richness_modified=estimate_richness(sub_modified_data, measures = "Observed")
    richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
    
    ratio[j]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
    
    }
    
    
    
    
    
    
    # need to get the plot ID of these
   
     all_land_cover%>%filter(broad_type%in%historical_type)%>%
       filter(broad_type==unique_historical_land[2])%>%
       pull(type)
    
    subset_samples(modified_data,type==)
    
    
    
    
    
    
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
