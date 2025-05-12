library(geosphere)
library(dplyr)

# to compare fungal diversity among natural and human-dominated land-systems
# the natural communities are selected to represent the pre-modified land use type for the current croplands
# get the distance between sites

plot_coordinates= readRDS("plot_coordinates.rds")
human_dominated_plots=readRDS("human_dominated_plots.rds")
rare_all_guild_biome=readRDS("rare_all_guild_biome.rds")


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


#to assign the current land cover type a broad types

all_land_cover=sample_data(rare_all_guild_biome)%>%data.frame()%>%distinct(type)%>%
  mutate(broad_type=case_when(type%in%c("evergreenForest", "deciduousForest","mixedForest")~"forest", 
                              type%in%c("emergentHerbaceousWetlands","woodyWetlands")~"wetland", 
                              type%in%c("grasslandHerbaceous","sedgeHerbaceous")~"grassland", 
                              type==c("shrubScrub","dwarfScrub")~"shrub", 
                              TRUE~"other"))

# to see if we can find pair plot for each of the croplands

plotid=human_dominated_plots$plotIDM
site=human_dominated_plots$Site

match_land=list()
for (i in 1:45)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  subset_samples(get(data[m]),type!="cultivatedCrops")->natural_data
  
  Site=human_dominated_plots$Site
  #when no sites are available, we can find the plots based on the most close sites
  nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site
  
  matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)# for these plots, they do not have shared sites
  
  
  if(i%in%c(matching_elements,24))
  {
    subset_samples(natural_data,Site==near_site)->natural_data_sub
    
  }
  else if(i==12)#can not find analogous plots
  {
    subset_samples(natural_data,Site==site[i])->natural_data_sub
  }
  
  else{
    subset_samples(natural_data,Site==site[i])->natural_data_sub
    
  }
  match_land[[i]]=natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type,historical_type)%>%distinct()
}
# to see if all the sites have analogous plotss for the comparision
# all have pair sites except for i=12

compare=list()
for (i in 1:45)
{
  compare[[i]]=list(match_land[[i]],human_dominated_plots$historical_type[i])
}


site=human_dominated_plots$Site

plotid=human_dominated_plots$plotIDM


sample_data(rare_all_guild_biome)%>%data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%left_join(human_dominated_plots,by="Site")%>%filter(!is.na(plotIDM))%>%head(45)->temp_data


biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

biome_group=list()
for (i in 1:4){
  human_dominated_plots%>%left_join(temp_data%>%dplyr::select(plotIDM,LABEL),by="plotIDM")%>%mutate(code=1:45)%>%
    group_by(LABEL)%>%
    filter(LABEL==biomes[i])%>%pull(code)->biome_group[[i]]
}



# updated selection of the paired plots 

data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")
data_EM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild_biome, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "epiphyte")


species_com_guild_adjust_natural=list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = ""))
  
  
  species_com=list()
  for (i in c(1:45))
  {
    #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
    
    subset_samples(get(data[m]),plotIDM%in%human_dominated_plots$plotIDM)%>%
      subset_samples(plotIDM==plotid[i])->modified_data
    # defined the historical land cover type for the croplands
    # for comparison, we also need to select current land cover that have the same historical land cover
    human_dominated_plots%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type0
    
    # the land cover type to be selected from the adjacent natural plots
    all_land_cover%>%filter(broad_type%in%historical_type0)%>%pull(type)->select_type
    unique(human_dominated_plots$Site)->modified_site
    
    subset_samples(get(data[m]),type!="cultivatedCrops")->natural_data
    
    
    sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
    
    #to see if the two data sets have shared sites
    #setdiff(modified_site,natural_site)->sites_no_in_nature
    
    # find the plot within the same site
    # all sites in the 
    Site=human_dominated_plots$Site
    #when no sites are available, we can find the plots based on the most close sites
    nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site
    
    matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)# for these plots, they do not have shared sites
    
    if(i%in%c(matching_elements,24))
    {
      subset_samples(natural_data,Site==near_site&type%in%select_type&historical_type%in%historical_type0)->natural_data_sub
      
    }
    else if(i==12)#can not find analogous plots
    {
      subset_samples(natural_data,Site==site[i])->natural_data_sub
    }
    
    else{
      subset_samples(natural_data,Site==site[i]&type%in%select_type&historical_type%in%historical_type0)->natural_data_sub
      
    }
    #when i=12, the historical land cover type was grassland but all current plots are deciduous forest
    # the nearest plots are also forest, so we used forest as the analogous land use type
    #when i=24 the historical land cover was forest but currently all are woodyWetlands
    #and we get the adjacent sites
    #for the LAJA site we do not have historical data and we assumed the historical land use type for all the plots was forest
    #n_natural_sample=nsamples(natural_data_sub)
    #n_modified_sample=nsamples(modified_data)
    # need to get the land use type for the natural data
    natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type)->natural_land_use
    modified_data%>%sample_data()%>%data.frame()%>%dplyr::select(type)->modified_land_use
    #combined_land_use=rbind(modified_land_use,natural_land_use)
    #richness_modified=estimate_richness(modified_data, measures = "Observed")
    #richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
    
    #otu_table(modified_data)%>%data.frame()->otu_modified
    #otu_table(natural_data_sub)%>%data.frame()->otu_natural
    
    merged_physeq <- merge_phyloseq(modified_data, natural_data_sub)
    
    #sample_data_df <- sample_data(merged_physeq )
    
    # Create a new column (e.g., 'new_column') with your data
    # Let's assume 'new_data' is a vector containing the new data, matching the order of samples
    
    #species_com[[i]]=bind_rows(otu_modified,otu_natural)%>%bind_cols(combined_land_use)
    species_com[[i]]=merged_physeq 
    #compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)
    
  }
  
  species_com_guild_adjust_natural[[m]]=species_com
}


saveRDS(species_com_guild_adjust_natural,file="species_com_guild_adjust_natural.rds")

species_com_guild_adjust_natural=readRDS("species_com_guild_adjust_natural.rds")


#

##
sample_data(rare_all_guild_biome)%>%data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%left_join(human_dominated_plots,by="Site")%>%filter(!is.na(plotIDM))%>%head(45)->temp_data

biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

biome_group=list()
for (i in 1:4){
  human_dominated_plots%>%left_join(temp_data%>%dplyr::select(plotIDM,LABEL),by="plotIDM")%>%mutate(code=1:45)%>%
    group_by(LABEL)%>%
    filter(LABEL==biomes[i])%>%pull(code)->biome_group[[i]]
}

plot_number=numeric()
for(i in 1:4)
{
  plot_number[i]=length(biome_group[[i]])
}


species_com_guild=readRDS("species_com_guild.rds")


richness_compare_crop_nature_guild=list()
for(m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = ""))
  
  
  richness_ratio_with_rarefaction=list()
  for (j in 1:4)
  {
    #cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = ""))
    
    richness_pair=matrix(ncol=4,nrow=plot_number[j])
    for (i in 1:plot_number[j])
    {
      sub= species_com_guild_adjust_natural[[m]][biome_group[[j]]][[i]] #bio me first plot with croplands for the first biome
      
      #need to convert the binary category of the pairwise
      sample_data_df <- data.frame(sample_data(sub))
      sample_data_df$type <- ifelse(sample_data_df$type != "cultivatedCrops", "Natural", sample_data_df$type)
      sample_data(sub) <- sample_data(sample_data_df)
      
      # get the number of the sample of the modified and the natural communities
      
      data_natural=subset_samples(sub,type=="Natural")
      data_crop=subset_samples(sub,type!="Natural")
      
      data_natural =transform_sample_counts(data_natural, function(x) ifelse(x>0, 1, 0))
      data_crop =transform_sample_counts(data_crop, function(x) ifelse(x>0, 1, 0))
      
      n_sample_natural= nsamples(data_natural)
      n_sample_crop= nsamples(data_crop)
      
      # richness in the croplands
      otu_tb <- otu_table(data_crop)%>%as.matrix()
      # Calculate total unique richness across samples
      # This involves summing the unique non-zero values across all samples
      total_occurrences <- colSums(otu_tb)
      # Count the number of species with more than one occurrence
      richness_crop <- table(total_occurrences >0)[2]%>%as.numeric()
      
      # richness in the natural communities
      # select the same number of the sample from the natural communities to get the richness
      ot=t(otu_table(data_natural))%>%data.frame()
      
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)
      
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, n_sample_crop, length.out=n_sample_crop)), nboot=100,se=TRUE)
      
      richness_nature=out3$iNextEst$size_based%>%filter(t==n_sample_crop)%>%dplyr::select(qD,   qD.LCL,   qD.UCL )
      
      richness_pair[i,1] =richness_nature[,2]#lower
      richness_pair[i,2] =richness_nature[,1]#the mean
      richness_pair[i,3] =richness_nature[,3]#the upper
      
      richness_pair[i,4] =richness_crop
    }  
    
    richness_ratio_with_rarefaction[[j]]=richness_pair
  }   
  #   
  
  richness_compare_crop_nature_guild[[m]]=richness_ratio_with_rarefaction
}




saveRDS(richness_compare_crop_nature_guild,file="richness_compare_crop_nature_guild.rds")


human_dominated_plots%>%mutate(plotid=1:45)%>%dplyr::select(Site,plotid,plotIDM)->temp


biome_select=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
               "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")


do.call(rbind,richness_ratio_with_rarefaction  )%>%data.frame()%>%
  mutate(plotid=unlist(biome_group))%>%mutate(biomes=rep(biome_select,times=plot_number))%>%
  left_join(temp,by="plotid")%>%rename_all(~paste0(c("low","mean_nature","up","mean_crop","plotid","biome","site","plotIDM")))->data_mean_richness_biome


data_mean_richness_biome%>%dplyr::select(mean_nature,biome, site,plotIDM)%>%
  rename(mean_rich=mean_nature)%>%
  rbind(data_mean_richness_biome%>%dplyr::select(mean_crop,biome, site,plotIDM)%>%rename(mean_rich=mean_crop))%>%
  mutate(type=rep(c("nature","crop"),each=45))->data_mean_richness_biome_model


###

