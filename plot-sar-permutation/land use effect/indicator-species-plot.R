#create the plot showing the indicator species 

biome_group=readRDS("biome_group.rds")
species_com_guild_adjust_natural=readRDS("species_com_guild_adjust_natural.rds")

# add a new column 

species_com_guild_adjust_natural%>%sample_data()%>%data.frame()%>%dplyr::select(type)->broad_type

species_com_guild_adjust_natural[[1]][biome_group[[i]]]->tem_select_data

merged_phyloseq <- do.call(merge_phyloseq, tem_select_data)
#need to have the same sample efforts

merged_phyloseq%>%sample_data()%>%data.frame()%>%dplyr::select(type)->broad_type

broad_type%>%mutate(broad_type = if_else(broad_type != "cultivatedCrops", "nature", type))%>%
  dplyr::select(broad_type)->broad_type

row.names(broad_type)=rownames(sample_data(merged_phyloseq))
broad_type=sample_data(broad_type)
  
merged_phyloseq=merge_phyloseq(merged_phyloseq,broad_type)

sample_meta <- sample_data(merged_phyloseq )

group_col <- "broad_type"

min_size <- min(table(sample_meta[[group_col]]))

subset_samples_list <- lapply(unique(sample_meta[[group_col]]), function(broad_type) {
  samples_in_group <- rownames(sample_meta[sample_meta[[group_col]] == broad_type, ])
  sampled <- sample(samples_in_group, size = min_size, replace = FALSE)
  return(sampled)
})

selected_samples <- unlist(subset_samples_list)
ps_subset <- prune_samples(selected_samples, merged_phyloseq)

#based on the sample size to subset different types

ps_subset%>%sample_data()%>%data.frame()%>%
  dplyr::select(broad_type)%>%pull(broad_type)->select_type

ps_subset%>%otu_table()%>%as.matrix()->
  species_presence
species_presence[species_presence>0]=1

isa_result[[i]] <- multipatt(species_presence, select_type, func = "IndVal.g", control = how(nperm = 999))

#get the life style for each otu

tax_table(rare_all_assign_type)%>%data.frame()->temp

temp%>%select(ta2)%>%bind_cols(rownames(temp))%>%
 rename_all(~paste0(c("primary_lifestyle","sp")))->species_life_style
  



#read in the data with the same sampling effort
#different land cover types were evaluated
isa_result=readRDS("isa_result.RData")

#to see how many habitat types are there for each biome
# for the first biome, there were six groups
habitat_type=list()
result_biome=list()
for (i in 1:4)
  {
 isa_result[[i]]$sign[which(isa_result[[i]]$sign$p.value < 0.05), ]%>%
  mutate(sp=rownames(.))%>%
    left_join(species_life_style,by="sp")%>%
    filter(primary_lifestyle%in%c("arbuscular_mycorrhizal",
                                  "ectomycorrhizal",
                                  "epiphyte",
                                  "litter_saprotroph",
                                  "mycoparasite",
                                  "plant_pathogen",
                                  "soil_saprotroph",
                                  "wood_saprotroph"))->result_biome[[i]]
  unique(result_biome[[i]]$index)%>%sort()->habitat_type[[i]]
}

# for each biome, get the habitat name rather than a value

land_cover=list()
for (m in 1:4){
  
  cover_temp=list()
for (i in 1:length(habitat_type[[m]]))
{
  isa_result[[m]]$sign[which(isa_result[[m]]$sign$p.value < 0.05), ]%>%filter(index==habitat_type[[m]][i])%>%
    select(-index)%>%
    select(where(~ all(. == 1)))%>%colnames()->cover_temp[[i]]
#bind all the cover types for the first biome
}
  land_cover[[m]]=do.call(rbind,cover_temp)%>%data.frame()
}
  
# to determine the land cover types for each biome by pasting the combination
  land_cover_biome=list()
  for (m in 1:4)
    {
    apply(land_cover[[m]], 1, function(row) paste(row, collapse = "_"))->land_cover_biome[[m]]
    
  }
  


#for each type and each land cover type, to get the number of species of different fungal guilds
  total_richness_biome=list()
  pie_data_biome=list()
for (m in 1:4)
{
  total_richness=numeric()
  pie_data=list()
  
  for (i in 1:length(habitat_type[[m]]))
  {
    result_biome[[m]]%>%filter(index==habitat_type[[m]][i])%>%count(primary_lifestyle)%>%
      mutate(percent=n/sum(n))->pie_data[[i]]
    result_biome[[m]]%>%filter(index==habitat_type[[m]][i])%>%count(primary_lifestyle)%>%
      mutate(percent=n/sum(n))%>%pull(n)%>%sum()->total_richness[i]
  }
  pie_data_biome[[m]]=pie_data
  total_richness_biome[[m]]=total_richness
}


#to create the plots and need to bind all the data based on different land cover
#to see how many guilds are associated with each type
nrow_biome=list()
for (m in 1:4)
  {
  nrows <- sapply(pie_data_biome[[m]], nrow)# the number of each land cover
  nrow_biome[[m]]=nrows
}

#combinations of different land cover types
#need to know what does each value mean 
comb_data_biome=list()
for (m in 1:4)
  {
  comb_data=do.call(rbind,pie_data_biome[[m]])%>%
    mutate(cover=rep(land_cover_biome[[m]],times=nrow_biome[[m]]))
  comb_data_biome[[m]] =comb_data
}



colors <- c("#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557", "#F4A261", "#2A9D8F", "#264653")

#change the names
#to see the unique color
do.call(rbind,comb_data_biome)->dk

custom_labels <- c(s.cultivatedCrops_s.cultivatedCrops = "Crop",
                   s.deciduousForest_s.deciduousForest = "DF", 
                   s.evergreenForest_s.evergreenForest="EF",
                   s.cultivatedCrops_s.deciduousForest= "C:DF",
                   s.cultivatedCrops_s.evergreenForest="C:EF",
                   s.deciduousForest_s.evergreenForest="DF:EF",
                   s.cultivatedCrops_s.cultivatedCrops_s.cultivatedCrops_s.cultivatedCrops="Crop",
                   s.deciduousForest_s.deciduousForest_s.deciduousForest_s.deciduousForest="DF",
                   s.evergreenForest_s.evergreenForest_s.evergreenForest_s.evergreenForest="EF",
                   s.mixedForest_s.mixedForest_s.mixedForest_s.mixedForest="MF",
                   s.woodyWetlands_s.woodyWetlands_s.woodyWetlands_s.woodyWetlands="WL",
                   s.cultivatedCrops_s.evergreenForest_s.cultivatedCrops_s.evergreenForest="Crop:EF",
                   s.cultivatedCrops_s.woodyWetlands_s.cultivatedCrops_s.woodyWetlands="Crop:WL",
                   s.deciduousForest_s.evergreenForest_s.deciduousForest_s.evergreenForest="DF:EF",
                   s.deciduousForest_s.mixedForest_s.deciduousForest_s.mixedForest="DF:MF",
                   s.evergreenForest_s.mixedForest_s.evergreenForest_s.mixedForest="EF:MF",
                   s.evergreenForest_s.woodyWetlands_s.evergreenForest_s.woodyWetlands="EF:WL",
                   s.deciduousForest_s.evergreenForest_s.mixedForest_s.deciduousForest="DF:EF:MF",
                   s.evergreenForest_s.mixedForest_s.woodyWetlands_s.evergreenForest="EF:MF:WL",
                   s.deciduousForest_s.evergreenForest_s.mixedForest_s.woodyWetlands="DF:EF:MF:WL",
                   s.grasslandHerbaceous_s.grasslandHerbaceous="GH",
                   s.cultivatedCrops_s.grasslandHerbaceous="Crop:GH",
                   s.deciduousForest_s.grasslandHerbaceous="DF:GH",
                   s.cultivatedCrops="Crop",
                   s.evergreenForest="EF"
                   )


number_row=c(3,5,3,2)

pp_indicator_result=list()
for (m in 1:4){

  pp_indicator_result[[m]]=ggplot(comb_data_biome[[m]], aes(x=1, y=percent, fill=primary_lifestyle)) +
  ## geom_col is geom_bar(stat = "identity")(bit shorter)
  ## use color = black for the outline
  geom_col(width=5,position="fill", color = "black",size=0.25)+
  coord_polar("y", start=0)+
  geom_text(aes(x = 4.5, label = paste0(round(percent*100), "%")), size=2.5, 
  position = position_stack(vjust = 0.5))+
  facet_wrap(~cover,ncol=number_row[m],labeller = labeller(cover = custom_labels))+
    xlab(1:6)+
    ylab("")+
  theme(legend.position ="bottom",
        legend.key.size = unit(0.3, "cm"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.background = element_blank(),
        panel.spacing.y  = unit(0.2, "lines"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5))+
  scale_fill_manual("Guild",breaks=c("arbuscular_mycorrhizal",
                                         "ectomycorrhizal",
                                         "epiphyte",
                                         "litter_saprotroph",
                                         "mycoparasite",
                                         "plant_pathogen",
                                         "soil_saprotroph",
                                         "wood_saprotroph"),
                         labels=c("AM","EM","Epiphyte","Litter saprotroph",
                                  "Parasite","Plant pathogen","Soil saprotroph","Wood saprotroph"),
                    values =colors )+
  theme(strip.text = element_text(size = 12, color = "black"),  # Customize facet titles
    strip.background = element_rect(fill = "white")  # Customize background color of facet titles
  )
}

  
plot_grid(pp_indicator_result[[1]],
          pp_indicator_result[[2]],
          pp_indicator_result[[3]],
          pp_indicator_result[[4]],ncol=2)


## for the broad classification of the types

isa_result=readRDS("isa_result_broad.RData")

#to see how many habitat types are there for each biome
# for the first biome, there were six groups
habitat_type=list()
for (i in 1:4)
{
  isa_result[[i]]$sign[which(isa_result[[i]]$sign$p.value < 0.05), ]%>%
    mutate(sp=rownames(.))%>%
    left_join(species_life_style,by="sp")%>%
    filter(primary_lifestyle%in%c("arbuscular_mycorrhizal",
                                  "ectomycorrhizal",
                                  
                                 
                                  "plant_pathogen",
                                  "soil_saprotroph"
                                  ))->result_biome[[i]]
  unique(result_biome[[i]]$index)%>%sort()->habitat_type[[i]]
}

# for each biome, get the habitat name rather than a value

land_cover=list()
for (m in 1:4){
  
  cover_temp=list()
  for (i in 1:length(habitat_type[[m]]))
  {
    isa_result[[m]]$sign[which(isa_result[[m]]$sign$p.value < 0.05), ]%>%filter(index==habitat_type[[m]][i])%>%
      select(-index)%>%
      select(where(~ all(. == 1)))%>%colnames()->cover_temp[[i]]
    #bind all the cover types for the first biome
  }
  land_cover[[m]]=do.call(rbind,cover_temp)%>%data.frame()
}

# to determine the land cover types for each biome by pasting the combination
land_cover_biome=list()
for (m in 1:4)
{
  apply(land_cover[[m]], 1, function(row) paste(row, collapse = "_"))->land_cover_biome[[m]]
  
}



#for each type and each land cover type, to get the number of species of different fungal guilds
total_richness_biome=list()
pie_data_biome=list()
for (m in 1:4)
{
  total_richness=numeric()
  pie_data=list()
  
  for (i in 1:length(habitat_type[[m]]))
  {
    result_biome[[m]]%>%filter(index==habitat_type[[m]][i])%>%count(primary_lifestyle)%>%
      mutate(percent=n/sum(n))->pie_data[[i]]
    result_biome[[m]]%>%filter(index==habitat_type[[m]][i])%>%count(primary_lifestyle)%>%
      mutate(percent=n/sum(n))%>%pull(n)%>%sum()->total_richness[i]
  }
  pie_data_biome[[m]]=pie_data
  total_richness_biome[[m]]=total_richness
}


#to create the plots and need to bind all the data based on different land cover
#to see how many guilds are associated with each type
nrow_biome=list()
for (m in 1:4)
{
  nrows <- sapply(pie_data_biome[[m]], nrow)# the number of each land cover
  nrow_biome[[m]]=nrows
}

#combinations of different land cover types
#need to know what does each value mean 
comb_data_biome=list()
for (m in 1:4)
{
  comb_data=do.call(rbind,pie_data_biome[[m]])%>%
    mutate(cover=rep(land_cover_biome[[m]],times=nrow_biome[[m]]))
  comb_data_biome[[m]] =comb_data
}



colors <- c("#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557", "#F4A261", "#2A9D8F", "#264653")

#change the names
#to see the unique color
do.call(rbind,comb_data_biome)->dk

custom_labels <- c(s.cultivatedCrops_s.cultivatedCrops = "Crop",
                   s.deciduousForest_s.deciduousForest = "DF", 
                   s.evergreenForest_s.evergreenForest="EF",
                   s.cultivatedCrops_s.deciduousForest= "C:DF",
                   s.cultivatedCrops_s.evergreenForest="C:EF",
                   s.deciduousForest_s.evergreenForest="DF:EF",
                   s.cultivatedCrops_s.cultivatedCrops_s.cultivatedCrops_s.cultivatedCrops="Crop",
                   s.deciduousForest_s.deciduousForest_s.deciduousForest_s.deciduousForest="DF",
                   s.evergreenForest_s.evergreenForest_s.evergreenForest_s.evergreenForest="EF",
                   s.mixedForest_s.mixedForest_s.mixedForest_s.mixedForest="MF",
                   s.woodyWetlands_s.woodyWetlands_s.woodyWetlands_s.woodyWetlands="WL",
                   s.cultivatedCrops_s.evergreenForest_s.cultivatedCrops_s.evergreenForest="Crop:EF",
                   s.cultivatedCrops_s.woodyWetlands_s.cultivatedCrops_s.woodyWetlands="Crop:WL",
                   s.deciduousForest_s.evergreenForest_s.deciduousForest_s.evergreenForest="DF:EF",
                   s.deciduousForest_s.mixedForest_s.deciduousForest_s.mixedForest="DF:MF",
                   s.evergreenForest_s.mixedForest_s.evergreenForest_s.mixedForest="EF:MF",
                   s.evergreenForest_s.woodyWetlands_s.evergreenForest_s.woodyWetlands="EF:WL",
                   s.deciduousForest_s.evergreenForest_s.mixedForest_s.deciduousForest="DF:EF:MF",
                   s.evergreenForest_s.mixedForest_s.woodyWetlands_s.evergreenForest="EF:MF:WL",
                   s.deciduousForest_s.evergreenForest_s.mixedForest_s.woodyWetlands="DF:EF:MF:WL",
                   s.grasslandHerbaceous_s.grasslandHerbaceous="GH",
                   s.cultivatedCrops_s.grasslandHerbaceous="Crop:GH",
                   s.deciduousForest_s.grasslandHerbaceous="DF:GH",
                   s.cultivatedCrops="Crop",
                   s.evergreenForest="EF"
)

custom_labels=c(s.cultivatedCrops="Modified",
                s.nature="Natural")

number_row=c(3,5,3,2)

number_row=c(2,2,2,2)

pp_indicator_result=list()
for (m in 1:4){
  
  pp_indicator_result[[m]]=ggplot(comb_data_biome[[m]], aes(x=1, y=percent, fill=primary_lifestyle)) +
    ## geom_col is geom_bar(stat = "identity")(bit shorter)
    ## use color = black for the outline
    geom_col(width=5,position="fill", color = "black",size=0.25)+
    coord_polar("y", start=0)+
    geom_text(aes(x = 5, label = paste0(round(percent*100), "%")), size=2.5, 
              position = position_stack(vjust = 0.5))+
    facet_wrap(~cover,ncol=number_row[m],labeller = labeller(cover = custom_labels))+
    xlab("")+
    ylab("")+
    theme(legend.position ="bottom",
          legend.key.size = unit(0.3, "cm"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0, -1, 0), "cm"),
          strip.background = element_blank(),
          panel.spacing.y  = unit(0.2, "lines"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_blank(),
          plot.title = element_text(size = 15, hjust = 0.5))+
    scale_fill_manual("Guild",breaks=c("arbuscular_mycorrhizal",
                                       "ectomycorrhizal",
                                       
                                       
                                       "plant_pathogen",
                                       "soil_saprotroph"
                                       ),
                      labels=c("AM","EM",
                               "Plant pathogen","Soil saprotroph"),
                      values =colors )+
    theme(strip.text = element_text(size = 12, color = "black"),  # Customize facet titles
          strip.background = element_rect(fill = "white")  # Customize background color of facet titles
    )
}


plot_grid(pp_indicator_result[[1]],
          pp_indicator_result[[2]],
          pp_indicator_result[[3]],
          pp_indicator_result[[4]],ncol=2)






