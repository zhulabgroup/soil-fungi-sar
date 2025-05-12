#
setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/land-use-effect")
library(dplyr)
library(phyloseq)
library(vegan)

species_com_guild_adjust_natural=readRDS("species_com_guild_adjust_natural.rds")

# to see the overall fungi beta diversity among the human modified and natural plots
# 
beta=numeric()
for (i in 1:45){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  #to see the number of the natural and modified plots
  species_com_guild_adjust_natural[[1]][[i]]%>%sample_data()%>%count(type)->dk
  
  #convert the otu table into 0-1 data
  species_com_guild_adjust_natural[[1]][[i]]%>%otu_table()%>%data.frame()%>%
    mutate_all(~ ifelse(. > 0, 1, 0))->d
  
  # determine the beta diversity of the data
  
  beta_jaccard <- vegdist(d, method = "bray")# this quantifys the dissimilarity based on soren simiarity
  
  jaccard_matrix <- as.matrix(beta_jaccard)
  
  # select the distance between natural and modified plots
  jaccard_matrix[(1:dk$n[1]),(dk$n[1]+1):sum(dk$n)]->natural_modified_pair
  # determine the mean of the dissimilarity 
  beta[i]= as.vector(natural_modified_pair)%>%mean()
  
}

beta_diversity_natural_crop=beta
saveRDS(beta_diversity_natural_crop,file="beta_diversity_natural_crop.rds")

