###
library(indicspecies)
library(phyloseq)
library(dplyr)

biome_group=readRDS("biome_group.rds")
species_com_guild_adjust_natural=readRDS("species_com_guild_adjust_natural.rds")

isa_result=list()
for(i in 1:4){
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
  
  species_com_guild_adjust_natural[[1]][biome_group[[i]]]->tem_select_data
  
  merged_phyloseq <- do.call(merge_phyloseq, tem_select_data)
  
  merged_phyloseq%>%sample_data()%>%data.frame()%>%
    dplyr::select(type)%>%pull(type)->select_type
  
  merged_phyloseq%>%otu_table()%>%as.matrix()->
    species_presence
  species_presence[species_presence>0]=1
  
  isa_result[[i]] <- multipatt(species_presence, select_type, func = "IndVal.g", control = how(nperm = 999))
}


saveRDS(isa_result,file="isa_result.RData")

set.seed(123)
sensitivity_broad=list()
for(j in 1:15){
  cat('\r',paste(paste0(rep("*", round(j/ 1, 0)), collapse = ''), j, collapse = ''))# informs the processing
  
  
  isa_result_broad=list()
  for(i in 1:4){
    
    species_com_guild_adjust_natural[[1]][biome_group[[i]]]->tem_select_data
    
    merged_phyloseq <- do.call(merge_phyloseq, tem_select_data)
    
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
    
    isa_result_broad[[i]] <- multipatt(species_presence, select_type, func = "IndVal.g", control = how(nperm = 999))
  }
  sensitivity_broad[[j]]=isa_result_broad
}


saveRDS(sensitivity_broad,file="sensitivity_broad.rds")




