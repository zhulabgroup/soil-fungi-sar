###

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


