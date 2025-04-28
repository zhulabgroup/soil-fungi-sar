# to compare the richness and composition of fungal communities among land use types
# the natural communities are selected to represent the pre-modifed land use type for the current croplands

data=c("rare_all_guild_biome","data_AM","data_EM","data_plapat","data_soilsap","data_littersap","data_woodsap","data_epiphy","data_para")
data_EM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild_biome, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild_biome, primary_lifestyle == "epiphyte")




 df4=readRDS("df4.rds")

  species_com=list()
  for (i in c(1:45))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
    
    subset_samples(get(data[m]),plotIDM%in%df4$plotIDM)%>%
      subset_samples(plotIDM==plotid[i])->modified_data
    
    modified_data <- prune_samples(sample_sums(modified_data) > 0, modified_data)
    
    df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type
    # the land cover to be selected from the adjacent natural plots
    all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type
    unique(df4$Site)->modified_site
    subset_samples(get(data[m]),type!="cultivatedCrops")->natural_data
    
    natural_data <- prune_samples(sample_sums(natural_data) > 0, natural_data)
    
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
    #n_natural_sample=nsamples(natural_data_sub)
    #n_modified_sample=nsamples(modified_data)
    # need to get the land use type for the natural data
    #natural_data_sub%>%sample_data()%>%data.frame()%>%dplyr::select(type)->natural_land_use
    #modified_data%>%sample_data()%>%data.frame()%>%dplyr::select(type)->modified_land_use
    #combined_land_use=rbind(modified_land_use,natural_land_use)
    #richness_modified=estimate_richness(modified_data, measures = "Observed")
    #richness_natural=estimate_richness(natural_data_sub, measures = "Observed")
    
    otu_table(modified_data)%>%data.frame()->otu_modified
    otu_table(natural_data_sub)%>%data.frame()->otu_natural
    
    species_com[[i]]=bind_rows(otu_modified,otu_natural)
    #compare_richness[[i]]=bind_rows(richness_modified,richness_natural)%>%bind_cols(combined_land_use)
  }
 
  
  
  #need to bind all the list
  
  do.call(rbind,compare_richness)->d
  # we don't have data on the shrublands for the comparision
  # write a function for a given the guilds
  
  #
  species_com_guild=readRDS("species_com_guild.rds")
  #species_com_guild and species_com_guild_adjust_natural should be identical
  
  
  
  # to compared the richness by rarefying the richness
  

  
  sample_data(rare_all_guild_biome)%>%data.frame()%>%
    dplyr::select(Site,LABEL)%>%distinct()%>%left_join(df4,by="Site")%>%filter(!is.na(plotIDM))%>%head(45)->temp_data
  
  biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")
  
  biome_group=list()
  for (i in 1:4){
    df4%>%left_join(temp_data%>%dplyr::select(plotIDM,LABEL),by="plotIDM")%>%mutate(code=1:45)%>%
      group_by(LABEL)%>%
      filter(LABEL==biomes[i])%>%pull(code)->biome_group[[i]]
  }
  
  
  
 # compare the richness between the natural and the modified plots 
  
  plot_number=numeric()
    for(i in 1:4)
      {
      plot_number[i]=length(biome_group[[i]])
    }
    
  dd=list()
  for(m in 2:9)
  {
    cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = ""))
    
  
  richness_ratio_with_rarefaction=list()
  for (j in 1:4)
  {
    #cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = ""))
    
    richness_pair=matrix(ncol=4,nrow=plot_number[j])
      for (i in 1:plot_number[j])
      {
       sub=species_com_guild[[m]][biome_group[[j]]][[i]]#the first plot with croplands for the first biome
     
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
   dd[[m]] =richness_ratio_with_rarefaction
  }
  
 
  
 dd[[1]]=richness_ratio_with_rarefaction  
 
 richness_ratio_with_rarefaction_guild[[1]]=richness_ratio_with_rarefaction  

 saveRDS(richness_ratio_with_rarefaction,file="richness_ratio_with_rarefaction.rds")
 
 richness_ratio_with_rarefaction =readRDS("richness_ratio_with_rarefaction.rds")
    
  #to bind all the data
 
 df4%>%mutate(plotid=1:45)%>%dplyr::select(Site,plotid,plotIDM)->temp
 
 
 biome_select=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
                "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")
 

 
 do.call(rbind,richness_ratio_with_rarefaction  )%>%data.frame()%>%
   mutate(plotid=unlist(biome_group))%>%mutate(biomes=rep(biome_select,times=plot_number))%>%
   left_join(temp,by="plotid")%>%rename_all(~paste0(c("low","mean_nature","up","mean_crop","plotid","biome","site","plotIDM")))->data_mean_richness_biome
 
 saveRDS(data_mean_richness_biome,file="mean_richness_biome.rds")
 
 data_mean_richness_biome=readRDS("mean_richness_biome.rds")
 
 
 #construct a model to look at the effect
 
 
 data_mean_richness_biome%>%dplyr::select(mean_nature,biome, site,plotIDM)%>%
   rename(mean_rich=mean_nature)%>%
   rbind(data_mean_richness_biome%>%dplyr::select(mean_crop,biome, site,plotIDM)%>%rename(mean_rich=mean_crop))%>%
   mutate(type=rep(c("nature","crop"),each=45))->data_mean_richness_biome_model
 
 
 mod=lm(mean_rich~type,data=data_mean_richness_biome_model%>%filter(biome==biome_select[4]))
   
 p=numeric()
 for(i in 1:4)
   {
 d=data_mean_richness_biome_model%>%filter(biome==biome_select[i])
 mod=aov(mean_rich~type,data=d)
 result=summary(mod)
 p[i]=result[[1]][["Pr(>F)"]][1]
 }
  
 data_mean_richness_biome_model$type=factor(data_mean_richness_biome_model$type,levels = c("nature","crop"))
 
 #lack the codes that produce the the guild-specific richness
 data_mean_richness_biome_model_guild=readRDS("data_mean_richness_biome_model_guild.rds")
 
 # construct the model for each guild
 
 data_mean_richness_biome_model_guild[[1]]$type=factor(data_mean_richness_biome_model_guild[[1]]$type,levels = c("crop","nature"))
 #the order of the factor will affect the out put
 
 dd=data_mean_richness_biome_model_guild[[2]]
 
 dd$mean_rich[is.na(dd$mean_rich)] <- 0
 #replace the initial data after converting na to 0
 dd=data_mean_richness_biome_model_guild[[2]]
 
 
 convert_to_integer <- function(df) {
   # Apply as.integer to numeric columns only
   df[] <- lapply(df, function(col) {
     if (is.numeric(col)) {
       as.integer(col)
     } else {
       col
     }
   })
   return(df)
 }
 
 df_list <- lapply(df_list, function(df) {
   df[is.na(df)] <- 0
   return(df)
 })
 
 my_list_integer <- lapply(data_mean_richness_biome_model_guild, convert_to_integer)
 
 data_mean_richness_biome_model_guild=my_list_integer
 #conver the na to 0, no species found 
 data_mean_richness_biome_model_guild <- lapply(data_mean_richness_biome_model_guild, function(df) {
   df[is.na(df)] <- 0
   return(df)
 })
 
 
 # test the significance of the results based on glm models
 result=list()
 for (m in 1:9)
 {
   mod=glmer(mean_rich~type+biome+(1 | site),family = poisson(link = "log"),data=data_mean_richness_biome_model_guild[[m]])
   
   result[[m]]=summary(mod) 
 }
 
 #get the species response ratio for different guilds
 # need to convert the na to zero for the richness
 # for each biome, we get a single value
 
 richness_ratio_rarefy_guild=list()
 for(m in 1:9)
   {
 data_mean_richness_biome_model_guild[[m]]%>%group_by(biome,type)%>%summarise(mean_richness=mean(mean_rich),sd_richness=sd(mean_rich))->temp
 
richness_ratio_rarefy=numeric()
 for (i in 1:4)
 {
   temp%>%filter(biome==biome_select[i])%>%data.frame()->df1
    richness_ratio_rarefy[i]=df1[1,3]/df1[2,3]
 }
   
 richness_ratio_rarefy_guild[[m]]=richness_ratio_rarefy
 }
 
 unlist(richness_ratio_rarefy_guild)%>%data.frame()%>%mutate(guild=rep(guild,times=4))%>%
   mutate(biome=rep(biome_select,times=9))%>%rename_all(~paste0(c("ratio_biome","guild","biome")))%>%
   mutate(guild_biome=paste(guild,"_",biome))->biome_level_richness_ratio
 
 # if we get the individual ratio and then get the mean
 
 
 richness_ratio_rarefy_guild_site=list()
 for(m in 1:9)
 {
   data_mean_richness_biome_model_guild[[m]]->temp
   
   
   site_level_ratio=list()
   for (i in 1:4)
   {
   temp%>%filter(biome==biome_select[i])%>%dplyr::select(mean_rich,type,plotIDM)->d
   
     site_level_ratio[[i]]= reshape(d, idvar = "plotIDM", timevar = "type", direction = "wide")%>%mutate(ratio=mean_rich.crop/mean_rich.nature)%>%pull(ratio)
   
   }
   
   richness_ratio_rarefy_guild_site[[m]]=unlist(site_level_ratio)
 }
 
 
 guild=c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para")
 
 do.call(cbind,richness_ratio_rarefy_guild_site)%>%data.frame()%>%rename_all(~paste0(guild))%>%
   mutate(biome=rep(biome_select,times=c(plot_number)))%>%melt()%>%group_by(biome,variable)%>%
   summarise(mean_ratio=mean(value),sd_ratio=sd(value))%>%mutate(guild_biome=paste(variable,"_",biome))%>%
   rename(guild=variable)->biome_site_level_richness_ratio
 
 # we can choose this for further modeling work
 # to create the plots
 
 biome_site_level_richness_ratio=readRDS("biome_site_level_richness_ratio.rds")
 
 pp=list()
 for (i in 1:4)
 {
   pp[[i]]=ggplot(biome_site_level_richness_ratio%>%filter( biome==biome_select[i]), aes(x=1:9,y = mean_ratio)) +
     geom_bar(stat = "identity",position = "dodge",width =0.5) +
     guides(fill="none")+
     geom_errorbar(aes(ymin = mean_ratio-sd_ratio, ymax = mean_ratio+sd_ratio), width = 0.1)+
     geom_hline(yintercept = 1,color="red",linetype="dashed")+
     scale_x_continuous (breaks=1:9,labels = c("All", "AM","EM","Pla. patho.", "Soil sapro.","Litter sapro.","Wood sapro.","Epiphyte","Parasite"))+
     ggtitle("Response ratio")+
     xlab("")+
     ylab("")+
     theme(legend.position = c(0.8,0.75),
           legend.text = element_text(size=8),
           legend.title  = element_text(size=10),
           text = element_text(size = 18),
           plot.title = element_text(size = 15, hjust = 0.5), 
           axis.text.y = element_text(size = 12), 
           axis.text.x = element_text(size = 12,angle=90,hjust=0), 
           axis.title.y = element_text(size = 18), 
           axis.title.x = element_text(size = 18), 
           legend.key.size = unit(0.3, "cm"),
           plot.margin = unit(c(0.3, 0.5, -0.5, 0.1), "cm"),
           panel.background = element_rect(fill = "NA"),
           panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
     geom_text(aes(x=1:9,y=(mean_ratio+sd_ratio)*1.085),label =df_significance[[i]],size=5)+
     coord_flip()
     
 }
 
 

   
 
 
 
 
 plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]])
 
 
 comp_ratio=biome_level_richness_ratio%>%bind_cols(biome_site_level_richness_ratio,by="guild_biome")
 

 comp_ratio%>%dplyr::select(ratio_biome,mean_ratio,guild...2 , biome...3)%>%melt()->df3
 
 pp=list()
 for (i in 1:9)
 {
   pp[[i]]=ggplot(df3%>%filter( guild...2 ==guild[i]), aes(x =  biome...3 , y = value, fill = variable)) +
     geom_bar(stat = "identity",position = "dodge") +
     guides(fill="none")+
     geom_hline(yintercept = 1,color="red",linetype="dashed")+
     scale_x_discrete(labels=c("temperate forest","conifer","grassland","dry forest"))
 }
 
 plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],ncol=3)
 
 mod=lmer(mean_rich~type+biome+(1 | site),data=data_mean_richness_biome_model)
 
 p=numeric()
 for(i in 1:4)
 {
   d=data_mean_richness_biome_model%>%filter(biome==biome_select[i])
   mod=aov(mean_rich~type,data=d)
   result=summary(mod)
   p[i]=result[[1]][["Pr(>F)"]][1]
 }
 
 
 
 
 
 
 
 
 
 
 
     # Print the total unique richness
      
      ##
 data_natural <- prune_samples(sample_names(data_natural)[1:20], data_natural)
 
    sample_names_vector <- sample_names(data_natural)
    
   
    
    
  
    pairwise_combinations <- combn(sample_names_vector, 3, simplify = FALSE)
    
    #based on each to get the samples
    
    two_sample_richness=numeric()
    for (i in 1:length( pairwise_combinations))
    {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
    ps_subset <- prune_samples(sample_names(data_natural) %in%pairwise_combinations[[i]], data_natural)
    # get the total richness in the two samples
    
    ps_subset =transform_sample_counts(ps_subset, function(x) ifelse(x>0, 1, 0))
    
    otu_tb <- otu_table(ps_subset)%>%as.matrix()
    
    # Calculate total unique richness across samples
    # This involves summing the unique non-zero values across all samples
    total_occurrences <- colSums(otu_tb)
    # Count the number of species with more than one occurrence
    two_sample_richness[i] <-table(total_occurrences >0)[2]%>%as.numeric()

    }
        
      
    ot=t(otu_table(data_natural))%>%data.frame()
    
    n=dim(ot)[2]# the number of "sites" for each plot
    ot1=rowSums(ot)
    
    out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, n_sample_crop, length.out=n_sample_crop)), nboot=100,se=TRUE)
    
    richness_nature=out3$iNextEst$size_based%>%filter(t==n_sample_crop)%>%dplyr::select(qD,   qD.LCL,   qD.UCL )
    
    
    
    
    
    
    
       
       merged_phyloseq <- Reduce(function(x, y) merge_phyloseq(x, y), phylo_list)
      # select the samples with richness
      phylo_temp <- prune_samples(sample_sums(merged_phyloseq ) > 0, merged_phyloseq )
      phylo_temp%>%sample_data()%>%data.frame()%>%dplyr::select(type)->land_cover_type
      vege_com=otu_table(phylo_temp)%>%data.frame()%>%bind_cols(land_cover_type)
      ordination <- metaMDS(vege_com[, -ncol(vege_com)], distance = "bray")
      ordination_data[[i]]=ordination$points%>%data.frame()%>%mutate(land_cover_type)
      
    }
    
    
    
    
  
  
  
  
