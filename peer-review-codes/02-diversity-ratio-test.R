# for the updated species response ratio
# considering the historical land use type for natural plots
setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-climate-historical")

biome_select=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
               "Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

species_com_guild_adjust_natural=readRDS("species_com_guild_adjust_natural.rds")#each element represent a phyloseq

#same data species_com_guild
#the data set for each guild and biome
# the data was used for comparing species composition

rare_all_guild_biome=readRDS("rare_all_guild_biome.rds")
#load in the human-modified plots
human_dominated_plots=readRDS("human_dominated_plots.rds")

sample_data(rare_all_guild_biome)%>%data.frame()%>%
  dplyr::select(Site,LABEL)%>%distinct()%>%left_join(human_dominated_plots,by="Site")%>%filter(!is.na(plotIDM))%>%head(45)->temp_data

biomes=c("Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands","Tropical & Subtropical Moist Broadleaf Forests")

biome_group=list()
for (i in 1:4){
  human_dominated_plots%>%left_join(temp_data%>%dplyr::select(plotIDM,LABEL),by="plotIDM")%>%mutate(code=1:45)%>%
    group_by(LABEL)%>%
    filter(LABEL==biomes[i])%>%pull(code)->biome_group[[i]]
}

# the plots belonging to the same biome
plot_number=numeric()
for(i in 1:4)
{
  plot_number[i]=length(biome_group[[i]])
}


# this is computed on great lakes


richness_compare_crop_nature_guild=readRDS("richness_compare_crop_nature_guild.rds")

human_dominated_plots%>%mutate(plotid=1:45)%>%dplyr::select(Site,plotid,plotIDM)->temp


# for the whole-community data

do.call(rbind,richness_compare_crop_nature_guild[[1]] )%>%data.frame()%>%
  mutate(plotid=unlist(biome_group))%>%mutate(biomes=rep(biome_select,times=plot_number))%>%
  left_join(temp,by="plotid")%>%rename_all(~paste0(c("low","mean_nature","up","mean_crop","plotid","biome","site","plotIDM")))->data_mean_richness_biome


# for different guilds
# the estimated diversity among the natural and human-dominated land systems
data_mean_richness_guild_consider_nature_history=list()
for(m in 1:9){
  
  data_mean_richness_guild_consider_nature_history[[m]]=do.call(rbind,richness_compare_crop_nature_guild[[m]] )%>%data.frame()%>%
    mutate(plotid=unlist(biome_group))%>%mutate(biomes=rep(biome_select,times=plot_number))%>%
    left_join(temp,by="plotid")%>%rename_all(~paste0(c("low","mean_nature","up","mean_crop","plotid","biome","site","plotIDM")))
}


saveRDS(data_mean_richness_guild_consider_nature_history,file="data_mean_richness_guild_consider_nature_history.rds")

# convert the data into integer

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



my_list_integer <- lapply(data_mean_richness_guild_consider_nature_history, convert_to_integer)

data_mean_richness_guild_consider_nature_history=my_list_integer

#conver the na to 0, no species found 
data_mean_richness_guild_consider_nature_history<- lapply(data_mean_richness_guild_consider_nature_history, function(df) {
  df[is.na(df)] <- 0
  return(df)
})

# reorganize the data
# so that for each guild, we have mean gamma diversity among natural and modified plots for all the 45 plots

data_mean_richness_guild_consider_nature_history_model=list()
for (m in 1:9)
{
  data_mean_richness_guild_consider_nature_history_model[[m]]=data_mean_richness_guild_consider_nature_history[[m]]%>%dplyr::select(mean_nature,biome, site,plotIDM)%>%
    rename(mean_rich=mean_nature)%>%
    rbind(data_mean_richness_guild_consider_nature_history[[m]]%>%dplyr::select(mean_crop,biome, site,plotIDM)%>%rename(mean_rich=mean_crop))%>%
    mutate(type=rep(c("nature","crop"),each=45))
}


saveRDS(data_mean_richness_guild_consider_nature_history_model,file="data_mean_richness_guild_consider_nature_history_model.rds")

data_mean_richness_guild_consider_nature_history_model=readRDS(file="data_mean_richness_guild_consider_nature_history_model.rds")
# test the significance of the results based on generalized linear mixed-effect models
# testing if gamma diversity is affected both by land cover type and biome
# all are significant except for the whole fungal community


result=list()
for (m in 1:9)
{
  mod=glmer(mean_rich~type+biome+(1 | site),family = poisson(link = "log"),data=data_mean_richness_guild_consider_nature_history_model[[m]])
  
  result[[m]]=summary(mod) 
}

#get the overall effect across biomes
biomes_effect=list()
for (i in 1:9)
{
  coef(result[[i]])%>%data.frame()%>%slice(2)->biomes_effect[[i]]
  
}

do.call(rbind,biomes_effect)%>%bind_cols(guild)%>%round(digits = 2)

#we need to test the difference for each biome
#note changes in the random effect term in the fourth biome
df_significance=list()
for (m in 1:9)
{
  result=numeric()
  for (j in 1:4)
  {
    if(j<4)
    {
      mod=glmer(mean_rich~type+(1 | site),family = poisson(link = "log"),
                data=data_mean_richness_guild_consider_nature_history_model[[m]]%>%filter(biome==biome_select[j]))
    }
    else
    {
      mod=glm(mean_rich~type,family = poisson(link = "log"),
              data=data_mean_richness_guild_consider_nature_history_model[[m]]%>%filter(biome==biome_select[j]))
      
    }
    
    result[[j]]=coef(summary(mod))[2,4]
  }
  df_significance[[m]]=result
  
}

# to bind all the results

guild=guild_select

do.call(rbind,df_significance)%>%data.frame()%>%
  bind_cols(guild)%>%rename_all(~paste0(c(biome_select,"guild")))%>%
  mutate(across(where(is.numeric), ~ round(., digits = 2)))->df_significance_temp


#change the values into significance levels

df_significance_temp %>%
  mutate(across(everything(), ~ case_when(
    is.numeric(.) | !is.na(suppressWarnings(as.numeric(.))) & !is.character(.) ~ case_when(
      as.numeric(.) < 0.001 ~ "***",
      as.numeric(.) < 0.01  ~ "**",
      as.numeric(.) <= 0.1  ~ "*",
      as.numeric(.) > 0.1  ~ "ns",
      TRUE~ as.character(.)
    ),
    TRUE ~ as.character(.) # Keep non-numeric values unchanged
  )))->df_significance_plot


### get the gamma diversity ratio for different guilds
# for each biomes, we get the mean


# if we get the site-level response ratio

# if we get the individual ratio and then get the mean


richness_ratio_rarefy_guild_site_consider_nature_history=list()
for(m in 1:9)
{
  data_mean_richness_guild_consider_nature_history_model[[m]]->temp
  
  site_level_ratio=list()
  for (i in 1:4)
  {
    temp%>%filter(biome==biome_select[i])%>%dplyr::select(mean_rich,type,plotIDM)->d
    
    site_level_ratio[[i]]= reshape(d, idvar = "plotIDM", timevar = "type", direction = "wide")%>%
      mutate(ratio=mean_rich.crop/mean_rich.nature)%>%pull(ratio)
    
  }
  
  richness_ratio_rarefy_guild_site_consider_nature_history[[m]]=unlist(site_level_ratio)
}


guild=c("all","AM","EM","plapat","soilsap","littersap","woodsap","epiphy","para")

do.call(cbind,richness_ratio_rarefy_guild_site_consider_nature_history)%>%data.frame()%>%rename_all(~paste0(guild))%>%
  mutate(biome=rep(biome_select,times=c(plot_number)))%>%melt()%>%group_by(biome,variable)%>%
  summarise(mean_ratio=mean(value),sd_ratio=sd(value))%>%mutate(guild_biome=paste(variable,"_",biome))%>%
  dplyr::rename(guild=variable)->biome_site_level_richness_ratio_consider_nature_history

saveRDS(biome_site_level_richness_ratio_consider_nature_history,file="biome_site_level_richness_ratio_consider_nature_history.rds")

# the quantified response ratio will be used to determine species relative habitat affinity

## create figures for the updated data
#

biome_site_level_richness_ratio_consider_nature_history=readRDS(file="biome_site_level_richness_ratio_consider_nature_history.rds")
# for the updated, we exclude some of the guilds

guild=c("all","AM","EM","plapat","soilsap")

df_significance_plot=df_significance_plot%>%filter(guild%in%c("all","AM" ,"EM","plapat","soilsap"))


biome_site_level_richness_ratio_consider_nature_history_sub=biome_site_level_richness_ratio_consider_nature_history%>%
  filter(guild%in%c("all","AM" ,"EM","plapat","soilsap"))

