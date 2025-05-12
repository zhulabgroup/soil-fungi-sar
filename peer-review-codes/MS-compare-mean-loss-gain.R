#to compare the mean diversity change rate among different guilds
# for the species loss rate and gain rate

setwd("/Volumes/seas-zhukai/proj-soil-fungi/land-use-climate-historical")

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")

species_change_climate_rcp245=species_change_climate_rcp245[,-c(1,2)]

species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")

species_change_climate_rcp585=species_change_climate_rcp585[,-c(1,2)]

species_change_land_rcp245=readRDS("species_change_land_rcp245.rds")
species_change_land_rcp585=readRDS("species_change_land_rcp585.rds")

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

biomes_four=c( "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests" ,                    
               "Temperate Grasslands, Savannas & Shrublands" ,"Tropical & Subtropical Dry Broadleaf Forests")


data=c("species_change_land_rcp245","species_change_land_rcp585","species_change_climate_rcp245","species_change_climate_rcp585")


# for the climate change impact, 41394 cells were included

observations=c(41725,41725,41394,41394)


guild_type=c("AM"  , "EM"  , "soilsap" ,  "littersap", "woodsap",   "plapat" ,   "para" ,     "epiphy" ,   "all" )

com_data_climate_guild=list()
for (m in 1:9){
  com_data_climate=list()
  
  for (j in 1:4){
    
    get(data[j])%>%
      bind_cols(rep(grid_level_biomes$LABEL,9))%>%
      rename_all(~paste0(c("guild","value","LABEL")))%>%
      filter(LABEL%in%biomes_four&!is.na(value))%>%
      mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered
    
    all_biome=dim(df1)[1]
    
    # for all fungal guilds
    
    df1%>%mutate(LABEL_all=rep("all",all_biome))%>%
      dplyr::select(guild,value,LABEL_all,plotid)%>%
      rename(LABEL=LABEL_all)%>%bind_rows(df1)->df1
    
    # for the overall effect
    df1%>%filter(guild==guild_type[m]&LABEL!="all")->df_overall
    
    # to compare the net change 
    
    #to give a plot id for each cell
    
    # to test difference for species gains for the overall fungal diversity
    
    df1%>%filter(value>0&guild%in%c(guild_type[m]))->df_gain
    
    #need to transform the resposne variable
    
    df_gain$log_response <- log(df_gain$value)
    
    
    model_all <- gls(log_response ~ LABEL, data = df_gain)# performe the model and get the mean for all groups
    
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="all")%>%mutate(.group="")%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))->temp_mean
    
    # performe the model just for four biomes
    
    model_four <- gls(log_response ~ LABEL, data = df_gain%>%filter(LABEL!="all"))
    
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    # to bind the letters with the mean and add a column specifying the different letters
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))%>%
      bind_rows(temp_mean)->data_gain # the estimated mean for all the groups
    # the data set for all the groups
    # the data used for species gains
    
    # for species loss rate
    
    df1%>%filter(value<0&guild%in%c(guild_type[m]))->df_loss
    
    df_loss$sqrt_response <- sqrt(abs(df_loss$value))
    
    model_all <- gls(sqrt_response  ~ LABEL, data = df_loss)# perform the model and get the mean for all groups
    
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    
    
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="all")%>%mutate(.group="")%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)->temp_mean
    
    model_four <- gls(sqrt_response ~ LABEL, data = df_loss%>%filter(LABEL!="all"))
    
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    #back transforming the rate
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)%>%
      bind_rows(temp_mean) ->data_loss
    #combined data
    # need to consider the overall effect
    
    
    #model_overal <- gls(sqrt_response ~ LABEL, data = df1%>%filter(LABEL!="all"))
    
    
    com_data=rbind(data_gain,data_loss)
    
    com_data%>%mutate(change=ifelse(origin_mean > 0, "Positive", "Negative"))->com_data
    
    com_data$LABEL=factor(com_data$LABEL,
                          levels=c("all",
                                   "Temperate Broadleaf & Mixed Forests",
                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                   "Tropical & Subtropical Dry Broadleaf Forests")%>%rev())
    
    
    com_data_climate[[j]]=com_data
  }
  com_data_climate_guild[[m]]=com_data_climate
}


## to combined the net effect based on the raw data
data_biome_change_rate_guild=readRDS("data_biome_change_rate_guild.rds")

biome_mean_guild=list()
for (m in 1:9){
  biome_mean=list()
  for (i in 1:4)
  {
    com_data_climate_guild[[m]][[i]]%>%
      rename(biome=LABEL)->df1
    
    data_biome_change_rate_guild[[m]][[i]]%>%
      select(biome,overal_mean)->df2
    
    df1%>%left_join(df2,by="biome")%>%distinct()->biome_mean[[i]]
    
  }
  biome_mean_guild[[m]]=biome_mean
}

com_data_climate_guild=biome_mean_guild#data used to compared the mean among biomes

saveRDS(com_data_climate_guild,file="com_data_climate_guild.rds")





