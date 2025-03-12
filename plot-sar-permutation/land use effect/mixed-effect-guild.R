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

# if we consider species loss

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

 ##################################################################################
######################codes below this line shoud be disregarded##################################


# the below codes could be disregarded
# for individual guilds, we focused on testing difference in species loss and gains within each biome 

#based on the difference to label letters

# test the difference among biomes with one-way anova

## create the plots

# need to bind the data with the net effect data
# for each guild data, the bind data shows difference in the mean species loss and gain among biomes but not guilds

com_data_guild_letter=list()

for(m in 1:9){
  
  com_temp=list()
  for(i in 1:4)
  {
    
    com_data_climate_guild[[m]][[i]]%>%mutate(new=paste(LABEL,"*",change))%>%
      dplyr::select(new,.group)->df1
    
    data_biome_change_rate_guild[[m]][[i]]%>%mutate(new=paste(biome,"*",type))%>%left_join(df1,by="new")%>%
      rename(letter=.group)%>%select(-new)%>%mutate(letter = str_replace_all(letter, " ", ""))->com_temp[[i]]
    
  }
  com_data_guild_letter[[m]]=com_temp
}

data_biome_change_rate_guild=com_data_guild_letter


#compare the mean among guilds within guilds
# bind the data with the other data


data_temp_final=list()
for(j in 1:4)
{
  data_temp=list()
  for(i in 1:9){
    
    
    data_biome_change_rate_guild[[i]][[j]]%>%
      mutate(new=paste(biome,"*",type))->df1
    
    data_compare_LETTER[[j]]%>%filter(guild==guild_type[i])%>%
      mutate(new=paste(LABEL,"*",type))%>%
      dplyr::select(new,LETTER)->df2
    
    df1%>%left_join(df2,by="new")%>%
      mutate(note=paste0(letter,"",LETTER))->data_final
    data_temp[[i]]=data_final
  }
  data_temp_final[[j]]=data_temp
}
## for each data, only four guils have data,

data_temp_final[[1]][[1]]
data_temp_final[[1]][[2]]
data_temp_final[[1]][[3]]
data_temp_final[[1]][[6]]

## create the plot that compare the means among guilds

# use the estimated data to create the plots


##### need to add the data that compared the means within guilds

observations=c(41725,41725,41394,41394)
####################################################################################
model_mean_guild=list()
for(m in 1:9)
  
{
  model_mean=list()
  for (j in 1:4){
    get(data[j])%>%
      bind_cols(rep(grid_level_biomes$LABEL,9))%>%
      rename_all(~paste0(c("guild","value","LABEL")))%>%
      filter(LABEL%in%biomes_four&!is.na(value))%>%
      mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered
    
    all_biome=dim(df1)[1]
    
    # for all fungal guilds
    
    df1%>%mutate(LABEL_all=rep("All",all_biome))%>%
      dplyr::select(guild,value,LABEL_all,plotid)%>%
      rename(LABEL=LABEL_all)%>%bind_rows(df1)->df1
    
    
    #to give a plot id for each cell
    
    # to test difference for species gains for the overall fungal diversity
    
    df1%>%filter(value>0&guild%in%c(guild_type[m]))->df_gain
    
    #need to transform the response variable
    
    df_gain$log_response <- log(df_gain$value)
    
    
    model_all <- gls(log_response ~ LABEL, data = df_gain)# performe the model and get the mean for all groups
    
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    
    
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="All")%>%mutate(.group="")%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))->temp_mean
    
    # performe the model just for four biomes
    
    model_four <- gls(log_response ~ LABEL, data = df_gain%>%filter(LABEL!="All"))
    
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    # to bind the letters with the mean and add a column specifying the different letters
    
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=exp(emmean),ori_low=exp(lower.CL),ori_up=exp(upper.CL))%>%
      bind_rows(temp_mean)->data_gain # the estimated mean for all the groups
    # the data set for all the groups
    # the data used for species gains
    
    # for species loss
    
    df1%>%filter(value<0&guild%in%c(guild_type[m]))->df_loss
    
    df_loss$sqrt_response <- sqrt(abs(df_loss$value))
    
    model_all <- gls(sqrt_response  ~ LABEL, data = df_loss)# performe the model and get the mean for all groups
    
    pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
    
    
    pairwise_all%>%data.frame()%>%
      filter(LABEL=="All")%>%mutate(.group="")%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)->temp_mean
    
    model_four <- gls(sqrt_response ~ LABEL, data = df_loss%>%filter(LABEL!="All"))
    pairwise_four <- emmeans(model_four, pairwise ~ LABEL)
    
    #back transforming the rate
    
    multcomp::cld(object = pairwise_four$emmeans,
                  Letters = letters)%>%data.frame()%>%
      mutate(origin_mean=-1*emmean^2,ori_low=-1*lower.CL^2,ori_up=-1*upper.CL^2)%>%
      bind_rows(temp_mean) ->data_loss
    #combined data
    
    com_data=rbind(data_gain,data_loss)
    
    
    com_data%>%mutate(change=ifelse(origin_mean > 0, "Positive", "Negative"))%>%
      mutate(new=paste(LABEL,"*",change))->com_data
    
    com_data$LABEL=factor(com_data$LABEL,
                          levels=c("All",
                                   "Temperate Broadleaf & Mixed Forests",
                                   "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                   "Tropical & Subtropical Dry Broadleaf Forests"))
    model_mean[[j]]=com_data
  }
  model_mean_guild[[m]]=model_mean
}


# need to bind with the net effect and the among-guild comparisons
# to get the full data



data_model_mean_guild=list()
for(m in c(1,2,3,6))
{
  data_model_mean=list()
  for(j in 1:4){
    if(j>1){
      data_compare_LETTER[[j]]%>%filter(guild==guild_type[m])%>%
        dplyr::select(LETTER,new,ori_mean,ori_lower,ori_upper )->df1
      
      model_mean_guild[[m]][[j]]%>%left_join(df1,by="new")%>%
        mutate(letter=paste(.group,LETTER))->df2
      
      df2$letter= str_replace_all(df2$letter, " ", "")
      
      #need to combined the net effect
      do.call(rbind, change_rate_mean_guild_scenario[[j]][[m]])%>%
        mutate(LABLE=rep(c(biomes_four,"All"),each=2))%>%
        dplyr::select(LABLE,overal_mean)%>%distinct()%>%
        rename(LABEL=LABLE)%>%left_join(df2,by="LABEL")%>%
        rename(type=change)->df3
      
      data_model_mean[[j]]=df3
      
    }
    else{
      data_compare_LETTER[[j]]%>%filter(guild==guild_type[m])%>%
        dplyr::select(LETTER,new,ori_mean,ori_lower,ori_upper )->df1
      
      model_mean_guild[[m]][[j]]%>%left_join(df1,by="new")%>%
        mutate(letter=paste(.group,LETTER))->df2
      
      df2$letter= str_replace_all(df2$letter, " ", "")
      
      #need to combined the net effect
      do.call(rbind, change_rate_mean_guild_scenario[[j]][[m]])%>%
        mutate(LABLE=rep(c(biomes_four,"All"),times=c(2,2,2,1,2)))%>%
        dplyr::select(LABLE,overal_mean)%>%distinct()%>%
        rename(LABEL=LABLE)%>%left_join(df2,by="LABEL")%>%
        rename(type=change)->df3
      
      
      df3$LABEL=factor(df3$LABEL,
                       levels=rev(c("All",
                                    "Temperate Broadleaf & Mixed Forests",
                                    "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                    "Tropical & Subtropical Dry Broadleaf Forests")))
      
      data_model_mean[[j]]=df3
    }
    
  }
  data_model_mean_guild[[m]]=data_model_mean
}

  
  
  

  