##get species change rate under two scenarios
##both climate and land cover were considered with four data

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")
species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")


data=c("species_change_land_rcp245","species_change_land_rcp585","species_change_climate_rcp245","species_change_climate_rcp585")

scenario=c("SSP2-4.5","SSP5-8.5","SSP2-4.5","SSP5-8.5")

#(1) get the proportions of pixels that show either species gain or loss

change_rate_guild_scenarios=list()# the j corresponds to four scenarios for both land cover and climate effect
for (j in 1:4)
  {
change_rate_guild=list()
for (m in 1:9)
{
  change_rate=matrix(ncol=3,nrow=5)
  
  for (i in 1:5)
  {
    if(i==5)#all biomes were included
    {
      get(data[j])%>%filter(variable==guild_type[m])%>%
        bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
        filter(LABEL%in%biomes_four&!is.na(value))->df1
      
      df1%>%filter(value>0)%>%dim()%>%head(1)->gain
      df1%>%filter(value==0)%>%dim()%>%head(1)->no_change
      df1%>%filter(value<0)%>%dim()%>%head(1)->loss
    }
    else{# select just one biome
      get(data[j])%>%filter(variable==guild_type[m])%>%
        bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
        filter(LABEL== biomes_four[i]&!is.na(value))->df1
      
      df1%>%filter(value>0)%>%dim()%>%head(1)->gain
      df1%>%filter(value==0)%>%dim()%>%head(1)->no_change
      df1%>%filter(value<0)%>%dim()%>%head(1)->loss
    }
    
    change_rate[i,1]=gain
    change_rate[i,2]=no_change
    change_rate[i,3]=loss
  }
  change_rate_guild[[m]]=change_rate
}
change_rate_guild_scenarios[[j]]=change_rate_guild

}

change_rate_guild_scenarios[[1]]#means land cover effect in the low-emission scenario for all guilds

saveRDS(change_rate_guild_scenarios,file="change_rate_guild_scenarios.rds")

#for each guild, to get the data to create the pie plot
pie_data_guild=list()
for (m in 1:9){
  pie_data=list()
  for (i in 1:4)# for both climate and land cover
  {
    change_rate_guild_scenarios[[i]][[m]]%>%
      data.frame()%>%rename_all(~paste0(c("gain","no","loss")))%>%
      mutate(biome=c("Temperate","Conifer","Grassland","Dry","All"))%>%
      melt()%>%group_by(biome)%>%mutate(percent = value / sum(value))->temp_pie_data
    temp_pie_data$biome=factor( levels=c("All" ,"Temperate", "Conifer" ,  "Grassland", "Dry" ),temp_pie_data$biome)
    
    pie_data[[i]]=temp_pie_data
  }
  pie_data_guild[[m]] =pie_data
}
#pie_data_guild[[1]], means the first fungal guild
#with four list, each correspond to land and cliamte effects under two scenarios

saveRDS(pie_data_guild,file="pie_data_guild.rds")

#to get the overall effect of both factors
# it is based on the species change rate

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

biomes_four=c( "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests" ,                    
               "Temperate Grasslands, Savannas & Shrublands" ,"Tropical & Subtropical Dry Broadleaf Forests")


change_rate_mean_guild_scenario=list()
for (j in 1:4)
  {
change_rate_mean_guild=list()
for (m in 1:9)
{
change_rate=list()
for (i in 1:5)
{
  if(i ==5)
  {
    get(data[j])%>%filter(variable==guild_type[m])%>%
      bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
      filter(LABEL%in%biomes_four&!is.na(value))->df1# na cells were filtered
    
    df1%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
      summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->data_temp
    
    df1%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
      overall_change
    data_temp%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")%>%
      data.frame()->change_rate[[i]]
  }
  else
  {
    get(data[j])%>%filter(variable==guild_type[m])%>%
      bind_cols(grid_level_biomes%>%dplyr::select(LABEL))%>%
      filter(LABEL== biomes_four[i]&!is.na(value))->df1
    
    df1%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
      summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->data_temp
    
    df1%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
      overall_change
    data_temp%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")%>%
      data.frame()->change_rate[[i]]
    
  }
  
}
change_rate_mean_guild[[m]]=change_rate#just for one global change scenario
}
change_rate_mean_guild_scenario[[j]]=change_rate_mean_guild
}

change_rate_mean_guild_scenario[[1]]#means for land cover effects in the low-emission scenarios

#to combine the list into a dataframe to create the plots
#for the two scenarios to get the replicates for the mean change
# we only select the 9 th guild when all the guilds were combined
data_biome_change_rate_guild=list()
for (m in 1:9)
  {
  data_biome_change_rate=list()
  for (i in 1:4)#i corresponds to the four data 2 scenarios x 2 effects(land and climate)
  {
    dimensions_matrix <- sapply(change_rate_mean_guild_scenario[[i]][[m]], dim)[1,]
    # the number of rows selected 
    do.call(rbind,change_rate_mean_guild_scenario[[i]][[m]])%>%
      mutate(biome=rep(c(biomes_four,"all"),time=dimensions_matrix))->data_mean_change
    data_mean_change$biome=factor(levels = rev(c("all",biomes_four)),data_mean_change$biome)
    data_biome_change_rate[[i]]=data_mean_change
  }
  data_biome_change_rate_guild[[m]]=data_biome_change_rate
}

saveRDS(data_biome_change_rate_guild,file="data_biome_change_rate_guild.rds")
#each list has 9 lists with found data sets


pie_data_guild=readRDS("pie_data_guild.rds")



# load in the data
data_biome_change_rate_guild=readRDS("data_biome_change_rate_guild.rds")


full_name=c("AM","EM","Soil saprotroph","Litter saprotroph","Wood saprotroph", "Plant pathogen", "Parasite","Epiphyte","All")

#set different
xlim_list=list(c(-30,30),c(-30,30),c(-60,60),c(-60,60))

xlim_list_am=list(c(-100,100),c(-100,100),c(-60,60),c(-60,60))

#set different space for the text
x_space=c(0.27,0.27,5.7,5.7)


###for the overal fungal diversity


# to combined the range and the magnitude of both effects
# each i corresponds to land and climate effect with two climate change scenarios



# for the individual guilds

pp_combine_effect_guild[[m]]

plot_grid(pp_combine_effect_guild[[1]][[3]],
          pp_combine_effect_guild[[2]][[3]],
          pp_combine_effect_guild[[3]][[3]],
          pp_combine_effect_guild[[6]][[3]],
          
          pp_combine_effect_guild[[1]][[1]],
          pp_combine_effect_guild[[2]][[1]],
          pp_combine_effect_guild[[3]][[1]],
          pp_combine_effect_guild[[6]][[1]],
          
          pp_combine_effect_guild[[1]][[4]],
          pp_combine_effect_guild[[2]][[4]],
          pp_combine_effect_guild[[3]][[4]],
          pp_combine_effect_guild[[6]][[4]],
          
          pp_combine_effect_guild[[1]][[2]],
          pp_combine_effect_guild[[2]][[2]],
          pp_combine_effect_guild[[3]][[2]],
          pp_combine_effect_guild[[6]][[2]],
          
          ncol=4,label_x = -0.1,label_y = 1.08,
          label_size = 10,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:16])


#for 




pp_mean_effect_guild[[1]][[1]]$widths=pp_mean_effect_guild[[2]][[1]]$widths
pp_mean_effect_guild[[1]][[2]]$widths=pp_mean_effect_guild[[2]][[2]]$widths

#to align the panels
pp_mean_effect_guild[[1]][[1]]$widths=pp_mean_effect_guild[[1]][[3]]$widths

pp_mean_effect_guild[[2]][[1]]$widths=pp_mean_effect_guild[[2]][[3]]$widths
pp_mean_effect_guild[[2]][[2]]$widths=pp_mean_effect_guild[[2]][[3]]$widths

pp_mean_effect_guild[[3]][[1]]$widths=pp_mean_effect_guild[[3]][[3]]$widths
pp_mean_effect_guild[[4]][[1]]$widths=pp_mean_effect_guild[[4]][[3]]$widths

pp_mean_effect_guild[[2]][[4]]$widths=pp_mean_effect_guild[[2]][[3]]$widths

pp_mean_effect_guild[[3]][[2]]$widths=pp_mean_effect_guild[[3]][[1]]$widths
pp_mean_effect_guild[[3]][[4]]$widths=pp_mean_effect_guild[[3]][[1]]$widths

pp_mean_effect_guild[[5]][[3]]$widths=pp_mean_effect_guild[[5]][[2]]$widths
pp_mean_effect_guild[[5]][[1]]$widths=pp_mean_effect_guild[[5]][[2]]$widths
pp_mean_effect_guild[[5]][[4]]$widths=pp_mean_effect_guild[[5]][[2]]$widths

pp_mean_effect_guild[[6]][[2]]$widths=pp_mean_effect_guild[[6]][[1]]$widths
pp_mean_effect_guild[[7]][[2]]$widths=pp_mean_effect_guild[[7]][[1]]$widths

pp_mean_effect_guild[[2]][[2]]$widths=pp_mean_effect_guild[[2]][[1]]$widths

plot_grid(pp_mean_effect_guild[[1]][[1]],
          pp_mean_effect_guild[[2]][[1]],
          pp_mean_effect_guild[[3]][[1]],
          pp_mean_effect_guild[[4]][[1]],
          pp_mean_effect_guild[[5]][[1]],
          pp_mean_effect_guild[[6]][[1]],
          pp_mean_effect_guild[[7]][[1]],
          pp_mean_effect_guild[[8]][[1]],
          
          pp_mean_effect_guild[[1]][[3]],
          pp_mean_effect_guild[[2]][[3]],
          pp_mean_effect_guild[[3]][[3]],
          pp_mean_effect_guild[[4]][[3]],
          pp_mean_effect_guild[[5]][[3]],
          pp_mean_effect_guild[[6]][[3]],
          pp_mean_effect_guild[[7]][[3]],
          pp_mean_effect_guild[[8]][[3]],
          
          pp_mean_effect_guild[[1]][[2]],
          pp_mean_effect_guild[[2]][[2]],
          pp_mean_effect_guild[[3]][[2]],
          pp_mean_effect_guild[[4]][[2]],
          pp_mean_effect_guild[[5]][[2]],
          pp_mean_effect_guild[[6]][[2]],
          pp_mean_effect_guild[[7]][[2]],
          pp_mean_effect_guild[[8]][[2]],
          
          
          
          pp_mean_effect_guild[[1]][[4]],
          pp_mean_effect_guild[[2]][[4]],
          pp_mean_effect_guild[[3]][[4]],
          pp_mean_effect_guild[[4]][[4]],
          pp_mean_effect_guild[[5]][[4]],
          pp_mean_effect_guild[[6]][[4]],
          pp_mean_effect_guild[[7]][[4]],
          pp_mean_effect_guild[[8]][[4]],
          ncol=8,label_x = -0.1,label_y = 1.08,
          label_size = 10,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:32])

##12x16



plot_grid(pp_pie_guild[[1]][[1]],
          pp_pie_guild[[2]][[1]],
          pp_pie_guild[[3]][[1]],
          pp_pie_guild[[4]][[1]],
          pp_pie_guild[[5]][[1]],
          pp_pie_guild[[6]][[1]],
          pp_pie_guild[[7]][[1]],
          pp_pie_guild[[8]][[1]],
          
          pp_pie_guild[[1]][[3]],
          pp_pie_guild[[2]][[3]],
          pp_pie_guild[[3]][[3]],
          pp_pie_guild[[4]][[3]],
          pp_pie_guild[[5]][[3]],
          pp_pie_guild[[6]][[3]],
          pp_pie_guild[[7]][[3]],
          pp_pie_guild[[8]][[3]],
          
          pp_pie_guild[[1]][[2]],
          pp_pie_guild[[2]][[2]],
          pp_pie_guild[[3]][[2]],
          pp_pie_guild[[4]][[2]],
          pp_pie_guild[[5]][[2]],
          pp_pie_guild[[6]][[2]],
          pp_pie_guild[[7]][[2]],
          pp_pie_guild[[8]][[2]],
          
          pp_pie_guild[[1]][[4]],
          pp_pie_guild[[2]][[4]],
          pp_pie_guild[[3]][[4]],
          pp_pie_guild[[4]][[4]],
          pp_pie_guild[[5]][[4]],
          pp_pie_guild[[6]][[4]],
          pp_pie_guild[[7]][[4]],
          pp_pie_guild[[8]][[4]],
          ncol=8,label_x = -0.01251,label_y = 1.07,
          label_size = 10,
          labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:32]
      )

