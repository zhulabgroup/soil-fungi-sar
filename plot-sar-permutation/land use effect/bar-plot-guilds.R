
#compare the mean values among guilds

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

observations=c(41725,41725,41394,41394)



data_compare_LETTER=list()

for(j in 1:4)
{
  get(data[j])%>%
    bind_cols(rep(grid_level_biomes$LABEL,9))%>%
    rename_all(~paste0(c("guild","value","LABEL")))%>%
    filter(LABEL%in%biomes_four&!is.na(value))%>%
    mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered
  
  #when all the biomes were included
  
  all_biome=dim(df1)[1]
  
  # add a column to show all biomes
  
  df1%>%mutate(LABEL_all=rep("All",all_biome))%>%
    dplyr::select(guild,value,LABEL_all,plotid)%>%
    rename(LABEL=LABEL_all)%>%bind_rows(df1)->df1
  
  # select the guild of your interest
  # for species gain for across all the biomes
  # compare the mean among guilds
  # for species gains
  df1%>%filter(value>0&guild%in%c("EM","soilsap","AM","plapat")&LABEL=="All")->df2
  
  # need to transform the data 
  
  df2$log_response=log(df2$value)
  
  mod=lmer(log_response~guild+(1 | plotid),data=df2)# compare the means among guilds across all biomes
  
  anova(mod, type = "III")
  
  compare_mean=emmeans(mod, pairwise ~ guild ,pbkrtest.limit = 93160)#the last parameter could change
  
  # get the values associated with the letters
  # the letetrs among guilds when all the guilds were combined
  
  multcomp::cld(object = compare_mean$emmeans,
                Letters = letters)->compare_data_gain_all
  
  # need to get the overall mean and when all the guilds were combined
  
  #compare the mean values among guilds for each biome
  
  df1%>%filter(value>0&guild%in%c("EM","soilsap","AM","plapat")&LABEL!="All")->df2
  
  df2$log_response=log(df2$value)
  
  mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df2)
  
  compare_mean=emmeans(mod, pairwise ~ guild| LABEL,pbkrtest.limit = 93160)#the last parameter could change
  
  #add letters for each biome
  # this can be add to the horizontal dimension of the figure
  
  multcomp::cld(object = compare_mean$emmeans,
                Letters = letters)%>%data.frame()%>%
    mutate(ori_mean =exp(emmean),ori_lower =  exp(lower.CL), ori_upper = exp(upper.CL))->compare_gain_guild_biome
  
  
  compare_gain_guild_biome%>%mutate(type=rep("Positive",length.out=nrow(compare_gain_guild_biome)))->compare_gain_guild_biome
  
  ## for species loss rate
  # when all the biomes were considered
  
  df1%>%filter(value<0&guild%in%c("EM","soilsap","AM","plapat")&LABEL=="All")->df3
  
  # need to transform the data 
  
  df3$sqrt_response=sqrt(abs(df3$value))
  
  mod=lmer(sqrt_response~guild+(1 | plotid),data=df3)# across biomes
  
  anova(mod, type = "III")
  
  compare_mean=emmeans(mod, pairwise ~ guild ,pbkrtest.limit = 125497)#the last parameter could change
  
  # get the values associated with the letters
  # the letetrs among guilds when all the guilds were combined
  
  multcomp::cld(object = compare_mean$emmeans,
                Letters = letters)->compare_data_loss_all
  # need to get the overall mean and when all the guilds were combined
  #compare the mean among the four biomes
  
  df1%>%filter(value<0&guild%in%c("EM","soilsap","AM","plapat")&LABEL!="All")->df4
  
  df4$sqrt_response=sqrt(abs(df4$value))
  mod=lmer(sqrt_response~guild*LABEL+(1 | plotid),data=df4)
  compare_mean=emmeans(mod, pairwise ~ guild| LABEL,pbkrtest.limit = 125497)#the last parameter could change
  
  #add letters for each biome
  # this can be add to the horizontal dimension of the figure
  
  multcomp::cld(object = compare_mean$emmeans,
                Letters = letters)%>%data.frame()%>%
    mutate(ori_mean =-1*emmean^2,ori_lower =  -1*lower.CL^2, ori_upper = -1*upper.CL^2)->compare_loss_guild_biome
  
  compare_loss_guild_biome%>%mutate(type=rep("Negative",length.out=nrow(compare_loss_guild_biome)))->compare_loss_guild_biome
  
  #for the overall mean among guilds
  
  compare_data_gain_all%>%data.frame()%>%bind_rows(
    compare_data_loss_all%>%data.frame())%>%
    mutate(type=rep(c("Positive","Negative"),each=4))%>%
    mutate(biome=rep("All",8))%>%
    rename(LETTER=.group,LABEL=biome)->data_compare_guild
  data_compare_guild$LETTER <- str_replace_all(data_compare_guild$LETTER, " ", "")#to replace the space
  
  #based on the estimated values, we back transformed the values when all the biomes were included
  
  data_compare_guild%>%mutate(ori_mean = ifelse(type =="Positive", exp(emmean), -1*emmean^2),
                              ori_lower = ifelse(type =="Positive", exp(lower.CL), -1*lower.CL^2),
                              ori_upper = ifelse(type =="Positive", exp(upper.CL), -1*upper.CL^2))%>%
    dplyr::select(guild, LABEL, emmean,   SE , df ,lower.CL,upper.CL, LETTER , type, LABEL,ori_mean,ori_lower,ori_upper)->data_compare_guild
  # to bind the biome-specific values and when all the biomes were considered
  compare_gain_guild_biome%>%bind_rows(compare_loss_guild_biome)%>%
    rename(LETTER=.group)%>%bind_rows(data_compare_guild)%>%
    mutate(LETTER = str_replace_all(LETTER, " ", ""))%>%
    mutate(new=paste0(LABEL," * ",type))->data_among_guilds_biomes
  
  data_compare_LETTER[[j]]=data_among_guilds_biomes
  
}

saveRDS(data_compare_LETTER,file="data_compare_LETTER.RDS")






impact=c("Land use impact","Land use impact","Climate impact","Climate impact")

scenario=c("SSP2-4.5","SSP5-8.5","SSP2-4.5","SSP5-8.5")


pp=list()
for (j in 1:4)
  {
  data_compare_LETTER[[j]]$LABEL=factor(data_compare_LETTER[[j]]$LABEL,
                                        levels=rev(c("All",
                                                     "Temperate Broadleaf & Mixed Forests",
                                                     "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                                     "Tropical & Subtropical Dry Broadleaf Forests")))
pp[[j]]=ggplot(data=data_compare_LETTER[[j]],aes(fill=guild,y=LABEL ,x=ori_mean*100,group=guild))+
          geom_col(width = 0.8,position = "dodge",color="black",size=0.1)+
      geom_errorbar(aes(xmin = ori_lower*100, xmax = ori_upper*100), 
                    position = position_dodge(width = 0.8), 
                    width = 0.15,size=0.3) +
      geom_text(aes(label=LETTER),position = position_dodge(width = 0.8),
                hjust = ifelse(data_compare_LETTER[[j]]$ori_mean> 0, -7.5, 10),size=3)+
      geom_text(aes(label= c(paste0("(",sprintf("%.2f",ori_mean*100,")"),")"))),
                position = position_dodge(width = 0.8),
                hjust = ifelse(data_compare_LETTER[[j]]$ori_mean > 0, -.25, 1.25),size=3)+
      xlim(-30,30)+
      scale_y_discrete(breaks=unique(data_compare_LETTER[[j]]$LABEL),position="right",
                       labels=paste0(c("",""," ","","")))+
      xlab("Diversity change rate (%)")+
      ylab("")+
      ggtitle(scenario[j])+
    theme(legend.position = c(0.8,0.25),
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10), 
          axis.title.y = element_text(size = 12), 
          axis.title.x = element_text(size = 12), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0, 0, 0.1, 0.1), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
    geom_vline(xintercept =0,color="gray",linetype="dashed")+
    scale_fill_manual("",breaks=c("AM","EM","plapat","soilsap"),
                      labels=c("AM","EM","Plant pathogens","Soil saprotrophs"),
                      values=c("#c94e65","#037f77","royalblue","chocolate1"))
    
} 


#to arrange the four plots           
  plot_grid(pp[[3]],pp[[1]],pp[[4]],pp[[2]],
            ncol=2,
            label_x = 0,label_y = 1.03,
            labels = paste0("(", c(letters, outer(letters, letters, paste0)), ")") [1:4])           
      
     
    
       
      
     
     