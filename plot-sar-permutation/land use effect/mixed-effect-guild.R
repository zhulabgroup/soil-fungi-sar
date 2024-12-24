#to compare the mean change rate among different guilds
# for the species loss rate and gain rate

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")

species_change_climate_rcp245=species_change_climate_rcp245[,-c(1,2)]

species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")

species_change_climate_rcp585=species_change_climate_rcp585[,-c(1,2)]
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

df1%>%filter(value>0&guild%in%c(guild_type[i]))->df_gain

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



ggplot(data=com_data_climate[[j]],aes(fill=change,y=LABEL ,x=origin_mean))+
  geom_col(width = 0.3)+
  geom_segment(data=com_data_climate[[j]],
               aes(y=LABEL,yend=LABEL,x= ori_low, xend = ori_up),
               arrow = arrow(length = unit(0.1, "cm"),angle=90))+
  
  geom_text(aes(y=LABEL ,x=origin_mean*1.25,label=.group))



  
  
  # for individual guilds, we focused on testing difference in species loss and gains within each biome 
  
  
  
  
  
  #based on the difference to label letters
  
  
  
  # test the difference among biomes with one-way anova
  
oneway.test(log_response~LABEL,data=df1, var.equal = FALSE)

library(DescTools)

dk=oneway.test(value ~ LABEL, data = df1, var.equal = FALSE)

TukeyHSD(dk, method = "hsd")

PostHocTest(aov(value ~ LABEL, data =df1), method = "hsd")

tukey(df1$value,df1$guild,method="G")

df1%>%group_by(LABEL)%>%summarise(m=mean(value))

df1%>%group_by(LABEL)%>%summarise(m=mean(log_response))




summary(model)


df1%>%filter(value>0&!(guild%in%c("all","para","epiphy","littersap","woodsap")))->df1

df1%>%filter(value<0&guild%in%c("all"))->df1

df1$log_response <- log(df1$value)


mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df1)

mod=lmer(log_response~LABEL,data=df1)

mod=lmer(log_response~guild+(1 | plotid),data=df1)# compare the mean among

emmeans(mod, pairwise ~ guild,pbkrtest.limit = 60649)

emm <- emmeans(mod, ~ guild,pbkrtest.limit = 37377)



comp_mean= emmeans(mod, pairwise ~ guild | LABEL,pbkrtest.limit = 60649)


### for species loss

df1%>%filter(value<0&!(guild%in%c("all","para","epiphy","littersap","woodsap")))->df1

df1%>%filter(value<0&!(variable%in%c("all")))->df1

df1$sqrt_response <- sqrt(abs(df1$value))

mod=lmer(sqrt_response~guild*LABEL+(1 | plotid),data=df1)

mod=lmer(sqrt_response~LABEL+(1 | plotid),data=df1)

emm=emmeans(mod, ~ LABEL ,pbkrtest.limit = 250350)
##

mod=lm(sqrt_response~LABEL,data=df1)

emmeans(mod, pairwise ~ LABEL)

## create the plots

emm_df <- as.data.frame(emm)
ggplot(emm_df, aes(x = LABEL, y = emmean)) +
  geom_bar(stat = "identity", position = "dodge", fill = "skyblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(x = "Group", y = "Adjusted Mean Richness") +
  theme_minimal()



## when we look at the net effect

df1%>%filter(!(guild%in%c("all","para","epiphy","woodsap","littersap")))->df1


df1$log_response <- sign(df1$value) * log(abs(df1$value) + 1)


mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df1)

comp_mean= emmeans(mod, pairwise ~ guild | LABEL,pbkrtest.limit = 250350)

# if we do not consider the biomes

mod=lmer(log_response~guild+(1 | plotid),data=df1)




#### for the species gain and loss for individual guilds
## this analysis focuses on the difference in the mean values among different guilds

species_change_climate_rcp245=species_change_climate_rcp245[,-c(1,2)]
species_change_climate_rcp585=species_change_climate_rcp585[,-c(1,2)]

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

df1%>%mutate(LABEL_all=rep("all",all_biome))%>%
  dplyr::select(guild,value,LABEL_all,plotid)%>%
  rename(LABEL=LABEL_all)%>%bind_rows(df1)->df1

# select the guild of your interest
# for species gain for across all the biomes

df1%>%filter(value>0&guild%in%c("EM","soilsap","AM","plapat")&LABEL=="all")->df2

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

df1%>%filter(value>0&guild%in%c("EM","soilsap","AM","plapat")&LABEL!="all")->df2

df2$log_response=log(df2$value)

mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df2)

compare_mean=emmeans(mod, pairwise ~ guild| LABEL,pbkrtest.limit = 93160)#the last parameter could change

#add letters for each biome
# this can be add to the horizontal dimension of the figure

multcomp::cld(object = compare_mean$emmeans,
              Letters = letters)%>%data.frame()%>%
  dplyr::select(guild,LABEL,.group)->compare_gain_guild_biome

compare_gain_guild_biome%>%mutate(type=rep("Positive",length.out=nrow(compare_gain_guild_biome)))->compare_gain_guild_biome

## for species loss

df1%>%filter(value<0&guild%in%c("EM","soilsap","AM","plapat")&LABEL=="all")->df3

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

df1%>%filter(value<0&guild%in%c("EM","soilsap","AM","plapat")&LABEL!="all")->df4

df4$sqrt_response=sqrt(abs(df4$value))
mod=lmer(sqrt_response~guild*LABEL+(1 | plotid),data=df4)
compare_mean=emmeans(mod, pairwise ~ guild| LABEL,pbkrtest.limit = 125497)#the last parameter could change

#add letters for each biome
# this can be add to the horizontal dimension of the figure

multcomp::cld(object = compare_mean$emmeans,
              Letters = letters)%>%data.frame()%>%
  dplyr::select(guild,LABEL,.group)->compare_loss_guild_biome

compare_loss_guild_biome%>%mutate(type=rep("Negative",length.out=nrow(compare_loss_guild_biome)))->compare_loss_guild_biome

#for the overall mean among guilds

compare_data_gain_all%>%data.frame()%>%bind_rows(
compare_data_loss_all%>%data.frame())%>%
  mutate(type=rep(c("Positive","Negative"),each=4))%>%
  mutate(biome=rep("all",8))%>%
  mutate(.group = toupper(.group))%>%
  rename(LETTER=.group,LABEL=biome)->data_compare_guild
data_compare_guild$LETTER <- str_replace_all(data_compare_guild$LETTER, " ", "")

data_compare_guild%>%dplyr::select(guild,LABEL,LETTER,type)->data_compare_guild
compare_gain_guild_biome%>%bind_rows(
  compare_loss_guild_biome)%>%
  mutate(.group = toupper(.group))%>%
  rename(LETTER=.group)->data_compare_guild_biome
data_compare_guild_biome$LETTER <- str_replace_all(data_compare_guild_biome$LETTER, " ", "")

data_compare_guild_biome%>%
  bind_rows(data_compare_guild)->data_compare_guild_biome

data_compare_LETTER[[j]]=data_compare_guild_biome
}


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
guild_type

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
  
  ## for species loss
  
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
    mutate(.group = toupper(.group))%>%
    rename(LETTER=.group,LABEL=biome)->data_compare_guild
  data_compare_guild$LETTER <- str_replace_all(data_compare_guild$LETTER, " ", "")
  #based on the estimated values, we back transforme the values
  
  data_compare_guild%>%mutate(ori_mean = ifelse(type =="Positive", exp(emmean), -1*emmean^2),
                              ori_lower = ifelse(type =="Positive", exp(lower.CL), -1*lower.CL^2),
                              ori_upper = ifelse(type =="Positive", exp(upper.CL), -1*upper.CL^2))%>%
    dplyr::select(guild, LABEL, emmean,   SE , df ,lower.CL,upper.CL, LETTER , type, LABEL,ori_mean,ori_lower,ori_upper)->data_compare_guild
  
  compare_gain_guild_biome%>%bind_rows(compare_loss_guild_biome)%>%
    rename(LETTER=.group)%>%bind_rows(data_compare_guild)%>%mutate(LETTER = toupper(LETTER))%>%
    mutate(LETTER = str_replace_all(LETTER, " ", ""))%>%
    mutate(new=paste0(LABEL," * ",type))->data_among_guilds_biomes
  
  data_compare_LETTER[[j]]=data_among_guilds_biomes
  }

  
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

  
  
  
  ## create the plots
  
  ###
  
  
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
    
    ## for species loss
    
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
      mutate(.group = toupper(.group))%>%
      rename(LETTER=.group,LABEL=biome)->data_compare_guild
    data_compare_guild$LETTER <- str_replace_all(data_compare_guild$LETTER, " ", "")
    #based on the estimated values, we back transforme the values
    
    data_compare_guild%>%mutate(ori_mean = ifelse(type =="Positive", exp(emmean), -1*emmean^2),
                                ori_lower = ifelse(type =="Positive", exp(lower.CL), -1*lower.CL^2),
                                ori_upper = ifelse(type =="Positive", exp(upper.CL), -1*upper.CL^2))%>%
      dplyr::select(guild, LABEL, emmean,   SE , df ,lower.CL,upper.CL, LETTER , type, LABEL,ori_mean,ori_lower,ori_upper)->data_compare_guild
    
    compare_gain_guild_biome%>%bind_rows(compare_loss_guild_biome)%>%
      rename(LETTER=.group)%>%bind_rows(data_compare_guild)%>%mutate(LETTER = toupper(LETTER))%>%
      mutate(LETTER = str_replace_all(LETTER, " ", ""))%>%
      mutate(new=paste0(LABEL," * ",type))->data_among_guilds_biomes
    
    data_compare_LETTER[[j]]=data_among_guilds_biomes
  }
  
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ggplot(data=data_model_mean_guild[[m]][[i]],aes(fill=type,y=biome ,x=100*mean_value))+
    geom_col(width = 0.3)+
    geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"), 
                 aes(y=biome,yend=biome,xend= 100*mean_value- 100*sd_value, x = 100*mean_value),
                 arrow = arrow(length = unit(0.1, "cm"),angle=90))+
    geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                 aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                 arrow = arrow(length = unit(0.1, "cm"),angle=90))+
    scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#d95f02", "#1b9e77"))+
    theme(legend.position = "none",
          legend.text = element_text(size=8),
          legend.title  = element_text(size=10),
          text = element_text(size = 18),
          plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10), 
          axis.title.x = element_text(size = 10), 
          legend.key.size = unit(0.3, "cm"),
          plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.6, fill = NA))+
    geom_vline(xintercept =0,color="gray",linetype="dashed")+
    ylab("")+
    #xlab(paste(full_name[m]))+
    xlim(xlim_list_am[[i]])+
    xlab("")+
    scale_y_discrete(breaks=unique(data_biome_change_rate_guild[[m]][[i]]$biome),position="right",
                     labels=paste0(rev(c("","","","","")),"(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%distinct( overal_mean )%>%pull()%>%round(3)*100,")"),")"))+
    geom_hline(yintercept = 9,color="gray48",size=10,alpha=0.2)+
    #ggtitle(paste0(scenario[i]))+
    geom_point(data=data_biome_change_rate_guild[[m]][[i]],aes(y=biome,x=100*overal_mean),pch=23,color="black",size=2,fill="#7570b3")+
    
    geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
              aes(y=biome,x = -50),size=3,vjust=-1.6,
              label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
    geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
              aes(y=biome,x = 50),size=3,vjust=-1.6,
              label=c(paste0("(",sprintf("%.1f",data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull( mean_value )%>%round(3)*100,")"),")")))+
    geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative"),
              aes(y=biome,x = -50),size=3,vjust=-3,
              label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Negative")%>%pull(letter))+
    geom_text(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"),
              aes(y=biome,x = 50),size=3,vjust=-3,
              label=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive")%>%pull(letter))
  

