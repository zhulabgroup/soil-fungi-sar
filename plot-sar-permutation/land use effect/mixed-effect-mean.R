#to compare the mean change rate among different guilds
# for the species loss rate and gain rate

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")

species_change_climate_rcp245=species_change_climate_rcp245[,-c(1,2)]

species_change_climate_rcp585=readRDS("species_change_climate_rcp585.rds")

species_change_climate_rcp585=species_change_climate_rcp585[,-c(1,2)]
# for the climate change impact, 41394 cells were included

observations=c(41725,41725,41394,41394)

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
        
        df1%>%filter(value>0&guild%in%c("all"))->df_gain
        
        #need to transform the response variable
        
        df_gain$log_response <- log(df_gain$value)
        
          
        model_all <- lmer(log_response ~ LABEL*guild+(1|plotid), data = df_gain)# performe the model and get the mean for all groups
        
        pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
        
        pairwise_all%>%data.frame()%>%
        filter(LABEL=="All")%>%mutate(.group="hh")%>%
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
        
   # if we consider species loss
        
        df1%>%filter(value<0&guild%in%c("all"))->df_loss
        
        df_loss$sqrt_response <- sqrt(abs(df_loss$value))
        
        model_all <- gls(sqrt_response  ~ LABEL, data = df_loss)# performe the model and get the mean for all groups
        
        pairwise_all=emmeans(model_all, ~ LABEL)# get the mean for each biome
        
        
        pairwise_all%>%data.frame()%>%
          filter(LABEL=="All")%>%mutate(.group="hk")%>%
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
        
        
  com_data%>%mutate(change=ifelse(origin_mean > 0, "positive", "negative"))->com_data
        
        
  
  com_data$LABEL=factor(com_data$LABEL,
                       levels=c("All",
                                "Temperate Broadleaf & Mixed Forests",
                                "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                "Tropical & Subtropical Dry Broadleaf Forests"))
  
        
        
          
        
         
    # for individual guilds, we focused on difference response within each biome 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        #based on the difference to label letters
        
        
          
         
         
        
         
      #to compare the mean among the guilds
         
         observations=c(41725,41725,41394,41394)
         
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
         
         
        
         
        
        df1%>%filter(value>0&!(guild%in%c("all","para","epiphy","littersap","woodsap")))->df1
        
        #df1%>%filter(value<0&guild%in%c("all"))->df1
        
        df1$log_response <- log(df1$value)
        
        
        mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df1)
        
        mod=lmer(log_response~LABEL,data=df1)
        
        mod=lmer(log_response~guild+(1 | plotid),data=df1)# compare the mean among
        
        emmeans(mod, pairwise ~ guild,pbkrtest.limit = 60649)
        
        emm <- emmeans(mod, ~ guild,pbkrtest.limit = 37377)
        
        
        
    comp_mean= emmeans(mod, pairwise ~ guild | LABEL,pbkrtest.limit = 60649)
        
    
    # when we look at species loss
    
    
    
    
    
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
     
     
    
     ## when we look at the net effect
    
    df1%>%filter(!(guild%in%c("all","para","epiphy","woodsap","littersap")))->df1
    
    
    df1$log_response <- sign(df1$value) * log(abs(df1$value) + 1)
    
   
    mod=lmer(log_response~guild*LABEL+(1 | plotid),data=df1)
    
    comp_mean= emmeans(mod, pairwise ~ guild | LABEL,pbkrtest.limit = 250350)
    
    # if we do not consider the biomes
    
    mod=lmer(log_response~guild+(1 | plotid),data=df1)
    
    
    
    
#### for the species gain


get(data[j])%>%
  bind_cols(rep(grid_level_biomes$LABEL,9))%>%
  rename_all(~paste0(c("guild","value","LABEL")))%>%
  filter(LABEL%in%biomes_four&!is.na(value))%>%
  mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered



# select the guild of your interest
# for species gain

df1%>%filter(value>0&guild%in%c("EM","soilsap","AM","plapat"))->df2
# need to transform the data 

df2$log_response=log(df2$value)

mod=lmer(log_response~LABEL*guild+(1 | plotid),data=df2)

anova(mod, type = "III")

compare_mean=emmeans(mod, pairwise ~ guild | LABEL,pbkrtest.limit = 71384)#the last parameter could change


# get the values associated with the letters


multcomp::cld(object = compare_mean$emmeans,
              Letters = letters)->compare_data


# need to get the overall mean and when all the guilds were combined

get(data[j])%>%
  bind_cols(rep(grid_level_biomes$LABEL,9))%>%
  rename_all(~paste0(c("guild","value","LABEL")))%>%
  filter(LABEL%in%biomes_four&!is.na(value))%>%
  mutate(plotid=rep(1:observations[j],9))->df1# na cells were filtered




