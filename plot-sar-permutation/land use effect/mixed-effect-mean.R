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
      
      df$LABEL=factor(df$LABEL,
                           levels=c("All",
                                    "Temperate Broadleaf & Mixed Forests",
                                    "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                    "Tropical & Subtropical Dry Broadleaf Forests"))
      
      
          #to give a plot id for each cell
    
         # to test difference for species gains for the overall fungal diversity
        
        df1%>%filter(value>0&guild%in%c("all"))->df_gain
        
        #need to transform the resposne variable
        
        df_gain$log_response <- log(df_gain$value)
        
        # reorder the biomes
        df_gain$LABEL=factor(df_gain$LABEL,
                         levels=c("All",
                                  "Temperate Broadleaf & Mixed Forests",
                                                                    "Temperate Conifer Forests","Temperate Grasslands, Savannas & Shrublands",
                                  "Tropical & Subtropical Dry Broadleaf Forests"))
        
        model_all <- gls(log_response ~ LABEL, data = df_gain)# performe the model and get the mean for all groups
        
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
  
        
        
        ggplot(data=com_data,aes(fill=change,y=LABEL ,x=origin_mean))+
          geom_col(width = 0.3)+
          geom_segment(data=com_data,
                       aes(y=LABEL,yend=LABEL,x= ori_low, xend = ori_up),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
        
          geom_text(aes(y=LABEL ,x=origin_mean*1.25,label=.group))+
          xlim(-0.1,0.1)
          
          
          geom_segment(data=data_biome_change_rate_guild[[m]][[i]]%>%filter(type=="Positive"), 
                       aes(y=biome,yend=biome,x = 100*mean_value, xend = 100*mean_value +100*sd_value),
                       arrow = arrow(length = unit(0.1, "cm"),angle=90))+
          scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#ee8c11","#1173EE"))+
          
        
         
    # for individual guilds, we focused on difference response within each biome 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
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
    
    
    
    
    
        
        
        #because the richness is based on the biomes we have, so we do not need to exclud the those biomes we did not consider
        
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
    change_rate_mean_guild[[m]]=change_rate
  }
  change_rate_mean_guild_scenario[[j]]=change_rate_mean_guild
}
