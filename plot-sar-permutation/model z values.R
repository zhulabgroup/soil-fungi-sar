library(dplyr)
library(glmmLasso)
library(ggplot2)

range01 <- function(x) 
{
  return((x - min(x)) / (max(x) - min(x))) 
}



## for the full data
model_data_SAR$logc=2.71828^model_data_SAR$logc

model_data_SAR=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc, zvalue, soilMoisture ,soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)

model_data_SAR=model_data_SAR[complete.cases(model_data_SAR),]

data_all=model_data_SAR%>%filter(guild=="all")

ggcorrplot(cor(data_all%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#

##nitrogen and soil carbon were correlated 
ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[5]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[6]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[7]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[8]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
ggcorrplot(cor(data[[9]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon
# bio1 and bio4 are correlated and i retained the former


# Find pairs with high correlation

data_all[,c(6:21)]=apply(data_all[,c(6:21)],2,range01)


# when use a lambda of 5, several key climate predictors were selected, and when use 8, only 1-2 were selected
# when using 0.01 all were selected
#when using 0.1
# it looks increasing lambda does not lead to less variables selected by identities of the selected variables

#mod <- cv.glmmLasso(form.fix=zvalue ~ logc + soilInCaClpH +organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
form.rnd=list(siteIDD=~1),
dat = data_all,family = gaussian(link = "identity"))#（0.01-10）

# select the lambda for the linear mixed effect model

model_data_SAR$siteIDD=as.factor(model_data_SAR$siteIDD)
guild_select=unique(model_data_SAR$guild)

data=list()
for (i in 1:9)
{
  d=model_data_SAR%>%filter(guild==guild_select[i]) 
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}



set.seed(235)
sel_vab=list()# the variables selected based on the optimal lambda
lambb=numeric()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  if(i<9)
  {
    mod1=cv.glmmLasso(fix=zvalue ~ logc + soilInCaClpH+organicCPercent+soilMoisture+ cec+sand+bio1+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                      rnd=list(siteIDD=~1), 
                      data = data[[i]],
                      family = gaussian(link = "identity"),
                      kfold = 5)#（0.01-10） 
    lambb[i]=mod1$lambda.min
    kk=mod1$glmmLasso.final$coefficients%>%data.frame()%>%filter(.!=0)
    sel_vab[[i]]=rownames(kk)[2:dim(kk)[1]]
    
  }
  else# here bio1 was excluded
  {
    mod1=cv.glmmLasso(fix=zvalue ~ logc + soilInCaClpH+ nitrogenPercent+organicCPercent+soilMoisture+ cec+sand+ bio2 +bio4+ bio8  + bio12 + bio15 + bio18 + richness,
                      rnd=list(siteIDD=~1), 
                      data = data[[i]],
                      family = gaussian(link = "identity"),
                      kfold = 5)#（0.01-10） 
    lambb[i]=mod1$lambda.min
    kk=mod1$glmmLasso.final$coefficients%>%data.frame()%>%filter(.!=0)
    sel_vab[[i]]=rownames(kk)[2:dim(kk)[1]]# the intercept was excluded
  }
}

# based on the select lambda to select the variables with linear mixed effect model 
# extract the effect size of different variables

effect=list()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  random_effect <- "(1|siteIDD)"
  response="zvalue"
  fixed_effects=sel_vab[[i]]
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")# the selected variables
  formula_string <- paste(response, "~", fixed_effects_formula, "+", random_effect)%>%as.formula()
  mod <- lmer(formula_string, data = data[[i]])
  mod2=summary(mod)
  effect[[i]]=mod2$coefficients[,c(1,2,5)]%>%data.frame()%>%mutate(var=rownames(mod2$coefficients))
}
## to plot the effect size
## to merge with the full variables
##get the modified effect size

modified_effect=list()
for (i in 1:9)
{
  modified_effect[[i]]=effect[[2]]%>%dplyr::select(var)%>%left_join(effect[[i]],by="var")%>%filter(!var%in%c("logc","(Intercept)"))%>%replace_na(list(Estimate=0))
}




ggplot()+
  geom_point(data=modified_effect[[1]],pch=21,color="black",aes(x= Estimate,y=1:dim(modified_effect[[1]])[1]),size=4,
             fill=rev(c("seagreen1", "royalblue","royalblue","royalblue","royalblue","royalblue", "royalblue","royalblue","peru","peru","peru", "peru","gray")))+
  geom_segment(data=modified_effect[[1]],size=.8,
               aes(x=modified_effect[[1]]$Estimate-1.96*modified_effect[[1]]$Std..Error,
                   y=1:dim(modified_effect[[1]])[1],
                   xend=modified_effect[[1]]$Estimate+1.96*modified_effect[[1]]$Std..Error,
                   yend=1:dim(modified_effect[[1]])[1]),
               color=rev(c("seagreen1", "royalblue","royalblue","royalblue","royalblue","royalblue", "royalblue","royalblue","peru","peru","peru", "peru","peru")))+
  geom_vline(xintercept = 0,color="red",linetype="dashed",size=.8)+
  scale_y_continuous(breaks=1:13,labels=rev(c("Pla.rich", "Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.", "MDR","MAT","Sand","CEC","Moisture", "SoilC","pH")))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlab("Effect size")+
  ylab("")


### when only soil and climate variables were considered, more plots were included


effect_guild_specific=bind_rows(modified_effect[[2]][,c(1,2,4)],
                                modified_effect[[3]][,c(1,2,4)],
                                modified_effect[[4]][,c(1,2,4)],
                                modified_effect[[5]][,c(1,2,4)],
                                modified_effect[[6]][,c(1,2,4)],
                                modified_effect[[7]][,c(1,2,4)],
                                modified_effect[[8]][,c(1,2,4)],
                                modified_effect[[9]][,c(1,2,4)])%>%mutate(guild=rep(c(guild_select[2:9]),each=13))%>%rename_all(~paste0(c("var","effect","pval","guild")))


effect_guild_specific$var=factor(effect_guild_specific$var,levels=c("soilInCaClpH", "organicCPercent", "soilMoisture" ,   "cec", "sand" , "bio1" ,    "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18" , "richness"))        
effect_guild_specific$guild=factor(effect_guild_specific$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           


effect_guild_specific$pval[effect_guild_specific$pval<0.01]="**"
effect_guild_specific$pval[effect_guild_specific$pval>0.01&effect_guild_specific$pval<0.05]="*"
effect_guild_specific$pval[effect_guild_specific$pval>0.05]=""



ggplot(effect_guild_specific, aes(x =var , y = guild, fill = effect)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  scale_x_discrete(breaks=as.character(unique(effect_guild_specific$var)),
                   labels = c("pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Tem.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" , "Pla.rich."))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15),
        plot.margin = margin(b=-0.5, unit="cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = effect_guild_specific$pval),size=6)

save(effect_guild_specific,file="effect_guild_specific.RData")# need to be corrected





# use the more conventional way for model selection 
# select the data for each fungal guild

data=list()
for (i in 1:9)
{
  d=model_data_SAR%>%filter(guild==guild_select[i]) # this data set included both richness and climate
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}


sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2 +bio4 +bio8+ bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  
  else{
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2 +bio8 + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
}

### if we don't do model selection and just to see all the significant variables

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2+bio8  + bio12 + bio15 + bio18 + richness + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}



### we just focus on the soil and climate variables, load in the initial data

load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")

na_ind <- which(is.na(model_data_SAR$cec))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["cec"]], xy), is.na(r_present_northam[["cec"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$cec[i] <- values["cec"]
}

na_ind <- which(is.na(model_data_SAR$sand))
for(i in na_ind) {
  xy <- cbind(model_data_SAR$lon[i], model_data_SAR$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["sand"]], xy), is.na(r_present_northam[["sand"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  model_data_SAR$sand[i] <- values["sand"]
}


load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")
model_data_SAR$logc=2.71828^model_data_SAR$logc
# do not include plotID
# use this data for modeling
model_data_SAR_climate=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)
model_data_SAR_climate=model_data_SAR_climate[complete.cases(model_data_SAR_climate),]

data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate%>%filter(guild==guild_select[i]) 
  
  d[,c(6:20)]=apply(d[,c(6:20)],2,range01)
  
  data[[i]]=d
}

# construct the model when plant richness was not included in the model

sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture +cec+ sand +bio1+ bio2 +bio4 +bio8+ bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+cec+ sand +bio1+ bio2 +bio8 + bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])#bio4 removed
      mod_sel=step(mod)
      kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}

### if we do not do model selection

### if we don't do model selection and just to see all the significant variables

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}





## extract the effect size

effect=list()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  random_effect <- "(1|siteIDD)"
  response="zvalue"
  fixed_effects=sel_vab_step[[i]]
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")# the selected variables
  formula_string <- paste(response, "~", fixed_effects_formula, "+", random_effect)%>%as.formula()
  mod <- lmer(formula_string, data = data[[i]])
  mod2=summary(mod)
  effect[[i]]=mod2$coefficients[,c(1,2,5)]%>%data.frame()%>%mutate(var=rownames(mod2$coefficients))
}

## to create a data.frame with all variables shown

all_variables=c("(Intercept)", "logc" , "soilInCaClpH" , "organicCPercent", "soilMoisture" ,"cec" , "sand" , "bio1" ,  "bio2" , "bio4" ,  "bio8" , "bio12" , "bio15" ,  "bio18" )%>%data.frame()
names(all_variables)="var"

modified_effect_no_plant=list()
for (i in 1:9)
{
  modified_effect_no_plant[[i]]=all_variables%>%dplyr::select(var)%>%left_join(effect[[i]],by="var")%>%filter(!var%in%c("(Intercept)"))%>%replace_na(list(Estimate=0))
}
# we may consider adding the logc term


# to combine all the effects for different guilds

effect_non_plant=bind_rows(modified_effect[[1]],
         modified_effect[[2]],
         modified_effect[[3]],
         modified_effect[[4]],
         modified_effect[[5]],
         modified_effect[[6]],
         modified_effect[[7]],
         modified_effect[[8]],
         modified_effect[[9]])%>%data.frame()%>%mutate(guild=rep(guild_select,each=12))%>%rename_all(~paste0(c("var","estimate","sd","pva","guild")))

save(effect_non_plant,file="effect_non_plant.RData")

## when all the guilds were included



effect_non_plant%>%filter(guild!="all")%>%select(var,estimate,guild,pva)->guild_effect_no_plant

guild_effect_no_plant$var=factor(guild_effect_no_plant$var,levels=c("soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18"))        
guild_effect_no_plant$guild=factor(guild_effect_no_plant$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           


guild_effect_no_plant$pva[guild_effect_no_plant$pva<0.01]="**"
guild_effect_no_plant$pva[guild_effect_no_plant$pva>0.01&guild_effect_no_plant$pva<0.05]="*"
guild_effect_no_plant$pva[guild_effect_no_plant$pva>0.05]=""



p2=ggplot(guild_effect_no_plant, aes(y =guild , x = var, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+

  scale_x_discrete(breaks=as.character(unique(guild_effect_no_plant$var)),
                   labels = c("pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Tem.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  scale_y_discrete(breaks=as.character(unique(guild_effect_no_plant$guild)),
                   labels = c("ECM(N=304)","Litter saprotroph(N=304)","Parasite(N=301)",
                   "Soil saprotroph(N=304)","Wood saprotroph(N=304)","Plant pathogen(N=303)","ACM(N=268)", "Epiphyte(N=270)"))+

  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15,hjust = 1),
        plot.margin = margin(b=-0.5, unit="cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = guild_effect_no_plant$pva),size=6)

plot_grid(p1,p3,ncol=2,labels=c("(a)","(b)"),label_x = 0.25)
  

##if we do not do model selection for the data that include 

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ soilInCaClpH +soilorganicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ soilInCaClpH +organicCPercent +cec+ sand +bio1+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}



## when richness was included in the model

load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")

model_data_SAR_rich=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18,richness)

model_data_SAR_rich=model_data_SAR_rich[complete.cases(model_data_SAR_rich),]

data=list()
for (i in 1:9)
{
  d=model_data_SAR_rich%>%filter(guild==guild_select[i]) 
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}

# construct the model

sel_vab_step=list()
for (i in 1:9)
  
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +soilMoisture+cec+ sand +bio1+ bio2 +bio4 + bio12 + bio15 + bio18 + richness+(1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+cec+ sand +bio1+ bio2  + bio12 + bio15 + bio18  + richness+(1 |siteIDD), data = data[[i]])#bio4 removed
      mod_sel=step(mod)
      kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}

###

effect=list()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  random_effect <- "(1|siteIDD)"
  response="zvalue"
  fixed_effects=sel_vab_step[[i]]
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")# the selected variables
  formula_string <- paste(response, "~", fixed_effects_formula, "+", random_effect)%>%as.formula()
  mod <- lmer(formula_string, data = data[[i]])
  mod2=summary(mod)
  effect[[i]]=mod2$coefficients[,c(1,2,5)]%>%data.frame()%>%mutate(var=rownames(mod2$coefficients))
}






all_variables=c("(Intercept)", "logc" , "soilInCaClpH" , "organicCPercent", "soilMoisture" ,"cec" , "sand" , "bio1" ,  "bio2" , "bio4" ,  "bio8" , "bio12" , "bio15" ,  "bio18","richness" )%>%data.frame()
names(all_variables)="var"

modified_effect=list()
for (i in 1:9)
{
  modified_effect[[i]]=all_variables%>%dplyr::select(var)%>%left_join(effect[[i]],by="var")%>%filter(!var%in%c("(Intercept)"))%>%replace_na(list(Estimate=0))
}

# to combine all the effects for different guilds

effect_with_plant=bind_rows(modified_effect[[1]],
                           modified_effect[[2]],
                           modified_effect[[3]],
                           modified_effect[[4]],
                           modified_effect[[5]],
                           modified_effect[[6]],
                           modified_effect[[7]],
                           modified_effect[[8]],
                           modified_effect[[9]])%>%data.frame()%>%mutate(guild=rep(guild_select,each=13))%>%rename_all(~paste0(c("var","estimate","sd","pva","guild")))
###create the plots



## for individual guilds



effect_with_plant%>%filter(guild!="all")%>%select(var,estimate,guild,pva)->guild_effect_with_plant

guild_effect_no_plant$var=factor(guild_effect_no_plant$var,levels=c("soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18"))        
guild_effect_no_plant$guild=factor(guild_effect_no_plant$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           


guild_effect_no_plant$pva[guild_effect_no_plant$pva<0.01]="**"
guild_effect_no_plant$pva[guild_effect_no_plant$pva>0.01&guild_effect_no_plant$pva<0.05]="*"
guild_effect_no_plant$pva[guild_effect_no_plant$pva>0.05]=""



p2=ggplot(guild_effect_no_plant, aes(x =var , y = guild, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  
  scale_x_discrete(breaks=as.character(unique(guild_effect_no_plant$var)),
                   labels = c("pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Tem.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  scale_y_discrete(breaks=as.character(unique(guild_effect_no_plant$guild)),
                   labels = c("ECM(N=304)","Litter saprotroph(N=304)","Parasite(N=301)",
                              "Soil saprotroph(N=304)","Wood saprotroph(N=304)","Plant pathogen(N=303)","ACM(N=268)", "Epiphyte(N=270)"))+
  
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15,hjust = 0),
        plot.margin = margin(b=-0.5, unit="cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = guild_effect_no_plant$pva),size=6)



####

guild_effect_with_plant$var=factor(guild_effect_with_plant$var,levels=c("soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18","richness"))        
guild_effect_with_plant$guild=factor(guild_effect_with_plant$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           


guild_effect_with_plant$pva[guild_effect_with_plant$pva<0.01]="**"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.01&guild_effect_with_plant$pva<0.05]="*"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.05]=""



p4=ggplot(guild_effect_with_plant, aes(x =var , y = guild, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  
  scale_x_discrete(breaks=as.character(unique(guild_effect_with_plant$var)),
                   labels = c("pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Tem.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ","Pla.rich" ))+
  scale_y_discrete(breaks=as.character(unique(guild_effect_with_plant$guild)),
                   labels = c("ECM(N=279)","Litter saprotroph(N=279)","Parasite(N=279)",
                              "Soil saprotroph(N=276)","Wood saprotroph(N=279)","Plant pathogen(N=278)","ACM(N=245)", "Epiphyte(N=245)"))+
  
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15,hjust = 1),
        plot.margin = margin(b=-0.5, unit="cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = guild_effect_with_plant$pva),size=6)

a=numeric()
for (i in 1:9){
  a[i]=dim(data[[i]])[1]
}


# variance partitioning for the relative impact of different variables
# when plants were taken into account


effect=list()
for (i in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  random_effect <- "(1|siteIDD)"
  response="zvalue"
  fixed_effects=sel_vab_step[[i]]
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")# the selected variables
  formula_string <- paste(response, "~", fixed_effects_formula, "+", random_effect)%>%as.formula()
  mod <- lmer(formula_string, data = data[[i]])
  mod2=summary(mod)
  effect[[i]]=mod2$coefficients[,c(1,2,5)]%>%data.frame()%>%mutate(var=rownames(mod2$coefficients))
}


# for the data when plant richness was included

soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")
plant=("richness")

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_n=which(sel_vab_step[[i]]%in%soil)%>%length()
  climate_n=which(sel_vab_step[[i]]%in%climate)%>%length()
  plant_n=which(sel_vab_step[[i]]%in%plant)%>%length()
  
  if (soil_n>0&&climate_n>0&&plant_n>0)
  {
    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],soil_variable,climate_variable,plant_variable)%>%rename(soil=product...22,climate=product...23,plant=product...24)
    mod=lmer(zvalue~logc+soil+climate+plant+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when climate variables are lacking
  else if (soil_n>0&&climate_n<1&&plant_n>0)
  {
    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    
    plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],soil_variable,plant_variable)%>%rename(soil=product...22,plant=product...23)
    mod=lmer(zvalue~logc+soil+plant+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when soil variables are lacking
  #else if (soil_n<1&&climate_n>0&&plant_n>0)
  else
  {
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],climate_variable,plant_variable)%>%rename(climate=product...22,plant=product...23)
    mod=lmer(zvalue~logc+climate+plant+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
}

# to combine all the effects

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))

kk$type=gsub("\\d", "", kk$type)

kk$type=gsub("\\.{3}", "", kk$type)

var_part_data=bind_cols(k,kk)

var_part_data$type=gsub(" ","",var_part_data$type)

var_part_data%>%data.frame()%>%rename_all(~paste0(c("fraction","total","guild","type")))# total 16 types


# to see the total variance for each guild

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var
total_var%>%select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp
temp%>%select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp
var_part_data%>%select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new

varp_with_plant=varp_new

save(varp_with_plant,file="varp_with_plant.RData")



varp_with_plant$guild=factor(varp_with_plant$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           

ggplot(data=varp_with_plant%>%filter(Fractions>0), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Component",breaks=unique(varp_with_plant$type),
                    labels=c("Logc","Soil","Climate","Plant","Logc+Soil","Logc+Climate",
                             "Soil+Climate","Logc+Plant","Soil+Plant","Climate+Plant",
                             "Logc+Soil+Plant","logc+Soil+Plant","Logc+Climate+Plant","Soil+Climate+Plant",
                             "Logc+Soil+Climate+Plant","Random effect","Residuls"),
                    values=c("tan", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple","lavender", "orange","black","yellow","red","blue","royalblue","mediumpurple","pink","#FFFFCC","gray"))+
  scale_x_discrete(breaks=unique(varp_new$guild),
                   labels = c("All(N=279)","ECM(N=279)","Litter sapro.(N=279)","Parasite(N=279)",
                              "Soil sapro.(N=276)","Wood sapro.(N=279)","Plant patho.(N=278)","AM(N=245)", "Epiphyte(N=245)"))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 1),
        plot.margin = margin(b=-0.5, unit="cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow = 17, byrow = TRUE))+
  ylab("Variance explained")+
  coord_flip()+
  xlab("")

# for the plot when all guilds were included

# order the different components
varp_with_plant$type=factor(varp_with_plant$type,levels=rev(unique(varp_with_plant$type)))

p2=ggplot(data=varp_with_plant%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Components",breaks=varp_with_plant%>%filter(Fractions>0&guild=="all")%>%distinct(type)%>%pull(type),
                    labels=c(expression(italic(C)),"Climate","Plant","Logc+Soil","Soil+Plant"
                             ,"Climate+Plant","Logc+Soil+Climate","Soil+Climate+Plant","Random effect","Residuls"),
                    values=c("tan", "royalblue", "seagreen1", "greenyellow", "forestgreen", "purple","yellow", "orange","white","lavender"))+
  theme(legend.position = c(0.5,0.750),
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "NA"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size=8,hjust=0))+
  guides(fill = guide_legend(nrow = 10, byrow = TRUE))+
  ylab("")+
  ylim(0,2)+
ggtitle(expression(italic(R["Fixed"]^2)*" =65.4%"))



p22=ggplot()+
  geom_point(data=modified_effect[[1]],pch=21,color="black",aes(x= Estimate,y=1:dim(modified_effect[[1]])[1]),size=4,
             fill=rev(c("seagreen1","white","white","royalblue","white","royalblue", "white","royalblue","white","white","darkgoldenrod1", "white","white","tan")))+
  geom_segment(data=modified_effect[[1]],size=.8,
               aes(x=modified_effect[[1]]$Estimate-1.96*modified_effect[[1]]$Std..Error,
               y=1:dim(modified_effect[[1]])[1],
                xend=modified_effect[[1]]$Estimate+1.96*modified_effect[[1]]$Std..Error,
                yend=1:dim(modified_effect[[1]])[1]),
               color=rev(c("seagreen1","gray","royalblue","royalblue","gray","royalblue", "gray","royalblue","gray","gray","darkgoldenrod1", "gray","gray","tan")))+
  geom_vline(xintercept = 0,color="red",linetype="dashed",size=.8)+
  scale_y_continuous(breaks=1:14,labels=rev(c("Pla.rich.","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.", "MDR","MAT","Sand","Cec","Moisture", "SoilC","pH",expression(italic(C)))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(-0.10, "cm"),
        axis.text.y  = element_text(size=12,hjust=0,margin = margin(r=-60),
color=c("tan","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","seagreen1")),
        plot.title=element_text(hjust=0.5,size=15),
        axis.text.x  = element_text(size=15))+
  ggtitle("Climate+soil+plant (N=279)")+
  xlab("Effect size ± 95% CI")+
  ylab("")+
  xlim(-1.5,0.6)

p2=ggplotGrob(p2)

p22=ggplotGrob(p22)

p2$heights=p22$heights

p6=plot_grid(p2,p22,rel_widths = c(1,2))

# when plant richness was not included


soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_n=which(sel_vab_step[[i]]%in%soil)%>%length()
  climate_n=which(sel_vab_step[[i]]%in%climate)%>%length()

  if (soil_n>0&&climate_n>0)
  {
    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
     data_temp=bind_cols(data[[i]],soil_variable,climate_variable)%>%rename(soil=product...21,climate=product...22)
    mod=lmer(zvalue~logc+soil+climate+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when climate variables are lacking
  else if (soil_n>0&&climate_n<1)
  {
    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    
    data_temp=bind_cols(data[[i]],soil_variable)%>%rename(soil=product)
    mod=lmer(zvalue~logc+soil+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when soil variables are lacking
  else if (soil_n<1&&climate_n>0)
  {
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],climate_variable)%>%rename(climate=product)
    mod=lmer(zvalue~logc+climate+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  # when both climate and soil variables are lacking
  else 
  {
    mod=lmer(zvalue~logc+(1 |siteIDD),data=data[[i]])
    #a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=r.squaredGLMM(mod)
    var_par[[i]]=r.squaredGLMM(mod)
  }
}

# to combine all the effects

var_climate=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

var_climate_name=rownames(var_climate)%>%data.frame()%>%rename_all(~paste0("type"))


var_climate_data=bind_cols(var_climate_name,var_climate)


var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))%>%
select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()%>%
select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

var_climate_data%>%select(Fractions,guild,type)%>%filter(!str_detect(type,"Total"))%>%bind_rows(temp)%>%filter(type!="fix")->tem_climate

tem_climate$type=gsub("\\d", "", tem_climate$type)

tem_climate$type=gsub("\\.{3}", "", tem_climate$type)
tem_climate$type=gsub(" ", "", tem_climate$type)

# create a data.frame for the two

df=data.frame(Fractions=c(0.66806127,0.73490360),guild=c("littersap","para"),type=c("Uniquetologc","Uniquetologc"))

tem_climate%>%bind_rows(df)%>%filter(!is.na(Fractions))->var_part_climate_data

######

save(var_part_climate_data,file="var_part_climate_data.RData")

var_part_climate_data$guild=factor(var_part_climate_data$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           

ggplot(data=var_part_climate_data%>%filter(Fractions>0), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+

  scale_fill_manual("Component",breaks=unique(var_part_climate_data$type),
                    labels=c("Logc","Soil","Climate","Logc+Soil","Logc+Climate",
                             "Soil+Climate","Logc+Soil+Climate","Random effect","Residuls"),
                    values=c("tan", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple","lavender", "#FFFFCC","gray"))+
  scale_x_discrete(breaks=unique(varp_new$guild),
                   labels = c("All(N=304)","ECM(N=304)","Litter sapro.(N=304)","Parasite(N=301)",
                              "Soil sapro.(N=304)","Wood sapro.(N=304)","Plant patho.(N=303)","AM(N=268)", "Epiphyte(N=270)"))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 1),
        plot.margin = margin(b=-0.5, unit="cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow = 17, byrow = TRUE))+
  ylab("Variance explained")+
  coord_flip()+
  xlab("")


# when plant richness was not included
# display the component with order

p1=ggplot()+
  geom_point(data=modified_effect_no_plant[[1]],pch=21,color="black",aes(x= Estimate,y=1:dim(modified_effect_no_plant[[1]])[1]),size=4,
             fill=rev(c("white","blue","royalblue","white","royalblue", "white","royalblue","white","white","darkgoldenrod1", "white","white","tan")))+
  geom_segment(data=modified_effect_no_plant[[1]],size=.8,
               aes(x=modified_effect_no_plant[[1]]$Estimate-1.96*modified_effect_no_plant[[1]]$Std..Error,
                   y=1:dim(modified_effect_no_plant[[1]])[1],
                   xend=modified_effect_no_plant[[1]]$Estimate+1.96*modified_effect_no_plant[[1]]$Std..Error,
                   yend=1:dim(modified_effect_no_plant[[1]])[1]),
               color=rev(c("gray","royalblue","royalblue","gray","royalblue", "gray","royalblue","gray","gray","darkgoldenrod1", "gray","gray","tan")))+
  geom_vline(xintercept = 0,color="red",linetype="dashed",size=.8)+
  scale_y_continuous(breaks=1:13,labels=rev(c("Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.", "MDR","MAT","Sand","Cec","Moisture", "SoilC","pH",expression(italic(C)))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(-0.10, "cm"),
        axis.text.y  = element_text(size=12,
        color=c("tan","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue","royalblue"),
        hjust=0,
       margin = margin(r=-60)),
        plot.title=element_text(hjust=0.5,size=15),
        axis.text.x  = element_text(size=15))+
  ggtitle("Climate+soil (N=304)")+
  xlab("Effect size ± 95% CI")+
  xlim(-1.2,0.5)


var_part_climate_data$type=factor(var_part_climate_data$type,levels = rev(unique(var_part_climate_data$type)))

p11=ggplot(data=var_part_climate_data%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  ggtitle(expression(italic(R["Fixed"]^2)*" =63.8%"))+
  scale_fill_manual("Components",breaks=var_part_climate_data%>%filter(Fractions>0&guild=="all")%>%distinct(type)%>%pull(type),
                    labels=c(expression(italic(C)),"Soil","Climate","Logc+Soil",
                             "Soil+Climate","Random effect","Residuals"),
                    values=c("tan", "darkgoldenrod1", "royalblue", "greenyellow", "seagreen1", "white","lavender"))+
  theme(legend.position = c(0.5,0.7802),
    panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=15),
        axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "NA"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size=8))+
  guides(fill = guide_legend(nrow = 7, byrow = TRUE))+
  ylab("")+
  xlab("")+
  ylim(0,2)+
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  xlab("")
  
  
## for individual guilds
p1=ggplotGrob(p1)
p11=ggplotGrob(p11)
p2=ggplotGrob(p2)
p22=ggplotGrob(p22)



p1$heights=p11$heights

p1$heights=p2$heights

p2$heights=p22$heights

p5=plot_grid(p11,p1,rel_widths = c(1,1.5))

p5=ggplotGrob(p5)
p6=ggplotGrob(p6)

p5$heights=p6$heights

plot_grid(p5,p6,ncol=2)


p1$heights=p2$heights

p11$heights=p22$heights

plot_grid(p11,p1,p2,p22,ncol=4)

plot_grid(p11,p2,ncol=2)

plot_grid(p11,p1,p2,p22,ncol=4)







