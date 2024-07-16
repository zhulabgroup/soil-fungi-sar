
# model the z values
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


effect=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()
    effect[[i]]=kk[2:14,][,c(1,2,5)]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+ bio2+bio4+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()# for i=9, bio4 was excluded to avoid colinearity
      effect[[i]]=kk[2:14,][,c(1,2,5)]
    }
  }
}

effect[[1]]%>%mutate(var=rownames(effect[[1]]))%>%
  bind_rows(effect[[2]]%>%mutate(var=rownames(effect[[2]])))%>%
  bind_rows(effect[[3]]%>%mutate(var=rownames(effect[[3]])))%>%
  bind_rows(effect[[4]]%>%mutate(var=rownames(effect[[4]])))%>%
  bind_rows(effect[[5]]%>%mutate(var=rownames(effect[[5]])))%>%
  bind_rows(effect[[6]]%>%mutate(var=rownames(effect[[6]])))%>%
  bind_rows(effect[[7]]%>%mutate(var=rownames(effect[[7]])))%>%
  bind_rows(effect[[8]]%>%mutate(var=rownames(effect[[8]])))%>%
  bind_rows(effect[[9]]%>%mutate(var=rownames(effect[[9]])))%>%
  mutate(guild=rep(guild_select,each=13))%>%
  rename_all(~paste0(c("estimate","sd","pva","var","guild")))-> effect_no_plant
###


effect_no_plant$var=factor(effect_no_plant$var,levels=c("logc","soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18"))        
effect_no_plant$guild=factor(effect_no_plant$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           
effect_no_plant$pva[effect_no_plant$pva<0.01]="**"
effect_no_plant$pva[effect_no_plant$pva>0.01&effect_no_plant$pva<0.05]="*"
effect_no_plant$pva[effect_no_plant$pva>0.05]=""


## when plant richness were included

effect=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + richness+(1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()
    effect[[i]]=kk[2:15,][,c(1,2,5)]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+ bio2+bio4+bio8  + bio12 + bio15 + bio18 + richness+(1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()# for i=9, bio4 was excluded to avoid co-linearity
      effect[[i]]=kk[2:15,][,c(1,2,5)]
    }
  }
}

# combine all the effects

effect[[1]]%>%mutate(var=rownames(effect[[1]]))%>%
  bind_rows(effect[[2]]%>%mutate(var=rownames(effect[[2]])))%>%
  bind_rows(effect[[3]]%>%mutate(var=rownames(effect[[3]])))%>%
  bind_rows(effect[[4]]%>%mutate(var=rownames(effect[[4]])))%>%
  bind_rows(effect[[5]]%>%mutate(var=rownames(effect[[5]])))%>%
  bind_rows(effect[[6]]%>%mutate(var=rownames(effect[[6]])))%>%
  bind_rows(effect[[7]]%>%mutate(var=rownames(effect[[7]])))%>%
  bind_rows(effect[[8]]%>%mutate(var=rownames(effect[[8]])))%>%
  bind_rows(effect[[9]]%>%mutate(var=rownames(effect[[9]])))%>%
  mutate(guild=rep(guild_select,each=14))%>%
  rename_all(~paste0(c("estimate","sd", "pva","var","guild")))-> effect_with_plant

effect_with_plant%>%select(var,estimate,sd,guild,pva)->guild_effect_with_plant
guild_effect_with_plant$var=factor(guild_effect_with_plant$var,levels=c("logc","soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18","richness"))        
guild_effect_with_plant$guild=factor(guild_effect_with_plant$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap","all"))                           
guild_effect_with_plant$pva[guild_effect_with_plant$pva<0.01]="**"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.01&guild_effect_with_plant$pva<0.05]="*"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.05]=""


p1=ggplot(guild_effect_no_plant, aes(y =guild , x = var, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "olivedrab1", mid = "#FFFFCC", high = "magenta",midpoint = -0.2)+
  
  scale_x_discrete(breaks=as.character(unique(guild_effect_no_plant$var)),
                   labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  scale_y_discrete(breaks=as.character(unique(guild_effect_no_plant$guild)),
                   labels = c("ECM(N=304)","Litter saprotroph(N=304)","Parasite(N=301)",
                              "Soil saprotroph(N=304)","Wood saprotroph(N=304)","Plant pathogen(N=303)","ACM(N=268)", "Epiphyte(N=270)"))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15,hjust = 1),
        plot.margin = unit(c(1, 0.5, -0.5, 1.5), "cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(aes(x = var, y = guild, label = guild_effect_no_plant$pva),size=6)


ggplot(data=guild_effect_with_plant%>%filter(guild!="all"), aes(x =var , y = guild, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Eeffect size",low = "olivedrab1", mid = "#FFFFCC", high = "magenta")+
  
  scale_x_discrete(breaks=as.character(unique(guild_effect_with_plant$var)),
                   labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ","Pla.rich" ))+
  
  scale_y_discrete(breaks=c("EM","littersap", "para" , "soilsap" ,"woodsap" ,  "plapat" ,   "AM" ,"epiphy"),
                   labels = c("ECM(N=279)","Litter saprotroph(N=279)","Parasite(N=279)",
                              "Soil saprotroph(N=276)","Wood saprotroph(N=279)","Plant pathogen(N=278)","ACM(N=245)", "Epiphyte(N=245)"))+
  
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=90,size=15,hjust = 1),
        plot.margin = unit(c(0, 0.5, 1, 1.5), "cm"))+
  #ggtitle("Guild-specific effect")+
  xlab("")+
  ylab("")+
  geom_text(data=guild_effect_with_plant%>%filter(guild!="all"),aes(x = var, y = guild, label = guild_effect_with_plant%>%filter(guild!="all")%>%
                                                                      dplyr::select(pva)%>%pull()),size=6)



plot_grid(p1,p2,ncol=1,labels=c("(a)","(b)"),label_x = 0.25,label_y = c(1,1,5))

###












# when plant was not included
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

## to exclude the variables that show high correlations
ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[5]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[6]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
ggcorrplot(cor(data[[7]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
ggcorrplot(cor(data[[8]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
ggcorrplot(cor(data[[9]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
#bio1 and bio 4 were correlated

## get the variables included in the model for each guild

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +soilMoisture+cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18+(1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()
    sel_vab_step[[i]]=rownames(kk)[2:dim(kk)[1]]# the first row is the intercept, which is not a targeted variable
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+cec+ sand +bio1+ bio2+ bio8 + bio12 + bio15 + bio18  +(1 |siteIDD), data = data[[i]])#bio4 removed
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()
      sel_vab_step[[i]]=rownames(kk)[2:dim(kk)[1]]
    }
  }
}

## for the variance partitioning part

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  #plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
  #mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],soil_variable*10,climate_variable)%>%rename(soil=product...21,climate=product...22)
  mod=lmer(zvalue~logc+soil+climate+(1|siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

## to get the random effect part and the residuals of the model

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

# to get the fixed effect part of the model

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed
var_fixed%>%select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_noplant

## when the plant diversity was included


soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")
plant=("richness")

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  
  data_temp=bind_cols(data[[i]],soil_variable*10,climate_variable,plant_variable/10)%>%rename(soil=product...22,climate=product...23,plant=product...24)
  mod=lmer(zvalue~logc+soil+climate+plant+(1 |siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}
## when the plant richness was nont included in the model

load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")
model_data_SAR_rich=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18,richness)
model_data_SAR_rich=model_data_SAR_rich[complete.cases(model_data_SAR_rich),]
ggcorrplot(cor(model_data_SAR_rich%>%filter(guild=="all")%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#

# soil c and nitrogen were correlated and need to remove one of them
data=list()
for (i in 1:9)
{
  d=model_data_SAR_rich%>%filter(guild==guild_select[i]) 
  
  d[,c(6:21)]=apply(d[,c(6:21)],2,range01)
  
  data[[i]]=d
}
ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)),hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[5]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[6]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[7]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
ggcorrplot(cor(data[[8]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen

ggcorrplot(cor(data[[9]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand, bio1, bio2,  bio4, bio8, bio12, bio15, bio18, richness)), hc.order = TRUE, type = "lower", lab = TRUE)#
#soil carbon and nitrogen
#bio4 and bio1

### select the variables

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +soilMoisture+cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18+richness+(1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()
    sel_vab_step[[i]]=rownames(kk)[2:dim(kk)[1]]# the first row is the intercept, which is not a targeted variable
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent+soilMoisture+cec+ sand +bio1+ bio2+ bio8 + bio12 + bio15 + bio18  +richness+(1 |siteIDD), data = data[[i]])#bio4 removed
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()
      sel_vab_step[[i]]=rownames(kk)[2:dim(kk)[1]]
    }
  }
}

## for the variance partitioning part
var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_n=which(sel_vab_step[[i]]%in%soil)%>%length()
  climate_n=which(sel_vab_step[[i]]%in%climate)%>%length()
  plant_n=which(sel_vab_step[[i]]%in%plant)%>%length()

    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],soil_variable*100,climate_variable*10,plant_variable)%>%rename(soil=product...22,climate=product...23,plant=product...24)
    mod=lmer(zvalue~logc+soil+climate+plant+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }


## to get the random effect part and the residuals of the model

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

# to get the fixed effect part of the model

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed
var_fixed%>%select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_noplant

# creat the variance partitioning plots

p1=ggplot(data=varp_new_noplant%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Component",breaks=unique(varp_new_noplant$type),
                    labels=c(expression(italic(C)),"S","Cli.",expression(italic(C)*"+S"),expression(italic(C)*"+Cli."),
                             "S+Cli.",expression(italic(C)*"+S+Cli."),"Random","Residuals"),
                    values=c("#F3A332", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple","lavender", "#FFFFCC","gray"))+
  theme(legend.position = c(0.5,0.75),
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.key.size = unit(0.25, "cm") ,
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.5,0, 0.5,0.5), "cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow =7, byrow = TRUE))+
  ylab("Variance explained")+
  xlab("")+
  ggtitle(expression(italic(R["Fixed"]^2)*" =60.6%"))+
  ylim(0,2)

p11=ggplot(data=effect_no_plant%>%filter(guild=="all"),aes(x=estimate,y=1:13))+
  geom_point(size=4,pch=21,color="black",fill=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  
  geom_errorbar(aes(xmin = estimate-1.96*sd, xmax = 1.96*sd+estimate,width=0.2),
                color=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  geom_vline(xintercept = 0,linetype="dashed")+
  scale_y_continuous (breaks=1:13,labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                                             "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=15),
        axis.ticks.length = unit(-0.150, "cm"),
        axis.text.y  = element_text(size=12,hjust=0,margin = margin(r=-70)),
        plot.title=element_text(hjust=0.5,size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 0),
        panel.background = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))+
  ylab("")+
  xlab("Effect size ± 95%CI")+
  xlim(-1.1,0.6)+
  ggtitle("Climate+Soil(N=304)")





# when plants were included

p2=ggplot(data=varp_new%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Component",breaks=varp_new%>%filter(Fractions>0&guild=="all")%>%select(type)%>%pull(),
                    labels=c(expression(italic(C)),"Cli.","P", expression(italic(C)*"+S"),
                             "S+Cli.","S+P", "Cli.+P",expression(italic(C)*"+Cli.+P"),
                             expression(italic(C)*"+S+Cli.+P"),
                             "Random","Residuals"),
                    values=c("#F3A332", "royalblue", "seagreen1", "greenyellow", "forestgreen", "purple","lavender", "white","pink","#FFFFCC","gray"))+
  
  theme(legend.position = c(0.5,0.75),
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=18),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow =11, byrow = TRUE))+
  ylab("")+
  xlab("")+
  ggtitle(expression(italic(R["Fixed"]^2)*" =62.2%"))+
  ylim(0,2)

p22=ggplot(data=guild_effect_with_plant%>%filter(guild=="all"),aes(x=estimate,y=1:14))+
  geom_point(size=4,pch=21,color="black",fill=rep(c("#DE582E","#F3A332","#1868b2","#018a67"),times=c(1,5,7,1)))+
  geom_errorbar(aes(xmin = estimate-1.96*sd, xmax = 1.96*sd+estimate,width=0.2),
                color=rep(c("#DE582E","#F3A332","#1868b2","#018a67"),times=c(1,5,7,1)))+
  geom_vline(xintercept = 0,linetype="dashed")+
  scale_y_continuous (breaks=1:14,labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                                             "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ","Pla.rich." ))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=15),
        axis.ticks.length = unit(-0.150, "cm"),
        plot.title=element_text(hjust=0.5,size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 1),
        axis.text.y  = element_text(size=12,hjust=0,margin = margin(r=-70)),
        panel.background = element_blank(),
        plot.margin = unit(c(0.5, 0.2, 0.5, 0.5), "cm"))+
  ylab("")+
  xlim(-1.5,0.8)+
  xlab("Effect size ± 95%CI")+
  ggtitle("Climate+Soil(N=279)")



plot_grid(p11,p1,p22,p2,ncol=4)

p1=ggplotGrob(p1)
p11=ggplotGrob(p11)

p2=ggplotGrob(p2)
p22=ggplotGrob(p22)


p1$heights=p11$heights
p2$heights=p22$heights
p1$widths=p2$widths


plot_grid(p1,p11,ncol=2,rel_widths = c(0.5,1))

plot_grid(p2,p22,ncol=2,rel_widths = c(0.4,1))

plot_grid(p1,p11,p2,p22,ncol=4,rel_widths = c(0.4,1,0.4,1),labels = c("","(a)","","(b)"),label_size = 14)

