
# model the z values

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standard sar neon")

full_parameter_data=readRDS("full_parameter_data.rds")

# need to combine the complete climate and soil variables
#
full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%
  left_join(model_data%>%dplyr::select(plotID,siteIDD,bio1,bio2,bio4,bio8,bio12,bio15,bio18,richness)%>%distinct(),by="plotID")%>%
  left_join(soil_mean%>%dplyr::select(-lon,-lat)%>%distinct(),by="plotID")->model_data_SAR_rarefaction



model_data_SAR_rarefaction$logc=2.71828^model_data_SAR_rarefaction$logc

# do not include plotID
# use this data for modeling when plant richness was not included
model_data_SAR_climate=model_data_SAR_rarefaction%>%dplyr::select(siteIDD,plotID,guild,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)
model_data_SAR_climate=model_data_SAR_climate[complete.cases(model_data_SAR_climate),]

guild_select=unique(model_data_SAR_climate$guild)

data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:18)]=apply(d[,c(4,6:18)],2,range01)
  
  data[[i]]=d
}

a=numeric()
for(i in 1:9)
  {
  a[i]=dim(data[[i]])[1]
}


# to check the linearity of the variables for each data set
high_cor_variables=list()
for (i in 1:9)
  {
  cor_matrix <- cor(data[[i]][,6:18])
  # Find columns with correlation coefficients higher than 0.8
  high_cor_variables[[i]] <- colnames(data[[i]][,6:18])[apply(cor_matrix, 2, function(x) any(abs(x) > 0.70 & abs(x) < 1))]
  
}
# we select all the variables with r <0.7

# when plant richness was considered

model_data_SAR_climate_plant=model_data_SAR_rarefaction%>%dplyr::select(siteIDD,plotID,guild,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18,richness)

model_data_SAR_climate_plant=model_data_SAR_climate_plant[complete.cases(model_data_SAR_climate_plant),]

guild_select=unique(model_data_SAR_climate_plant$guild)

data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate_plant%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:19)]=apply(d[,c(4,6:19)],2,range01)
  
  data[[i]]=d
}

a=numeric()
for(i in 1:9)
{
  a[i]=dim(data[[i]])[1]
}

high_cor_variables=list()
for (i in 1:9)
{
  cor_matrix <- cor(data[[i]][,6:19])
  # Find columns with correlation coefficients higher than 0.8
  high_cor_variables[[i]] <- colnames(data[[i]][,6:19])[apply(cor_matrix, 2, function(x) any(abs(x) > 0.70 & abs(x) < 1))]
  
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
      mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+bio4+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
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
effect_no_plant$pva[effect_no_plant$pva<=0.001]="***"
effect_no_plant$pva[effect_no_plant$pva>0.001&effect_no_plant$pva<=0.01]="**"
effect_no_plant$pva[effect_no_plant$pva>0.01&effect_no_plant$pva<=0.05]="*"
effect_no_plant$pva[effect_no_plant$pva>0.05]=""

guild_effect_no_plant=effect_no_plant%>%filter(guild!="all")


## when plant richness were included



high_cor_variables=list()
for (i in 1:9)
{
  cor_matrix <- cor(data[[i]][,6:18])
  # Find columns with correlation coefficients higher than 0.8
  high_cor_variables[[i]] <- colnames(data[[i]][,6:18])[apply(cor_matrix, 2, function(x) any(abs(x) > 0.70 & abs(x) < 1))]
  
}

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

effect_with_plant%>%dplyr::select(var,estimate,sd,guild,pva)%>%filter(guild!="all")->guild_effect_with_plant

guild_effect_with_plant$var=factor(guild_effect_with_plant$var,levels=c("logc","soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18","richness"))        
guild_effect_with_plant$guild=factor(guild_effect_with_plant$guild,levels=c("AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           
guild_effect_with_plant$pva[guild_effect_with_plant$pva<=0.001]="***"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.001&guild_effect_with_plant$pva<=0.01]="**"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.01&guild_effect_with_plant$pva<=0.05]="*"
guild_effect_with_plant$pva[guild_effect_with_plant$pva>0.05]=""






p1=ggplot(guild_effect_no_plant, aes(y =guild , x = var, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Effect size",low = "darkseagreen2", mid = "#FFFFCC", high = "violetred1",limits=c(-1,1))+
  geom_text(aes(x = var, y = guild, label = guild_effect_no_plant$pva),size=6)+

  
  scale_x_discrete(breaks=levels(guild_effect_no_plant$var),
                   labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  scale_y_discrete(breaks=as.character(unique(guild_effect_no_plant$guild)),
                   labels = c("ACM(N=329)","ECM(N=447)","Soil saprotroph(N=448)","Litter saprotroph(N=445)","Wood saprotroph(N=430)","Plant pathogen(N=430)","Parasite(N=416)",
                              "Epiphyte(N=332)"))+
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


p2=ggplot(data=guild_effect_with_plant, aes(x =var , y = guild, fill = estimate)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  scale_fill_gradient2("Effect size",low = "darkseagreen2", mid = "#FFFFCC", high = "violetred1",limits=c(-1,1))+
  geom_text(data=guild_effect_with_plant,aes(x = var, y = guild, label = guild_effect_with_plant%>%
                                               dplyr::select(pva)%>%pull()),size=6)+
  scale_x_discrete(breaks=levels(guild_effect_with_plant$var),
                   labels = c(expression(italic(C)),"pH", "SoilC", "Moisture", "CEC", "Sand" , "MAT" ,   
                              "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ","Pla.rich" ))+
  
  scale_y_discrete(breaks=c("EM","littersap", "para" , "soilsap" ,"woodsap" ,  "plapat" ,   "AM" ,"epiphy"),
                   labels = c("ECM(N=410)","Litter saprotroph(N=409)","Parasite(N=382)",
                              "Soil saprotroph(N=412)","Wood saprotroph(N=394)","Plant pathogen(N=395)","ACM(N=302)", "Epiphyte(N=300)"))+
  
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
  geom_text(data=guild_effect_with_plant,aes(x = var, y = guild, label = guild_effect_with_plant%>%
                                                                      dplyr::select(pva)%>%pull()),size=6)
plot_grid(p1,p2,ncol=1,label_x = 0.25,label_y = c(1,1,5))

###



# when plant was not included


# do not include plotID
# use this data for modeling

model_data_SAR_rarefaction$logc=2.71828^model_data_SAR_rarefaction$logc
# do not include plotID
# use this data for modeling
model_data_SAR_climate=model_data_SAR_rarefaction%>%dplyr::select(siteIDD,plotID,guild,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)
model_data_SAR_climate=model_data_SAR_climate[complete.cases(model_data_SAR_climate),]

guild_select=unique(model_data_SAR_climate$guild)



data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:18)]=apply(d[,c(4,6:18)],2,range01)
  
  data[[i]]=d
}


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
    mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
  climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
  #plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
  #mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],soil_variable*10,climate_variable)%>%rename(soil=product...19,climate=product...20)
  mod=lmer(zvalue~logc+soil+climate+(1|siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

## to get the random effect part and the residuals of the model

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%dplyr::select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%dplyr::select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

# to get the fixed effect part of the model

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed
var_fixed%>%dplyr::select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_noplant

## when the plant diversity was included


data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate_plant%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:19)]=apply(d[,c(4,6:19)],2,range01)
  
  data[[i]]=d
}


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


soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")
plant=("richness")


## when the plant richness was included in the model






## for the variance partitioning part
var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
    soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],soil_variable*100,climate_variable*10,plant_variable)%>%rename(soil=product...20,climate=product...21,plant=product...22)
    mod=lmer(zvalue~logc+soil+climate+plant+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }


## to get the random effect part and the residuals of the model

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%dplyr::select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%dplyr::select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

# to get the fixed effect part of the model

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed
var_fixed%>%dplyr::select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_withplant

# create the variance partitioning plots


ggplot(data=varp_new_withplant%>%filter(Fractions>0&guild!="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.5)+

  scale_fill_manual("Component",breaks=unique(varp_new_withplant$type),
                    labels=c(expression(italic(C)),"S","Cli.","P", expression(italic(C)*"+S"),expression(italic(C)*"+Cli."),
                             "S+Cli.",  expression(italic(C)*"+P"),"S+P","Cli.+P",
                             
                             expression(italic(C)*"+S+Cli."),
                             expression(italic(C)*"+S+P"),
                             expression(italic(C)*"+P+Cli."),
                             "Cli.+P+S",expression(italic(C)*"+P+Cli.+S"),  
                             "Random","Residuals"),
                    values=c("#F3A332", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple","lavender", "red","tan","chocolate1", "#037f77", "royalblue", "yellow", "forestgreen", "#7c1a97","#FFFFCC", "gray"))+

  theme(legend.position = "right",
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.key.size = unit(0.45, "cm") ,
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.5,0, 0.5,0.5), "cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow =17, byrow = TRUE))+
  ylab("Variance explained")+
  xlab("")+
    coord_flip()+
  scale_x_discrete(labels=c("Wood saprotroph(58.6%)","Soil saprotroph(38.5%)","Plant pathogen(49.0%)","Parasite(40.7%)","Litter saprotroph(63.0%)","Epiphyte(22.0%)","EM(62.5%)","AM(10.2%)"),
                   breaks=c("woodsap","soilsap","plapat","para","littersap","epiphy","EM","AM"))


#when all the guilds were considered

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
  guides(fill = guide_legend(nrow =17, byrow = TRUE))+
  ylab("Variance explained")+
  xlab("")+
  ggtitle(expression(italic(R["Fixed"]^2)*" =56.4%"))+
  ylim(0,2)

p11=ggplot(data=effect_no_plant%>%filter(guild=="all"),aes(x=estimate,y=1:13))+
  geom_point(size=4,pch=21,color="black",fill=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  
  geom_errorbar(aes(xmin = estimate-1.96*sd, xmax = 1.96*sd+estimate,width=0.2),
                color=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  geom_vline(xintercept = 0,linetype="dashed")+
  scale_y_continuous (breaks=1:13,labels = c(expression(italic(C)),"pH",  "Moisture","SoilC", "CEC", "Sand" , "MAT" ,   
                                             "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=15),
        axis.ticks.length = unit(-0.150, "cm"),
        axis.text.y  = element_text(size=12,hjust=0,margin = margin(r=-70)),
        plot.title=element_text(hjust=0.5,size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 0),
        panel.background = element_blank(),
        plot.margin = unit(c(0.5, 0.2, 0.5, 0.5), "cm"))+
  ylab("")+
  xlab("Effect size ± 95%CI")+
  xlim(-0.6,0.6)+
  ggtitle("Climate+Soil(N=453)")+
  annotate("text",x=-0.31,y=1,label="***",size=8)+
annotate("text",x=0.265,y=7,label="*",size=8)+
annotate("text",x=0.265,y=9,label="*",size=8)+
  annotate("text",x=0.24583,y=11,label="*",size=8)+
annotate("text",x=-0.16820583,y=3,label="*",size=8)+
  xlim(-0.5,0.5)




# when plants were included

p2=ggplot(data=varp_new_withplant%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Component",breaks=varp_new_withplant%>%filter(Fractions>0&guild=="all")%>%dplyr::select(type)%>%pull(),
                    labels=c(expression(italic(C)),"S","Cli.","P", expression(italic(C)*"+S"),
                             "S+Cli.", "Cli.+P","random","residuals"),
                    values=c("#F3A332", "violetred1", "royalblue", "greenyellow", "forestgreen", "purple","lavender","#FFFFCC","gray"))+
  
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

p22=ggplot(data=effect_with_plant%>%filter(guild=="all"),aes(x=estimate,y=1:14))+
  geom_point(size=4,pch=21,color="black",fill=rep(c("#DE582E","#F3A332","#1868b2","#018a67"),times=c(1,5,7,1)))+
  geom_errorbar(aes(xmin = estimate-1.96*sd, xmax = 1.96*sd+estimate,width=0.2),
                color=rep(c("#DE582E","#F3A332","#1868b2","#018a67"),times=c(1,5,7,1)))+
  geom_vline(xintercept = 0,linetype="dashed")+
  scale_y_continuous (breaks=1:14,labels = c(expression(italic(C)),"pH",  "Moisture","SoilC", "CEC", "Sand" , "MAT" ,   
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
  xlim(-0.5,0.5)+
  xlab("Effect size ± 95%CI")+
  ggtitle("Climate+Soil+Plant(N=416)")+
  annotate("text",x=0.202583,y=11,label="*",size=8)+
  
annotate("text",x=-0.1583,y=3,label="*",size=8)+
annotate("text",x=0.1767583,y=14,label="***",size=8)+
annotate("text",x=-0.3202850,y=1,label="***",size=8)


  




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

