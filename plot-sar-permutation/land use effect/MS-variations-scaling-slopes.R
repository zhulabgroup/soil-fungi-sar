################################################################################
#compare the mean z values among fungal guilds

load("~/soil-sar/plot-sar-permutation/model_comp.RData")

model_comp%>%filter(guild=="all")%>%group_by(type)%>%dplyr::summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%group_by(type)%>%summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%filter(guild=="all")%>%group_by(type)%>%summarise(count = n())

model_comp%>%filter(guild!="all")->guild_mean

guild_mean%>%filter(!guild%in%c("para","epiphy","littersap", "woodsap" ))->guild_mean

# we updated the z value analysis by excluding the parasitic and epiphytic fungi
# based on the model to test significant difference

mod=lmer(zvalue~guild+(1 | plotID),data=guild_mean)
# to compare the means
library(emmeans)
emmeans_model <- emmeans(mod, ~ guild)
pairwise_comparisons <- contrast(emmeans_model, method = "pairwise")
summary(pairwise_comparisons, adjust = "tukey") # Adjust for multiple testing

k=aggregate(zvalue~guild,data=guild_mean,FUN=mean)
od <- k[order(k$zvalue, decreasing = TRUE),] 
guild_mean$guild <- factor(guild_mean$guild, levels = od$guild)


#################################################################################
#the impact of climate variables on the variability of the estimated z values

load("~/soil-sar/plot-sar-permutation/soil_mean.RData")

load("~/soil-sar/plot-sar-permutation/model_data.RData")

full_parameter_data=readRDS("full_parameter_data.rds")

# need to combine the complete climate and soil variables

full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%
  left_join(model_data%>%dplyr::select(plotID,siteIDD,bio1,bio2,bio4,bio8,bio12,bio15,bio18,richness)%>%distinct(),by="plotID")%>%
  left_join(soil_mean%>%dplyr::select(-lon,-lat)%>%distinct(),by="plotID")->model_data_SAR_rarefaction

model_data_SAR_rarefaction$logc=2.71828^model_data_SAR_rarefaction$logc

model_data_SAR_climate=model_data_SAR_rarefaction%>%dplyr::select(siteIDD,plotID,guild,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)
model_data_SAR_climate=model_data_SAR_climate[complete.cases(model_data_SAR_climate),]

guild_select=unique(model_data_SAR_climate$guild)

# to standardize the data

data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:18)]=apply(d[,c(4,6:18)],2,range01)
  
  data[[i]]=d
}

# the number of plots included in the data
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

# the effect size of the variables for each of the guilds
# do not do a model selection
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

##

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


effect_no_plant$var=factor(effect_no_plant$var,levels=c("logc","soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18"))        
effect_no_plant$guild=factor(effect_no_plant$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           
effect_no_plant$pva[effect_no_plant$pva<=0.001]="***"
effect_no_plant$pva[effect_no_plant$pva>0.001&effect_no_plant$pva<=0.01]="**"
effect_no_plant$pva[effect_no_plant$pva>0.01&effect_no_plant$pva<=0.05]="*"
effect_no_plant$pva[effect_no_plant$pva>0.05]=""



##################################################################################
## the variance partitioning analysis 
## only significant variables were selected for the variance partitioning. 

sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +soilMoisture +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +soilMoisture +organicCPercent +cec+ sand +bio1+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}


# grouping the selected variables as climate and soil effects
# the production of the variables were treated as a combined variable

soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")

# the out put is based on the below codes
# when soil and climate variable become constants, the model will drop these terms automatcally.
# whether or not multiplying 10 does not affect the results

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_variable=data[[i]][,intersect( sel_vab_step[[i]],soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
  climate_variable=data[[i]][,intersect( sel_vab_step[[i]],climate)]%>%data.frame()%>%rowwise() %>%
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

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),
            var_par[[3]]%>%data.frame(),
            var_par[[4]]%>%data.frame(),
            var_par[[5]]%>%data.frame(),
            var_par[[6]]%>%data.frame(),
            var_par[[7]]%>%data.frame(),
            var_par[[8]]%>%data.frame(),
            var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed# some do have the data
var_fixed%>%dplyr::select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_noplant

