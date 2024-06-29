
#1. ## load in the new climate data 

climate=read.csv("climate.csv",sep=",",header=T)

save(rare_all,file="rare_all.Rdata")# save the rarefied data

d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site correspondes to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(rare_all, plotIDM)# merge the new plotid with the initial data 
# select an unique plot and build a SAR within the plot
a1= sample_data(d)# the unique plotID, we have 476 plots
a1=unique(a1$plotIDM)

location=sample_data(a1)%>%data.frame()%>%dplyr::select(plotIDM)


location%>%mutate(id=rownames(location))->location

location=location%>%rename(ID1=id)

climate=left_join(location%>%dplyr::select(ID1,plotIDM),climate,by="ID1")%>%dplyr::select(-ID2,-Latitude,-Longitude)%>%group_by(plotIDM)%>%summarise(across(c(MAT, MWMT,  MCMT,   TD, MAP, MSP,  AHM , SHM, DD_0,  DD5,DD_18, DD18, NFFD, bFFP, eFFP, FFP, PAS,   EMT,  EXT,  MAR, Eref, CMD, RH,   CMI, DD1040), mean, na.rm = TRUE))

save(climate,file="climate.RData")

#2. to see the variables with high correlations


ggcorrplot(cor(climate%>%dplyr::select(MAT, MWMT,  MCMT,   TD, MAP, MSP,  AHM , SHM, DD_0,  DD5,DD_18, DD18, NFFD, bFFP, eFFP, FFP, PAS,   EMT,  EXT,  MAR, Eref, CMD, RH,   CMI, DD1040)), hc.order = TRUE, type = "lower", lab = TRUE)#


df=climate%>%dplyr::select(MAT, MWMT,  MCMT,   TD, MAP, MSP,  AHM , SHM, DD_0,  DD5,DD_18, DD18, NFFD, bFFP, eFFP, FFP, PAS,   EMT,  EXT,  MAR, Eref, CMD, RH,   CMI, DD1040)%>%filter(MAT!=-9999.0)

cor_matrix <- cor(df, use = "pairwise.complete.obs")

# Convert the correlation matrix to a long format
cor_long <- as.data.frame(as.table(cor_matrix))
# Rename columns for clarity
names(cor_long) <- c("Var1", "Var2", "Correlation")

# Filter for pairs with high correlation (e.g., |correlation| > 0.8) and avoid self-correlation
high_cor_pairs <- cor_long %>%
  filter(abs(Correlation) > 0.75 & Var1 != Var2) %>%
  arrange(desc(abs(Correlation)))

df%>%dplyr::select(-DD_18,-Eref,-MCMT,-EMT,-DD5,-eFFP,-MWMT, -NFFD,-FFP, -DD1040,-bFFP,-DD_0,-DD18,-EXT,-CMI,-CMD)->dd

# the variables selected in the new climate data
ggcorrplot(cor(dd%>%dplyr::select(  MAT,    TD,   MAP ,  MSP  ,  SHM ,  PAS  , MAR ,AHM  ,  RH)), hc.order = TRUE, type = "lower", lab = TRUE)#


## to merge the climate data with the z value

climate%>%dplyr::select(plotIDM,MAT,    TD,   MAP ,  MSP  ,  SHM ,  PAS  , MAR ,AHM  ,  RH)%>%filter(MAT!=-9999.0)%>%rename(plotID=plotIDM)->temp_climate


model_data_SAR_NAclimate=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand)
model_data_SAR_NAclimate%>%left_join(temp_climate,by="plotID")->model_data_SAR_NAclimate_soil
model_data_SAR_NAclimate_soil=model_data_SAR_NAclimate_soil[complete.cases(model_data_SAR_NAclimate_soil),]

## get the different data

data=list()
for (i in 1:9)
{
  d=model_data_SAR_NAclimate_soil%>%filter(guild==guild_select[i]) 
  
  d[,c(6:22)]=apply(d[,c(6:22)],2,range01)
  
  data[[i]]=d
}

## to see colinearity of the variables

ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon



sel_vab_step=list()
for (i in 1:9)
  
{
    mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +MAT+   TD+ MAP+ MSP+  SHM+ PAS + MAR+  AHM +RH + (1 |siteIDD), data = data[[i]])
    mod_sel=step(mod)
    kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
}

#5. when plant richness was included in the model

model_data_SAR_NAclimate_richness=model_data_SAR%>%dplyr::select(siteIDD,plotID,guild,lon, lat,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,richness)
model_data_SAR_NAclimate_richness%>%left_join(temp_climate,by="plotID")->model_data_SAR_NAclimate_richness
model_data_SAR_NAclimate_richness=model_data_SAR_NAclimate_richness[complete.cases(model_data_SAR_NAclimate_richness),]

## get the different data

data=list()
for (i in 1:9)
{
  d=model_data_SAR_NAclimate_richness%>%filter(guild==guild_select[i]) 
  
  d[,c(6:22)]=apply(d[,c(6:22)],2,range01)
  
  data[[i]]=d
}

## to see colinearity of the variables

ggcorrplot(cor(data[[1]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[2]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[3]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon

ggcorrplot(cor(data[[4]]%>%dplyr::select(logc,soilInCaClpH, nitrogenPercent,organicCPercent, soilMoisture, cec,sand,MAT,   TD, MAP,MSP,  SHM, PAS,  MAR,  AHM ,RH )), hc.order = TRUE, type = "lower", lab = TRUE)#
#nitrogen and soil carbon



sel_vab_step=list()
for (i in 1:9)
  
{
  mod <- lmer(zvalue ~ logc+soilInCaClpH +organicCPercent +cec+ sand +MAT+  TD+ MAP+ MSP+  SHM+ PAS + MAR+  AHM +RH + richness+(1 |siteIDD), data = data[[i]])
  mod_sel=step(mod)
  kk=mod_sel$fixed%>%data.frame()%>%filter(!Pr..F.>0.05)
  sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
}






