library(ggcorrplot)
library(lme4)
library(lmerTest)
library(tidyverse)

# sensitive analysis for the climate and soil model
# how the results change when differ number of samples were included in the analysis
## for the z values
d=list()
for (i in 1:515){
  
  d1=subset(all_z_ran1000,a1==a[i])
  d[[i]]=t(d1[,3:1002])
}

z_sen=d[[1]]
for(i in 2:515){
  z_sen=rbind(z_sen,d[[i]])
}
## for the c values
d=list()
for (i in 1:515){
  
  d1=subset(all_c_ran1000,a1==a[i])
  d[[i]]=t(d1[,3:1002])
}

c_sen=d[[1]]
for(i in 2:515){
  c_sen=rbind(c_sen,d[[i]])
}

# sensitive analysis data
plotID=rep(a,each=1000)
data_sen=data.frame(cbind(plotID,c_sen,z_sen))
data_sen$V2=as.numeric(data_sen$V2)
data_sen$X1=as.numeric(data_sen$X1)
names(data_sen)=c("plotID","logc","z")
# to merge the data with the environmental variables
data_sen_mod=merge(data_sen,model_data[,c(2:3,6:28)],by="plotID")
data_sen_mod=subset(data_sen_mod,siteIDD != "GUAN"&z<10)# excluded some plots that have less than 3 cores
# select the data
range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

data_sen_mod1=data_sen_mod[,c(1,2,3,5:17,19,27)]# select only climate,fungal diversity and soil variables 

data_sen_mod1[,c(2,4:18)]=apply(data_sen_mod1[,c(2,4:18)],2,range01)%>% data.frame()

#. construct the model

mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = data_sen_mod1)


# 

write.csv(data_sen_mod1,"data_sen_mod1.csv")

pv=list()
ef=list()
for (j in 20:100)
  {
  cat('\r',paste(paste0(rep("*", round(j/ 1, 0)), collapse = ''), j, collapse = ''))# informs the processing
  a=unique(data_sen_mod1$plotID)
  ddf=list()
  for (i in 1:length(a))
  {
    a1=subset(data_sen_mod1,plotID==a[i])
    
    ddf[[i]]=sample_n(a1,j,replace = FALSE)#we selected 20 rows for one site 
  }
  
  ddf1=ddf[[1]]
  for(i in 2:length(a)){
    ddf1=rbind(ddf1,ddf[[i]])
  }
  
  mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = ddf1)
  op=summary(mod)
  pv[[j]]=op$coefficients[1:15,5]# the pvalue for each selected variables
  ef[[j]]=op$coefficients[1:15,1]# the pvalue for each selected variables
}
  }

pvv=pv[[20]]
for(j in 21:30){
  pvv=rbind(pvv,pv[[j]])
}
pvv=data.frame(pvv)
q=tidyr::gather(pvv)
boxplot(value~key,data=q)# to see the 

## for the effect size,

