# get the variation of the core-level variation of soil variables
library(ggcorrplot)
library(phyloseq)


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


# we do not distinguish between the O and M horizons

soil_mean=matrix(nrow=length(a1),ncol=6)
soil_sd=matrix(nrow=length(a1),ncol=6)
for (i in 1:length(a1))
     {
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
data_sub <- subset_samples(d, plotIDM==a1[i])
data_sub =sample_data(data_sub)
# the mean and the variation of the soil variables
soil_var=data_sub[,c("soilInCaClpH","nitrogenPercent","organicCPercent","soilMoisture","Perc.Soil.Moisture","CNRatio")]%>%data.frame()
soil_var$Perc.Soil.Moisture=as.numeric(soil_var$Perc.Soil.Moisture)

soil_var$CNRatio=as.numeric(soil_var$CNRatio)

soil_var_mean=apply(soil_var,2,FUN=mean,na.rm=TRUE)
soil_var_sd=apply(soil_var,2,FUN=sd,na.rm=TRUE)

soil_mean[i,]=soil_var_mean
soil_sd[i,]=soil_var_sd
}
# the ph and soil variables were more completed with more data available
# the unite of the soil variables
# the 
soil_mean=soil_mean[,-c(5,6)]
soil_sd=soil_sd[,-c(5,6)]

soil_mean=data.frame(soil_mean)
soil_sd=data.frame(soil_sd)

names(soil_mean)=c("plotID","ph","organicCPercent","nitrogen","soilmoisture")
names(soil_sd)=c("phsd","organicCPercentsd","nitrogensd","soilmoisturesd")

soil_mean=cbind(a1,soil_mean)
names(soil_mean)[1]="plotID"

model_var_unique=unique(model_var)
soil_sd=cbind(a1,soil_sd)

names(soil_sd)[1]="plotID"
# ph and soil moisture have the highest completeness

# add the sd of the soil variables to the model variables

model_var_unique=merge(model_var_unique,soil_sd[,c(1,2,5)],by="plotID",all.x = TRUE)
# to see the raltionthip betwen the observed and the rastered

soil_com=merge(model_var_unique,soil_mean,by="plotID")

## combind with the z data
zc_nest_mean_var=merge(plot_level_zc_nest_mean,model_var_unique,by="plotID")

zc_nest_mean_var$meanz=as.numeric(zc_nest_mean_var$meanz)

zc_nest_mean_var$meanc=as.numeric(zc_nest_mean_var$meanc)

# get the saved climate and soil data

plot_loca_all_soil_climate_mean=plot_loca_all_soil_climate_mean[,-1]

# based on the 30 permutations to get the 
env=list()
a2=unique(model_data$plotID)

for (i in 1:length(a2))
  {
  a3=subset(model_data,plotID==a2[i])[,c(1,2,5:dim(model_data)[2])]
  
  env[[i]]=unique(a3)
}

kk=env[[1]]
for (i in 2:515)
  {
 kk=rbind(kk,env[[i]]) 
}
# all the variables for the 515 plots

zc_nest_mean_var=merge(plot_level_zc_nest_mean,kk,by="plotID")

# the coordinats of the observed 

kk=sample_data(d)[,c("plotIDM", "lon","lat")]%>%data.frame()

names(kk)[1]="plotID"

kk=aggregate(kk[,2:3],by=list(kk$plotID),FUN=mean)
names(kk)[1]="plotID"

obs=merge(kk,zc_nest_mean_var[,c(1,2,5,6)],by="plotID")


# standardized the data

1. # climate and soil model:both dob and neon sites were included, here only the plotID was treated as a random effect
mode.data1 <- zc_nest_mean_var[, c(2:4,5:19)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)

mode.data1 <- subset(mode.data1, siteIDD != "GUAN"&meanz>0)

# check colinearity among variables



# mapping the correlation

ggcorrplot(cor(mode.data1[, c(1,2,4:16)]), hc.order = TRUE, type = "lower", lab = TRUE)#

mode.data1$meanz=as.numeric(mode.data1$meanz)
mode.data1$meanc=as.numeric(mode.data1$meanc)

# standardize the data

mode.data1[, c(1,2,4:16)]=apply(mode.data1[, c(1,2,4:16)],2,range01)

# bold was related with soil C, and MAT;cec was related with soil C; I decided to exclude cec and bold


range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}
# for the linear regression

mod=lm(meanz ~ organicCPercent+  ph  + nitrogen  + sand + bio1 +     bio2  +    bio4+bio8 +    bio12  +   bio15  +   bio18,data=mode.data1)

# for the lasson regression model

set.seed(66771)
train_set <- sample(c(TRUE, FALSE), size=dim(mode.data1)[1],
                    prob=c(0.70, 0.30), replace=TRUE)

test_set=!train_set

train_set[!train_set]=NA
mm=train_set

which(!is.na(mm))

train_data=mode.data1[which(!is.na(mm)),]

## for the test data
test_set[!test_set]=NA
mmt=test_set

which(!is.na(mmt))
test_data=mode.data1[which(!is.na(mmt)),]


set.seed(56885)
dd=matrix(ncol=2,nrow=1000)
for (i in 1:1000){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  mod_cv <- cv.glmnet( as.matrix(train_data[,c(4,5,7,9:16)]), train_data[,"meanz"], family='gaussian', intercept = F, alpha=1)
  
  dd[i,1]=mod_cv$lambda.min
  dd[i,2]=mod_cv$lambda.1se
}
dd=data.frame(dd)

head(dd)


lam=unique(dd$X1)
R=numeric()
term=list()
for (i in 1:length(lam)){
  
  la.eq <- glmnet( train_data[,c(4,5,7,9:16)], train_data[,"meanz"],lambda=lam[i], family='poisson', intercept = F, alpha=1)
  matrix(coef(la.eq ))[2:12,]
  terms=matrix(nrow=11,ncol=2)
  terms[,1]=colnames(train_data[,c(4,5,7,9:16)])
  terms[,2]=coef(la.eq )[2:12,]
  terms=data.frame(terms)
  best_model <- glmnet( train_data[,c(subset(terms,X2!=0)[,"X1"])], train_data[,"meanz"], 
                        family='gaussian', intercept = F, alpha=1,lambda = lam[i],final.re=TRUE)
  pre_zvalue=predict(best_model,lambda = lam[i],newx=as.matrix(test_data[,c(subset(terms,X2!=0)[,"X1"])]))
  
  a4=summary(lm(pre_zvalue~test_data[,"meanz"]))
  R[i]=a4$adj.r.squared
  term[[i]]=terms
}

# the best 

## add the coordinates to the data

coordinate=sample_data(d)[,c("plotIDM","lon","lat")]%>%data.frame()

coordinate=aggregate(coordinate[,2:3],by=list(coordinate$plotIDM),FUN=mean)

names(coordinate)[1]="plotID"

obs=merge(obs,coordinate,by="plotID")
