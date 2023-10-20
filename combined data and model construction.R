# combing data for modeling
plot_plant_rich#plant data#p
core_mass_type#root mass data
root.chemi.mean#root trait data
rich_plot#fungal richness
plot_loca_all_soil_climate_mean#soil variables
# mergeing these datasets generates an obhect of 'd5'
#computing the z and the log(c) value 
a=list()
for (i in 1:dim(all_z)[1])
  {
  a[[i]]=t(all_z[i,2:31])
}

b=a[[1]]
for(i in 2:dim(all_z)[1])
  {
  b=rbind(b,a[[i]])
}

plotID=rep(all_z$a1,each=30)

all_z_30=cbind(plotID,b)# for each plot,with 30 estimated z values
all_z_30=data.frame(all_z_30)
all_z_30$X1=as.numeric(all_z_30$X1)
# for the c values

a=list()
for (i in 1:dim(all_c)[1])
{
  a[[i]]=t(all_c[i,2:31])
}

b=a[[1]]
for(i in 2:dim(all_c)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(all_c$a1,each=30)
all_c_30=cbind(plotID,b)# for each plot,with 30 estimated c values, which can be used to estimated species per unit sample(estimated=log(c))

all_c_30=data.frame(all_c_30)
all_c_30$X1=as.numeric(all_c_30$X1)
# combind the z and c
com_z_c=cbind(all_z_30,all_c_30[,2])
names(com_z_c)=c("plotID","z","log(c)")
com_z_c=cbind(com_z_c,c=2.71828^com_z_c$`log(c)`)# the estimated c and z are negatively correlated
model_data=merge(com_z_c[,c(1,2,4)],d5,by="plotID",all.x = TRUE)
# adding the site 
siteIDD=substr(model_data$plotID,1,4)
model_data=cbind(siteIDD,model_data)

write.csv(model_data,"model_data.csv")
model_data=model_data[,-5]
model_data$rich=as.numeric(model_data$rich)

write.csv(model_data,"model_data.csv")

# climate and soil model:both dob and neon sites were included, here only the plotID was treated as a random effect
mode.data1=model_data[,c(1:17,19,27)]# GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data1=subset(mode.data1,siteIDD!="GUAN"&z<10)
# check co linearity among variables
# climate and soil model:both dob and neon sites were included, here only the plotID was treated as a random effect
mode.data1=model_data[,c(1:17,19,27)]# GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data1=subset(mode.data1,siteIDD!="GUAN"&z<10)
# check colinearity among variables

cor(mode.data1[,3:19])

# mapping the correlation
ggcorrplot(cor(mode.data1[,3:19]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil C, and MAT;cec was related with soil C; I decided to exclude cec and bold
# build model with plot included as the random effect
# standardized the data.

range01=function(x)## to
{
  return((x-min(x))/(max(x)-min(x)) )
}

mode.data1[,4:19]=apply(mode.data1[,4:19], 2, range01)%>%data.frame

# site and plot are nested
mod=lmer(z ~ c+organicCPercent  + ph+ nitrogen+   sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |siteIDD/plotID),data=mode.data1)

step(mod)
#best model
mod=lmer(z ~ c + nitrogen + sand + bio2 + bio12 + bio15 + funrich + (1 | siteIDD/plotID),data=mode.data1)
summary(mod)


# plotid as the only random effect
mod=lmer(z ~ c+organicCPercent  + ph+ nitrogen+   sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |plotID),data=mode.data1)
step(mod)
#best model
mod=lmer(z ~  c + sand + bio2 + bio4 + bio12 + bio15 + funrich + (1 | plotID),data=mode.data1)
summary(mod)
# site as the only random effect

mod=lmer(z ~ c+organicCPercent  + ph+ nitrogen+   sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |siteIDD),data=mode.data1)
step(mod)
#best model
mod=lmer(z ~ c + organicCPercent + nitrogen + sand + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD),data=mode.data1)
summary(mod)


## climate, plant diversity and soil model: only neon sites were included, dob sites do not have plant data
mode.data2=model_data[,c(1:19,27)]# GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data2=subset(mode.data2,siteIDD!="GUAN"&z<10&rich>0)# some sites do have plant data
# check colinearity among variables
cor(mode.data2[,3:20])
ggcorrplot(cor(mode.data2[,3:20]), hc.order = TRUE, type = "lower", lab = TRUE)#
mode.data2[,3:20]=apply(mode.data2[,3:20], 2, range01)%>%data.frame

# site and plot are nested,  c was included

mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |siteIDD/plotID),data=mode.data2)
step(mod)
#best model
mod=lmer(z ~ organicCPercent + c + sand + bio12 + bio15 + funrich + bio1 + (1 | siteIDD/plotID),data=mode.data2)
## only site as random
mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |siteIDD),data=mode.data2)

step(mod)
# only plot as random
mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +    spei+ funrich +bio1+ (1 |plotID),data=mode.data2)

# climate, soil, plant and root model# will have 104 plots

mode.data3=model_data[,c(1:27)]# GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data3=subset(mode.data3,siteIDD!="GUAN"&z<10&rich>0&fine>0&rootc>0)# some sites do have plant data
# check colinearity among variables
cor(mode.data3[,3:27])
#coarse and fine roots are correlated;bio1 and bio4, bio1 and bold; organic and bold; organic and cec; root n and root cn;
# i therefor excluded coarse, bio4,bold,cec and root n
ggcorrplot(cor(mode.data3[,3:27]), hc.order = TRUE, type = "lower", lab = TRUE)#
mode.data3[,3:27]=apply(mode.data3[,3:27], 2, range01)%>%data.frame
#site and plot are nested
mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  +bio12+ bio15  +    spei+ funrich +bio1+fine+d13C+rootc  +  rootcn + (1 |siteIDD/plotID),data=mode.data3)
step(mod)
# best model
mode=lmer(z ~ c + rich + funrich + bio1 + (1 | siteIDD/plotID),data=mode.data3)
# site as the only random effect
mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  +bio12+ bio15  +    spei+ funrich +bio1+fine+d13C+rootc  +  rootcn + (1 |siteIDD),data=mode.data3)
step(mod)
#best model
mod=lmer(z ~ c + rich + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD),data=mode.data3)

# plot as the random
# site as the only random effect
mod=lmer(z ~ organicCPercent + c+ ph+ nitrogen+ rich+  sand +bio2 +bio8+ bio18+  +bio12+ bio15  +    spei+ funrich +bio1+fine+d13C+rootc  +  rootcn + (1 |plotID),data=mode.data3)
step(mod)
#best model
mod=lmer(z ~ c + rich + sand + bio8 + bio15 + funrich + rootc + (1 | plotID),data=mode.data3)
#


