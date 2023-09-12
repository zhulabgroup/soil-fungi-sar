# fungi
hope things go well

this is very helpful for me but

can you help me with that.

d=c(1:4)
a=rep(1:3,each=5)
setwd()
getwd()

setwd("/Users/luowenqi/desktop")

a=b+c
summary(mod_)
# the combined data 
head(com.neon.dob)
# get the root data,11 tables in total

root.data <- loadByProduct(dpID="DP1.10067.001")
root.data[[1]]

root.data[[2]]# data contains root mass and root chemical traits?(what id somDryMass)

root.data[[3]]# about sampling location and time

root.data[[4]]# root chemical traits



root.data[[5]]#dry mass with different root orders

root.data[[6]]# site description?

root.data[[7]]#describtion
root.data[[8]]# site description?
root.data[[9]]# site description?
root.data[[10]]# site description?


root.chemi=root.data[[4]]

head(root.chemi)

colnames(root.chemi)

root.chemi=root.chemi[,c(5,8,15:19)]

class(root.chemi)
head(root.chemi)
table(root.chemi$siteID)
# to see the samples based on the years
year=substr(root.chemi$collectDate,1,4)
root.chemi=cbind(root.chemi,year)# seems that each year, there are many samples for a site, should we get a mean for that?

head(root.chemi)# any NAS?

is.na(root.chemi)
aggregate(root.chemi[,c(3:7)],by=list(root.chemi$siteID),FUN=mean,na.omit=TRUE)# a number of sites have NAS

d=subset(root.chemi,d15N!="NA"&d13C!="NA"&nitrogenPercent!="NA"&carbonPercent!="NA"&CNratio!="NA")

root.mean=aggregate(d[,c(3:7)],by=list(d$siteID),FUN=mean,na.omit=TRUE)# a number of sites have NAS
head(root.mean)
cor(root.mean[,2:6])

dim(root.mean)

names(root.mean)[1]="site"
neon.root=merge(com.neon.dob,root.mean,by="site")

head(neon.root)

dim(neon.root)

cor(neon.root[,2:15])

mod=lm(ranz30~ soilInCaClpH+ nitrogenPercent.x+ organicCPercent+ soilMoisture+ mat_celsius+  map_mm+ temp_seasonality+prec_seasonality +        d15N +     d13C+ nitrogenPercent.y +carbonPercent+  CNratio,data=neon.root)

step(mod)
 summary(lm(ranz30 ~ soilInCaClpH + soilMoisture + mat_celsius + 
   temp_seasonality + d13C, data = neon.root))# best model for the NEON sites