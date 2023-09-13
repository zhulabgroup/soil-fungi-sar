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
library(neon)

library(neonUtilities)

root.data <- loadByProduct(dpID="DP1.10067.001")
root.data[[1]]

root.data[[2]]# data contains root mass and root chemical traits?(what is id somDryMass)

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


## for the megapit data
root.mega.data <- loadByProduct(dpID="DP1.10066.001")

root.mega.data[[1]]
root.mega.data[[2]]
root.mega.data[[3]]
root.mega.data[[4]]   
root.mega.data=root.mega.data[[5]] # total root biomass along the profile
root.mega.data[[6]]   
root.mega.data[[7]]   
root.mega.data[[8]]   

 abby=subset(root.mega.data,siteID=="ABBY")[,c(12,14:19)]
subset(abby,pitProfileID=="ABBY.1"&sizeCategory=="<=4mm") # the first profile was sampled to a maximum depth of 180 mm  
root.mega.data.full <- loadByProduct(dpID="DP1.10066.001")
y   
root.mega.data.full[[1]]

root.mega.data.full[[2]]#issues and corrections
class(root.mega.data.full[[2]])   
 mega.root.dep= root.mega.data.full[[3]] 
colnames(root.mega.data.full[[3]] )# more detail classification of root traits data
 mega.root.dep= root.mega.data.full[[3]] [,c("siteID" ,"plotID" ,"collectDate" , "cnSampleID" , "cnSampleCode"  , "sampleType" , "d15N",   "d13C" ,  "nitrogenPercent" ,   
"carbonPercent"    ,   "CNratio")]
head(mega.root.dep)# each site has only one megapit with three profiles   
   subset()
 mega.root.dep=k
 
 #selecting the living roots 
 d=subset(mega.root.dep,live=="LIVE"&siteID=="BART"&file==1&thick=="<2MM")
 
 dep=c("0","10","20","30","40","50","60","70","80","90","100","110","120","130","140","150","160","170","180","190","200")
 # check the megapit data for each site
 
 