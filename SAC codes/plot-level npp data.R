# extract the plot-level NPP data for all plots
#https://github.com/claraqin/fungal-climate-niche/blob/master/R/1-data.R

# the function for standardization 

range01 <- function(x)
{
  return((x - min(x)) / (max(x) - min(x)))
}


library(neonUtilities)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(patchwork)

############

## extract the the NPP data based on the coordinates of the plots
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
# get the plotid and the coordinates of each plot
a1= sample_data(d)[,c("lon","lat","plotIDM")]
a1=data.frame(a1)

# get the mean value of npp across 19 years from 2001-2019

pa=paste("modis-250-npp-",2001:2019,".tif")# creat a name for each tif object
pa=gsub(" ","",pa)
a=list()
for(i in 1:19)
{
  str_name=pa[i]
  b<- raster(str_name)
  a[[i]]=raster::extract(b[["annualNPP"]], a1[, 1:2])
}

anpp=matrix(ncol = 19,nrow=6378)
for (i in 1:19)
{
  anpp[,i]=a[[i]]
}

# to check whetehr some plots have none npp data for all years
names(anpp)[1]="npp1"
anpp=anpp[,-c(23:25)]

anpp_na=anpp %>% filter(is.na(npp1))# yes! these sites do not have npp data for all years

# determine the mean of the 19 years

anpp=cbind(anpp,a1)

anpp_mean=apply(anpp[,1:19],1,mean)

anpp_mean=cbind(anpp_mean,a1["plotIDM"])
anpp_mean=aggregate(anpp_mean~plotIDM,data=anpp_mean,FUN=mean)# leads to 374 plots with npp data

write.csv(anpp_mean,"anpp_mean.csv")

# bind the npp data with the initial model data

anpp_mean=anpp_mean[,-1]
names(anpp_mean)=c("plotID","npp")

d=merge(model_data[,2:28],anpp_mean,by="plotID",all.x=TRUE)

write.csv(d,"model_data.csv")






