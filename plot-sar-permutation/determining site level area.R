library(phyloseq)
library(ggplot2)

d=sample_data(rare_all)
# select the NEON SITE
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

a1= subset_samples(d,Project=="NEON")# 

d <- sample_data(a1) #
d=data.frame(d)
plotID=substr(d$geneticSampleID,6,8)
d=cbind(d,plotID)
a=unique(d$Site)

site=numeric()# look at the number of plots per site
for (i in 1:length(a))
{
  a1=subset(d,Site==a[i])
  site[i]=length(unique(a1$plotID))# the number of plots in one site
}



pp <- vector('list', 45)# location of the plots within a site
for (i in seq_along(pp))
{
  a1=subset(d,Site==a[i]) [,c("lon","lat")]
  a1=unique(a1)
  pp[[i]]=ggplot()+geom_point(data=a1,aes(x=lon,y=lat),color="red",size=3)# the number of 
}

plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],
          pp[[10]],pp[[11]],pp[[12]],pp[[13]],pp[[14]],pp[[15]],pp[[16]],pp[[17]],pp[[18]],pp[[19]],
          ncol=9)


dtlon=numeric()# the lon and lat range of the plots within a site
for (i in 1:length(a))
{
  a1=subset(d,Site==a[i])[,c("lon","lat")]
  dtlon[i]=max(a1$lon)-min(a1$lon)
}

dtlat=numeric()
for (i in 1:length(a))
{
  a1=subset(d,Site==a[i])[,c("lon","lat")]
  dtlat[i]=max(a1$lat)-min(a1$lat)
}

site_area=data.frame(cbind(a,dtlon,dtlat))
site_area$dtlon=as.numeric(site_area$dtlon)
site_area$dtlat=as.numeric(site_area$dtlat)

#estimate an area that encompass all the plots in a site
# 1 degree of lon==54.6 miles
# 1 degree of lat=69 miles
#https://www.usgs.gov/faqs/how-much-distance-does-a-degree-minute-and-second-cover-your-maps#:~:text=One%2Ddegree%20of%20longitude%20equals,one%20second%20equals%2080%20feet.

dit1=site_area$dtlon*54.6
dit2=site_area$dtlat*69
site_area=cbind(site_area,dit1,dit2)#
# the site for the STEI site covers a quite large area, the maximum sampling area is 27.61029*22.67657(28*23 miles)

# need to find a area that more samples would be included in the range

site_area%>%filter(dit1>10)
# three sites were non included in the data, including STEI,TOOL and YELL

site_area%>%mutate(area=dit1* dit2/0.3861019,core=con_number)%>%filter(area<120&core>119)->site_temp
# we will have 21 sites

site0=unique(site_temp$a)# focused on sites with 120 cores and then the area is about 120 km^2
set.seed(567)
site_rich=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:21)
  {
    dk=subset_samples(rare_all,Site==site0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 120,replace = FALSE)
    sampled_physeq <- prune_samples(sampled_names, dk)
    otu=otu_table(sampled_physeq )
    richness[i]=table(apply(otu,2,sum)>0)[[2]]#
  } 
  site_rich[[j]]=richness
}








# assume that we are taking samples for a 10 by 10 mile (16.09344*16.09344 km)square for the 42 sites

# in that case, we need to sample  2589988 cores in toal to make the sampling effort comparable, which is unrealistic 

# to see how much cores are there per site

d <- sample_data(a1) #
d=data.frame(d)
plotID=substr(d$geneticSampleID,6,8)
d=cbind(d,plotID)
a=unique(d$Site)

con_number=numeric()# look at the number of plots per site
for (i in 1:length(a))
{
  a1=rare_all%>%sample_data()%>%data.frame()%>%filter(Site==a[i])
  
  con_number[i]=dim(a1)[1]# the number of plots in one site
}


# for each site we select a



