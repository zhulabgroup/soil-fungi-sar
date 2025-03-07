library(phyloseq)
library(ggplot2)

load("~/soil-sar/plot-sar-permutation/rare_all.Rdata")
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

# to see the number of soil cores per site

core_number=numeric()# look at the number of plots per site
for (i in 1:length(a))
{
  a1=subset(d,Site==a[i])
  core_number[i]=dim(a1)[1]# the number of plots in one site
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
# creat a plot to show the distribution of the sites
#i=8

a1=subset(d,Site==a[i]) [,c("lon","lat")]

p_plot_dis=ggplot()+
  geom_point(data=a1,pch=22,size=10,fill="seagreen2",alpha=0.5,aes(x=lon,y=lat))+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=10), 
        axis.text.y = element_text(hjust = 1,size=10), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA))+
  xlab("Longitude")+
  ylab("Latitude")




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

# need to find an area that more samples would be included in the range

site_area%>%filter(dit1>10)
# three sites were non included in the data, including STEI,TOOL and YELL

site_area%>%mutate(area=dit1* dit2/0.3861019,core=core_number)%>%filter(area<120&core>119)->site_temp
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

richness%>%bind_cols(site0)%>%mutate(area=rep(1.2e+08,21))%>%rename_all(~paste0(c("mean_value","site","area")))->site_richness_sample_120


# if we use the extrapolation approach to estimate the richness in the site

site0=unique(site_temp$a)# focused on sites with 120 cores and then the area is about 120 km^2
set.seed(567)
  richness=numeric()
  for (i in 1:21)
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # 
    dk=subset_samples(rare_all,Site==site0[i])
    data_sub =transform_sample_counts(dk, function(x) ifelse(x>0, 1, 0))
    ot=t(otu_table(data_sub))%>%data.frame()
    n=dim(ot)[2]
 
    ot1=rowSums(ot)
    out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(20, 4000, length.out=600)), se=FALSE)
    richness[i]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
  } 
  
  
  richness%>%bind_cols(site0)%>%mutate(area=rep(1.2e+08,21))%>%rename_all(~paste0(c("mean_value","site","area")))->site_richness
  
# to merge with the plot-scale diversity data
  setwd("/Volumes/seas-zhukai/proj-soil-fungi/SAR-30m")
  
  full_dob_neon_richness_data=readRDS("full_dob_neon_richness_data.rds")

  full_dob_neon_richness_data%>%filter(guild=="all")%>%mutate(site=substr(plotid,1,4))%>%
  group_by(site,area)%>%summarise(mean_value=mean(mean_value,na.rm=TRUE))%>%
    dplyr::select(site,area,mean_value)%>%bind_rows(site_richness)->site_level_sar

  site_level_sar%>%filter(site=="KONZ")
  
  site_level_sar%>%filter(site%in%site_richness$site)%>%
    mutate(type=rep("landscape",n()))%>%
    bind_rows(site_level_sar%>%filter(area<2000)%>%mutate(type=rep("local",n())))->local_landscape_richness
  
  site_id=unique(site_richness$site)
  
  pp=list()
  for (i in 1:21){
  data=local_landscape_richness%>%filter(site==site_id[i])
 
   pp[[i]]=ggplot(data=data,aes(x=log(area),y=log(mean_value),color=type))+
     geom_point(size=3)+
     geom_smooth(data=data,aes(x=log(area),y=log(mean_value),color=type),method="lm")+
     theme(legend.position =c(0.7,0.3),
           text = element_text(size = 13),
           plot.title = element_text(size = 12, hjust = 0.5), 
           axis.text.x = element_text(hjust = 1,size=8), 
           axis.title.y = element_text(size = 15), 
           axis.title.x = element_text(size = 15), 
           axis.ticks.x = element_line(color="black"), 
           panel.grid = element_blank(),
           panel.background = element_rect(fill = "NA"),
           panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
     xlab(expression(ln(Area) *" m"^2))+
     ylab("ln(Richness)")+
     ylim(0,11)+
     #ggtitle(paste0(site_id[i]))+
     scale_color_manual("",breaks=c("landscape","local"),labels=c("landscape","local"),values=c("#037f77","#7c1a97"))
  }
  
  plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],
            pp[[5]],pp[[6]],pp[[7]],pp[[8]],
            pp[[9]],pp[[10]],pp[[11]],
            pp[[12]],pp[[13]],pp[[14]],pp[[15]],
            pp[[16]],pp[[17]],pp[[18]],pp[[19]],pp[[20]],pp[[21]],ncol=5)
  
# to see the landscape z values
  zva=numeric()
  for (i in 1:21){
    data1=local_landscape_richness%>%filter(type=="landscape"&site==site_id[i])
    
    d=lm(log(mean_value)~log(area),data=data1)%>%summary()
    zva[i]=d$coefficients[2,1]
  }
    
  # to see how the two set of z values are related
  
 bind_cols (site_id,zva)%>%data.frame()%>%rename_all(~paste0(c("site","value_site")))%>%
   left_join(full_parameter_data%>%filter(guild=="all")%>%mutate(site=substr(plotID,1,4))
      %>%group_by(site)
      %>%summarise(mean_value=mean(zvalue,na.rm=TRUE),sd=sd(zvalue,na.rm=TRUE)),by="site")->compare_zvalue
 
 p_zvalue_corr=ggplot()+
   geom_point(data=compare_zvalue,aes(x=value_site,y= mean_value),size=3)+
   geom_segment(data=compare_zvalue, aes(x=value_site,xend=value_site,y= mean_value-sd,yend=mean_value+sd))+
   theme(legend.position =c(0.7,0.3),
         text = element_text(size = 13),
         plot.title = element_text(size = 12, hjust = 0.5), 
         axis.text.x = element_text(hjust = 1), 
         axis.title.y = element_text(size = 15), 
         axis.title.x = element_text(size = 15), 
         axis.ticks.x = element_line(color="black"), 
         panel.grid = element_blank(),
         panel.background = element_rect(fill = "NA"),
         panel.border = element_rect(color = "black", size = 0.8, fill = NA))+
   geom_smooth(data=compare_zvalue,aes(x=value_site,y= mean_value),method="lm")+
   xlab(expression("Landscape-level "*italic(z)*" value"))+
   ylab(expression("Plot-level "*italic(z)*" value"))+
   annotate("text",x=0.25,y=0.9,size=5,label=expression(italic(R^2)*"=0.57***"))
   
   
 # to combine the three plots
 p_plot_dis=ggplotGrob(p_plot_dis)
 pp[[4]]=ggplotGrob(pp[[4]])
 p_zvalue_corr=ggplotGrob(p_zvalue_corr)
 
 p_plot_dis$heights=pp[[4]]$heights
 p_plot_dis$heights=p_zvalue_corr$heights
 
 plot_grid(p_plot_dis,pp[[4]],p_zvalue_corr,ncol=3,labels = c("(a)","(b)","(c)"))
 
 
  
# assume that we are taking samples for a 10 by 10 mile (16.09344*16.09344 km)square for the 42 sites

# in that case, we need to sample  2589988 cores in total to make the sampling effort comparable, which is unrealistically too much 

# to see how much cores are there per site

d <- sample_data(a1) #
d=data.frame(d)
plotID=substr(d$geneticSampleID,6,8)
d=cbind(d,plotID)
a=unique(d$Site)



# if we look at the diversity at a finer resolution of 0.05 by 0.05 degree, which corresponds to an area of about 30 km^2
# we have 24 sites

site_area%>%mutate(area=dit1* dit2/0.3861019,core=core_number)%>%filter(area<30)->fine_data

a3=unique(fine_data$a)
# get the richness in these sites
site_rich=data.frame(nrow=length(a3),ncol=3)
for (i in 1:length(a3))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, Site==a3[i])
  data_sub =transform_sample_counts(data_sub, function(x) ifelse(x>0, 1, 0))
  ot=t(otu_table(data_sub))%>%data.frame()
  
  n=dim(ot)[2]# the number of "sites" for each plot
  
  if(n>2)
  {
    ot1=rowSums(ot)
    out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(10, 600, length.out=10)), se=FALSE)
    
    site_rich[i,1]=a3[i]
    site_rich[i,2]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    site_rich[i,3]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    
  }
  
  else{
    
    site_rich[i,1]=a3[i]
    site_rich[i,2]=NA
    site_rich[i,3]=NA
  }
  
}

# to get the richness data for the finer spatial scales
load("~/soil-sar/plot-sar-permutation/standardize sar/data_full_neon.RData")


data_full_neon%>%mutate(site=substr(plotid,1,4))%>%group_by(site,area)%>%summarise(rich=mean(mean_value,na.rm=TRUE))%>%

bind_rows(site_rich%>%mutate(area=20000000)%>%select(nrow,area,ncol)%>%rename_all(~paste0(c("site","area","rich"))))->sar_grid



sar_grid%>%filter(!is.na(rich))%>%group_by(site)%>%mutate(number=n())%>%filter(number>4)->sar_grid_five_points

a5=sar_grid_five_points%>%distinct(site)%>%pull()

pp <- vector('list', length(a5))
for (i in 1:24){
sar_grid%>%filter(site==a5[i])->op

pp[[i]]=ggplot()+
geom_point(data=op%>%filter(area<2000),aes(y=log(rich),x=log(area)),size=3)+
  geom_smooth(data=op%>%filter(area<2000),aes(y=log(rich),x=log(area)),method="lm",color="seagreen1",fill="seagreen1")+
  geom_point(data=op%>%filter(area>=2000),aes(y=log(rich),x=log(area)),size=3)+
  geom_smooth(data=op,aes(y=log(rich),x=log(area)),method="lm",linetype="dashed",color="purple",fill="purple")+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=8), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
  xlab(expression(Area *" m"^2))+
  ylab("Richness")+
  ylim(0,11)
}

plot_grid(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],ncol=3)

## create a plot to show the site level area
pp <- vector('list', length(a5))
for (i in 1:19){
sub=subset_samples(a1,Site==a5[i])

sample_data(sub)%>%data.frame()%>%select(lon,lat,plotIDM)%>%group_by(plotIDM)%>%
  summarise(mean_lon=mean(lon),mean_lat=mean(lat))->temp_plot
pp[[i]]=ggplot()+
  geom_point(data=temp_plot,pch=22,size=10,aes(x=mean_lon,y=mean_lat))
}


#i=1
p1=pp[[1]]

p2=ggplot()+
  geom_point(data=temp_plot,pch=22,size=10,fill="seagreen2",alpha=0.5,aes(x=mean_lon,y=mean_lat))+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=10), 
        axis.text.y = element_text(hjust = 1,size=10), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA))+
  xlab("Longitude")+
  ylab("Latitude")


## to see the slope between the two data sets
zval=numeric()
for (i in 1:length(a5))
{
  mod=lm(log(rich)~log(area),data=sar_grid_five_points%>%filter(site==a5[i]))
  zval[i]=coef(mod)[2]
}
# with a mean of 0.43
zval=numeric()
for (i in 1:length(a5))
{
  mod=lm(log(rich)~log(area),data=sar_grid_five_points%>%filter(site==a5[i]&area<2000))
  zval[i]=coef(mod)[2]
}

# with a mean of 0.69

## to show some plots for the NEON sites

load("~/soil-sar/plot-sar-permutation/standardize sar/data_full_neon.RData")

data_full_neon%>%filter(!is.na(mean_value)&guild=="all")%>%group_by(plotid)%>%summarise(n=n())%>%filter(n>2)->n_neon

plotid=unique(n_neon$plotid)

data_full_neon%>%filter(plotid==plotid[i]&guild=="all")->df_temp

data_full_neon%>%filter(guild=="all"&!is.na(mean_value))%>%filter(plotid%in%n_neon$plotid)->df_temp

# to show several 


ggplot()+
geom_point(data=df_temp,aes(x=log(area),y=log(mean_value),color=plotid))+
 guides(color="none")+
  geom_smooth(data=df_temp,size=0.15,aes(x=log(area),y=log(mean_value),color=plotid),method="lm",se=FALSE)+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=12), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
  xlab(expression(ln(Area) *" m"^2))+
  ylab("log(Richness)")+
  ylim(0,8)+
  scale_color_manual("",breaks=unique(df_temp$plotid),labels=unique(df_temp$plotid),values=rep("black",278))
  
data_full_neon%>%filter(guild=="soilsap"&!is.na(mean_value))%>%filter(plotid%in%n_neon$plotid)->df_temp_soilsap

p2=ggplot()+
  geom_point(data=df_temp_soilsap,aes(x=log(area),y=log(mean_value),color=plotid))+
  guides(color="none")+
  geom_smooth(data=df_temp_soilsap,size=0.15,aes(x=log(area),y=log(mean_value),color=plotid),method="lm",se=FALSE)+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=12), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
  xlab(expression(ln(Area) *" m"^2))+
  ylab("log(Richness)")+
  ylim(0,8)+
  scale_color_manual("",breaks=unique(df_temp_soilsap$plotid),labels=unique(df_temp_soilsap$plotid),values=rep("black",278))



# create a plot with sd for the neon sites
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar.RData")

richness_subplot40_neon_standar$mean_value=as.numeric(richness_subplot40_neon_standar$mean_value)
richness_subplot40_neon_standar$sd_value=as.numeric(richness_subplot40_neon_standar$sd_value)
richness_subplot40_neon_standar$area=as.numeric(richness_subplot40_neon_standar$area)


data_with_sd=bind_rows(richness_subplot10_neon_standar,richness_subplot20_neon_standar,richness_subplot30_neon_standar,richness_subplot40_neon_standar)

data_with_sd%>%filter(!is.na(mean_value))%>%group_by(plotid)%>%summarise(n=n())%>%filter(n>2)->n_plot

plot_ID=n_plot%>%pull(plotid)

ggplot(data=data_with_sd%>%filter(!is.na(mean_value))%>%filter(plotid==plot_ID[5]),aes(x=log(area),y=log(mean_value)))+
geom_point()+
geom_errorbar(data=data_with_sd%>%filter(!is.na(mean_value))%>%filter(plotid==plot_ID[5]),aes(ymin =log(mean_value)-log(sd_value), ymax =log(mean_value)+log(sd_value),width=0.2))
                
data_with_sd%>%filter(!is.na(sd_value ))%>%group_by(plotid)%>%summarise(n=n())%>%filter(n>2)->data_with_sd_select

 


data_with_sd%>%filter(plotid%in%data_with_sd_select$plotid)->data_example 

a6=unique(data_example$plotid)

pp <- vector('list', length(a6))
for (i in 1:length(a6)){
  pp[[i]]=ggplot(data=data_example%>%filter(plotid==a6[i]),aes(x=area,y=mean_value,color=plotid))+
    geom_point(size=3,alpha=0.5)+
    geom_errorbar(aes(ymin =mean_value-sd_value, ymax =mean_value+sd_value,width=0.2))+
    theme(legend.position =c(0.83,0.26),
          text = element_text(size = 13),
          plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.x = element_text(hjust = 1,size=12), 
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          axis.ticks.x = element_line(color="black"), 
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
    xlab(expression(Area*" m"^2))+
    ylab("Richness")+
    guides(color="none")
}

  



ggplot(data=data_example,aes(x=log(area),y=log(mean_value)))+
  geom_point(size=3)+
  theme(legend.position =c(0.83,0.26),
        text = element_text(size = 13),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.x = element_text(hjust = 1,size=12), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_line(color="black"), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 0.8, fill = NA)) +
  xlab(expression(ln(Area)*" m"^2))+
  ylab("ln(Richness)")+
  geom_smooth(method="lm")



  p3=ggplotGrob(p3)
  p4=ggplotGrob(p4)
  
  p3$widths=p4$widths

