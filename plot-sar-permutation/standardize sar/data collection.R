### determination of the z and c values for the full data 

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar.RData")

richness_subplot40_neon_standar$area=as.numeric(richness_subplot40_neon_standar$area)
richness_subplot40_neon_standar$mean_value=as.numeric(richness_subplot40_neon_standar$mean_value)
richness_subplot40_neon_standar$sd_value=as.numeric(richness_subplot40_neon_standar$sd_value)

## for the dob sites
setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot10_dob_standar.RData")
load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot20_dob_standar.RData")
load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot30_dob_standar.RData")
load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot40_dob_standar.RData")

richness_subplot30_dob_standar=richness_subplot30_dob_standar%>%rename_all(~paste0(c("plotid","mean_value","area")))
richness_subplot40_dob_standar=richness_subplot40_dob_standar%>%rename_all(~paste0(c("plotid","mean_value","area")))

richness_subplot40_dob_standar$mean_value=as.numeric(richness_subplot40_dob_standar$mean_value)
richness_subplot40_dob_standar$area=as.numeric(richness_subplot40_dob_standar$area)


 data_dob_all=bind_rows(richness_subplot10_dob_standar%>%select( plotid, mean_value,area),richness_subplot20_dob_standar%>%select( plotid, mean_value,area),richness_subplot30_dob_standar%>%select( plotid, mean_value,area),richness_subplot40_dob_standar%>%select( plotid, mean_value,area))


###for different guilds
# for the ECM fungi with two datasets

 load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot10_dob_standar_EM.RData")
 # to modify the data
 
 richness_subplot10_dob_standar_EM%>%mutate(area=rep(100,length.out=n))->richness_subplot10_dob_standar_EM
 richness_subplot20_dob_EM%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_EM
 richness_subplot30_dob_standar_EM%>%rename(plotid=plotID)->richness_subplot30_dob_standar_EM
 richness_subplot40_dob_standar_EM=richness_subplot40_dob_standar_EM%>%as.matrix()%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value","area")))
 richness_subplot40_dob_standar_EM$mean_value=as.numeric( richness_subplot40_dob_standar_EM$mean_value)
 richness_subplot40_dob_standar_EM$area=as.numeric( richness_subplot40_dob_standar_EM$area)
 

save(richness_subplot10_dob_standar_EM,file="richness_subplot10_dob_standar_EM.RData")
save(richness_subplot20_dob_standar_EM,file="richness_subplot20_dob_standar_EM.RData")
save(richness_subplot30_dob_standar_EM,file="richness_subplot30_dob_standar_EM.RData")
save(richness_subplot40_dob_standar_EM,file="richness_subplot40_dob_standar_EM.RData")
 
 
 load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot10_dob_standar_EM.RData")
 load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot20_dob_EM.RData")
 load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot30_dob_standar_EM.RData")
 load("~/soil-sar/plot-sar-permutation/standardize sar/richness_subplot40_dob_standar_EM.RData")

data_dob_EM=bind_rows(richness_subplot10_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_EM%>%select( plotid, mean_value,area))
 

## for the soil saprophytic fungi
 richness_subplot10_dob_soilsap%>%mutate(area=rep(100,length.out=n))-> richness_subplot10_dob_soilsap
 richness_subplot10_dob_standar_soilsap=richness_subplot10_dob_soilsap
 save(richness_subplot10_dob_standar_soilsap,file="richness_subplot10_dob_standar_soilsap.RData") 
 richness_subplot20_dob_soilsap%>%mutate(area=rep(400,length.out=n))-> richness_subplot20_dob_soilsap
 richness_subplot20_dob_standar_soilsap=richness_subplot20_dob_soilsap
 save(richness_subplot20_dob_soilsap,file="richness_subplot20_dob_soilsap.RData") 
 richness_subplot30_dob_standar_soilsap%>%rename(plotid=plotID)->richness_subplot30_dob_standar_soilsap
 save(richness_subplot30_dob_standar_soilsap,file="richness_subplot30_dob_standar_soilsap.RData")
 
 data_dob_soilsap=bind_rows(richness_subplot10_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_EM%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_EM%>%select( plotid, mean_value,area))

 ###for the wood saprotrophic
 
 richness_subplot10_dob_standar_woodsap=richness_subplot10_dob_standar_woodsap%>%mutate(area=rep(100,length.out=n))
 richness_subplot20_dob_standar_woodsap=richness_subplot20_dob_woodsap
 richness_subplot20_dob_standar_woodsap%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_woodsap
 save(richness_subplot20_dob_standar_woodsap,file="richness_subplot20_dob_standar_woodsap.RData")
 richness_subplot30_dob_standar_woodsap=richness_subplot30_dob_standar_woodsap%>%rename(plotid=plotID)
 save(richness_subplot30_dob_standar_woodsap,file="richness_subplot30_dob_standar_woodsap.RData")
data_dob_woodsap=bind_rows(richness_subplot10_dob_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_woodsap%>%select( plotid, mean_value,area))

## for the litter saprotroph
richness_subplot10_dob_standar_littersap=richness_subplot10_dob_standar_littersap%>%mutate(area=rep(100,length.out=n))
richness_subplot20_dob_standar_littersap=richness_subplot20_dob_littersap
richness_subplot20_dob_standar_littersap%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_littersap
save(richness_subplot20_dob_standar_littersap,file="richness_subplot20_dob_standar_littersap.RData")
richness_subplot30_dob_standar_littersap=richness_subplot30_dob_standar_littersap%>%rename(plotid=plotID)
save(richness_subplot30_dob_standar_littersap,file="richness_subplot30_dob_standar_littersap.RData")
data_dob_littersap=bind_rows(richness_subplot10_dob_standar_littersap%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_littersap%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_littersap%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_littersap%>%select( plotid, mean_value,area))

### for the plant pathogens
richness_subplot10_dob_standar_plapat=richness_subplot10_dob_standar_plapat%>%mutate(area=rep(100,length.out=n))
richness_subplot20_dob_standar_plapat=richness_subplot20_dob_plapat
richness_subplot20_dob_standar_plapat%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_plapat
save(richness_subplot20_dob_standar_plapat,file="richness_subplot20_dob_standar_plapat.RData")
richness_subplot30_dob_standar_plapat=richness_subplot30_dob_standar_plapat%>%rename(plotid=plotID)
save(richness_subplot30_dob_standar_plapat,file="richness_subplot30_dob_standar_plapat.RData")
data_dob_plapat=bind_rows(richness_subplot10_dob_standar_plapat%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_plapat%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_plapat%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_plapat%>%select( plotid, mean_value,area))

#### for the parasitic fungi

richness_subplot10_dob_standar_para=richness_subplot10_dob_standar_para%>%mutate(area=rep(100,length.out=n))
richness_subplot20_dob_standar_para=richness_subplot20_dob_para
richness_subplot20_dob_standar_para%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_para
save(richness_subplot20_dob_standar_para,file="richness_subplot20_dob_standar_para.RData")
richness_subplot30_dob_standar_para=richness_subplot30_dob_standar_para%>%rename(plotid=plotID)
save(richness_subplot30_dob_standar_para,file="richness_subplot30_dob_standar_para.RData")
data_dob_para=bind_rows(richness_subplot10_dob_standar_para%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_para%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_para%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_para%>%select( plotid, mean_value,area))

### for the AM fungi

richness_subplot10_dob_standar_AM=richness_subplot10_dob_standar_AM%>%mutate(area=rep(100,length.out=n))
richness_subplot20_dob_standar_AM=richness_subplot20_dob_AM
richness_subplot20_dob_standar_AM%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_AM
save(richness_subplot20_dob_standar_AM,file="richness_subplot20_dob_standar_AM.RData")
richness_subplot30_dob_standar_AM=richness_subplot30_dob_standar_AM%>%rename(plotid=plotID)
save(richness_subplot30_dob_standar_AM,file="richness_subplot30_dob_standar_AM.RData")

data_dob_AM=bind_rows(richness_subplot10_dob_standar_AM%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_AM%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_AM%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_AM%>%select( plotid, mean_value,area))

## for the epiphytes

richness_subplot10_dob_standar_epiphy=richness_subplot10_dob_standar_epiphy%>%mutate(area=rep(100,length.out=n))
richness_subplot20_dob_standar_epiphy=richness_subplot20_dob_epiphy
richness_subplot20_dob_standar_epiphy%>%mutate(area=rep(400,length.out=n))->richness_subplot20_dob_standar_epiphy
save(richness_subplot20_dob_standar_epiphy,file="richness_subplot20_dob_standar_epiphy.RData")
richness_subplot30_dob_standar_epiphy=richness_subplot30_dob_standar_epiphy%>%rename(plotid=plotID)
save(richness_subplot30_dob_standar_epiphy,file="richness_subplot30_dob_standar_epiphy.RData")

data_dob_epiphy=bind_rows(richness_subplot10_dob_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot20_dob_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot30_dob_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot40_dob_standar_epiphy%>%select( plotid, mean_value,area))




data_full_dob=bind_rows(data_dob_all%>%mutate(guild=rep("all",length.out=n)),data_dob_AM%>%mutate(guild=rep("AM",length.out=n)),data_dob_EM%>%mutate(guild=rep("EM",length.out=n)),data_dob_soilsap%>%mutate(guild=rep("soilsap",length.out=n)),data_dob_littersap%>%mutate(guild=rep("littersap",length.out=n)),data_dob_woodsap%>%mutate(guild=rep("woodsap",length.out=n)),data_dob_plapat%>%mutate(guild=rep("plapat",length.out=n)),data_dob_para%>%mutate(guild=rep("para",length.out=n)),data_dob_epiphy%>%mutate(guild=rep("epiphy",length.out=n)))

save(data_full_dob,file="data_full_dob.RData")
 ###3

data_neon_all=bind_rows(richness_subplot10_neon_standar%>%select(plotid, mean_value,area),richness_subplot20_neon_standar%>%select( plotid, mean_value,area),richness_subplot30_neon_standar%>%select( plotid, mean_value,area),richness_subplot40_neon_standar%>%select( plotid, mean_value,area))

### for neon sites for the EM
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_EM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_EM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_EM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_EM.RData")

richness_subplot40_neon_standar_EM$mean_value=as.numeric(richness_subplot40_neon_standar_EM$mean_value)
richness_subplot40_neon_standar_EM$area=as.numeric(richness_subplot40_neon_standar_EM$area)


data_neon_EM=bind_rows(richness_subplot10_neon_standar_EM%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_EM%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_EM%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_EM%>%select( plotid, mean_value,area))

######
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_AM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_AM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_AM.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_AM.RData")

richness_subplot40_neon_standar_EM$mean_value=as.numeric(richness_subplot40_neon_standar_EM$mean_value)
richness_subplot40_neon_standar_EM$area=as.numeric(richness_subplot40_neon_standar_EM$area)


data_neon_AM=bind_rows(richness_subplot10_neon_standar_AM%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_AM%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_AM%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_AM%>%select( plotid, mean_value,area))

##########

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_soilsap.RData")

richness_subplot40_neon_standar_soilsap$mean_value=as.numeric(richness_subplot40_neon_standar_soilsap$mean_value)
richness_subplot40_neon_standar_soilsap$area=as.numeric(richness_subplot40_neon_standar_soilsap$area)


data_neon_soilsap=bind_rows(richness_subplot10_neon_standar_soilsap%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_soilsap%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_soilsap%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_soilsap%>%select( plotid, mean_value,area))

####

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_littersap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_littersap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_littersap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_littersap.RData")

richness_subplot40_neon_standar_littersap$mean_value=as.numeric(richness_subplot40_neon_standar_littersap$mean_value)
richness_subplot40_neon_standar_littersap$area=as.numeric(richness_subplot40_neon_standar_littersap$area)


data_neon_littersap=bind_rows(richness_subplot10_neon_standar_littersap%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_littersap%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_littersap%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_littersap%>%select( plotid, mean_value,area))

#######3

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_woodsap.RData")

richness_subplot40_neon_standar_woodsap$mean_value=as.numeric(richness_subplot40_neon_standar_woodsap$mean_value)
richness_subplot40_neon_standar_woodsap$area=as.numeric(richness_subplot40_neon_standar_woodsap$area)


data_neon_woodsap=bind_rows(richness_subplot10_neon_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_woodsap%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_woodsap%>%select( plotid, mean_value,area))

########


load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_plapat.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_plapat.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_plapat.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_plapat.RData")

richness_subplot40_neon_standar_plapat$mean_value=as.numeric(richness_subplot40_neon_standar_plapat$mean_value)
richness_subplot40_neon_standar_plapat$area=as.numeric(richness_subplot40_neon_standar_plapat$area)


data_neon_plapat=bind_rows(richness_subplot10_neon_standar_plapat%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_plapat%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_plapat%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_plapat%>%select( plotid, mean_value,area))

######

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_para.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_para.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_para.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_para.RData")

richness_subplot40_neon_standar_para$mean_value=as.numeric(richness_subplot40_neon_standar_para$mean_value)
richness_subplot40_neon_standar_para$area=as.numeric(richness_subplot40_neon_standar_para$area)


data_neon_para=bind_rows(richness_subplot10_neon_standar_para%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_para%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_para%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_para%>%select( plotid, mean_value,area))

#####


load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar_epiphy.RData")

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot20_neon_standar_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot30_neon_standar_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot40_neon_standar_epiphy.RData")

richness_subplot40_neon_standar_epiphy$mean_value=as.numeric(richness_subplot40_neon_standar_epiphy$mean_value)
richness_subplot40_neon_standar_epiphy$area=as.numeric(richness_subplot40_neon_standar_epiphy$area)


data_neon_epiphy=bind_rows(richness_subplot10_neon_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot20_neon_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot30_neon_standar_epiphy%>%select( plotid, mean_value,area),richness_subplot40_neon_standar_epiphy%>%select( plotid, mean_value,area))

### all the data for the neon sites


data_full_neon=bind_rows(data_neon_all%>%mutate(guild=rep("all",length.out=n)),data_neon_AM%>%mutate(guild=rep("AM",length.out=n)),data_neon_EM%>%mutate(guild=rep("EM",length.out=n)),data_neon_soilsap%>%mutate(guild=rep("soilsap",length.out=n)),data_neon_littersap%>%mutate(guild=rep("littersap",length.out=n)),data_neon_woodsap%>%mutate(guild=rep("woodsap",length.out=n)),data_neon_plapat%>%mutate(guild=rep("plapat",length.out=n)),data_neon_para%>%mutate(guild=rep("para",length.out=n)),data_neon_epiphy%>%mutate(guild=rep("epiphy",length.out=n)))

save(data_full_neon,file="data_full_neon.RData")

# bind the two data sets

data_neon_dob=bind_rows(data_full_neon,data_full_dob)

save(data_neon_dob,file="data_neon_dob.RData")

######

data_neon_dob%>%mutate(cb_guild=paste(plotid,"_",guild))%>%filter(mean_value!=0)->temp_data

## to estimate the z and the c value at the plot scale
## only select the plots with at least three cores available to construct a sar

table(temp_data$cb_guild)%>%data.frame()%>%filter(Freq>2)->temp_plot

plot_id=unique(temp_plot$Var1)
sar_para=matrix(ncol=2,nrow=length(plot_id))
for (i in 1:length(plot_id))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  mod_data=temp_data%>%filter(cb_guild==plot_id[i])
  
  mod=lm(log(mean_value)~log(area),data=mod_data)
  
  sar_para[i,1]=coef(mod)[1]%>%as.numeric()
  sar_para[i,2]=coef(mod)[2]%>%as.numeric()
}

sar_para%>%data.frame()%>%mutate(plotid=plot_id)%>%rename_all(~paste0(c("logc","zvalue","plotid")))->temp_zvalue

#add a plot id to the data
# to detect if there is 0 in the plot id

mm=str_detect(temp_zvalue$plotid,"0")

plot_ID_z=matrix(ncol=1,nrow=length(plot_id))

for (i in 1:length(plot_id))
  {
  
  if(str_detect(temp_zvalue[i,3],"0")){
    plot_ID_z[i,1] =substr(temp_zvalue[i,3],1,8)
  }
  else{
    plot_ID_z[i,1] =substr(temp_zvalue[i,3],1,3)
  }
}

temp_zvalue%>%mutate(plot_ID_z=plot_ID_z)%>%head()



# to spilit the plotid based on the

plot_ID_z=matrix(ncol=2,nrow=length(plot_id))

for (i in 1:length(plot_id))
{
  plot_ID_z[i,]=strsplit(temp_zvalue$plotid%>%as.character()," _ ")[[i]]
}

plot_ID_z=plot_ID_z%>%as.matrix()%>%data.frame()

data_zvalue=bind_cols(temp_zvalue,plot_ID_z)%>%rename_all(~paste0(c("logc","zvalue","plot_guild","plotid","guild")))

data_zvalue%>%rename(plotID=plotid)->data_zvalue
  
save(data_zvalue,file="data_zvalue.RData")


# add the environmental variables to the data

load("~/soil-sar/plot-sar-permutation/model_data.RData")

plot_env=model_data%>%select(plotID,siteIDD,organicCPercent,ph,bdod,nitrogen,cec,sand,bio1,bio2,bio4,bio8,bio12,bio15,bio18,richness)%>%distinct()%>%rename()

#standardize the environmental data


k=left_join(data_zvalue,plot_env,by="plotID")

data=k%>%filter(guild=="all")

## to see the core level richness


data_neon_dob%>%filter(area==100)->temp

ggplot(temp%>%filter(guild!="all"),aes(x=guild,y=mean_value,fill=guild),alpha=0.5)+
  sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(width = 0.1, color = "black", size=0.1,outlier.size = 1) 


df=data.frame(site=rep(c("A","B"),times=100),value=rnorm(100, mean = 2))


ggplot(df,aes(x=site,y=value,fill=site),alpha=0.5)+
  sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(width = 0.1, color = "black", size=0.1,outlier.size = 1) 


  

