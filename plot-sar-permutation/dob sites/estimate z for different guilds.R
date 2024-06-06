## construct the sar based on different guilds 

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_EM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_EM.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale

  
richness_subplot5_dob_EM=richness_subplot5_dob_EM%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_EM=richness_subplot10_dob_EM%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_EM=richness_subplot20_dob_EM%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_EM=bind_rows(richness_subplot30A_dob_EM,richness_subplot30B_dob_EM,richness_subplot30C_dob_EM,richness_subplot30D_dob_EM)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_EM_scale=bind_rows(richness_subplot5_dob_EM,richness_subplot10_dob_EM,richness_subplot20_dob_EM,richness_mean_subplot30_dob_EM,richness_mean_subplot40_dob_EM%>%mutate(across(c(richness,area),as.numeric)))

save(data_EM_scale,file="data_EM_scale.RData")

###

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_soilsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_soilsap.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_soilsap=richness_subplot5_dob_soilsap%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_soilsap=richness_subplot10_dob_soilsap%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_soilsap=richness_subplot20_dob_soilsap%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_soilsap=bind_rows(richness_subplot30A_dob_soilsap,richness_subplot30B_dob_soilsap,richness_subplot30C_dob_soilsap,richness_subplot30D_dob_soilsap)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_soilsap_scale=bind_rows(richness_subplot5_dob_soilsap,richness_subplot10_dob_soilsap,richness_subplot20_dob_soilsap,richness_mean_subplot30_dob_soilsap,richness_mean_subplot40_dob_soilsap%>%mutate(across(c(richness,area),as.numeric)))

save(data_soilsap_scale,file="data_soilsap_scale.RData")

####

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_AM.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_AM.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_AM=richness_subplot5_dob_AM%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_AM=richness_subplot10_dob_AM%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_AM=richness_subplot20_dob_AM%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_AM=bind_rows(richness_subplot30A_dob_AM,richness_subplot30B_dob_AM,richness_subplot30C_dob_AM,richness_subplot30D_dob_AM)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_AM_scale=bind_rows(richness_subplot5_dob_AM,richness_subplot10_dob_AM,richness_subplot20_dob_AM,richness_mean_subplot30_dob_AM,richness_mean_subplot40_dob_AM%>%mutate(across(c(richness,area),as.numeric)))

save(data_AM_scale,file="data_AM_scale.RData")
####

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_littersap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_littersap.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_littersap=richness_subplot5_dob_littersap%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_littersap=richness_subplot10_dob_littersap%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_littersap=richness_subplot20_dob_littersap%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_littersap=bind_rows(richness_subplot30A_dob_littersap,richness_subplot30B_dob_littersap,richness_subplot30C_dob_littersap,richness_subplot30D_dob_littersap)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_littersap_scale=bind_rows(richness_subplot5_dob_littersap,richness_subplot10_dob_littersap,richness_subplot20_dob_littersap,richness_mean_subplot30_dob_littersap,richness_mean_subplot40_dob_littersap%>%mutate(across(c(richness,area),as.numeric)))

save(data_littersap_scale,file="data_littersap_scale.RData")


####

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_plapat.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_plapat.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_plapat=richness_subplot5_dob_plapat%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_plapat=richness_subplot10_dob_plapat%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_plapat=richness_subplot20_dob_plapat%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_plapat=bind_rows(richness_subplot30A_dob_plapat,richness_subplot30B_dob_plapat,richness_subplot30C_dob_plapat,richness_subplot30D_dob_plapat)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_plapat_scale=bind_rows(richness_subplot5_dob_plapat,richness_subplot10_dob_plapat,richness_subplot20_dob_plapat,richness_mean_subplot30_dob_plapat,richness_mean_subplot40_dob_plapat%>%mutate(across(c(richness,area),as.numeric)))

save(data_plapat_scale,file="data_plapat_scale.RData")


####

## construct the sar based on different guilds 

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_woodsap.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_woodsap.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_woodsap=richness_subplot5_dob_woodsap%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_woodsap=richness_subplot10_dob_woodsap%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_woodsap=richness_subplot20_dob_woodsap%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_woodsap=bind_rows(richness_subplot30A_dob_woodsap,richness_subplot30B_dob_woodsap,richness_subplot30C_dob_woodsap,richness_subplot30D_dob_woodsap)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_woodsap_scale=bind_rows(richness_subplot5_dob_woodsap,richness_subplot10_dob_woodsap,richness_subplot20_dob_woodsap,richness_mean_subplot30_dob_woodsap,richness_mean_subplot40_dob_woodsap%>%mutate(across(c(richness,area),as.numeric)))

save(data_woodsap_scale,file="data_woodsap_scale.RData")

##

## construct the sar based on different guilds 

load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_epiphy.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_epiphy.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_epiphy=richness_subplot5_dob_epiphy%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_epiphy=richness_subplot10_dob_epiphy%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_epiphy=richness_subplot20_dob_epiphy%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_epiphy=bind_rows(richness_subplot30A_dob_epiphy,richness_subplot30B_dob_epiphy,richness_subplot30C_dob_epiphy,richness_subplot30D_dob_epiphy)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_epiphy_scale=bind_rows(richness_subplot5_dob_epiphy,richness_subplot10_dob_epiphy,richness_subplot20_dob_epiphy,richness_mean_subplot30_dob_epiphy,richness_mean_subplot40_dob_epiphy%>%mutate(across(c(richness,area),as.numeric)))

save(data_epiphy_scale,file="data_epiphy_scale.RData")


load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot5_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot10_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot20_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30A_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30B_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30C_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_subplot30D_dob_para.RData")
load("~/soil-sar/plot-sar-permutation/dob sites/richness_mean_subplot40_dob_para.RData")

# to see the richness for each spatial scale
# get the mean at the 30 by 30 scale


richness_subplot5_dob_para=richness_subplot5_dob_para%>%mutate(area=rep(25,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot10_dob_para=richness_subplot10_dob_para%>%mutate(area=rep(100,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_subplot20_dob_para=richness_subplot20_dob_para%>%mutate(area=rep(400,length.out=n))%>%filter(!is.na(mean_value))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

richness_mean_subplot30_dob_para=bind_rows(richness_subplot30A_dob_para,richness_subplot30B_dob_para,richness_subplot30C_dob_para,richness_subplot30D_dob_para)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out=n))%>%select(plotid, mean_value, area)%>%rename_all(~paste0(c("plotid","richness","area")))

data_para_scale=bind_rows(richness_subplot5_dob_para,richness_subplot10_dob_para,richness_subplot20_dob_para,richness_mean_subplot30_dob_para,richness_mean_subplot40_dob_para%>%mutate(across(c(richness,area),as.numeric)))

save(data_para_scale,file="data_para_scale.RData")



