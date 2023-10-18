# get the plot-level root traits#
library(neonUtilities)
root.data <- loadByProduct(dpID="DP1.10067.001")#
root.chemi=root.data[[4]]
head(root.chemi)
root.chemi=root.chemi[,c("siteID","plotID","d15N","d13C","nitrogenPercent","carbonPercent","CNratio")]
head(root.chemi)
table(root.chemi$plotID)
root.chemi.mean=aggregate(root.chemi[,3:7], by=list(root.chemi$plotID), mean)
names(root.chemi.mean)[1]="plotID"
<<<<<<< HEAD
write.csv(root.chemi.mean,"root.chemi.mean.csv")

## get the simulated values for each as the response variable
a=list()
for(i in 1:473)
  {
  a[[i]]=t(rare_neon_z_simu[i,3:32])
}
b=a[[1]]
for (i in 2:473)
  {
  b=rbind(b,a[[i]])
}
b=data.frame(b)
pid=rep(rare_neon_z_simu$a1,each=30)

b=cbind(b,pid)
names(b)[1]="z"
plot_simu_z=b
## get the soil variables for the neon site

neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)# 

d=sample_data(neon)
plotID=substr(d$geneticSampleID,1,8)
d=cbind(d,plotID)
# soil pH and moisture data were derived from the Qin's data while for soil C and N, i got them from the soil product

d=d[,c("plotID","soilInCaClpH","soilMoisture")]

soil <- loadByProduct(dpID="DP1.10047.001")
soil.chemi=soil[[5]]

soil.chemi=soil.chemi[,c("phCacl2","carbonTot","nitrogenTot","ctonRatio","horizonID")]

k=soil$spc_biogeochem %>%
  # Get soil horizons that are completely less than 30 cm in depth
  dplyr::filter(horizonID %in% soil$spc_perhorizon$horizonID[which(soil$spc_perhorizon$horizonBottomDepth < 30)]) %>%
  left_join(dplyr::select(soil$spc_perplot, plotID, decimalLongitude, decimalLatitude), by="plotID") 

k=k[,c("plotID","carbonTot","nitrogenTot")]

k=subset(k,carbonTot!="NA"&nitrogenTot!="NA")
# the mean for soil C and N
soil.mean=aggregate(k[,2:3],by=list(k$plotID),mean)# there are 691 plot-level soil measurements
names(soil.mean)[1]="plotID"
# for both soil C and N and pH, the data are not complete, suggesting that some data may still need to be filled.
# the mean for soil pH
soil.mean.ph=aggregate(d[,2:3],by=list(d$plotID),mean)
names(soil.mean.ph)[1]="plotID"
# combine all the soil variables at the plot level
neon_soil=merge(soil.mean,soil.mean.ph,by="plotID")# the pH data would be more completed
p=merge(plot_simu_z,neon_climate_aridity,by="plotID",all.x = TRUE)# merge z with climate variables

core_mass_type_mean=aggregate(core_mass_type[,2:3],by=list(core_mass_type$plotID),mean)
names(core_mass_type_mean)[1]="plotID"
p=merge(p,core_mass_type_mean,by="plotID")# merge the root mass data but many plots do not have root mass data


write.csv(root.chemi.mean,"root.chemi.mean.csv")
