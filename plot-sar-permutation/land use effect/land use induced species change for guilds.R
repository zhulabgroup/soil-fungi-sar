
setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

library(phyloseq)
library(reshape2)
library(stringr)
library(ggplot2)
library(terra)
library(raster)


#(1)# determining the habitat affinity for each fungal guild and land use type

#add the fungal guild data to the full sample data
ft<- read_csv("~/soil-sar/plot-sar-permutation/FungalTraits_1.2_ver_16Dec_2020.csv")
ft%>%data.frame()%>%rename(genus=GENUS)%>%select(genus,primary_lifestyle)->ft_temp
tax_table(rare_all)%>%data.frame()%>%left_join(ft_temp,by="genus")->guild_temp
taxa_matrix <- as.matrix(guild_temp)
new_taxa_table <- tax_table(taxa_matrix)
row.names(new_taxa_table )=row.names(tax_table(rare_all))
tax_table(rare_all) <- new_taxa_table
rare_all_guild=rare_all

save(rare_all_guild,file="rare_all_guild.RData")

#(2) add the land use type data to the whole data set

load("~/soil-sar/plot-sar-permutation/plot_diversity_env_land.RData")

d=sample_data(rare_all_guild)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site corresponds to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(rare_all_guild, plotIDM)# merge the new plotid with the initial data 

sample_data(d)%>%data.frame()%>%dplyr::select(plotIDM)%>%rename(plotID=plotIDM)->temp

plot_diversity_env_land%>%dplyr::select(plotID,type)%>%distinct()%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))%>%
filter(!is.na(type))->land#(if we do not filter out the NAs, the join data will be longer)

temp%>%left_join(land,by="plotID")%>%dplyr::select(plotID,type)%>%dplyr::select(type)->temp


#assign NA rows with "evergreen forest"
row.names(temp)=row.names(sample_data(d))
temp=sample_data(temp)
d<- merge_phyloseq(rare_all_guild, temp)#
# to see the number of different soil cores among different land use types
sample_data(d)%>%data.frame()%>%dplyr::select(type)%>%group_by(type)%>%summarize(count=n())

# cultivatedCrops              397
# deciduousForest             1445
# dwarfScrub                    12
# emergentHerbaceousWetlands    46
# evergreenForest             2066
# grasslandHerbaceous         1018
# mixedForest                  342
# sedgeHerbaceous               14
# shrubScrub                   556
# woodyWetlands                312
# NA                           170#some plots still do not have land use data

rare_all_guild=d
save(rare_all_guild,file="rare_all_guild.RData")

plot_diversity_env_land%>%dplyr::select(plotID,type)->temp

# combined the pasture hay and cultivated croplands

temp%>%mutate(type = ifelse(type == "pastureHay", "cultivatedCrops" , type))->temp

#need to exclude the na
temp=temp%>%filter(!is.na(type))

#(3) get the guild-specific c and z values 

full_parameter_data%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%filter(!is.na(type))%>%
  group_by(guild)%>%summarise(mean_cvalue=mean(logc,na.rm=TRUE))%>%mutate(c=2.71828^mean_cvalue)->c_temp
#AM             -3.15  0.0427
#EM             -0.807 0.446 
#all             1.84  6.32  
#epiphy         -2.96  0.0520
#littersap      -1.07  0.344 
#para           -1.36  0.257 
#plapat         -1.51  0.221 
#soilsap         0.367 1.44  
#woodsap        -1.88  0.152 


full_parameter_data%>%dplyr::select(plotID,logc,zvalue,guild)%>%left_join(temp,by="plotID")%>%filter(!is.na(type))%>%
  group_by(guild)%>%summarise(mean_zvalue=mean(zvalue,na.rm=TRUE))->z_temp

#guild     mean_zvalue
#AM              0.775
#EM              0.771
#all             0.714
#epiphy          0.745
#littersap       0.717
#para            0.657
#plapat          0.691
#soilsap         0.653
#woodsap         0.751


#(4)# get the guild-specific habitat affinity with the formula of affinity=(S1/S2)^1/z
# where S1 and S2 indicate the core-level of fungal richness in the human-modified and the natural habitats, respectively.
# get the same number of cores from each land use type and determined the mean value of core-level richness
#####
# split the data into different guilds

data_EM <- subset_taxa(rare_all_guild, primary_lifestyle == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_guild, primary_lifestyle == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_guild, primary_lifestyle == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_guild, primary_lifestyle == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_guild, primary_lifestyle == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_guild, primary_lifestyle == "wood_saprotroph")
data_para <- subset_taxa(rare_all_guild, primary_lifestyle%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_guild, primary_lifestyle == "epiphyte")

### write a function to manipulate all the data sets


###
type0=c("cultivatedCrops","deciduousForest","evergreenForest","grasslandHerbaceous","mixedForest","shrubScrub","woodyWetlands")

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(rare_all_guild,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_all_updated=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_all_updated[i,]=kk[[i]]
}

land_rich_all_updated%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_all_updated

save(land_rich_all_updated,file="land_rich_all_updated.RData")
####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_plapat,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_plapat=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_plapat[i,]=kk[[i]]
}

land_rich_plapat%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_plapat_updated

save(land_rich_plapat_updated,file="land_rich_plapat_updated.RData")

###
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_soilsap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_soilsap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_soilsap[i,]=kk[[i]]
}

land_rich_soilsap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_soilsap_updated

save(land_rich_soilsap_updated,file="land_rich_soilsap_updated.RData")
#####
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_littersap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_littersap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_littersap[i,]=kk[[i]]
}

land_rich_littersap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_littersap_updated

save(land_rich_littersap_updated,file="land_rich_littersap_updated.RData")

##

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_woodsap,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_woodsap=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_woodsap[i,]=kk[[i]]
}

land_rich_woodsap%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_woodsap_updated

save(land_rich_woodsap_updated,file="land_rich_woodsap_updated.RData")
###

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_para,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_para=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_para[i,]=kk[[i]]
}

land_rich_para%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_para_updated

save(land_rich_para_updated,file="land_rich_para_updated.RData")
####
set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_epiphy,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_epiphy=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_epiphy[i,]=kk[[i]]
}

land_rich_epiphy%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_epiphy_updated

save(land_rich_epiphy_updated,file="land_rich_epiphy_updated.RData")




####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_EM,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_EM=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_EM[i,]=kk[[i]]
}

land_rich_EM%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_EM_updated

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

save(land_rich_EM_updated,file="land_rich_EM_updated.RData")


####

set.seed(123)
kk=list()
for(j in 1:500)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # 
  richness=numeric()
  for (i in 1:7)
  {
    dk=subset_samples(data_AM,type==type0[i])
    sample_names=sample_names(dk)
    sampled_names <- sample(sample_names, 300)
    sampled_physeq <- prune_samples(sampled_names, dk)
    richness[i]=estimate_richness(sampled_physeq, measures = "Observed")%>%summarize(mean_A = mean(Observed, na.rm = TRUE))%>%as.numeric()
  } 
  kk[[j]]=richness
}

land_rich_AM=matrix(ncol=7,nrow=500)
for (i in 1:500)
{
  land_rich_AM[i,]=kk[[i]]
}

land_rich_AM%>%data.frame()%>%rename_all(~paste0(c(type0)))%>%melt()->land_rich_AM_updated

save(land_rich_AM_updated,file="land_rich_AM_updated.RData")


####

land_rich_AM_updated=richness_among_land(data_AM)

land_rich_soilsap_updated=richness_among_land(data_soilsap)
land_rich_woodsap_updated=richness_among_land(data_woodsap)
land_rich_littersap_updated=richness_among_land(data_littersap)
land_rich_plapat_updated=richness_among_land(data_plapat)
land_rich_para_updated=richness_among_land(data_para)
land_rich_epiphy_updated=richness_among_land(data_epiphy)
##save all the results
save(land_rich_AM_updated,file="land_rich_AM_updated.RData")
save(land_rich_soilsap_updated,file="land_rich_soilsap_updated.RData")
save(land_rich_woodsap_updated,file="land_rich_woodsap_updated.RData")
save(land_rich_littersap_updated,file="land_rich_littersap_updated.RData")
save(land_rich_plapat_updated,file="land_rich_plapat_updated.RData")
save(land_rich_para_updated,file="land_rich_para_updated.RData")
save(land_rich_epiphy_updated,file="land_rich_epiphy_updated.RData")



mean_richness_guild=bind_rows(land_rich_AM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_EM_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_soilsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_littersap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_woodsap_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_plapat_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_para_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_epiphy_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)),
                              land_rich_all_updated%>%group_by(variable)%>%summarize(mean=mean(value,na.rm=TRUE)))%>%
  mutate(guild=rep(c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all"),each=7))%>%data.frame()

guild_type=c("AM","EM","soilsap","littersap","woodsap","plapat","para","epiphy","all")

# determining the sensitivity for different guilds
sensitivity=numeric()
for (i in 1:9)
{
  df=mean_richness_guild%>%filter(guild==guild_type[i])
  sensitivity[i]=df[1,2]/df[2:7,2] %>%mean()
}



sensitivity%>%data.frame()%>%bind_cols(guild_type)%>%data.frame()%>%rename_all(~paste0(c("sensitivity","guild")))%>%
  left_join(z_temp,by="guild")%>%mutate(affinity=sensitivity^(1/mean_zvalue))%>%left_join(c_temp,by="guild")->parameters

#sensitivity     guild mean_zvalue  affinity mean_cvalue          c
#1   1.2114232        AM   0.7747318 1.2809016  -3.1526586 0.04273844
#2   0.7160444        EM   0.7713955 0.6485615  -0.8067486 0.44630710
#3   0.7202833   soilsap   0.6533710 0.6052092   0.3667578 1.44304802
#4   1.1929157 littersap   0.7174158 1.2787501  -1.0679487 0.34371309
#5   0.9605541   woodsap   0.7507297 0.9478037  -1.8844469 0.15191325
#6   1.2578671    plapat   0.6912454 1.3935997  -1.5114489 0.22059035
#7   0.8884922      para   0.6571817 0.8353507  -1.3579871 0.25717818
#8   0.7177464    epiphy   0.7449075 0.6406908  -2.9560207 0.05202563
#9   0.9266827       all   0.7143864 0.8988971   1.8441516 6.32272547

save(parameters,file="parameters.RData")
# replace the prior parameters

# guild specific richness change caused by land use change
# use the 2015 land use change data as the baseline

#(1) load in the land use change data

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/land use effect")

raster1 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2015.nc")# use the 2015 data as the base line

ext(raster1) <- c(-90, 90, -180, 180)
crs(raster1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent (-72, -18, -170, -55), "SpatialPolygons")

coarser_raster <- aggregate(raster1, fact = 3, fun = mean) # convert it to a coarse resolution


PFT_2015 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster<- crop(coarser_raster[[i]], b)
  cropped_raster=flip(t(cropped_raster),direction = "horizontal")
  ext <- ext(cropped_raster)
  # Create a new extent by flipping the y-axis
  new_ext <- ext(-169.95, -55.05, 18, 72)
  new_raster <- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
                     xmin=-169.95, xmax=-55.05, 
                     ymin=18, ymax=72, 
                     crs=crs(cropped_raster))
  values(new_raster) <- values(cropped_raster)
  
  coords_present <- xyFromCell(new_raster, cell = 1:ncell(new_raster)) 
  
  # get the coordinates
  # the coordinates are not in the right format and should be converted for extracting othere data
  cell_values <- raster::extract(new_raster, coords_present) %>% as.matrix()
  PFT_2015[, i] <- cell_values
}



#raster_equal<- project(raster1, equal_area_crs)
#coarser_raster <- aggregate(raster1, fact = 3, fun = mean) # convert it to a coarse resolution
# make projection and then 
#equal_area_crs <- "EPSG:5070"
#raster_equal<- projectRaster(r, equal_area_crs)



#(2) determine the different proportions of land use type per grid

#PFT_2015 <- matrix(nrow = 275760, ncol = 33)
#for (i in 1:33) {
#cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
#cropped_raster<- crop(coarser_raster[[i]], b)
#coords_present <- xyFromCell( cropped_raster, cell = 1:ncell( cropped_raster)) 
   
#cell_values <- raster::extract( cropped_raster, coords_present) %>% as.matrix()
  
#PFT_2015[, i] <- cell_values
}#


PFT_2015 <- cbind(coords_present, PFT_2015) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster1)))) 

PFT_2015=PFT_2015%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2015[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2015%>%bind_cols(No_crop)->PFT_2015

#(3)# determine the total area for each grid cell
# Function to calculate the area of a 10-minute grid cell
calculate_grid_area <- function(lat, lon) {
  # Define the half-width and half-height of a 10-minute cell in degrees
  half_width <- 10 / 60 / 2
  half_height <- 10 / 60 / 2
  
  # Define the four corners of the grid cell
  top_left <- c(lat + half_height, lon - half_width)
  top_right <- c(lat + half_height, lon + half_width)
  bottom_left <- c(lat - half_height, lon - half_width)
  bottom_right <- c(lat - half_height, lon + half_width)
  
  # Calculate the area of the quadrilateral defined by the four corners
  area <- areaPolygon(matrix(c(top_left[2], top_left[1],
                               top_right[2], top_right[1],
                               bottom_right[2], bottom_right[1],
                               bottom_left[2], bottom_left[1],
                               top_left[2], top_left[1]), 
                             ncol = 2, byrow = TRUE))
  return(area)
}


## get the latitude of the data

coords_present%>%data.frame()%>%rename_all(~c("lon","lat"))%>%dplyr::select(lat,lon)->coords_present
coords_present$area <- mapply(calculate_grid_area, coords_present$lat, coords_present$lon)
coords_present=data.frame(coords_present)
area=coords_present$area
PFT_2015$area=coords_present$area



#(4)# determine the richness of different guilds based on the countryside SAR parameters and habitat area
# the affinity in the no_crop habitat is assumed to be 1
# for each guild, the z value is the mean across all land use types
#affinity=(sensitivity)1/z
# the sensitivity in the cropland is determined as the richness ratio between the cropland and the non-cropland areas

guild_specific_richness=matrix(ncol=9,nrow=dim(PFT_2015)[1])

for(i in 1:9)
{
  guild_specific_richness[,i]=parameters$c[i]*((PFT_2015$No_crop/100*PFT_2015$area)+(parameters$affinity[i]*PFT_2015$crop/100*PFT_2015$area))^parameters$mean_zvalue[i]
  
}

guild_specific_richness%>%data.frame()%>%rename_all(~paste0(guild_type))->current_richness_guild_land


## for future richness in the year of 2100 for the scenario of rcp245


raster2 <- rast("GCAM_Demeter_LU_ssp2_rcp45_hadgem_2100.nc")# use the 2015 data as the base line
ext(raster2) <- c(-90, 90, -180, 180)
crs(raster2) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent(-72, -18, -170, -55), "SpatialPolygons")

coarser_raster_2100 <- aggregate(raster2, fact = 3, fun = mean) # convert it to a coarse resolution



PFT_2100 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster<- crop(coarser_raster_2100[[i]], b)
  cropped_raster=flip(t(cropped_raster),direction = "horizontal")
  ext <- ext(cropped_raster)
  # Create a new extent by flipping the y-axis
  new_ext <- ext(-169.95, -55.05, 18, 72)
  new_raster <- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
                     xmin=-169.95, xmax=-55.05, 
                     ymin=18, ymax=72, 
                     crs=crs(cropped_raster))
  values(new_raster) <- values(cropped_raster)
  
  coords_present <- xyFromCell(new_raster, cell = 1:ncell(new_raster)) 
  
  # get the coordinates
  # the coordinates are not in the right format and should be converted for extracting othere data
  cell_values <- raster::extract(new_raster, coords_present) %>% as.matrix()
  PFT_2100[, i] <- cell_values
}



#PFT_2100 <- matrix(nrow = 275760, ncol = 33)
#for (i in 1:33) {
#  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
 # cropped_raster_2100<- crop(coarser_raster_2100[[i]], b)

 # coords_present <- xyFromCell(cropped_raster_2100, cell = 1:ncell(cropped_raster_2100)) 
  # get the coordinates
  # the coordinates are not in the right format and should be converted for extracting othere data
 # cell_values <- raster::extract(cropped_raster_2100, coords_present) %>% as.matrix()
  #PFT_2100[, i] <- cell_values
#}



PFT_2100 <- cbind(coords_present, PFT_2100) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster2)))) 

PFT_2100=PFT_2100%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2100[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2100%>%bind_cols(No_crop)->PFT_2100

PFT_2100%>%mutate(area=area)->PFT_2100


## determine the richness per cell for the year of 2100

guild_specific_richness=matrix(ncol=9,nrow=dim(PFT_2100)[1])
for(i in 1:9)
{
  guild_specific_richness[,i]=parameters$c[i]*((PFT_2100$No_crop/100*PFT_2100$area)+(parameters$affinity[i]*PFT_2100$crop/100*PFT_2100$area))^parameters$mean_zvalue[i]
  
}
guild_specific_richness%>%data.frame()%>%rename_all(~paste0(guild_type))->richness_guild_land_2100

# determine the changes in richness between the two time points caused by land use change

df_rcp245_land=(richness_guild_land_2100-current_richness_guild_land)/current_richness_guild_land



#df_rcp585=(future_richness_rcp585-current_richness_guild)/current_richness_guild

## to look at the present and future richness for each grid

present_future_richness0=list()
for (i in 1:9)
{
  present_future_richness0[[i]]=bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))
}

# bind all the data.frames

present_future_richness=present_future_richness0[[1]]
for (i in 2:9)
{
  present_future_richness=rbind(present_future_richness,present_future_richness0[[i]])
}

present_future_richness%>%bind_cols(df_rcp245_land%>%melt())->species_change_land_rcp245

## look at the data to make necessary conversions

### there are 86818  cases that no species are there for a grid
d=list()
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&is.na(x2))%>%dim()
  
}

# there are 0 cases that there are species presently but NA in the future

d=list()
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(!is.na(x1)&is.na(x2))%>%dim()
  
}

d=list()# there were 0 cases where we have NA species at present but some species in future
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&!is.na(x2))%>%dim()
  
}

d=list()# for two guilds,  we have two cases where some species present while 0 species in the future (AM-158,epiphy-667 grids))
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1!=0&x2==0)%>%dim()
  
}


d=list()# for two guilds, we 0 cases where  0 species but some in the future (AM-158,ephiphy-667 grids)
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2!=0)%>%dim()
  
}

d=list()#  we have 867 cases where 0 species present for both time points
for (i in 1: 9)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2==0)%>%dim()
  
}

# make some changes for the species change ratio


#present_future_richness%>%bind_cols(df_rcp585%>%melt())->species_change_temp

present_future_richness%>%bind_cols(df_rcp245_land%>%melt())->species_change_land_rcp245

species_change_land_rcp245%>%mutate(value=if_else(x1==0&x2==0,0,value))->species_change_land_rcp245



species_change_land_rcp245%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
  summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->tem_df_rcp245_land



## for future richness in 2100 in the scenario of RCP585


raster3 <- rast("GCAM_Demeter_LU_ssp5_rcp85_hadgem_2100.nc")
ext(raster3) <- c(-90, 90, -180, 180)
crs(raster3) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b <- as(extent(-72, -18, -170, -55), "SpatialPolygons")
coarser_raster_2100 <- aggregate(raster3, fact = 3, fun = mean) # convert it to a coarser resolution


PFT_2100 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster<- crop(coarser_raster_2100[[i]], b)
  cropped_raster=flip(t(cropped_raster),direction = "horizontal")
  ext <- ext(cropped_raster)
  # Create a new extent by flipping the y-axis
  new_ext <- ext(-169.95, -55.05, 18, 72)
  new_raster <- rast(nrows=nrow(cropped_raster), ncols=ncol(cropped_raster), 
                     xmin=-169.95, xmax=-55.05, 
                     ymin=18, ymax=72, 
                     crs=crs(cropped_raster))
  values(new_raster) <- values(cropped_raster)
  
  coords_present <- xyFromCell(new_raster, cell = 1:ncell(new_raster)) 
  
  # get the coordinates
  # the coordinates are not in the right format and should be converted for extracting other data
  cell_values <- raster::extract(new_raster, coords_present) %>% as.matrix()
  PFT_2100[, i] <- cell_values
}

PFT_2100 <- cbind(coords_present, PFT_2100) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster3)))) 





PFT_2100=PFT_2100%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2100[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2100%>%bind_cols(No_crop)->PFT_2100

PFT_2100=PFT_2100%>%mutate(area=area)


## determine the richness per cell for the year of 2100 in the scenario of RCP585

guild_specific_richness=matrix(ncol=9,nrow=dim(PFT_2100)[1])
for(i in 1:9)
{
  guild_specific_richness[,i]=parameters$c[i]*((PFT_2100$No_crop/100*PFT_2100$area)+(parameters$affinity[i]*PFT_2100$crop/100*PFT_2100$area))^parameters$mean_zvalue[i]
  
}
guild_specific_richness%>%data.frame()%>%rename_all(~paste0(guild_type))->richness_guild_land_2100

# determine the changes in richness between the two time points caused by land use change

df_rcp585_land=(richness_guild_land_2100-current_richness_guild_land)/current_richness_guild_land





## to look at the present and future richness for each grid

present_future_richness0=list()
for (i in 1:9)
{
  present_future_richness0[[i]]=bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))
}

# bind all the data.frames

present_future_richness=present_future_richness0[[1]]
for (i in 2:9)
{
  present_future_richness=rbind(present_future_richness,present_future_richness0[[i]])
}

present_future_richness%>%bind_cols(df_rcp585_land%>%melt())->species_change_land_rcp585

## look at the values to make necessary conversions

### there are 86818 cell with no species (NA)
d=list()
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&is.na(x2))%>%dim()
}

# there are 0 cases where there are some species presently but NA in the future

d=list()
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(!is.na(x1)&is.na(x2))%>%dim()
}

d=list()# there were 0 cases where we have NA for present but some species in future
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(is.na(x1)&!is.na(x2))%>%dim()
  
}

d=list()# there are two cases where we have some species presently while 0 species in the future
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1!=0&x2==0)%>%dim()
  
}


d=list()# there are 0 cases where  0 species but some in the future
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2!=0)%>%dim()
  
}

d=list()#  we have 867 cases where 0 species present for both time points, in this case, the change would be 0%
for (i in 1: 8)
{
  d[[i]]= bind_cols(current_richness_guild_land[,i],richness_guild_land_2100[,i])%>%data.frame()%>%rename_all(~paste0(c("x1","x2")))%>%filter(x1==0&x2==0)%>%dim()
  
}


present_future_richness%>%bind_cols(df_rcp585_land%>%melt())->species_change_land_rcp585

species_change_land_rcp585%>%mutate(value=if_else(x1==0&x2==0,0,value))->species_change_land_rcp585



species_change_land_rcp585%>%mutate(type=ifelse(value > 0, "Positive", "Negative"))%>%filter(!is.na(type))%>%group_by(variable,type)%>%
  summarise(mean_value = mean(value,na.rm=TRUE),sd_value = sd(value,na.rm=TRUE),count=n())->tem_df_rcp585_land


species_change_land_rcp245%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
  overall_change

tem_df_rcp245_land%>%left_join(overall_change%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")->tem_df_rcp245_land

species_change_land_rcp585%>%group_by(variable)%>%summarize(overal_mean=mean(value,na.rm=TRUE),overal_sd=sd(value,na.rm=TRUE),count0=n())%>%data.frame()->
  overall_change_585

tem_df_rcp585_land%>%left_join(overall_change_585%>%dplyr::select(variable,overal_mean,overal_sd,count0),by="variable")->tem_df_rcp585_land




###################################################
# create a map to show the spatial variation in the richness

species_change_land_rcp245%>%filter(variable=="all")%>%bind_cols(coords_present) ->change_richness_rcp245

species_change_land_rcp585%>%filter(variable=="all")%>%bind_cols(coords_present) ->change_richness_rcp585




ggplot(change_richness_rcp245) +
  geom_point(data = change_richness_rcp245, pch = 15, aes(x = x, y = y, color = 100*value), size = 0.275) +
  
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid = "yellow", high = "purple", na.value = "white")+ 
  xlab("Predicted species loss") +
  ylab("")+
  theme(legend.position = c(0.25,0.38),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("RCP4.5 & SSP2")
  
 





ggplot(data=tem_df_rcp245_land,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  geom_errorbar(data=tem_df_rcp245_land, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
  scale_fill_manual("RCP245",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#8fd1e1","#fedc5e"))

  theme(legend.position = c(0.8,0.25),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  geom_vline(xintercept =0,color="gray",linetype="dashed")+
  ylab("")+

  xlab("Species change rate")+
  scale_y_discrete(breaks=guild_type,labels=c("AM","EM","Soil sapro.","Litter sapro.","Wood sapro.","Plant patho.","Parasite","Epiphyte","All"))+
  geom_segment(data=tem_df_rcp245_land,size=0.25,color="black",aes(x=overal_mean-overal_sd,xend=overal_mean+overal_sd,y=variable,yend=variable))+
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
  geom_hline(yintercept = 8.5,color="red",size=1,alpha=0.3,linetype="dotted")+
  xlim(-0.15,0.15)+
  ggtitle("RCP4.5 & SSP2")

  





p3=ggplot(change_richness_rcp585) +
  geom_point(data = change_richness_rcp585, pch = 15, aes(x = x, y = y, color = 100*value), size = 0.275) +
  
  scale_color_gradient2(expression("Change %"), low = "seagreen", mid = "yellow", high = "purple",midpoint = 0, na.value = "white")+ 
 
  ylab("")+
  theme(legend.position = c(0.25,0.38),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  ggtitle("RCP8.5 & SSP5")+
  xlab("Longitude")+
  ylab("Latitude")


p4=ggplot(data=tem_df_rcp585_land,aes(fill=type,y=variable ,x=mean_value))+
  geom_col(width = 0.5)+
  #geom_errorbar(data=tem_df_rcp585_land, aes(xmin = mean_value- sd_value/sqrt(count), xmax = mean_value +sd_value/sqrt(count)),width=0.2)+
  scale_fill_manual("RCP585",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#8fd1e1","#fedc5e"))+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  geom_vline(xintercept =0,color="gray",linetype="dashed")+
  ylab("")+
  geom_segment(data=tem_df_rcp585_land,size=0.25,color="black",aes(x=overal_mean-overal_sd,xend=overal_mean+overal_sd,y=variable,yend=variable)
               )+
  geom_point(aes(y=variable,x=overal_mean),pch=23,color="black",size=2,fill="seagreen1")+
  xlab("Species change rate")+
  scale_y_discrete(breaks=guild_type,labels=c("AM","EM","Soil sapro.","Litter sapro.","Wood sapro.","Plant patho.","Parasite","Epiphyte","All"))+
  ggtitle("RCP 8.5 & SSP 5")+
  geom_hline(yintercept = 8.5,color="red",size=1,alpha=0.3,linetype="dotted")+
  xlim(-0.15,0.15)

p1=ggplotGrob(p1)
p2=ggplotGrob(p2)
p3=ggplotGrob(p3)
p4=ggplotGrob(p4)

p1$heights=p2$heights
p3$heights=p4$heights


plot_grid(p1,p2,p3,p4,ncol=2)

dk=st_as_sf(north_america_cropped)

dd1=rast(dd)

