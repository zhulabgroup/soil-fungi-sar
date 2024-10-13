## get the history map for the conus
land_use_data <- rast("nlcd_2001_land_cover_l48_20210604.img")

result_raster <- project(land_use_data, 'EPSG:4326')# select a CRS


plot_coordinates_convert=plot_coordinates[,2:3]




plot_coordinates_convert%>%
  st_as_sf(coords = c('lon', 'lat'), crs ="EPSG:4326" )%>%
  vect()%>%terra::project(land_use_data )->temp_coords

values=terra::extract(result_raster,temp_coords)

values=terra::extract(result_raster,plot_coordinates_convert)

## to select the cropland
rare_all_guild_biome=readRDS("rare_all_guild_biome.rds")

sample_data(rare_all_guild_biome)%>%data.frame()%>%
filter(type=="cultivatedCrops")%>%dplyr::select(Site,lon,lat,plotIDM)%>%distinct(plotIDM)


#BLAN_033-M-20.5-2.5-20151110-GEN-DNA1 BLAN
#SCBI_006-M-9-3.5-20151113-GEN-DNA2    SCBI
#STER_031-M-25-7.5-20160420-GEN-DNA1   STER
#SERC_007-M-6.5-28-20170412-GEN-DNA1   SERC
#DSNY_041-M-36-36-20170829-GEN-DNA1    DSNY
#JERC_004-M-0.5-25-20170724-GEN-DNA1   JERC
#KONA_003-M-12-2.5-20170807-GEN-DNA1   KONA
#ORNL_021-O-28.5-33-20161102-GEN-DNA1  ORNL
#LAJA_044-M-39-2.5-20170307-GEN-DNA1   LAJA

# 45 plots from 9 sites were included

sample_data(rare_all_guild_biome)%>%data.frame()%>%
  filter(type=="cultivatedCrops")%>%dplyr::select(Site,lon,lat,plotIDM)%>%distinct()->crop_plot


# extrat the value

values=terra::extract(result_raster,crop_plot[,c("lon","lat")])

#based on the 2001 data, some with NA while
#croopland 20

#pasture 11
# decious forest

bind_cols(values,crop_plot)

# to get all the year data

values=list()
for (i in 1630:2020)
{
cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
file_name=paste("conus_lulc_",i,".tif")
file_name=gsub(" ","",file_name)
land_history_data <- rast(file_name)
result_raster_history <- project(land_history_data, 'EPSG:4326',method="near")# select a CRS
values[[i]]=terra::extract(result_raster_history,crop_plot[,c("lon","lat")],method="simple")%>%
  data.frame()%>%dplyr::select(-ID)%>%rename_all(~paste0(c("year")))
}

# to see when the crop land stars to appear in the 400 years


# to see the variation in the land use type for each plot


## for each year and then bind them

subset_columns <- values[1630:2020]

do.call(cbind,subset_columns)%>%data.frame()->land_use_among_year

# for each plot to select the most often or year land use type

type=numeric()
for (i in 1:45)
  {
 type[i]=land_use_among_year[i,]%>%t()%>%data.frame()%>%
   rename_all(~paste0(c("name")))%>%
 group_by(name)%>%count()%>%data.frame()%>%
 slice_max(n,n=1)%>%pull(name)
  
}

# to find the latest version? or based on the nearest plot
crop_plot%>%bind_cols(type)%>%
  rename(history=...5)%>%
  mutate(historical_type=case_when(history%in%c("4","NA")~"forest", 
                                   history=="7"~"wetland", 
                                   history=="6"~"grassland", 
           TRUE~"forest"))->df4



# to get the pairwise plot to see how the plots are mathces
#there are cases that the initial land use type has disappeared so we need to find the neariest
site=df4$Site

d=list()
for(i in 1:9)
{
  d[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(type)
  
}
  #to check if one site has several biomes?

# one site has only one biome
biome_type=list()
for(i in 1:9)
{
  biome_type[[i]]=rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i])%>%distinct(LABEL)
  
}

  
# for the first site
# for the cropland plot
rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(plotIDM=="BLAN_033")

rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site=="BLAN")

crop_type=list()
for (i in 1:9)
  {
  crop_type[[i]]=df4%>%filter(Site==site[i])
}

# in most cases, there are only one type for their historical landuse type


rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i]&type=="deciduousForest")
# for the cropland data



# need to fins the plots with forest
#create a dataframe so that the historical and present land use type can be matched

all_land_cover=sample_data(rare_all_guild_biome)%>%data.frame()%>%distinct(type)%>%
  mutate(broad_type=case_when(type%in%c("evergreenForest", "deciduousForest","mixedForest")~"forest", 
                                   type%in%c("emergentHerbaceousWetlands","woodyWetlands")~"wetland", 
                                   type%in%c("grasslandHerbaceous","sedgeHerbaceous")~"grassland", 
                                  type==c("shrubScrub","dwarfScrub")~"shrub", 
                              
                                   TRUE~"other"))

# to select the same land use type to get the mean richness
# determe the richness based on each site

crop_data=subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)

#crop_data%>%sample_data()%>%data.frame()%>%filter(Site==site[i])->crop_sub_data
plotid=unique(df4$plotIDM)

response_ratio=numeric()
for (i in c(1:45))
  {
 cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = ""))
subset_samples(rare_all_guild_biome,plotIDM%in%df4$plotIDM)%>%
  subset_samples(plotIDM==plotid[i])->modified_data

unique(df4$Site)->modified_site
subset_samples(rare_all_guild_biome,type!="cultivatedCrops")->natural_data
  
sample_data(natural_data)%>%data.frame()%>%distinct(Site)%>%pull(Site)->natural_site
 
 # find the sites not included in the natural 
setdiff(modified_site,natural_site)->sites_no_in_nature
   

df4%>%filter(plotIDM==plotid[i])%>%distinct(historical_type)%>%pull(historical_type)->historical_type

all_land_cover%>%filter(broad_type%in%historical_type)%>%pull(type)->select_type

#rare_all_guild_biome%>%sample_data()%>%data.frame()%>%filter(Site==site[i]&type%in%select_type)%>%distinct(type)->natural_sub_data
# the plot within the same site
Site=df4$Site
#to see if have overlap for the two data sets
nearest_df_site%>%filter(site==Site[i])%>%pull(nearest_site)->near_site

matching_elements <- grep(paste("STER", "KONA",sep = "|"), plotid)

if(i%in%c(matching_elements,24))
  {
  subset_samples(natural_data,Site==near_site&type%in%select_type)->natural_data_sub
  
}
else if(i==12)
  {
  subset_samples(natural_data,Site==site[i])->natural_data_sub
}

else{
  #subset_samples(natural_data,Site==site[i])%>%sample_data()%>%data.frame()%>%distinct(type)%>%pull(type)->natural_type
  
  subset_samples(natural_data,Site==site[i]&type%in%select_type)->natural_data_sub
  
}
#when i =12, the history landtype was grassland but currently all are deciduous forest
# the neariest is also forest

#when i=24 the history landtype was forest but currently all are woodyWetlands

# for the LAJA site we do not have historical data and we assumed it was forest

#n_natural_sample=nsamples(natural_data_sub)
#n_modified_sample=nsamples(modified_data)

richness_modified=estimate_richness(modified_data, measures = "Observed")
richness_natural=estimate_richness(natural_data_sub, measures = "Observed")

response_ratio[i]=mean(richness_modified$Observed)/mean(richness_natural$Observed)
}
# for sites[3,4,7,8]


# only select the forest type and should not select the grassland
# it might be better to do that for each plot



