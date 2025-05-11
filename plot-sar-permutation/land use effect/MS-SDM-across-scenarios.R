# climate effects on diversity losses and gains with SDM

library(phyloseq)
library(dplyr)
library(reshape2)
library(SSDM)
library(rJava)

agg_data <- readRDS("neon_dob_prevalent_v4.1.Rds")

env <- sample_data(agg_data) %>%data.frame() %>%dplyr::select(lon, lat, mat_celsius, temp_seasonality, soilInCaClpH, organicCPercent, mat_celsius_2, map_mm_2)

occ <- otu_table(agg_data) %>%data.frame() %>%bind_cols(env[, 1:2])# the coordinates

occurrence_data <- melt(occ[, 1:8597])
df <- occ[, 8598:8599]
replicated_df <- do.call(rbind, replicate(8597, df, simplify = FALSE))
occurrence_data <- cbind(occurrence_data, replicated_df)
# check if all the species have occurrence data
d <- colSums(occ[, 1:8597])
sp <- colnames(occ[, 1:8597]) # model each species and save the output

#(1)# diversity projections based on current climate scenarios


cl <- makeCluster(10)
registerDoParallel(cl)

sdm_present=list()
for(i in 1:8597){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  tryCatch({sdm_present[[i]] <- ensemble_modelling(algorithms = c('MAXENT', 'RF',"GLM"), 
                                                   subset(occurrence_data,variable==sp[i]),
                                                   r_present_northam, Xcol = 'lon', Ycol = 'lat',Pcol="value", 
                                                   Spcol="variable",
                                                   
                                                   rep = 10,  # Number of repetitions
                                                   cv = 'holdout',  # Cross-validation method
                                                   cv.param = c(0.7, 1),  # Parameters for cross-validation
                                                   # Metric for evaluating models
                                                   ensemble.thresh = 0.75,  # Threshold for including models in the ensemble
                                                   weight = TRUE,  # Whether to weight models by their performance
                                                   
                                                   verbose = FALSE)
  saveRDS(sdm_present[[i]],file=paste0("result_present_glm_ensemble", i, ".rds"))},
  
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))
  }
  )
}

#(2)# future diversity patterns under the rcp2-4.5 scenarios

cl <- makeCluster(15)
registerDoParallel(cl)

sdm_future_rcp245_ensemble=list()
for(i in 2516:2524){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  tryCatch({sdm_future_rcp245_ensemble[[i]] <- ensemble_modelling(algorithms = c('MAXENT', 'RF',"GLM"), 
                                                                  subset(occurrence_data,variable==sp[i]),
                                                                  r_future_northam, Xcol = 'lon', Ycol = 'lat',Pcol="value", 
                                                                  Spcol="variable",
                                                                  
                                                                  rep = 10,  # Number of repetitions
                                                                  cv = 'holdout',  # Cross-validation method
                                                                  cv.param = c(0.7, 1),  # Parameters for cross-validation
                                                                  # Metric for evaluating models
                                                                  ensemble.thresh = 0.75,  # Threshold for including models in the ensemble
                                                                  weight = TRUE,  # Whether to weight models by their performance
                                                                  
                                                                  verbose = FALSE)
  saveRDS(sdm_future_rcp245_ensemble[[i]],file=paste0("result_future_rcp245—ensemble", i, ".rds"))},
  
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))
  }
  )
}

#(3)# future species distributions for the rcp5-8.5 scenarios
cl <- makeCluster(10)
registerDoParallel(cl)

sdm_future_rcp858_ensemble=list()
for(i in 3055:4000){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  tryCatch({sdm_future_rcp858_ensemble[[i]] <- ensemble_modelling(algorithms = c('MAXENT', 'RF',"GLM"), 
                                                                  subset(occurrence_data,variable==sp[i]),
                                                                  r_future_northam_rcp585, Xcol = 'lon', Ycol = 'lat',Pcol="value", 
                                                                  Spcol="variable",
                                                                  
                                                                  rep = 10,  # Number of repetitions
                                                                  cv = 'holdout',  # Cross-validation method
                                                                  cv.param = c(0.7, 1),  # Parameters for cross-validation
                                                                  # Metric for evaluating models
                                                                  ensemble.thresh = 0.75,  # Threshold for including models in the ensemble
                                                                  weight = TRUE,  # Whether to weight models by their performance
                                                                  
                                                                  verbose = FALSE)
  saveRDS(sdm_future_rcp858_ensemble[[i]],file=paste0("result_future_rcp858_ensemble", i, ".rds"))},
  
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))
  }
  )
}


### some species could not be modeled with all the three models
#in that case, we used a single model for the projectio

#(4)# check the code of the species that needs to be modeled with a single model

directory_path <- "/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP245"

## get the lacking number

files <- list.files(path = directory_path , pattern = "\\.rds$", full.names = TRUE)

dd=files %>%str_extract_all("\\d+") 

finished_number =numeric()
for (i in 1:length(files))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  finished_number[i]=dd[[i]][3]%>%as.numeric()
}
# Generate the complete sequence of numbers
complete_sequence <- seq(1:8597)
# Identify the missing numbers
missing_numbers <- setdiff(complete_sequence, finished_number)

# Get file information
file_info <- file.info(files)

# Filter files smaller than 40 bytes
small_files <- rownames(file_info[file_info$size < 100, ])

ddk=small_files %>%str_extract_all("\\d+") 

small_number =numeric()
for (i in 1:length(small_files))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  small_number[i]=ddk[[i]][3]%>%as.numeric()
}
all_numbers=unique(c(missing_numbers, small_number))%>%sort()

all_numbers_rcp245=all_numbers#

#(5)# for guild-level species projections

tax_table(rare_all_guild)%>%data.frame()%>%dplyr::select(primary_lifestyle)->df
df%>%mutate(sp=rownames(df))->species_guild
sp%>%data.frame()%>%rename_all(~paste0("sp"))%>%left_join(species_guild,by="sp")->op
op%>%mutate(ta2 = ifelse(str_detect(primary_lifestyle, "parasite"), "para", primary_lifestyle))->op
guild_model=c("soil_saprotroph","plant_pathogen","ectomycorrhizal","litter_saprotroph","wood_saprotroph","arbuscular_mycorrhizal","epiphyte","para")

model_guild=matrix(ncol=8,nrow=275760)
for (j in 1:8)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  ## future richness for the scenario of RCP 585 for total richness
  #cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  species_otu=op%>%mutate(number=1:8597)%>%filter(ta2==guild_model[j])%>%pull(number)%>%as.character()
  species_number=op%>%mutate(number=1:8597)%>%filter(ta2==guild_model[j])%>%pull(number)%>%as.character()%>%length()
  data <- matrix(ncol = species_number, nrow = 275760)
  
  for (i in seq_along(species_otu)) {
    
    if(i%in%all_numbers_rcp245)
    {
      setwd("/home/wenqil/results")
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      df <- readRDS(paste0("result_future", i, ".rds"))
      
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
    else {
      setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP245")
      
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      df <- readRDS(paste0("result_future_rcp245—ensemble", i, ".rds"))
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
  }
  
  data <- data.frame(data)
  d <- rowSums(data) # the total richness for each grid cell
  model_guild[,j] <- d
}

future_richness_rcp245_guild=model_guild

save(future_richness_rcp245_guild,file="future_richness_rcp245_guild.RData")

#(6)# when all the guilds were combined

data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  if(i%in%all_numbers_rcp245)
  {
    setwd("/home/wenqil/results")
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_future", i, ".rds"))
    
    data[, i] <- raster::extract(df@binary, coords_present)
  }
  
  else {
    setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP245")
    
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    df <- readRDS(paste0("result_future_rcp245—ensemble", i, ".rds"))
    data[, i] <- raster::extract(df@binary, coords_present)
  }
  
}


#(7)# to check the completeness of the molded species for diversity projections based on current climate scenarios

setwd("/Volumes/seas-zhukai/proj-soil-fungi")
directory_path <- "/Volumes/seas-zhukai/proj-soil-fungi/SDM-present"

## get the lacking number

files <- list.files(path = directory_path , pattern = "\\.rds$", full.names = TRUE)

numbers <- files %>%
  basename() %>%
  str_extract("\\d+") %>%
  as.numeric()

# Generate the complete sequence of numbers
complete_sequence <- seq(1:8597)

# Identify the missing numbers
missing_numbers <- setdiff(complete_sequence, numbers)

# Get file information
file_info <- file.info(files)

# Filter files smaller than100 bytes
small_files <- rownames(file_info[file_info$size < 100, ])

ddk=small_files %>%str_extract_all("\\d+") 

small_number =numeric()
for (i in 1:length(small_files))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  small_number[i]=ddk[[i]][3]%>%as.numeric()
}

all_numbers_present=unique(c(missing_numbers, small_number))%>%sort()


#(8)# get the grid-cell-level whole-community fungal diversity

load(file="coords_present_new.RData")

data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  if(i%in%all_numbers_present)
  {
    setwd("/Volumes/seas-zhukai/proj-soil-fungi/results")
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_", i, ".rds"))
    
    data[, i] <- raster::extract(df@binary, coords_present)
    
  }
  
  else {
    setwd("/Volumes/seas-zhukai/proj-soil-fungi/SDM-present")
    
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_present_glm_ensemble", i, ".rds"))
    if(is.null(df))
    {
      setwd("/Volumes/seas-zhukai/proj-soil-fungi/results")
      cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      df <- readRDS(paste0("result_", i, ".rds"))
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    else{
      data[, i] <- raster::extract(df@binary, coords_present)
    }
  }
}

data <- data.frame(data)
d <- rowSums(data) # the total richness for each grid cell for the current climate

present_richness_all=bind_cols(coords_present, d)

save(present_richness_all,file="present_richness_all.RData")


#(9)# get the variable importance for individual guilds

model_guild=list()
for (j in 1:8)
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  #cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  species_otu=op%>%mutate(number=1:8597)%>%filter(ta2==guild_model[j])%>%pull(number)%>%as.character()
  species_number=op%>%mutate(number=1:8597)%>%filter(ta2==guild_model[j])%>%pull(number)%>%as.character()%>%length()
  
  data <- matrix(ncol =7,nrow=species_number)
  
  for (i in seq_along(species_otu)) {
    #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    if(i%in%all_numbers_present)
    {
      setwd("/home/wenqil/results")
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      df <- readRDS(paste0("result_", i, ".rds"))
      
      data[i, ] <- df@variable.importance%>%as.numeric()
    }
    
    else {
      setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-present")
      
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      df <- readRDS(paste0("result_present_glm_ensemble", i, ".rds"))
      data[i, ] <- df@variable.importance%>%as.numeric()
    }
    
  }
  
  data <- data.frame(data,guild_model[j])
  
  # the total richness for each grid cell
  model_guild[[j]] <- data
}

#(10)# bind the importance for all the eight fungal guilds

guild_importance=bind_rows(model_guild[[1]],
                           model_guild[[2]],
                           model_guild[[3]],
                           model_guild[[4]],
                           model_guild[[5]],
                           model_guild[[6]],
                           model_guild[[7]],
                           model_guild[[8]])%>%rename_all(~paste0(c(df@variable.importance,"guild")))
saveRDS(guild_importance, file = "guild_importance.rds")

guild_importance%>%group_by(guild)%>%summarise(mean1=mean(mat_celsius, na.rm = TRUE),
                                mean2=mean(temp_seasonality, na.rm = TRUE),
                                mean3=mean(map_mm, na.rm = TRUE),
                                mean4=mean(soilInCaClpH, na.rm = TRUE),
                                mean5=mean(organicCPercent, na.rm = TRUE),
                                mean6=mean(mat_celsius_2, na.rm = TRUE),
                                mean7=mean(map_mm_2, na.rm = TRUE))%>%
  rename_all(~paste0(c("guild","mat_celsius",
                         "temp_seasonality",
                         "map_mm",
                         "soilInCaClpH",
                         "organicCPercent",
                         "mat_celsius_2",
                         "map_mm_2")))%>%melt()->d_importance

guild_importance%>%group_by(guild)%>%summarise(mean1=sd(mat_celsius, na.rm = TRUE),
                                mean2=sd(temp_seasonality, na.rm = TRUE),
                                mean3=sd(map_mm, na.rm = TRUE),
                                mean4=sd(soilInCaClpH, na.rm = TRUE),
                                mean5=sd(organicCPercent, na.rm = TRUE),
                                mean6=sd(mat_celsius_2, na.rm = TRUE),
                                mean7=sd(map_mm_2, na.rm = TRUE))%>%
  rename_all(~paste0(c("guild","mat_celsius",
                         "temp_seasonality",
                         "map_mm",
                         "soilInCaClpH",
                         "organicCPercent",
                         "mat_celsius_2",
                         "map_mm_2")))%>%melt()->d_importance_sd

d_importance=bind_cols(d_importance,d_importance_sd%>%
                         dplyr::select(value))%>%
  rename_all(~paste0(c("guild","variable","mean_value","sd_value")))

d_importance$variable=factor(d_importance$variable,levels = c("temp_seasonality","soilInCaClpH",
                                                              "mat_celsius","organicCPercent","mat_celsius_2","map_mm_2","map_mm"))

