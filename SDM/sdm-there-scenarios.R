##

ggplot()+
  geom_point(data=d1,aes(x=gx,y=gy))+
  geom_text(data=d1,aes(x=gx,y=gy),label=d1$subplotID5 ,size=2)+
  #geom_vline(xintercept = 10,color="red")+
  
  geom_vline(xintercept = 20,color="red")+
  #geom_vline(xintercept = 30,color="red")+
  #geom_vline(xintercept = 40,color="red")+
  #geom_hline(yintercept = 10,color="red")+
  geom_hline(yintercept = 20,color="red")+
  #geom_hline(yintercept = 30,color="red")+
  geom_hline(yintercept = 40,color="red")

a5=unique(d$subplotID20)# means 20 

a6=numeric()
for (i in 1:length(a5)){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
data_sub <- subset_samples(rare_all_assign, subplotID20==a5[i])
a6[i]=dim(sample_data(data_sub))[1]
}
### for the rcp245 scenario

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


###3

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
####



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




cl <- makeCluster(10)
registerDoParallel(cl)

sdm_present=list()
for(i in 6078:6500){
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


####
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



cl <- makeCluster(10)
registerDoParallel(cl)

sdm_future_rcp858_ensemble=list()
for(i in 4000:8597){
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


### to extract the data based 

directory_path <- "/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP245"
# the above codes do not work

directory_path <- "/Volumes/seas-zhukai/proj-soil-fungi/SDM-RCP245"

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
all_numbers_rcp245=all_numbers


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


###

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

## when all the guilds were combined

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



#### to check the completeness of the molded species

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


## when all the guilds were combined


data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  if(i%in%all_numbers_present)
  {
    setwd("/home/wenqil/results")
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_", i, ".rds"))
    
    data[, i] <- raster::extract(df@binary, coords_present)
    
  }
  
  else {
    setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-present")
    
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_present_glm_ensemble", i, ".rds"))
    if(is.null(df))
    {
      setwd("/home/wenqil/results")
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
d <- rowSums(data) # the total richness for each grid cell



present_richness_all=bind_cols(coords_present, d)

save(present_richness_all,file=".present_richness_allRData")

##
setwd("/home/wenqil/results")
data <- matrix(ncol = 8597, nrow = 275760)
for (i in 3600:8597) {
  
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  df <- readRDS(paste0("result_", i, ".rds"))
  
  data[, i] <- raster::extract(df@binary, coords_present)
}



## when all the guilds were combined
# Set the directory path
dir_path <- "/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-present"

# Get the list of files in the directory
files <- list.files(dir_path, full.names = TRUE)

# Initialize a vector to store filenames with class NULL
null_files <- c()

# Loop through each file
for (file in files) {
  # Try to read the file (assuming text file for this example)
  file_content <- tryCatch(readLines(file), error = function(e) NULL)
  
  # Check if the class of file_content is NULL
  if (is.null(file_content)) {
    # If class is NULL, add the filename to null_files
    null_files <- c(null_files, file)
  }
}

# Print the filenames with class NULL
print(null_files)




###


# for i =3666, an error occur


df<- ensemble_modelling(algorithms = c('MAXENT', 'RF',"GLM"), 
                        subset(occurrence_data,variable==sp[i]),
                        r_present_northam, Xcol = 'lon', Ycol = 'lat',Pcol="value", 
                        Spcol="variable", rep = 10,  # Number of repetitions
                        cv = 'holdout',  # Cross-validation method
                        cv.param = c(0.7, 1),  # Parameters for cross-validation
                        # Metric for evaluating models
                        ensemble.thresh = 0.75,  # Threshold for including models in the ensemble
                        weight = TRUE,  # Whether to weight models by their performance
                        
                        verbose = FALSE)

saveRDS(sdm_present[[i]],file=paste0("result_present_glm_ensemble", i, ".rds"))





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


###

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
    
    if(i%in%all_numbers_present)
    {
      setwd("/home/wenqil/results")
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      df <- readRDS(paste0("result_", i, ".rds"))
      
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
    else {
      setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-present")
      
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      df <- readRDS(paste0("result_present_glm_ensemble", i, ".rds"))
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
  }
  
  data <- data.frame(data)
  
  d <- rowSums(data) # the total richness for each grid cell
  model_guild[,j] <- d
}


## read in the data

variable_importance_present=load("variable_importance_present.RData")
#need to open this file in the repository
k=load("/Volumes/seas-zhukai/proj-soil-fungi/variable_importance_present.RData")

variable_importance_present%>%melt()%>%group_by(variable)%>%summarise(mean_value=mean(value),sd_value=sd(value))->
  variable_importance_present_data

variable_importance_present_data$variable=factor(variable_importance_present_data$variable,
                                                 levels=c( "temp_seasonality", "soilInCaClpH" , "organicCPercent" , "mat_celsius_2" ,   "mat_celsius"  ,   
                                                            "map_mm"  , "map_mm_2"  ))
ggplot(variable_importance_present_data, aes(x = variable, y = mean_value)) +
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.2)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10,angle=60), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  ylab("Importance")+
  xlab("")+
  scale_x_discrete(breaks=c("temp_seasonality", "soilInCaClpH" ,    "organicCPercent",  "mat_celsius_2",    "mat_celsius"  ,   
                           "map_mm", "map_mm_2" ),labels=c("Temp.seas.","pH","SoilC",expression(MAT^2),"MAT","MAP",expression(MAP^2)))
  

model_evaluation_present=readRDS("model_evaluation_present.rds")

model_evaluation_present%>%melt()%>%group_by(variable)%>%summarise(mean_value=mean(value),sd_value=sd(value))->
  model_evaluation_present_data


ggplot(model_evaluation_present_data%>%filter(variable%in%c("AUC","sensitivity","specificity","prop.correct")), aes(x = variable, y = mean_value)) +
  geom_bar(stat = "identity",width=0.3)+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.2)+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10,angle=60), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("")

##get the variable importance for individual guilds for the current presence

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
    
    if(i%in%all_numbers_present)
    {
      setwd("/home/wenqil/results")
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      df <- readRDS(paste0("result_", i, ".rds"))
      
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
    else {
      setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-present")
      
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      df <- readRDS(paste0("result_present_glm_ensemble", i, ".rds"))
      data[, i] <- raster::extract(df@binary, coords_present)
    }
    
  }
  
  data <- data.frame(data)
  
  d <- rowSums(data) # the total richness for each grid cell
  model_guild[,j] <- d
}




###
tax_table(rare_all_guild)%>%data.frame()%>%dplyr::select(primary_lifestyle)->df
df%>%mutate(sp=rownames(df))->species_guild
sp%>%data.frame()%>%rename_all(~paste0("sp"))%>%left_join(species_guild,by="sp")->op
op%>%mutate(ta2 = ifelse(str_detect(primary_lifestyle, "parasite"), "para", primary_lifestyle))->op
guild_model=c("soil_saprotroph","plant_pathogen","ectomycorrhizal","litter_saprotroph","wood_saprotroph","arbuscular_mycorrhizal","epiphyte","para")


## the importance of the variables


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

#bind the importance for all the eight fungal guilds

guild_importance=bind_rows(model_guild[[1]],model_guild[[2]],model_guild[[3]],model_guild[[4]],
                           model_guild[[5]],model_guild[[6]],model_guild[[7]],model_guild[[8]])%>%rename_all(~paste0(c(df@variable.importance,"guild")))
saveRDS(guild_importance, file = "guild_importance.rds")

## to see the impact for individual guilds

d=readRDS("guild_importance.rds")

d%>%group_by(guild)%>%summarise(mean1=mean(mat_celsius, na.rm = TRUE),
                                mean2=mean(temp_seasonality, na.rm = TRUE),
                                mean3=mean(map_mm, na.rm = TRUE),
                                mean4=mean(soilInCaClpH, na.rm = TRUE),
                                mean5=mean(organicCPercent, na.rm = TRUE),
                                mean6=mean(mat_celsius_2, na.rm = TRUE),
                                mean7=mean(map_mm_2, na.rm = TRUE)
                                )%>%rename_all(~paste0(c("guild","mat_celsius",
                                                         "temp_seasonality",
                                                         "map_mm",
                                                         "soilInCaClpH",
                                                         "organicCPercent",
                                                         "mat_celsius_2",
                                                         "map_mm_2")))%>%melt()->d_importance
d%>%group_by(guild)%>%summarise(mean1=sd(mat_celsius, na.rm = TRUE),
                                mean2=sd(temp_seasonality, na.rm = TRUE),
                                mean3=sd(map_mm, na.rm = TRUE),
                                mean4=sd(soilInCaClpH, na.rm = TRUE),
                                mean5=sd(organicCPercent, na.rm = TRUE),
                                mean6=sd(mat_celsius_2, na.rm = TRUE),
                                mean7=sd(map_mm_2, na.rm = TRUE)
)%>%rename_all(~paste0(c("guild","mat_celsius",
                         "temp_seasonality",
                         "map_mm",
                         "soilInCaClpH",
                         "organicCPercent",
                         "mat_celsius_2",
                         "map_mm_2")))%>%melt()->d_importance_sd

d_importance=bind_cols(d_importance,d_importance_sd%>%dplyr::select(value))%>%
  rename_all(~paste0(c("guild","variable","mean_value","sd_value")))
                                                  
d_importance$variable=factor(d_importance$variable,levels = c("temp_seasonality","soilInCaClpH",
                                                              "mat_celsius","organicCPercent","mat_celsius_2","map_mm_2","map_mm"))

ggplot(data=d_importance, aes(x = guild, y = mean_value, fill = variable)) +
  geom_bar(stat = "identity",color="black",width = 0.8,position = "dodge")+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.5,position = position_dodge(0.8))+
  scale_fill_manual("variable",breaks=c("temp_seasonality","soilInCaClpH",
                                        "mat_celsius","organicCPercent","mat_celsius_2","map_mm_2","map_mm"),
                    labels=c("Temp.seas.","pH", "MAT","SoilC",expression(MAT^2),expression(MAP^2),"MAP"),
                                        values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan"))+
  theme(legend.position ="right",
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10,angle=60), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("Importance(%)")+
  scale_x_discrete(breaks=unique(d_importance$guild),
                   labels = c("AM","ECM","Epiphyte","Litter sapro.","Parasite","Plant patho.",
                              "Soil sapro.","Wood sapro." ))
  
  #just select the four guilds included in the main text

d_importance%>%filter(guild%in%c("arbuscular_mycorrhizal", "ectomycorrhizal" , 
                                 "plant_pathogen" ,        "soil_saprotroph"   ))->d_importance_sub



ggplot(data=d_importance_sub, aes(x = guild, y = mean_value, fill = variable)) +
  geom_bar(stat = "identity",color="black",width = 0.8,position = "dodge")+
  geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.5,position = position_dodge(0.8))+
  scale_fill_manual("Variables",breaks=c("temp_seasonality","soilInCaClpH",
                                        "mat_celsius","organicCPercent","mat_celsius_2","map_mm_2","map_mm"),
                    labels=c("Temp.seas.","pH", "MAT","SoilC",expression(MAT^2),expression(MAP^2),"MAP"),
                    values=c("#c94e65","#037f77","royalblue","forestgreen","chocolate1","#7c1a97","tan"))+
  theme(legend.position ="right",
        legend.text = element_text(size=10),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 0.5,size=15,angle=0,color="black"), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  xlab("")+
  ylab("Importance(%)")+
  scale_x_discrete(breaks=unique(d_importance_sub$guild),
                   labels = c("AM","ECM","Plant pathogens",
                              "Soil saprotrophs" ))

#just select the four guilds included in the main text

d_importance%>%filter(guild%in%c("arbuscular_mycorrhizal", "ectomycorrhizal" , 
                                 "plant_pathogen" ,        "soil_saprotroph"   ))->d_importance_sub






d%>%filter(guild=="soil_saprotroph")%>%summarise(mean1=mean(mat_celsius, na.rm = TRUE),
                                mean2=mean(temp_seasonality, na.rm = TRUE),
                                mean3=mean(map_mm, na.rm = TRUE),
                                mean4=mean(soilInCaClpH, na.rm = TRUE),
                                mean5=mean(organicCPercent, na.rm = TRUE),
                                mean6=mean(mat_celsius_2, na.rm = TRUE),
                                mean7=mean(map_mm_2, na.rm = TRUE)
)

  