##


directory_path <- "/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP245"

directory_path <- "/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP585"


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


all_numbers_rcp585=missing_numbers
all_numbers_rcp245=all_numbers

## get the future species composition

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











data <- matrix(ncol = 8597, nrow = 275760)
for (i in 1:8597) {
  if(i%in%all_numbers_rcp585)
  {
    setwd("/home/wenqil/results")
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    df <- readRDS(paste0("result_future", i, ".rds"))
    
    data[, i] <- raster::extract(df@binary, coords_present)
  }
  
  else {
    setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/SDM-RCP585")
    
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    df <- readRDS(paste0("result_future_rcp585_ensemble", i, ".rds"))
    data[, i] <- raster::extract(df@binary, coords_present)
  }
  
}

###

setwd("/nfs/turbo/seas-zhukai/proj-soil-fungi/land-use-climate-historical")

grid_biome=readRDS("grid_level_biomes.rds")

species_composition_rcp585=data%>%data.frame()

species_composition_rcp585=species_composition_rcp585%>%
  bind_cols(grid_biome)%>%filter(LABEL%in%c("Tropical & Subtropical Moist Broadleaf Forests",
                                            "Temperate Conifer Forests",
                                            "Temperate Grasslands, Savannas & Shrublands",
                                            "Temperate Broadleaf & Mixed Forests" ))
# for the present diversity
species_composition_present=species_composition_present%>%
  bind_cols(grid_biome)%>%filter(LABEL%in%c("Tropical & Subtropical Moist Broadleaf Forests",
                                            "Temperate Conifer Forests",
                                            "Temperate Grasslands, Savannas & Shrublands",
                                            "Temperate Broadleaf & Mixed Forests" ))






bray_dist=numeric()
for (i in 1:42048)
{
com1=species_composition_present[i,1:8597]
com1=species_composition_rcp585[i,1:8597]
combined_comm=rbind(com1,com2)
bray_dist[i]=vegdist(combined_comm, method = "bray")
}


