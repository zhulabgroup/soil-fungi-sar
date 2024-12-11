##for the hotspot analysis

species_change_land_rcp585=readRDS("species_change_land_rcp585.rds")

species_change_climate_rcp245=readRDS("species_change_climate_rcp245.rds")
  load("coords_present_new.RData")
  
 #funtion to compute the moran index 
  my_function_moran=function(data)
  {
    data%>%filter(variable=="all")%>%
      bind_cols(coords_present)%>%filter(!is.na(value))->df
    sf_df <- st_as_sf(df, coords = c("x", "y"), crs = 4326)
    threshold <- 2
    neighbors <- dnearneigh(st_coordinates(sf_df), 0, threshold)
    nb_listw <- nb2listw(neighbors)
    
    local_moran<- localmoran(sf_df$value, nb_listw)
    
    sf_df$local_moran_I <- local_moran[, 1]
    sf_df$local_moran_pvalue <- local_moran[, 5]
    
    return(sf_df)
    
  }
  
  dk=my_function_moran(species_change_land_rcp245)
  
  
  species_change_climate_rcp245%>%filter(variable=="all")%>%
    bind_cols(coords_present)%>%filter(!is.na(value))->df
  
  # select part of the data
  df=df[1:2000,]
  
  sf_df_land <- st_as_sf(df_land, coords = c("x", "y"), crs = 4326)
  
  
  sf_df <- st_as_sf(df, coords = c("x", "y"), crs = 4326)
  
  threshold <- 2
  
  neighbors <- dnearneigh(st_coordinates(sf_df), 0, threshold)
  
  neighbors_land <- dnearneigh(st_coordinates(sf_df_land), 0, threshold)
  
  
  nb_listw_land <- nb2listw(neighbors_land)
  
  nb_listw <- nb2listw(neighbors)
  
  # Apply Local Moran's I test
  local_moran_land <- localmoran(sf_df_land$value, nb_listw_land)
  local_moran_land <- localmoran(sf_df_land$value, nb_listw_land)
  
  
  sf_df_land$local_moran_I <- local_moran_land[, 1]
  sf_df_land$local_moran_pvalue <- local_moran_land[, 5]
  
  plot(st_geometry(sf_df_land), col = ifelse(sf_df_land$local_moran_I > 0, "red", "blue"), pch = 16, cex = .5)
  
  
  
  # Add results to the dataframe
  sf_df$local_moran_I <- local_moran[, 1]
  sf_df$local_moran_pvalue <- local_moran[, 5]
  
  
  plot(st_geometry(dk), col = ifelse(dk$local_moran_I > 0, "red", "blue"), pch = 16, cex = .5)
  
  
  sf_df%>%filter(local_moran_pvalue<0.05&value<0)->dk
  # this mean that for species gain, areas showing high species gain
  
  plot(st_geometry(dk), col = ifelse(dk$local_moran_I > 0, "brown", "blue"), pch = 21, cex = .05)
  
  
  dk <- st_set_crs(dk, 4326) 
  
  my_dk_equal_area <- sf::st_transform(dk,  5070)
  
  ggplot() +
    geom_sf(data=my_dk_equal_area , size=0.5,aes(fill=local_moran_I))
  
  
  