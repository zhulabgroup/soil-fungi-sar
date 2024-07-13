## countryside SAR model

load("~/soil-sar/plot-sar-permutation/PFT_2100_ssp245.RData")

head(PFT_2100)

PFT_2100 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster_2100 <- crop(flip(coarser_raster_2100[[i]]), b)
  coords_present <- xyFromCell(cropped_raster_2100, cell = 1:ncell(cropped_raster_2100)) # get the coordinates
  cell_values <- raster::extract(cropped_raster_2100, coords_present) %>% as.matrix()
  PFT_2100[, i] <- cell_values
}
PFT_2100 <- cbind(coords_present, PFT_2100) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster2)))) 

PFT_2100=PFT_2100%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)

temp=PFT_2100[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2100%>%bind_cols(No_crop)->PFT_2100

# the affinity of the crops to be 0.93
# the species loss rate

PFT_2100%>%dplyr::select(crop,No_crop)%>%mutate(crop_induce_loss=1-((1*No_crop+0.93*crop)/(crop+No_crop))^0.3877299)->temp


# species loss rate

species_change_climate%>%bind_cols(temp$crop_induce_loss)%>%rename_all(~paste0(c("lon","lat","present_rich","future_rich","change","rate", "crop_induce_loss")))->countrysar



ggplot(countrysar) +
  geom_point(data = countrysar, pch = 15, aes(x = lon, y = lat, color = -1*crop_induce_loss* 100), size = 0.275) +

scale_color_gradient2(expression("%"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = -1, na.value = "white")+ 
  xlab("Predicted species loss") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = c(0.2, 0.35),
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size = 8),
        text = element_text(size = 18),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  )+geom_sf(data=st_as_sf(north_america_cropped),
            size=0.1, col="black", fill=alpha("white", 0))


## let's consider different land use types based on the land-use type for two types

calculate_pixel_area <- function(latitude) {
  # Convert latitude to radians
  lat_rad <- latitude * pi / 180
  
  # Earth radius in meters
  earth_radius <- 6371000  # in meters
  
  # Calculate the length of one degree of latitude and longitude at the given latitude
  lat_deg_length <- earth_radius * pi / 180
  lon_deg_length <- earth_radius * cos(lat_rad) * pi / 180
  
  # Calculate the area of a square defined by one degree of latitude and longitude
  square_area <- lat_deg_length * lon_deg_length
  
  # Convert square area to the area of a 10-minute pixel
  pixel_area <- square_area * (10/60)^2
  
  return(pixel_area)
}

## get the latitude of the data

coords_present=data.frame(coords_present)

latitude=coords_present$x*-1

latitude%>%data.frame()%>%rename_all(~paste0(c("latitude")))->latitude

Area_pixel=calculate_pixel_area (latitude)

## get the land use type for different years

PFT_2020 <- matrix(nrow = 275760, ncol = 33)
for (i in 1:33) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  cropped_raster_2020 <- crop(flip(coarser_raster[[i]]), b)
  coords_present <- xyFromCell(cropped_raster_2020, cell = 1:ncell(cropped_raster_2020)) # get the coordinates
  cell_values <- raster::extract(cropped_raster_2020, coords_present) %>% as.matrix()
  PFT_2020[, i] <- cell_values
}

# total cover of non-crop land use types in each cell

PFT_2020 <- cbind(coords_present, PFT_2020) %>%
  data.frame() %>%
  rename_all(~ paste0(c("lon", "lat", names(raster1))))

temp=PFT_2020[,c("PFT1","PFT2","PFT3","PFT4","PFT5","PFT6","PFT7","PFT8","PFT9","PFT10","PFT11","PFT12","PFT13","PFT14")]

No_crop=apply(temp, 1, sum)%>%data.frame()
names(No_crop)="No_crop"
PFT_2020%>%bind_cols(No_crop)->PFT_2020
###

PFT_2020=PFT_2020%>%mutate(tree=PFT1+PFT2+PFT3+PFT4+PFT5+PFT6+PFT7+PFT8,shrup=PFT9+PFT10+PFT11,grass=PFT12+PFT13+PFT14,
                           crop=PFT15+PFT16+PFT17+PFT18+PFT19+PFT20+PFT21+PFT22+PFT23+PFT24+PFT25+PFT26+PFT27+PFT28+PFT29+PFT30)


## include the richness

PFT_2020=PFT_2020%>%mutate(area=Area_pixel$latitude)

## get the  current richness for each cell

PFT_2020%>%mutate(richness=61.65*((No_crop/100*area)+(0.93*crop/100*area))^0.38)->PFT_2020

## future richness

PFT_2100=PFT_2100%>%mutate(area=Area_pixel$latitude)

## get the  current richness for each cell

PFT_2100%>%mutate(richness=61.65*((No_crop/100*area)+(0.93*crop/100*area))^0.38)->PFT_2100

change_richness_land=data.frame(rate=(PFT_2100$richness-PFT_2020$richness)/PFT_2020$richness,PFT_2020[,1:2])


###


ggplot(change_richness_land) +
  geom_point(data = change_richness_land, pch = 15, aes(y = -1*lon, x = lat, color = rate* 100), size = 0.275) +
  
  scale_color_gradient2(expression("%"), low = "seagreen3", mid = "yellow", high = "purple", midpoint = -1, na.value = "white")+ 
  xlab("Predicted species changes") +
  ylab("")+
  ggtitle("RCP 4.5 & SSP 2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = c(0.2, 0.35),
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size = 8),
        text = element_text(size = 18),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  )+geom_sf(data=st_as_sf(north_america_cropped),
            size=0.1, col="black", fill=alpha("white", 0))


