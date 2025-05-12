
length(unique(get_variable(neon_dob, "Site")))
sites_to_project <- cbind.data.frame(
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(Site) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_site <- merge_samples(neon_dob, group="Site")
sites_df <- as(sample_data(neon_dob_site), "data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_site), sites_to_project$Site)])





theme_custom <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

theme_custom_map <- theme_custom +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggplot(sites_df) +
  geom_sf(data=st_transform(st_as_sf(north_america_cropped),size=0.1, col="black", fill=alpha("white", 0)) )+ 
  #geom_point(aes(x=lon, y=lat, col=Project), shape=21, fill=alpha("white", 0.3), size=1)+ 
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  guides(col="none") +
  theme_custom_map+
  geom_point(data=subset(no_na_prediction[,c("lon","lat","OTU520")],OTU520>0),aes(x=lon,y=lat),color="blue",size=0.15)

# to take the loop to look at the effect



  length(unique(get_variable(neon_dob, "Site")))
sites_to_project <- cbind.data.frame(
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(Site) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_site <- merge_samples(neon_dob, group="Site")
sites_df <- as(sample_data(neon_dob_site), "data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_site), sites_to_project$Site)])



# load the species presence data

# remove all the NAS

predict_presence=allpresence_absence[129:dim(allpresence_absence)[1],]

com_pre=predict_presence[which(complete.cases(grid_coordinates_climate)),]
# check 
# select the no na columns

no_na_prediction <- com_pre[, colSums(is.na(com_pre)) == 0]
