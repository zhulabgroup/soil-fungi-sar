
#diversity data for both whole-community and individual guilds

load(file="full_dob_neon_richness_data.rds")

full_dob_neon_richness_data_rarefaction=full_dob_neon_richness_data

# estimate the z values for each guild and plot

guild_types=unique(full_dob_neon_richness_data_rarefaction$guild)

parameter_list=list()

for (i in 1:9)#each list stores one fungal guild
{
  
  data=full_dob_neon_richness_data_rarefaction%>%filter(guild==guild_types[i]&mean_value>0)# excluded the plots with zeros
  plot_id=unique(data$plotid)
  parameter=data.frame(ncol=6,nrow=length(plot_id))
  for (j in 1:length(plot_id))
  {
    cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
    data_select=data%>%filter(plotid==plot_id[j])# select a plot to calculate the z value
    
    point_number=dim(data_select)[1]
    if(point_number>3)# all plots have the same number of four points
    {
      mod=lm(log(mean_value)~log(area),data=data_select)
      model_summary=summary(mod)
      p_value <- model_summary$coefficients[2, 4]
      parameter[j,1]=coef(mod)[1]%>%as.numeric()
      parameter[j,2]=coef(mod)[2]%>%as.numeric()
      parameter[j,3]=p_value <- model_summary$coefficients[2, 4]
      parameter[j,4]=point_number
      parameter[j,5]=plot_id[j]
      parameter[j,6]=guild_types[i]
      
    }
    else 
    {
      parameter[j,1]=NA
      parameter[j,2]=NA
      parameter[j,3]=NA
      parameter[j,4]=point_number
      parameter[j,5]=plot_id[j]
      parameter[j,6]=guild_types[i]
    }
  }
  
  parameter_list[[i]]=parameter
}

#initial check suggests that the number of points included in the regression affects the estimated z parameters

p_value=numeric()
for (i in 1:9)
{
  mod=lm(nrow~V4,data=parameter_list[[i]]%>%filter(!is.na(nrow))) 
  model_summary=summary(mod)
  
  p_value[i]=model_summary$coefficients[2, 4]
}
# all are significant so we need to control for the number of point per plot

# combine the data for all the guilds

do.call(rbind,parameter_list)%>%
  rename_all(~paste0(c("logc","zvalue","pvale","point_number","plotID", "guild")))->full_parameter_data

saveRDS(full_parameter_data,file="full_parameter_data.rds")

