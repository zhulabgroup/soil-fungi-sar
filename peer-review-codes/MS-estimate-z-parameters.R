
saveRDS(full_dob_neon_richness_data_rarefaction,file="full_dob_neon_richness_data.rds")


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

# test if the number of points and the z values are correlated for different guilds and
p_value=numeric()
for (i in 1:9)
{
  mod=lm(nrow~V4,data=parameter_list[[i]]) 
  model_summary=summary(mod)
  
  p_value[i]=model_summary$coefficients[2, 4]
}
# all are significant so we need to control for the number of point per plot

# combine the data for all the guilds

bind_rows(parameter_list[[1]],
          parameter_list[[2]],
          parameter_list[[3]],
          parameter_list[[4]],
          parameter_list[[5]],
          parameter_list[[6]],
          parameter_list[[7]],
          parameter_list[[8]],# all the z and c values estimated for different guilds
          parameter_list[[9]])%>%
  rename_all(~paste0(c("logc","zvalue","pvale","point_number","plotID", "guild")))->full_parameter_data

saveRDS(full_parameter_data,file="full_parameter_data.rds")

## to see the correlation between the initial z and the rarefaction z for the dob sites

## 
### create several plots to show the linear regression between area and richness

guild_types=unique(full_dob_neon_richness_data_rarefaction$guild)
# when i=1

plot_id=unique(full_neon_data_with_sd$plotid)

pp=list()
for (j in 1:length(plot_id))
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  data_select=full_neon_data_with_sd%>%filter(plotid==plot_id[j])# select a plot to calculate the z value
  
  point_number=dim(data_select)[1]
  if(point_number>3)# all plots have the same number of four points
  {
    pp[[j]]=
      ggplot(data=data_select,aes(y=log(mean_value),x=log(area)))+
      geom_point(size=3)+
      
      xlab(expression("ln(A)"*m^2))+
      ylab("ln(S)")+
      theme(legend.position = c(0.75,0.28),
            legend.text = element_text(size=8),
            legend.title  = element_text(size=10),
            text = element_text(size = 18),
            plot.title = element_text(size = 12, hjust = 0.5), 
            axis.text.y = element_text(hjust = 0,size=10), 
            axis.text.x = element_text(hjust = 1,size=10), 
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15), 
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
            panel.background = element_rect(fill = "NA"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))+
      geom_smooth(method = "lm",color="blue")
    
    
  }
  else 
  {
    pp[[j]]=NULL
  }
}

for(j in 1:463){
  print(pp[[j]])
}



### to create the plots with sd for the neon sites

full_dob_neon_richness_data_rarefaction=full_dob_neon_richness_data
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar.RData")

full_dob_neon_richness_data_rarefaction%>%filter(area!=100&guild=="all"&!is.na(sd_value))%>%dplyr::select(plotid, mean_value, sd_value,  area)%>%
  bind_rows(richness_subplot10_neon_standar)->full_neon_data_with_sd
###

plot_id=unique(data_with_sd_select$plotid)

pp=list()
for (j in 1:length(plot_id))
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  data_select=full_neon_data_with_sd%>%filter(plotid==plot_id[j])# select a plot to calculate the z value
  
  point_number=dim(data_select)[1]
  if(point_number>3)# all plots have the same number of four points
  {
    pp[[j]]=
      ggplot(data=data_select,aes(y=mean_value,x=area))+
      geom_point(size=3,alpha=0.5)+
      geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.5)+
      xlab(expression("A"*m^2))+
      ylab("S")+
      theme(legend.position = c(0.75,0.28),
            legend.text = element_text(size=8),
            legend.title  = element_text(size=10),
            text = element_text(size = 18),
            plot.title = element_text(size = 12, hjust = 0.5), 
            axis.text.y = element_text(hjust = 0,size=10), 
            axis.text.x = element_text(hjust = 1,size=10), 
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15), 
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0.2, 0, 0, 0.2), "cm"),
            panel.background = element_rect(fill = "NA"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))+
      ylim(0,2000)
    
    
    
  }
  else 
  {
    pp[[j]]=NULL
  }
}



wrap_plots(pp[[10]],pp[[14]],pp[[15]],pp1[[10]],
           pp1[[14]],pp1[[15]],
           ncol=3)

# create a figure to show the relationship between area and richness in a log space
pp1=list()
for (j in 1:length(plot_id))
{
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  data_select=full_neon_data_with_sd%>%filter(plotid==plot_id[j])# select a plot to calculate the z value
  
  point_number=dim(data_select)[1]
  if(point_number>3)# all plots have the same number of four points
  {
    pp1[[j]]=
      ggplot(data=data_select,aes(y=log(mean_value),x=log(area)))+
      geom_point(size=3)+
      #geom_errorbar(aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value),width=0.5)+
      xlab(expression("ln(A)"*m^2))+
      ylab("ln(S)")+
      theme(legend.position = c(0.75,0.28),
            legend.text = element_text(size=8),
            legend.title  = element_text(size=10),
            text = element_text(size = 18),
            plot.title = element_text(size = 12, hjust = 0.5), 
            axis.text.y = element_text(hjust = 0,size=10), 
            axis.text.x = element_text(hjust = 1,size=10), 
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15), 
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0.2, 0, 0, 0.2), "cm"),
            panel.background = element_rect(fill = "NA"),
            panel.border = element_rect(color = "black", size = 1, fill = NA))+
      ylim(4,8)+
      geom_smooth(method="lm",color="blue")
  }
  else 
  {
    pp1[[j]]=NULL
  }
}
## to show the generall variaion of the z values

ggplot(full_parameter_data%>%filter(guild=="all"&!is.na(zvalue)), aes(x = zvalue)) +
  geom_histogram(binwidth = 0.02, fill = "#FFCB05", color = "black", alpha = 0.7) +
  ggtitle("") +
  xlab(expression(italic(z)*" value")) +
  ylab("Frequency") +
  geom_vline(xintercept = 0.713,linetype="dashed",color="red")+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  annotate("text",x=0.5,y=50,label="CV=11%",size=5)

## create different slopes

full_parameter_data%>%filter(point_number>3&guild=="all")%>%distinct(plotID)%>%pull(plotID) ->plotid_with_four_point


full_neon_data_with_sd%>%group_by(plotid)%>%summarize(n())%>%filter(`n()`>3)%>%pull(plotid)->plotid_with_four_point

set.seed(125)
rand_plot=sample(plotid_with_four_point,10,replace = FALSE)
full_neon_data_with_sd%>%filter(plotid%in%rand_plot&!is.na(mean_value) )->plot_data


ggplot(data=plot_data,aes(x=log(area),y=log(mean_value)))+
  geom_point(data=plot_data,aes(x=log(area),y=log(mean_value),color=plotid),size=3)+
  geom_smooth(data=plot_data,aes(x=log(area),y=log(mean_value),color=plotid),method="lm",se=FALSE)+
  theme(legend.position ="bottom",
        legend.text = element_text(size=10),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 2, 0, 1), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  scale_color_manual("PlotID",breaks=unique(plot_data$plotid),
                     labels=c(expression("Central Plains Experimental Range "*italic( z)*"=0.71"),
                              expression("Guanica Dry Forest Reserve "*italic(z )*"=0.77"),
                              expression("Harvard Forest "*italic(z)*"=0.70"),
                              expression("Niwot Ridge Mountain Research Station "*italic(z)*"=0.78"),
                              expression("Oak Ridge National Laboratory "*italic(z)*"=0.74"),
                              expression("Ordway-Swisher Biological Station "*italic(z)*"=0.82"),
                              expression("The San Joaquin Experimental range "*italic(z)*"=0.69"),
                              expression("The Talladega National Forest "*italic(z)*"=0.81"),
                              expression("University of Kansas Field Station "*italic(z)*"=0.76"),
                              expression("Woodworth "*italic(z)*"=0.67")),
                     values = c("chocolate1", "#037f77", "royalblue", "#f0a73a", "forestgreen", "#7c1a97","#c94e65", "tan","pink","gray"))+
  xlab(expression(ln(Area)*m^2))+
  ylab("ln(Richness)")+
  guides(color = guide_legend(nrow = 10, byrow = TRUE))


## with sd added

p1=ggplot(data=plot_data,aes(x=area,y=mean_value))+
  geom_point(data=plot_data,aes(x=area,y=mean_value,color=plotid),size=3)+
  geom_errorbar(data=plot_data,aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value,color=plotid),width=45)+
  geom_smooth(data=plot_data,aes(x=area,y=mean_value,color=plotid),method="lm",se=FALSE)+
  theme(legend.position ="right",
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=10), 
        axis.text.x = element_text(hjust = 1,size=10), 
        axis.title.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  scale_color_manual("PlotID",breaks=unique(plot_data$plotid),labels=unique(plot_data$plotid),
                     values = c("chocolate1", "#037f77", "royalblue", "#f0a73a", "forestgreen", "#7c1a97","#c94e65", "tan","pink","gray"))+
  xlab(expression("Area "*m^2))+
  ylab("Richness")+
  guides(color=FALSE)



p1=
  
  
  dd=numeric() 
for (i in seq_along(rand_plot)){
  d=summary(lm(log(mean_value)~log(area),data=plot_data%>%filter(plotid==rand_plot[i])))
  dd[i]=d$coefficients[2,1]
}


