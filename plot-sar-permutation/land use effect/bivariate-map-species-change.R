## to create a bivariate map for both species loss and gains for climate and land use effects


species_change_climate_rcp245=readRDS(file="species_change_climate_rcp245.rds")

species_change_climate_rcp585=readRDS(file="species_change_climate_rcp585.rds")

grid_level_biomes=readRDS(file="grid_level_biomes.rds")

species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value<0)%>%pull(value)%>%summary()->df1

species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value>=0)%>%pull(value)%>%summary()->df2


# assign the value to each cell based on the climate effect

species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df3

df3%>%mutate(type=case_when(value>=df1[1]&value<df1[2]~"1",
                            value>=df1[2]&value<df1[3]~"2",
                            value>=df1[3]&value<df1[4]~"3",
                            value>=df1[4]&value<df1[5]~"4",
                            value>=df1[5]&value<df1[6]~"5",
                            
                            value>=df2[1]&value<df2[2]~"6",
                            value>=df2[2]&value<df2[3]~"7",
                            value>=df2[3]&value<df2[4]~"8",
                            value>=df2[4]&value<df2[5]~"9",
                            value>=df2[5]&value<df2[6]~"10",
                            TRUE~"other" ))->df4


# for land use

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value<0)%>%pull(value)%>%summary()->df1

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value>=0)%>%pull(value)%>%summary()->df2


# assign the value to each cell based on the climate effect

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df3

df3%>%mutate(type=case_when(value>=df1[1]&value<df1[2]~"1",
                            value>=df1[2]&value<df1[3]~"2",
                            value>=df1[3]&value<df1[4]~"3",
                            value>=df1[4]&value<df1[5]~"4",
                            value>=df1[5]&value<df1[6]~"5",
                            
                            value>=df2[1]&value<df2[2]~"6",
                            value>=df2[2]&value<df2[3]~"7",
                            value>=df2[3]&value<df2[4]~"8",
                            value>=df2[4]&value<df2[5]~"9",
                            value>=df2[5]&value<df2[6]~"10",
                            TRUE~"other" ))->df5



# combined the climate and land use effects
# for the land use effect, there should be many NAS

bind_cols(df4,df5)->df6
df6%>%dplyr::select(type...5,type...10)%>%
  mutate(type=paste(type...5,"-",type...10))->df7
  
# with different types

coords_present%>%bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df8

# get the combined data
df8%>%bind_cols(df7)%>%
  dplyr::select(lon,lat,type...7)->df9
# need to convert the data 
# remove the cells with 
df9%>%filter(!grepl("other", type...7))->df10

### get different quantile with four
###

my_function_quantile=function(data)
  
  {
  
  data%>%
    filter(variable=="all")%>%
    bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%biome_select&value<0)%>%pull(value)%>%quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df1
  
  data%>%
    filter(variable=="all")%>%
    bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%biome_select&value>=0)%>%pull(value)%>%
    quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df2
  
  data%>%
    filter(variable=="all")%>%
    bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%biome_select)->df3
  
  df3%>%mutate(type=case_when(value>=df1[1]&value<df1[2]~"loss-high",
                              value>=df1[2]&value<df1[3]~"loss-medium",
                              value>=df1[3]&value<df1[4]~"loss-medium",
                              value>=df1[4]&value<df1[5]~"loss-low",
                              value>=df2[1]&value<df2[2]~"gain-low",
                              value>=df2[2]&value<df2[3]~"gain-medium",
                              value>=df2[3]&value<df2[4]~"gain-medium",
                              value>=df2[4]&value<df2[5]~"gain-high",
                              
                              TRUE~"other" ))->df4
return(df4)
}






species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value<0)%>%pull(value)%>%quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df1

species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value>=0)%>%pull(value)%>%
  quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df2


# assign the value to each cell based on the climate effect

species_change_climate_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df3

df3%>%mutate(type=case_when(value>=df1[1]&value<df1[2]~"loss-high",
                            value>=df1[2]&value<df1[3]~"loss-medium",
                            value>=df1[3]&value<df1[4]~"loss-medium",
                            value>=df1[4]&value<df1[5]~"loss-low",
                            
                            
                            value>=df2[1]&value<df2[2]~"gain-low",
                            value>=df2[2]&value<df2[3]~"gain-medium",
                            value>=df2[3]&value<df2[4]~"gain-medium",
                            value>=df2[4]&value<df2[5]~"gain-high",
                          
                            TRUE~"other" ))->df4


# for land use

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value<0)%>%pull(value)%>%quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df1

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select&value>=0)%>%pull(value)%>%quantile( probs = c(0, 0.25, 0.5, 0.75, 1))->df2


# assign the value to each cell based on the climate effect

species_change_land_rcp245%>%
  filter(variable=="all")%>%
  bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df3

df3%>%mutate(type=case_when(value>=df1[1]&value<df1[2]~"loss-high",
                            value>=df1[2]&value<df1[3]~"loss-medium",
                            value>=df1[3]&value<df1[4]~"loss-medium",
                            value>=df1[4]&value<df1[5]~"loss-low",
                            
                            
                            value>=df2[1]&value<=df2[2]~"gain-low",
                            value>=df2[2]&value<df2[3]~"gain-medium",
                            value>=df2[3]&value<df2[4]~"gain-medium",
                            value>=df2[4]&value<df2[5]~"gain-high",
                            
                            TRUE~"other" ))->df5

##bind the two columns

bind_cols(df4,df5)%>%dplyr::select(type...5,type...10)%>%
  mutate(type=paste(type...5,"-",type...10))->df7

# with different types

coords_present%>%bind_cols(grid_level_biomes)%>%
  filter(LABEL%in%biome_select)->df8

# get the combined data
df8%>%bind_cols(df7)%>%
  dplyr::select(lon,lat,type...7)->df9
# need to convert the data 
# remove the cells with 
df9%>%filter(!grepl("other", type...7))->df10


df9%>%mutate(type=case_when(grepl("other", type...7)~"other - other",
                            TRUE~type...7))->df11









d1=my_function_quantile(species_change_climate_rcp245)
d2=my_function_quantile(species_change_land_rcp245)

d3=my_function_quantile(species_change_climate_rcp585)
d4=my_function_quantile(species_change_land_rcp585)

data_quantile_rcp245=bind_cols(d1,d2)

data_quantile_rcp585=bind_cols(d3,d4)


# combined the climate and land use effects
# for the land use effect, there should be many NAS

bind_cols(d1,d2)->df6

data=c("data_quantile_rcp245","data_quantile_rcp585")

data_four_quantile=list()
for(i in 1:2)
{
 
  get(data[i])%>%dplyr::select(type...5,type...10)%>%
    mutate(type=paste(type...5,"-",type...10))->df7
  
  # with different types
  
  coords_present%>%bind_cols(grid_level_biomes)%>%
    filter(LABEL%in%biome_select)->df8
  
  # get the combined data
  df8%>%bind_cols(df7)%>%
    dplyr::select(lon,lat,type...7)->df9
  # need to convert the data 
  # remove the cells with 
  df9%>%filter(!grepl("other", type...7))->df10
  
  
  #df9%>%mutate(type=case_when(grepl("other", type...7)~"other - other",
                             # TRUE~type...7))->df11
  
  expand.grid(c("loss-high","loss-medium","loss-low","gain-low","gain-medium","gain-high"),
              c("loss-high","loss-medium","loss-low","gain-low","gain-medium","gain-high"))->d
  
  paste(d$Var1,"-",d$Var2)->d
  
  df10$type=factor(df10$type,levels =c(d,"other - other"))
  
  df5=my_function_project(df10)
  
  df5$type=factor(df5$type)
  data_four_quantile[[i]]=df5
  
}







c("#001f3f", "#FF7F50", "#008080", "#D3D3D3")


set.seed(42)
category=c("Loss-low","Loss-medium","Loss-high","Gain-low","Gain-medium","Gain-high")
data <- expand.grid(
Variable1 = factor(category,levels=c("Loss-high","Loss-medium","Loss-low","Gain-low","Gain-medium","Gain-high")),
  
Variable2 = factor(category,levels=c("Loss-high","Loss-medium","Loss-low","Gain-low","Gain-medium","Gain-high")))

data$Value <- runif(nrow(data))  #

data$Combination <- interaction(data$Variable1, data$Variable2)

data%>%mutate(Combination = gsub("\\.", " - ", Combination))->data

data$Combination=factor(data$Combination,levels=c("loss-low","loss-medium","loss-high","gain-low","gain-medium","gain-high"))

"#FFD700" "#FFAC00" "#FF8100" "#FF5500" "#FF2A00" "#FF0000"
"#CCDF32" "#CCB232" "#CC8532" "#CC5932" "#CC2C32" "#CC0032"



p2=ggplot(data, aes(x = Variable1, y = Variable2, fill = Combination)) +
  geom_tile(color = "gray")+
  scale_fill_manual("",breaks=data$Combination,
                    labels=data$Combination,values=col.matrix[2:7,2:7])+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90,hjust=0,vjust = 0.5),
        axis.text.y = element_text(hjust=0,vjust=0.5)
       )+
  
  ylab("Land use effect")+
  xlab("Climate effect")
  
  theme(legend.position = "none")+
  scale_y_discrete(breaks=1:8,labels=c("-0.22","-0.07","-0.03","-0.01","0","h","h","h"))+
scale_x_discrete(breaks=1:8,labels=c("-0.72--0.12","-0.12-0.06","-0.06-0.03","-0.03-0.00","0-0.03","0.03-0.08","0.08-0.16","0.16-1.12"))

smaller_panel_grob1 <- ggplotGrob(p2)

ggplot(data=df10,aes(x=lon,y=lat,color=type))+
  geom_point()+
  scale_color_manual("",breaks=unique(df10$type)%>%sort(),
                     labels=unique(df10$type)%>%sort(),
                     values=c(col.matrix[2:7,2:7],"white" ))+

  theme(legend.position = "none")
  
  #change the levels
  
  expand.grid(c("loss-high","loss-medium","loss-low","gain-low","gain-medium","gain-high"),
         c("loss-high","loss-medium","loss-low","gain-low","gain-medium","gain-high"))->d
  
  paste(d$Var1,"-",d$Var2)->d
  
  df10$type=factor(df10$type,levels =c(d,"other - other"))

# need to project the data


df5=my_function_project(df10)

df5$type=factor(df5$type)


pp_bivariate=list()
for(i in 1:2)
  {
  pp_bivariate[[i]]=ggplot()+
    geom_point(data=data_four_quantile[[i]],pch=15,aes(x=x,y=y,color=type),size=0.01)+
    scale_color_manual("",breaks=unique(data_four_quantile[[i]]$type)%>%sort()%>%as.factor(),
                       labels=unique(data_four_quantile[[i]]$type)%>%sort()%>%as.factor(),
                       values=c(col.matrix[2:7,2:7],"white" ))+
    geom_sf(data = us_projected, fill = NA, size=0.01,color = "black")+
    geom_sf(data = canada_clipped, fill=NA,size=0.01,color = "black")+
    geom_sf(data = rico_projected, fill = NA, size=0.01,color = "black")+
    geom_sf(data = cuba_projected, fill = NA, size=0.01,color = "black")+
    geom_sf(data = mexico_projected, fill = NA, size=0.01,color = "black")+
    geom_sf(data = haiti_projected, fill = NA, size=0.01,color = "black")+
    geom_sf(data = dominican_projected, fill = NA, size=0.01,color = "black")+
    coord_sf(xlim = c(-5000000 , 3000000), ylim = c(-252303 , 5980000))+
    common_theme+
    guides(color="none")+
    xlab("")+
    ylab("")+
    ggtitle(scenario[[i]])+
    annotation_custom(grob = smaller_panel_grob1, xmin = -5840567 , xmax = -2509970, ymin =-427265.6, ymax = 3012262.7)
  
}


  
####


  stevenpinkblue <- c(
    "#08306B", # Dark Blue
    "#2171B5", # Medium Blue
    "#6BAED6", # Light Blue
    "#FDD0A2", # Light Pink
    "#A50F15"  # Dark Pink
  )

  
  scale_fill_manual("",breaks=c("Negative","Positive"),labels=c("Loss","Gain"),values=c("#E1BE6A", "#40B0A6"

## need to change the color
col.matrix<-colmat(nquantiles=6, bottomleft="gold", 
                   bottomright="red",
                   upperleft="cyan", 
                   upperright="blue",
                   xlab="My x label", ylab="My y label")

# create the grids

col.matrix<-colmat(nquantiles=6, bottomleft="#FDD0A2", 
                   bottomright="#6BAED6",
                   upperright="#08306B",
                   upperleft="#A50F15", 
                   xlab="My x label", ylab="My y label")

## for the high emission scenarios

## to create a bivariate map for both species loss and gains for climate and land use effects

# for the map with data showing proportional cells







# for the bar plot


  
#bind all the plots
  
  pp_bivariate[[1]]=ggplotGrob(pp_bivariate[[1]])
  bar_plot1=ggplotGrob(bar_plot1)
  
  map_direction_overall[[1]]=ggplotGrob( map_direction_overall[[1]])
  
  
  pp_bivariate[[1]]$widths=bar_plot1$widths
  
  map_direction_overall[[1]]$widths=pp_bivariate[[1]]$widths

plot_grid(pp_bivariate[[1]],
          pp_bivariate[[2]],
          map_direction_overall[[1]],
          map_direction_overall[[2]],
          bar_plot1,
          bar_plot2,
          ncol=2,
          label_x = 0.1,label_y = 1.02,label_size = 16,
          labels = paste0("(", letters[1:6], ")"))






# create figures for the four main fungal guilds

plot_grid(map_direction_overall_guild[[1]][[1]],
          map_direction_overall_guild[[1]][[2]],
          map_direction_overall_guild[[2]][[1]],
          map_direction_overall_guild[[2]][[2]],
          map_direction_overall_guild[[3]][[1]],
          map_direction_overall_guild[[3]][[2]],
          map_direction_overall_guild[[6]][[1]],
          map_direction_overall_guild[[6]][[2]],
          label_x = 0.1,
          label_y = 0.85,
          label_size = 16,ncol=4,
          labels = paste0("(", letters[1:8], ")"))

# the final dimension for the output figure is 10 x 20 

plot_grid(biome_bar_plot_low[[1]],
          biome_bar_plot_high[[1]],
          biome_bar_plot_low[[2]],
          biome_bar_plot_high[[2]],
          biome_bar_plot_low[[3]],
          biome_bar_plot_high[[3]],
          biome_bar_plot_low[[6]],
          biome_bar_plot_high[[6]],
          label_x = 0.1,
          label_y = 1.02,
          label_size = 16,ncol=4,
          labels = paste0("(", letters[1:8], ")"))

# the final dimension for the output figure is 10 x 20 



