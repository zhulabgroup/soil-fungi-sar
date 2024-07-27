install.packages("rstatix")
install.packages("car")      # For type II and III ANOVA
install.packages("afex")     # For mixed ANOVA and other advanced options
install.packages("dplyr")

                 
library(rstatix)
library(car)
library(afex)
library(dplyr)     

# to see the effect of guild and land cover type on the estimated z values
# the land use type data for each plot

#load("~/soil-sar/plot-sar-permutation/model_data_SAR.RData")

load("full_parameter_data.RData")

load("~/soil-sar/data/comp_vege.RData")

land_use=comp_vege%>%dplyr::select(plotID,type)%>%distinct()

full_parameter_data%>%left_join(land_use,by="plotID")%>%mutate(type = ifelse(nchar(plotID)<4, "evergreenForest", type))->model_comp

model_comp%>%filter(is.na(type))%>%filter(guild=="all")%>%distinct()%>%dplyr::select(lon,lat,plotID)->coordinates_temp

## based on the coordinates to extract the land use information

raster_data <- rast("nlcd_2021_land_cover_l48_20230630.img")

coordinates <- coordinates_temp[,1:2]

coordinates%>%st_as_sf(coords = c('lon', 'lat'), crs ="EPSG:4326" )%>%vect()%>%project(raster_data )->temp_coords

values=extract(raster_data,temp_coords)


# there are still NA values for some points and we just find the nearest point to get the values

values=bind_cols(values,coordinates_temp)

na_ind <- which(is.na(values$`NLCD Land Cover Class`))
for(i in na_ind) {
  xy <- cbind(values$lon[i], values$lat[i])
  nearest_ind <- which.min(replace(raster::distanceFromPoints(raster_data, xy), is.na(raster_data), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values0 <- as.matrix(raster_data)[nearest_ind,]
  values$`NLCD Land Cover Class`[i] <- values0["NLCD Land Cover Class"]
}

# the code does not work



values=values%>%rename(type="NLCD Land Cover Class")%>%
  mutate(type = case_when(
      type == "Evergreen Forest" ~ "evergreenForest",      
      type == "Mixed Forest" ~ "mixedForest",     
      type == "Herbaceous" ~ "grasslandHerbaceous",  
      type == "Woody Wetlands" ~ "woodyWetlands",   
      type == "Deciduous Forest" ~ "deciduousForest",  
      type == "Deciduous Forest" ~ "deciduousForest",  
      type == "Shrub/Scrub" ~ "shrubScrub",  
      type=="Cultivated Crops" ~"cultivatedCrops",
      TRUE ~ type          
    )
  )


## to match the land use data with the full data

model_comp=model_comp%>%dplyr::left_join(values%>%dplyr::select(plotID,type),by="plotID")%>%dplyr::mutate(type.x=if_else(is.na(type.x),type.y,type.x))%>%dplyr::select(-type.y)%>%rename(type=type.x)



model_comp=model_comp[complete.cases(model_comp),]

mod=aov(zvalue~guild*type,data=model_comp)

shapiro.test(residuals(mod))

leveneTest(zvalue~guild*type,data=model_comp)

summary(mod)

model_comp$guild=as.factor(model_comp$guild)
model_comp$type=as.factor(model_comp$type)

k=aov(logc~guild*type,data=model_comp%>%filter(guild!="all"))

anova_type2 <- Anova(k, type = "II")

aov(logc~type,data=model_comp%>%filter(guild=="all"))%>%Anova(type = "II")




save(model_comp,file="model_comp.RData")

model_comp%>%filter(guild=="all")%>%group_by(type)%>%summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%group_by(type)%>%summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%filter(guild=="all")%>%group_by(type)%>%summarize(count = n())

model_comp%>%filter(guild!="all")->guild_mean



model_comp%>%filter(guild!="all")%>%group_by(guild)%>%summarize(mean_value=mean(zvalue,na.rm=TRUE))

model=lm(log(zvalue)~guild,data=guild_mean)

Anova(model,type = "III")
glht(model, linfct = mcp(guild = "Tukey"))%>%cld()

k=aggregate(zvalue~guild,data=guild_mean,FUN=mean)

od <- k[order(k$zvalue, decreasing = TRUE),] 

guild_mean$guild <- factor(guild_mean$guild, levels = od$guild)


p1=ggplot(guild_mean,aes(x=guild,y=zvalue,fill=guild),alpha=0.5)+
  sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(width = 0.2, color = "black", size=0.1,outlier.size = 1)+ 
  
  theme(legend.position =c(0.5,0.153) ,
        plot.title = element_text(hjust = 0.5, vjust = 2),
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(angle=90), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  ylab(expression(italic(Z)*" value")) +
  xlab("")+
  geom_hline(yintercept = 0.71,color="red",linetype="dashed",size=0.8)+
  annotate("text", x = 1, y = 1.1, label = "ab", size = 6) +
  annotate("text", x = 2, y = 1.15, label = "a", size = 6) +
  annotate("text", x = 3, y = 1.05, label = "b", size = 6) +
  annotate("text", x = 4, y = 1.15, label = "b", size = 6) +
  annotate("text", x = 5, y = 1.04, label = "c", size = 6) +
  annotate("text", x = 6, y = 1.07, label = "e", size = 6) +
  annotate("text", x = 7, y = 1.1, label = "d", size = 6) +
  annotate("text", x = 8, y = 0.99, label = "d", size = 6) +
  ylim(-0.12,1.3)+

  scale_fill_manual("", breaks = od$guild, values = c("chocolate1", "seagreen1", "cadetblue1", "greenyellow", "forestgreen", "purple","lavender", "tan"), 
                    labels = c( "EM(N=415)", "AM(N=326)", "Wood saprotroph(N=401)", 
                                "Epiphyte(N=337)" , 
                               "Litter saprotroph(N=414)", 
                               "Plant pathogen(N=403)", 
                               "Parasite(N=391)","Soil saprotroph(N=416)"))
                                                                                                                                                             )) 


## for the forest types
# just select some of types with more replicates
model_comp%>%filter(guild=="all")%>%group_by(type)%>%summarise(n())%>%filter( `n()`>=10)->many_type

model_comp%>%filter(guild=="all")%>%filter(!is.na(type)&type%in%many_type$type)->type_mean

type_mean$type=as.factor(type_mean$type)

model=lm(log(zvalue)~type,data=type_mean)
Anova(model,type = "III")

glht(model, linfct = mcp(type = "Tukey"))%>%cld()

k=aggregate(zvalue~type,data=type_mean,FUN=mean)

od <- k[order(k$zvalue, decreasing = TRUE),] 

type_mean$type <- factor(type_mean$type, levels = od$type)


p2=ggplot(data=type_mean,aes(x=type,y=zvalue,fill=type),alpha=0.5)+
  sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(width = 0.2, color = "black", size=0.1,outlier.size = 1)+ 
  theme(legend.position =c(0.5,0.174896203), 
        plot.title = element_text(hjust = 0.5, vjust = 2),
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(angle=90), 
        axis.title.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())+
  guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
  ylab(expression(italic(Z)*" value")) +
  xlab("")+
  geom_hline(yintercept = 0.706,color="red",linetype="dashed",size=0.8)+
  annotate("text", x = 1, y = 0.9, label = "d", size = 6) +
  annotate("text", x = 2, y = 0.935, label = "cd", size = 6) +
  annotate("text", x = 3, y = 1, label = "bd", size = 6) +
  annotate("text", x = 4, y = 0.89, label = "bd", size = 6) +
  annotate("text", x = 5, y = 0.85, label = "bc", size = 6) +
  annotate("text", x = 6, y = 0.85, label = "acd", size = 6) +
  annotate("text", x = 7, y = 0.95, label = "ab", size = 6) +
  annotate("text", x = 8, y =0.86307 , label = "a", size = 6)+ 
  
  ylim(0.2,1)+
  scale_fill_manual("", breaks = od$type, 
                    values = c(  "orange", "tan", "greenyellow", "mediumseagreen", "wheat", "pink","lavender","mediumpurple"), 
                    labels = c("MixedForest(N=23)", "DeciduousForest(N=100)","EvergreenForest(N=95)", "WoodyWetlands(N=30)","Grassland\nHerbaceous(N=71)","PastureHay(N=17)",  "CultivatedCrops(N=22)","ShrubScrub(N=56)"))
  
  
  
plot_grid(p1,p2,ncol=2,labels=c("(a)","(b)"))
  
  
leveneTest(zvalue ~ type, data = type_mean)

oneway.test(zvalue ~ type, data = type_mean, na.action = na.omit, var.equal = FALSE)

tukey(type_mean$zvalue, type_mean$type, method = "G")

###create the plots
k <- aggregate(core_rich ~ type, data = compare_richness, FUN = mean)
od <- k[order(k$core_rich), ] # with the increase trend to display the box plots
compare_richness$type <- factor(compare_richness$type, levels = od$type)

  
  





