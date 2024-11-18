library(lmerTest)
library(lme4)
library(smplot2)
library(glmm.hp)
library(reshape2)

# create figure 2
# an example to show the SAR for 10 randomly selected plots

full_dob_neon_richness_data <- readRDS("~/soil-sar/plot-sar-permutation/full_dob_neon_richness_data.rds")
full_dob_neon_richness_data_rarefaction=full_dob_neon_richness_data

load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar.RData")

full_dob_neon_richness_data_rarefaction%>%filter(area!=100&guild=="all"&!is.na(sd_value))%>%dplyr::select(plotid, mean_value, sd_value,  area)%>%
  bind_rows(richness_subplot10_neon_standar)->full_neon_data_with_sd


full_neon_data_with_sd%>%group_by(plotid)%>%summarize(n())%>%filter(`n()`>3)%>%pull(plotid)->plotid_with_four_point

set.seed(125)
rand_plot=sample(plotid_with_four_point,10,replace = FALSE)
full_neon_data_with_sd%>%filter(plotid%in%rand_plot&!is.na(mean_value) )->plot_data



  
 p_z_pva=ggplot(data=plot_data,aes(x=log(area),y=log(mean_value)))+
  geom_point(data=plot_data,aes(x=log(area),y=log(mean_value),color=plotid),size=3)+
  geom_smooth(data=plot_data,aes(x=log(area),y=log(mean_value),color=plotid),method="lm",se=FALSE)+
  theme(legend.position =c(0.5,0.26),
        legend.key.size = unit(0.5, "cm"),
       
        legend.text = element_text(size=10),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=15), 
        axis.text.x = element_text(hjust = 1,size=15), 
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
  guides(color = guide_legend(nrow = 10, byrow = TRUE))+
  ylim(2,8)

# the plot with the raw data
 
 # Fit a nonlinear least squares model for each group
 fit_power_law <- function(df) {
   nls(mean_value ~ c * area^z, data = df, start = list(c = 1, z = 0.5))
 }
 
 # Apply the NLS function to each group
 nls_fits <- plot_data%>%
   group_by(plotid) %>%
   do(model = fit_power_law(.))
 
 #to get the values
 
 nls_fits $model[[1]]
 
 
 # Predict values based on the model fits
 predicted_values <- nls_fits %>%
   rowwise() %>%
   mutate(Predicted = list(data.frame(
     area = seq(min(plot_data$area), max(plot_data$area), length.out = 100),
     mean_value= predict(model, newdata = data.frame(area = seq(min(plot_data$area), max(plot_data$area), length.out = 100)))
   ))) %>%
   unnest(cols = c(Predicted))
 
 # Plotting the data and the fitted NLS lines
 
 p_zvalue=ggplot(plot_data, aes(x = area, y = mean_value, color = plotid)) +
   geom_errorbar(data=plot_data,aes(ymin=mean_value-sd_value,ymax=mean_value+sd_value,color=plotid),width=45)+
   geom_point(size=3) +
   scale_color_manual("",breaks=unique(plot_data$plotid),
                      labels=c(expression("Central Plains Experimental Range "*italic( z)*"=0.65"),
                               expression("Guanica Dry Forest Reserve "*italic(z )*"=0.74"),
                               expression("Harvard Forest "*italic(z)*"=0.65"),
                               expression("Niwot Ridge Mountain Research Station "*italic(z)*"=0.70"),
                               expression("Oak Ridge National Laboratory "*italic(z)*"=0.69"),
                               expression("Ordway-Swisher Biological Station "*italic(z)*"=0.82"),
                               expression("The San Joaquin Experimental range "*italic(z)*"=0.65"),
                               expression("The Talladega National Forest "*italic(z)*"=0.79"),
                               expression("University of Kansas Field Station "*italic(z)*"=0.72"),
                               expression("Woodworth "*italic(z)*"=0.66")),
                      values = c("chocolate1", "#037f77", "royalblue", "#f0a73a", "forestgreen", "#7c1a97","#c94e65", "tan","pink","gray"))+
   geom_line(data = predicted_values, aes(x = area, y = mean_value, color = plotid), linetype = "solid",size=1.2)+ 
   theme(legend.position =c(0.37,0.780361205205),
         legend.text = element_text(size=10),
         legend.title  = element_text(size=10),
         legend.key.height = unit(0.48, "cm"),
         text = element_text(size = 18),
         plot.title = element_text(size = 12, hjust = 0.5), 
         axis.text.y = element_text(hjust = 0,size=15), 
         axis.text.x = element_text(hjust = 0.5,size=15), 
         axis.title.y = element_text(size = 18), 
         axis.title.x = element_text(size = 18), 
         plot.margin = unit(c(0.2, 2, 0, 1), "cm"),
         panel.background = element_rect(fill = "NA"),
         panel.border = element_rect(color = "black", size = 1, fill = NA))+
   xlab(expression(Area *" " * m^2))+
   ylab("Richness")+
   guides(color = guide_legend(nrow = 10))+
   ylim(0,3500)
   
 
 
 
 
 
# to compare the z values among land use types

load("~/soil-sar/plot-sar-permutation/model_comp.RData")
model_comp%>%filter(guild=="all")%>%group_by(type)%>%dplyr::summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%group_by(type)%>%summarize(mean_value=mean(zvalue,na.rm=TRUE))

model_comp%>%filter(guild=="all")%>%group_by(type)%>%summarise(count = n())

model_comp%>%filter(guild!="all")->guild_mean

k=aggregate(zvalue~guild,data=guild_mean,FUN=mean)
od <- k[order(k$zvalue, decreasing = TRUE),] 
guild_mean$guild <- factor(guild_mean$guild, levels = od$guild)

p_z_guild=ggplot(guild_mean,aes(x=guild,y=zvalue,fill=guild),alpha=0.5)+
  sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(width = 0.2, color = "black", size=0.1,outlier.size = 1)+
  theme(legend.position =c(0.5,0.15) ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, vjust = 2),
        legend.text = element_text(size = 10), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(angle=90,size=15,hjust=0.5), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color="white"),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),
        panel.border  = element_rect(size=1),
        )+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +

  ylab(expression(italic(z)*" value")) +
  xlab("Trophic guilds")+
  geom_hline(yintercept = 0.71,color="red",linetype="dashed",size=0.8)+
  annotate("text", x = 1, y = 1.3, label = "a", size = 6) +
  annotate("text", x = 2, y = 1.05, label = "ab", size = 6) +
  annotate("text", x = 3, y = 1.05, label = "b", size = 6) +
  annotate("text", x = 4, y = 1.15, label = "bc", size = 6) +
  annotate("text", x = 5, y = 1.04, label = "c", size = 6) +
  annotate("text", x = 6, y = 1.2, label = "d", size = 6) +
  annotate("text", x = 7, y = 1.014, label = "d", size = 6) +
  annotate("text", x = 8, y = 1.15, label = "d", size = 6) +
  ylim(0,1.3)+
  scale_fill_manual("", breaks = od$guild, values = c("chocolate1", "gray", "royalblue", "#f0a73a", "seagreen", "#7c1a97","#c94e65", "tan"), 
                    labels = c("EM (N=415)", "AM (N=326)", "Wood saprotroph (N=401)", 
                               "Epiphyte (N=337)" , 
                               "Litter saprotroph (N=414)", 
                               "Plant pathogen (N=403)", 
                               "Parasite (N=391)","Soil saprotroph (N=416)"))

# the impact of climate variables on the estimated z values
# plant richness was not included


setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standard sar neon")
#load("~/soil-sar/plot-sar-permutation/model_data.RData")
load("~/soil-sar/plot-sar-permutation/soil_mean.RData")

load("~/soil-sar/plot-sar-permutation/model_data.RData")

full_parameter_data=readRDS("full_parameter_data.rds")

# need to combine the complete climate and soil variables
#
full_parameter_data%>%dplyr::select(logc,zvalue,plotID,guild)%>%
  left_join(model_data%>%dplyr::select(plotID,siteIDD,bio1,bio2,bio4,bio8,bio12,bio15,bio18,richness)%>%distinct(),by="plotID")%>%
  left_join(soil_mean%>%dplyr::select(-lon,-lat)%>%distinct(),by="plotID")->model_data_SAR_rarefaction

model_data_SAR_rarefaction$logc=2.71828^model_data_SAR_rarefaction$logc


# use this data for modeling when plant richness was not included

model_data_SAR_climate=model_data_SAR_rarefaction%>%dplyr::select(siteIDD,plotID,guild,logc,  zvalue,  soilInCaClpH ,nitrogenPercent, organicCPercent, soilMoisture, cec, sand,  bio1, bio2,  bio4, bio8, bio12, bio15, bio18)
model_data_SAR_climate=model_data_SAR_climate[complete.cases(model_data_SAR_climate),]

guild_select=unique(model_data_SAR_climate$guild)

# to standardize the data

data=list()
for (i in 1:9)
{
  d=model_data_SAR_climate%>%filter(guild==guild_select[i]) 
  
  d[,c(4,6:18)]=apply(d[,c(4,6:18)],2,range01)
  
  data[[i]]=d
}

# the number of plots included in the data
a=numeric()
for(i in 1:9)
{
  a[i]=dim(data[[i]])[1]
}


# to check the linearity of the variables for each data set

high_cor_variables=list()
for (i in 1:9)
{
  cor_matrix <- cor(data[[i]][,6:18])
  # Find columns with correlation coefficients higher than 0.8
  high_cor_variables[[i]] <- colnames(data[[i]][,6:18])[apply(cor_matrix, 2, function(x) any(abs(x) > 0.70 & abs(x) < 1))]
  
}

# the effect size of the variables for each of the guilds
# do not do a model selection
effect=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()
    effect[[i]]=kk[2:14,][,c(1,2,5)]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH+soilMoisture +organicCPercent +cec+ sand +bio1+bio4+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()# for i=9, bio4 was excluded to avoid colinearity
      effect[[i]]=kk[2:14,][,c(1,2,5)]
    }
  }
}

##

effect[[1]]%>%mutate(var=rownames(effect[[1]]))%>%
  bind_rows(effect[[2]]%>%mutate(var=rownames(effect[[2]])))%>%
  bind_rows(effect[[3]]%>%mutate(var=rownames(effect[[3]])))%>%
  bind_rows(effect[[4]]%>%mutate(var=rownames(effect[[4]])))%>%
  bind_rows(effect[[5]]%>%mutate(var=rownames(effect[[5]])))%>%
  bind_rows(effect[[6]]%>%mutate(var=rownames(effect[[6]])))%>%
  bind_rows(effect[[7]]%>%mutate(var=rownames(effect[[7]])))%>%
  bind_rows(effect[[8]]%>%mutate(var=rownames(effect[[8]])))%>%
  bind_rows(effect[[9]]%>%mutate(var=rownames(effect[[9]])))%>%
  mutate(guild=rep(guild_select,each=13))%>%
  rename_all(~paste0(c("estimate","sd","pva","var","guild")))-> effect_no_plant


effect_no_plant$var=factor(effect_no_plant$var,levels=c("logc","soilInCaClpH", "organicCPercent", "soilMoisture" , "cec", "sand" , "bio1" , "bio2","bio4" , "bio8" ,"bio12" , "bio15" , "bio18"))        
effect_no_plant$guild=factor(effect_no_plant$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           
effect_no_plant$pva[effect_no_plant$pva<=0.001]="***"
effect_no_plant$pva[effect_no_plant$pva>0.001&effect_no_plant$pva<=0.01]="**"
effect_no_plant$pva[effect_no_plant$pva>0.01&effect_no_plant$pva<=0.05]="*"
effect_no_plant$pva[effect_no_plant$pva>0.05]=""




p11=ggplot(data=effect_no_plant%>%filter(guild=="all"),aes(x=estimate,y=1:13))+
  geom_point(size=4,pch=21,color="black",fill=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  
  geom_errorbar(aes(xmin = estimate-1.96*sd, xmax = 1.96*sd+estimate,width=0.2),
                color=rep(c("#DE582E","#F3A332","#1868b2"),times=c(1,5,7)))+
  geom_vline(xintercept = 0,linetype="dashed")+
  scale_y_continuous (breaks=1:13,labels = c(expression(italic(C)),"pH",  "Moisture","SoilC", "CEC", "Sand" , "MAT" ,   
                                             "MDR","Temp.seas." , "MTWQ" ,"MAP" , "Pre.seas." , "PWQ" ))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=15),
        
        axis.ticks.length.y  = unit(-0.150, "cm"),
        axis.text.y  = element_text(size=12,hjust=0,margin = margin(r=-70)),
        plot.title=element_text(hjust=0.5,size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 0.5),
        panel.background = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.3), "cm"))+
  ylab("")+
  xlab("Effect size ± 95%CI (N=453)")+
  xlim(-0.6,0.6)+
  #ggtitle("Climate+Soil(N=453)")+
  annotate("text",x=-0.31,y=1,label="***",size=8)+
  annotate("text",x=0.265,y=7,label="*",size=8)+
  annotate("text",x=0.265,y=9,label="*",size=8)+
  annotate("text",x=0.24583,y=11,label="*",size=8)+
  annotate("text",x=-0.16820583,y=3,label="*",size=8)+
  xlim(-0.5,0.5)

## the variance partition 
soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")
## do we need to select the variables



sel_vab_step=list()
for (i in 1:9)
{
  if(i<9)
  {
    mod <- lmer(zvalue ~ logc+soilInCaClpH +soilMoisture +organicCPercent +cec+ sand +bio1+ bio2 +bio4 + bio8+bio12 + bio15 + bio18  + (1 |siteIDD), data = data[[i]])
    mod=summary(mod)
    kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
    sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
  }
  else{
    {
      mod <- lmer(zvalue ~ logc+soilInCaClpH +soilMoisture +organicCPercent +cec+ sand +bio1+ bio2+bio8  + bio12 + bio15 + bio18 + (1 |siteIDD), data = data[[i]])
      mod=summary(mod)
      kk=mod$coefficients%>%data.frame()%>%filter(!Pr...t..>0.05)
      sel_vab_step[[i]]=rownames(kk)[1:dim(kk)[1]]
    }
  }
}


#
soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_n=which(sel_vab_step[[i]]%in%soil)%>%length()
  climate_n=which(sel_vab_step[[i]]%in%climate)%>%length()
  
  if (soil_n>0&&climate_n>0)
  {
    soil_variable=data[[i]][,intersect( sel_vab_step[[i]],soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    climate_variable=data[[i]][,intersect( sel_vab_step[[i]],climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],soil_variable,climate_variable)%>%rename(soil=product...19,climate=product...20)
    mod=lmer(zvalue~logc+soil+climate+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when climate variables are lacking
  else if (soil_n>0&&climate_n<1)
  {
    soil_variable=data[[i]][,intersect( sel_vab_step[[i]],soil)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    
    data_temp=bind_cols(data[[i]],soil_variable)%>%rename(soil=product)
    mod=lmer(zvalue~logc+soil+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  
  # when soil variables are lacking
  else if (soil_n<1&&climate_n>0)
  {
    climate_variable=data[[i]][,intersect( sel_vab_step[[i]],climate)]%>%data.frame()%>%rowwise() %>%
      mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
    data_temp=bind_cols(data[[i]],climate_variable)%>%rename(climate=product)
    mod=lmer(zvalue~logc+climate+(1 |siteIDD),data=data_temp)
    a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=a$r.squaredGLMM
    var_par[[i]]=a$commonality.analysis
  }
  # when both climate and soil variables are lacking
  else 
  {
    mod=lmer(zvalue~logc+(1 |siteIDD),data=data[[i]])
    #a=glmm.hp(mod,commonality=TRUE)
    var_diff[i,]=r.squaredGLMM(mod)
    var_par[[i]]=r.squaredGLMM(mod)
  }
}




var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
{
  soil_variable=data[[i]][,intersect( sel_vab_step[[i]],soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
  climate_variable=data[[i]][,intersect( sel_vab_step[[i]],climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%dplyr::select(product)%>%data.frame()
  #plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
  #mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],soil_variable*10,climate_variable)%>%rename(soil=product...19,climate=product...20)
  mod=lmer(zvalue~logc+soil+climate+(1|siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

## to get the random effect part and the residuals of the model

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%dplyr::select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%dplyr::select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

# to get the fixed effect part of the model

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),
            var_par[[3]]%>%data.frame(),
            var_par[[4]]%>%data.frame(),
            var_par[[5]]%>%data.frame(),
            var_par[[6]]%>%data.frame(),
            var_par[[7]]%>%data.frame(),
            var_par[[8]]%>%data.frame(),
            var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))
kk$type=gsub("\\d", "", kk$type)
kk$type=gsub("\\.{3}", "", kk$type)
var_part_data=bind_cols(k,kk)
var_part_data$type=gsub(" ","",var_part_data$type)
var_part_data%>%data.frame()%>%rename_all(~paste0(c("Fractions","total","guild","type")))->var_fixed# some do have the data
var_fixed%>%dplyr::select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new_noplant

p1=ggplot(data=varp_new_noplant%>%filter(Fractions>0&guild=="all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual(expression(italic(R["Fixed"]^2)*" =58.14%"),breaks=unique(varp_new_noplant$type),
                    labels=c(expression(italic(C)),"S","Cli.",expression(italic(C)*"+S"),expression(italic(C)*"+Cli."),
                             "S+Cli.",expression(italic(C)*"+S+Cli."),"Random","Residuals"),
                    values=c("#F3A332", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple","lavender", "#FFFFCC","gray"))+
  theme(legend.position = c(0.5,0.7),
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.key.size = unit(0.35, "cm") ,
        legend.key.height = unit(0.38, "cm"),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.5,0, 0.5,0.5), "cm"),
        panel.background = element_rect(fill = "NA"))+
  guides(fill = guide_legend(nrow =10, byrow = TRUE))+
  ylab("Variance explained")+
  xlab("")+
  #ggtitle(expression(italic(R["Fixed"]^2)*" =56.4%"))+
  ylim(0,2)

p_z_distribu=ggplot(full_parameter_data%>%filter(guild=="all"&!is.na(zvalue)), aes(x = zvalue)) +
  geom_histogram(binwidth = 0.02, fill = "#F3A332", color = "black", alpha = 0.7) +
  ggtitle("") +
  xlab(expression(italic(z)*" value")) +
  ylab("Frequency") +
  geom_vline(xintercept = 0.713,linetype="dashed",color="red")+
  theme(legend.position = c(0.75,0.28),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        text = element_text(size = 18),
        plot.title = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0,size=15), 
        axis.text.x = element_text(hjust = 0.5,size=15), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  annotate("text",x=0.5,y=50,label="CV=11%",size=5)



p_effects=plot_grid(p1,p11,ncol=2,rel_widths = c(1,2))

p1=ggplotGrob(p1)
p11=ggplotGrob(p11)

p1$heights=p11$heights


p_z_pva=ggplotGrob(p_zvalue)
p_effects=ggplotGrob(p_effects)
p_z_guild=ggplotGrob(p_z_guild)
p_z_distribu=ggplotGrob(p_z_distribu)

p_z_pva$heights=p_z_distribu$heights

p_z_pva$heights=p_z_guild$heights

p_z_pva$heights=p_z_distribu$heights
p_z_guild$heights=p_effects$heights



p_z_pva$widths=p_z_distribu$widths

p_z_pva$widths=p_z_guild$widths

p_z_pva$widths=p_z_distribu$widths
p_z_guild$widths=p_effects$widths
p_z_guild$widths=p_z_pva$widths
p_z_guild$heights=p1$heights

p_effects=plot_grid(p1,p11,ncol=2,rel_widths = c(1,2))




plot_grid(p_z_pva,p_z_distribu,p_z_guild,p_effects,ncol=2, label_size = 18, label_x = 0.1,
          label_y = c(1.01,1.01,1.03,1.03),
          labels = paste0("(", letters[1:4], ")"))

plot_grid(p_z_guild,p_effects,ncol=2)
