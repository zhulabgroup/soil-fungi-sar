

soil=c("soilInCaClpH", "nitrogenPercent", "organicCPercent","soilMoisture", "cec","sand")
climate=c("bio1", "bio2",  "bio4", "bio8", "bio12", "bio15", "bio18")
plant=("richness")

var_par=list()
var_diff=matrix(ncol=2,nrow=9)
for (i in 1:9)
  {
soil_n=which(sel_vab_step[[i]]%in%soil)%>%length()
climate_n=which(sel_vab_step[[i]]%in%climate)%>%length()
plant_n=which(sel_vab_step[[i]]%in%plant)%>%length()

if (soil_n>0&&climate_n>0&&plant_n>0)
{
  soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],soil_variable,climate_variable,plant_variable)%>%rename(soil=product...22,climate=product...23,plant=product...24)
  mod=lmer(zvalue~logc+soil+climate+plant+(1 |siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

# when climate variables are lacking
else if (soil_n>0&&climate_n<1&&plant_n>0)
{
  soil_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%soil)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  
  plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],soil_variable,plant_variable)%>%rename(soil=product...22,plant=product...23)
  mod=lmer(zvalue~logc+soil+plant+(1 |siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

# when soil variables are lacking
#else if (soil_n<1&&climate_n>0&&plant_n>0)
else
{
  climate_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%climate)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  plant_variable=data[[i]][,sel_vab_step[[i]]][,which(sel_vab_step[[i]]%in%plant)]%>%data.frame()%>%rowwise() %>%
    mutate(product = prod(c_across(everything())))%>%select(product)%>%data.frame()
  data_temp=bind_cols(data[[i]],climate_variable,plant_variable)%>%rename(climate=product...22,plant=product...23)
  mod=lmer(zvalue~logc+climate+plant+(1 |siteIDD),data=data_temp)
  a=glmm.hp(mod,commonality=TRUE)
  var_diff[i,]=a$r.squaredGLMM
  var_par[[i]]=a$commonality.analysis
}

}

# to combine all the effects

k=bind_rows(var_par[[1]]%>%data.frame(),var_par[[2]]%>%data.frame(),var_par[[3]]%>%data.frame(),var_par[[4]]%>%data.frame(),var_par[[5]]%>%data.frame(),var_par[[6]]%>%data.frame(),var_par[[7]]%>%data.frame(),var_par[[8]]%>%data.frame(),var_par[[9]]%>%data.frame())%>%
  mutate(guild=rep(guild_select,times=c(dim(var_par[[1]])[1],dim(var_par[[2]])[1],dim(var_par[[3]])[1],dim(var_par[[4]])[1],dim(var_par[[5]])[1],dim(var_par[[6]])[1],dim(var_par[[7]])[1],dim(var_par[[8]])[1],dim(var_par[[9]])[1])))

kk=rownames(k)%>%data.frame()%>%rename_all(~paste0("type"))

kk$type=gsub("\\d", "", kk$type)

kk$type=gsub("\\.{3}", "", kk$type)

var_part_data=bind_cols(k,kk)

var_part_data$type=gsub(" ","",var_part_data$type)

var_part_data%>%data.frame()%>%rename_all(~paste0(c("fraction","total","guild","type")))# total 16 types
## create a plot

 var_part_data%>%filter(!str_detect(type, "logc"))->var_part_data_no_logc

 
 varp_new$guild=factor(varp_new$guild,levels=c("all","AM","EM","epiphy","littersap","para","plapat","soilsap","woodsap"))                           
 

 ggplot(data=varp_new%>%filter(Fractions>0), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill",color="black",width = 0.6)+
  scale_fill_manual("Component",breaks=unique(varp_new$type),
                    labels=c("Logc","Soil","Climate","Plant","Logc+Soil","Logc+Climate",
                             "Soil+Climate","Logc+Plant","Soil+Plant","Climate+Plant",
                             "Logc+Soil+Plant","logc+Soil+Plant","Logc+Climate+Plant","Soil+Climate+Plant",
                             "Logc+Soil+Climate+Plant","Random effect","Residuls"),
                             values=c("tan", "seagreen1", "cadetblue1", "greenyellow", "forestgreen", "purple","lavender", "orange","black","yellow","red","blue","royalblue","mediumpurple","pink","#FFFFCC","gray"))+
  scale_x_discrete(breaks=unique(varp_new$guild),
                   labels = c("All(N=279)","ECM(N=279)","Litter sapro.(N=279)","Parasite(N=279)",
                              "Soil sapro.(N=276)","Wood sapro.(N=279)","Plant patho.(N=278)","AM(N=245)", "Epiphyte(N=245)"))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_text(angle=0,size=15,hjust = 1),
        plot.margin = margin(b=-0.5, unit="cm"),
        panel.background = element_rect(fill = "NA"))+
   guides(fill = guide_legend(nrow = 17, byrow = TRUE))+
  ylab("Variance explained")+
   coord_flip()+
   xlab("")
  



ggplot(data=ddm, aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity",position="fill")


# to see the total variance

var_diff%>%data.frame()%>%mutate(ree=1-X2,randd=X2-X1)%>%rename_all(~paste0(c("fix","ran","ree","randd")))->total_var

total_var%>%select(fix,randd,ree)%>%mutate(guild=guild_select)%>%melt()->temp

temp%>%select(value,guild,variable)%>%rename_all(~paste0(c("Fractions","guild","type")))->temp

var_part_data%>%select(Fractions,guild,type)%>%filter(type!="Total")%>%bind_rows(temp)%>%filter(type!="fix")->varp_new

varp_new%>%filter(guild=="all")

