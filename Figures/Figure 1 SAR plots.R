
library(ggplot2)
library(dplyr)
library(reshape)

set.seed(1010)
times=30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three cores
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
     
      return(ex)
    }
    stopCluster(cl)
  }
  
  else
  {
    power.z[[i]]=matrix(10,nrow=2,ncol = dim1[1])#the number 10 is randomly selected, to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}

# to rbind all the data

df=numeric()
for (i in 1:515){
  df[i]=dim(power.z[[i]])[1]
}
df=cbind(df,a1)%>%data.frame()
df=data.frame(df,or=1:515)
df$df=as.numeric(df$df)
df=subset(df,df>2)

# get all the 30 estimates for each plot

d30=list()
d30a=list()
dcb=list()
for (i in df$or[1:493])
  {
 d30[[i]]=melt(power.z[[i]][,c(seq(from=2,to=60,by=2))])["value"]
 d30a[[i]]=melt(power.z[[i]][,c(seq(from=1,to=59,by=2))])["value"]
 dcb[[i]]=cbind(d30a[[i]],d30[[i]])
}# some rows have NA because of < 3 cores


# cbind all the simulations 

dall=dcb[[1]]
for (i in df$or[2:493]){
  dall=rbind(dall,dcb[[i]])
}
names(dall)=c("A","spe")

# add the replication to each plot

ak=list()
for (i in 1:493)
{
  akk=df[i,]
  akkk=paste(akk$a1,"_",1:30)
  ak[[i]]=rep(akkk,each=df$df[i])%>%data.frame()
}

# rbind all the names
ak2=ak[[1]]
for (i in 2:493){
ak2=rbind(ak2,ak[[i]])
}

df30=cbind(ak2,dall)
names(df30)[1]="pt"
pid=substr(df30$pt,1,8)
df30=cbind(df30,pid)

#create plots



df1=power.z[[1]]
for (i in df$or[2:493])
  {
  df1=rbind(df1,power.z[[i]])
}

df2=data.frame(df1$A)

df1=df1[,c(seq(from=2,to=60,by=2))]

df3=apply(df1,1,mean)%>%data.frame()
df4=cbind(df2,df3)
names(df4)=c("A","spe")
pid=rep(df$a1,times=df$df)
df4=cbind(df4,pid)
names(df)[2]="pid"

df4=merge(df4,df,by="pid")
names(df4)[4]="cores"
names(df4)[1]="plotID"

df4=merge(df4,comp_vege[,c(2,4)],by="plotID")
df4=unique(df4)

ggplot(data=df30,aes(x=log(A),y=log(spe),color=pt))+
  #geom_point(size=2,alpha=0.25)+
  geom_smooth(data=df30,aes(x=log(A),y=log(spe),color=pt),method = "lm",se=FALSE,size=0.25,alpha=0.5) +
  guides(color="none")+
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_text(hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  xlab("Log(Area)") +
  ylab("Log(Species)")+
  scale_color_manual(breaks=unique(df30$pt),labels=unique(df30$pt),values=rep("pink",times=14790))

# the density plot
  
  
  # for one plot

  p1=ggplot(data=subset(df30,piddd=="WY1"),
           aes(x=log(A),y=log(spe),color=pid))+
    geom_point(size=3,alpha=0.1)+
    geom_smooth(data=subset(df30,piddd=="WY1"),
                aes(x=log(A),y=log(spe),color=pt),
                method = "lm",se=FALSE,size=0.45)+
    theme(legend.position = "bottom",
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(hjust = 0), 
          axis.text.x = element_text(hjust = 1), 
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18), 
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
    xlab("Log(Area)") +
    ylab("Log(Species)") +
    guides(color="none")
  
  ####
  p2=ggplot(data=subset(df30,pid%in%c("NIWO_043","TALL_044")),
         aes(x=log(A),y=log(spe),color=pt))+
    geom_point(size=3,alpha=0.1)+
    geom_smooth(data=subset(df30,pid%in%c("NIWO_043","TALL_044")),
                aes(x=log(A),y=log(spe),color=pt),
                method = "lm",se=FALSE,size=0.45)+
    scale_color_manual(breaks=subset(df30,pid%in%c("NIWO_043","TALL_044"))[,"pt"],labels=rep(c(1:2),each=750),values=rep(c("seagreen1","mediumpurple"),each=750))+
    guides(color="none")+
    theme(legend.position = "bottom",
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(hjust = 0), 
          axis.text.x = element_text(hjust = 1), 
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18), 
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
    xlab("Log(Area)") +
    ylab("Log(Species)") +
    geom_point(aes(x=2,y=5),size=3)+
    annotate("text", x = 2.5, y = 5, label = "NIWO_043", color="mediumpurple",size = 3) +
    annotate("text", x = 2.5, y = 4.8, label = "TALL_044", size = 3) +
    geom_point(aes(x=2,y=4.8),size=3,color="seagreen1")
  
  ####
  p3=ggplot(data = subset(df30mean,projd=="neon"&Freq>24), aes(x = log(A), y = spe, color = piddd)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_segment(data = subset(df30mean,projd=="neon"&Freq>24), size=0.3,
                 aes(x = log(A), xend=log(A),y = spe-sd,yend=spe+sd,
                     color = as.factor(piddd)))+
    guides(color = guide_legend(nrow = 6, byrow = TRUE))+
    theme(legend.position = c(0.7,0.213), 
          legend.title = element_text(size=10),
          text = element_text(size = 18), 
          legend.text = element_text(size=8),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(hjust = 0), 
          axis.text.x = element_text(hjust = 1), 
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18),
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "NA"), 
          panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
    xlab("Log(Area)") +
    ylab("Log(Species)") +
    scale_color_manual("plotID",breaks=unique(data$piddd),labels=unique(data$piddd),
                       values=c(c("purple", "gray", "cadetblue1", "tan1", "greenyellow", "mediumseagreen", "burlywood1", "tan", "wheat", "gold", "mediumpurple",
                                          "midnightblue","blue","cyan","black","gold3","deeppink","aquamarine2")))
                                          
  
  
  
  ggplot(data=subset(model_data,z<10),aes(x=z,color=projd))+
    geom_histogram(aes(y=..density..,fill=projd),color="gray")+
    geom_density(alpha=0.1)+
    theme(legend.position = c(0.5,0.7),
          text = element_text(size = 18),
          plot.title = element_text(size = 15, hjust = 0.5), 
          axis.text.y = element_text(hjust = 0), 
          axis.text.x = element_text(hjust = 1), 
          axis.title.y = element_text(size = 18), 
          axis.title.x = element_text(size = 18), 
          axis.ticks.x = element_blank(), 
          panel.background = element_rect(fill = "NA"),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
    xlab(expression(italic(z)))+
    ylab("Density")+
    geom_vline(xintercept  =0.7538,linetype="dashed",color="blue",size=1)+
    annotate("text",x=1.15,y=3.5,label="Mean=0.7538",size=6)+
  geom_vline(xintercept  =0.717,linetype="dashed",color="blue",size=1)
  
  scale_color_manual("Project",breaks=c("dob","neon"),labels=c("dob","neon"),values=c("mediumpurple","seagreen1"))+
  scale_fill_manual("Project",breaks=c("dob","neon"),labels=c("dob","neon"),values=c("mediumpurple","seagreen1"))
  
  
  
  

  

# add the plotID to the data
pid=vector("list",length=dim(df30)[1])
for (i in 1:dim(df30)[1])
  {
  pidd=df30[i,4]
  if (str_detect(pidd," "))
    {
   pid[i] =substr(pidd,1,3)
   
  }
  else{
    pid[i]=substr(pidd,1,8)
  }
}

piddd=matrix(ncol=1,nrow=dim(df30)[1])
for (i in 1:dim(df30)[1])
  {
  piddd[i,1]=pid[[i]]
}

df30=cbind(df30,piddd)

df30mean=aggregate(log(spe)~A*piddd,data=df30,FUN=mean)
df30sd=aggregate(log(spe)~A*piddd,data=df30,FUN=sd)

df30mean=cbind(df30mean,df30sd["log(spe)"])

names(df30mean)[4]="sd"

## add the project id to the data

proj=vector("list",length=dim(df30mean)[1])
for (i in 1:dim(df30mean)[1])
{
  pidd=df30mean[i,2]
  dim2=str_length(pidd)
  if (dim2<4)
  {
    proj[[i]] ="dob"
  }
  else{
    proj[[i]]="neon"
  }
}

projd=matrix(ncol=1,nrow=dim(df30mean)[1])
for (i in 1:dim(df30mean)[1])
{
  projd[i,1]=proj[[i]]
}

df30mean=cbind(df30mean,projd)

names(df30mean)[3]="spe"

names(df)[2]="piddd"

df30mean=merge(df30mean,df[,c(1,2)],by="piddd")
names(df30mean)[6]="Freq"



# head(df1)





# exmaine the patterns between area and the total richness

EP=numeric()
for (i in 1:120)
  {
  EP0=summary(lm(log(species)~log(A),data=power.z[[i]]))
  EP[i]=EP0[[9]]
}

# if we use a non-linear model

mod=nls(log(species)~c*log(A)^z,power.z[[1]],start = list(c=3,z=3))

fit1=numeric()
for(i in 1:120)
  {
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  d=summary(sar_power(power.z[[i]]))
  fit1[i]=d$AIC
}

fit2=numeric()
for(i in 1:120)
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  d=summary(sar_linear(power.z[[i]]))
  fit2[i]=d$AIC
}


all_z=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for all the plots
for(i in 1:length(a1)){
  all_z[i,]=power.z[[i]][1,]
}


##


op=aggregate(log(spe)~log(A),data=df30,FUN=mean)
opd=aggregate(log(spe)~log(A),data=df30,FUN=sd)

op=cbind(op,opd[,"log(spe)"])
names(op)=c("A","spe","sd")

ggplot()+
  geom_point(data=op,aes(x=A,y=spe))+
  geom_segment(data=op,aes(x=A,xend=A,y=spe-sd,yend=spe+sd))
  
