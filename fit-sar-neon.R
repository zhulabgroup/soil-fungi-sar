# to see the fit for the neon sites

table(sar_neon_permutation$plotid)

point_number=table(sar_neon_permutation$plotid)%>%data.frame()

point_number_five=subset(point_number,Freq>4)
point_number_four=subset(point_number,Freq==4)
point_number_three=subset(point_number,Freq==3)

## to see the fits

for (i in 1:dim(point_number_five)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_five$Var1)[i])
  
  fit <-sar_power(data_subb[,3:2] )
  plot=plot(fit)
  
  # Print the plot
  print(plot)
  
}

for (i in 1:dim(point_number_five)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_five$Var1)[i])
  
  fit <-plot(log(area)~log(richness),data=data_subb[,3:2] )
 
  
  # Print the plot
  print(fit)
  
}


for (i in 1:dim(point_number_four)[1])
{
  data_subb <- subset(sar_neon_permutation,plotid==unique(point_number_four$Var1)[i])
  
  fit <-plot(log(area)~log(richness),data=data_subb[,3:2] )
  
  
  # Print the plot
  print(fit)
  
}

## estimate the z value



zvalue=numeric()
pvalue=numeric()
a5=unique(sar_neon_permutation$plotid)
for (i in 1:472)
{
  df1=subset(sar_neon_permutation,plotid==a5[i])
  ft=lm(log(richness)~log(area),data=df1)%>%summary()
  zvalue[i]=ft$coefficients[2,1]
  pvalue[i]=ft$coefficients[1,4]
}

df=cbind(a5,zvalue)%>%data.frame()

df$zvalue=as.numeric(df$zvalue)
names(df)[1]="plotid"

names(point_number)[1]="plotid"

df=merge(df,point_number,by="plotid")
# when 4 and 5 points were included, there is no relation between core numbers and the estimated the z values
