
a=ggplot()+geom_point(data=pca_39_leaf_spe_8,alpha=0.8,size=4,aes(x=Comp.1  ,y=Comp.2 ,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_leaf_trait_8,size=0.7,color="green3",aes(x=0,y=0,xend=1.5*Comp.1  ,yend=1.5*Comp.2),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (55.5%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+
  geom_text(data=pca_39_leaf_trait_8,x=1.4*c(  0.15978201, -0.47839229, -0.29770778,  0.586344569681,  0.33634880,  0.02622771, -0.58449088118,  0.32179531),y=1.5*c(0.73169516274,  0.283900465,  0.382305429, -0.0531305629782, -0.2025663550,  0.3653031020 ,-0.03798437,  0.5126584300),size=5,label=c(expression("C"["leaf"]),expression("N"["leaf"]),expression("P"["leaf"]),expression("C:N"["leaf"]),expression("C:P"["leaf"]),"LA","SLA","LTD"))+
  theme(legend.position = c(0.27,0.14685291425),plot.title=element_text(hjust=0.5,face="bold",size=20),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+
  ylim(-1,1.1)+xlim(-1.1,1)+ylab("Second PCA axis (20.5%)")+
  labs(title="Standardized PCA",size=40)


ggplot()+geom_point(data=pca_39_leaf_spe_8_phy,alpha=0.8,size=4,aes(x=PC1  ,y=PC2 ,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_leaf_trait_8_phy,size=0.7,color="green",aes(x=0,y=0,xend=1.5*PC1  ,yend=1.5*PC2),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (61.8%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+geom_text(data=pca_39_leaf_trait_8,x=1.4*c(  0.294646425, -1.092079707, -0.82088597,  1.131096193153601,  1.088017087,  0.08361598, -1.0982323495, 0.78962697260),y=1.5*c(0.5969781505, -0.009170698,  0.3953720517,  0.087579436, -0.214805634 , 0.94618154352,-0.070564728,  0.10278209118),size=5,label=c(expression("C"["leaf"]),expression("N"["leaf"]),expression("P"["leaf"]),expression("C:N"["leaf"]),expression("C:P"["leaf"]),"LA","SLA","LTD"))+theme(legend.position = c(0.27,0.14685291425),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+ylab("Second PCA axis (14.2%)")+xlim(-1.6,1.7)

b=ggplot()+geom_point(data=pca_39_leaf_spe_8_phy,alpha=0.8,size=4,aes(x=PC1  ,y=PC2 ,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_leaf_trait_8_phy,size=0.7,color="green3",aes(x=0,y=0,xend=PC1  ,yend=PC2),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (61.8%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+
  geom_text(data=pca_39_leaf_trait_8,x=c(  0.294646425, -1.092079707, -0.82088597,  1.131096193153601,  1.088017087,  0.08361598, -1.0982323495, 0.78962697260),y=c(0.5969781505, -0.009170698,  0.3953720517,  0.087579436, -0.214805634 , 0.94618154352,-0.070564728,  0.140278209118),size=5,label=c(expression("C"["leaf"]),expression("N"["leaf"]),expression("P"["leaf"]),expression("C:N"["leaf"]),expression("C:P"["leaf"]),"LA","SLA","LTD"))+
  theme(legend.position = c(0.27,0.14685291425),plot.title=element_text(hjust=0.5,face="bold",size=20),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+ylab("Second PCA axis (14.2%)")+
  xlim(-1.2,1.3)+
  labs(title="Phylogenetic PCA")

c=ggplot()+geom_point(data=pca_39_root_spe_8,alpha=0.8,size=4,aes(x=Comp.1  ,y=Comp.2 ,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+geom_segment(data=pca_39_root_trait_8,size=0.7,color="orange",aes(x=0,y=0,xend=1.5*Comp.1  ,yend=1.5*Comp.2),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (50.8%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+geom_text(data=pca_39_leaf_trait_8,x=1.4*c(  0.0511284269,  0.48237202168,  0.62358290728451,  0.16087581,  0.01510978, -0.32926412, -0.533551851,-0.6051703797),y=1.5*c(0.41204451,  0.0710834744, -0.06863963 , 0.66892106380, -0.6671842613,  0.160019876,  0.02068744,   0.158605412),size=5,label=c(expression("C"["root"]),expression("N"["root"]),expression("P"["root"]),"AD","SRL","RTD",expression("C:N"["root"]),expression("C:P"["root"])))+theme(legend.position = c(0.27,0.14685291425),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+
  ylim(-1,1.1)+xlim(-1.1,1)+ylab("Second PCA axis (20.8%)")


d=ggplot()+geom_point(data=pca_39_root_spe_8_phy,alpha=0.8,size=4,aes(x=PC1  ,y=PC2 ,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_root_trait_8_phy,size=0.7,color="orange",aes(x=0,y=0,xend=1.5*PC1  ,yend=1.5*PC2),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (49.4%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+
  geom_text(data=pca_39_root_trait_8_phy,x=1.4*c(  -0.01057971936,  0.77680164,  1.052189006605, -0.13255315,  0.46821878, -0.7975879197, -1.034893892171,-1.0695188928215),y=1.5*c(-0.5542788104, -0.2119262746, -0.20892011, -0.91388428653,  0.850727107, -0.09277834,  0.03516724, 0.11818062),size=5,label=c(expression("C"["root"]),expression("N"["root"]),expression("P"["root"]),"AD","SRL","RTD",expression("C:N"["root"]),expression("C:P"["root"])))+
  theme(legend.position = c(0.27,0.14685291425),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+ylab("Second PCA axis (19.8%)")+
  xlim(-1.57,1.5)+
  ylim(-1.5,1.3)







e=ggplot()+geom_point(data=pca_39_ab_spe_8,alpha=0.8,size=4,aes(x=PC1 ,y=PC2,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_ab_trait_8,size=0.7,aes(x=0,y=0,xend=1.5*pc1 ,yend=1.5*pc2,color=k,linetype=k),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (34.2%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+
  theme(legend.position = c(0.77,0.16454685291425),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+
  ylab("Second PCA axis (23.6%)")+
  scale_color_manual("",breaks=c("ac","am","bc","bm"),labels=c("Leaf chemical trait","Leaf morphorlogical trait","Root chemical trait","Root morphorlogical trait"),values=c("green3","green3","orange","orange"))+
  scale_linetype_manual("",breaks=c("ac","am","bc","bm"),labels=c("Leaf chemical trait","Leaf morphorlogical trait","Root chemical trait","Root morphorlogical trait"),values=c("solid","dashed","solid","dashed"))+
  guides(fill=FALSE)+
  guides(shape=FALSE)+
  ylim(-1.2,1)+
  geom_text(data=pca_39_ab_trait_8,x=1.4*c( 0.12416621, -0.396229094920, -0.296319651747,  0.41529384221,  0.35293906320, -0.114118030, -0.29592831,  0.24627773, -0.04974902, -0.27622523, -0.36791911,  0.02401138, -0.13765149,  0.342035262798,  0.34889074,  0.5158559784),y=1.5*c(0.330185507, -0.21078677, -0.01941585,  0.23120650,  0.08792517,  0.2050252451, -0.317858961,  0.2794034040,  0.295708319,  0.3297777324,  0.296607814,  0.46664557, -0.362884091, -0.11687065, -0.263356287, -0.21317355),size=5,label=c(expression("C"["leaf"]),expression("N"["leaf"]),expression("P"["leaf"]),expression("C:N"["leaf"]),expression("C:P"["leaf"]),"LA","SLA","LTD",expression("C"["root"]),expression("N"["root"]),expression("P"["root"]),"AD","SRL","RTD",expression("C:N"["root"]),expression("C:P"["root"])),color=c("green3","green3","green3","green3","green3","green3","green3","green3","orange","orange","orange","orange","orange","orange","orange","orange"))



f=ggplot()+geom_point(data=pca_39_ab_spe_8_phy,alpha=0.8,size=4,aes(x=PC1 ,y=PC2,shape=m2,fill=m2))+
  geom_vline(xintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_hline(yintercept = 0,color="gray",linetype="dotted",size=1)+
  geom_segment(data=pca_39_ab_trait_8_phy,size=0.7,aes(x=0,y=0,xend=PC1 ,yend=PC2,color=k,linetype=k),arrow=arrow(length=unit(0.3,"cm")))+
  xlab("First PCA axis (36.4%)")+
  scale_shape_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c(22,23,24))+
  scale_fill_manual("",breaks=c("Asterids","Magnolids","Superrosids"),labels=c("Asterids","Magnolids","Superrosids"),values=c("red","orange","black"))+
  theme(legend.position = c(0.27,0.16685291425),panel.background = element_blank(), panel.border = element_rect(fill=NA,size=2,color="black"),legend.text = element_text(size=12),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),axis.text.x  = element_text(size=15))+
  ylab("Second PCA axis (19.4%)")+
  scale_color_manual("",breaks=c("ac","am","bc","bm"),labels=c("Leaf chemical trait","Leaf morphorlogical trait","Root chemical trait","Root morphorlogical trait"),values=c("green3","green3","orange","orange"))+
  scale_linetype_manual("",breaks=c("ac","am","bc","bm"),labels=c("Leaf chemical trait","Leaf morphorlogical trait","Root chemical trait","Root morphorlogical trait"),values=c("solid","dashed","solid","dashed"))+
  guides(fill=FALSE)+
  guides(shape=FALSE)+
  ylim(-1.3,0.8)+
  geom_text(data=pca_39_ab_trait_8,x=c(  0.39736604403 ,-0.92863781984278, -0.82793699749662,  0.848835120 , 1.02811977918,  0.001199331, -0.92890774932219,  0.68535953443359, -0.11075566401, -0.609900142, -0.583886146,  0.328616564506, -0.400817876,  0.5391814379,  0.83783676865202,  0.555492271),y=c(-0.34128118286, -0.5853894645, -0.3373091338 , 0.4651060899 , 0.29624632393164 , 0.1752765658, -0.46109390,  0.06526120, -0.282789438,  0.4185770103,  0.7530398446,  0.04794201,  0.32296595764, -0.5247736080, -0.4751555875, -0.79974864556),size=5,label=c(expression("C"["leaf"]),expression("N"["leaf"]),expression("P"["leaf"]),expression("C:N"["leaf"]),expression("C:P"["leaf"]),"LA","SLA","LTD",expression("C"["root"]),expression("N"["root"]),expression("P"["root"]),"AD","SRL","RTD",expression("C:N"["root"]),expression("C:P"["root"])),color=c("green3","green3","green3","green3","green3","green3","green3","green3","orange","orange","orange","orange","orange","orange","orange","orange"))