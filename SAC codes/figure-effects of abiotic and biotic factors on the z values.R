# effect_all_clima_soil

library(ggplot2)

a <- ggplot() +
  geom_segment(data = effect_all_clima_soil[3:15, ], color = "tan", size = 0.8, aes(y = 1:13, yend = 1:13, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_all_clima_soil[3:15, ], pch = 21, aes(x = Estimate, y = 1:13, fill = sig), size = 4) +
  scale_y_continuous(breaks = 1:13, label = c("SoilOC", "pH", "SoilN", "Sand", "Bio2", "Bio8", "Bio18", "Tem.seas.", "MAP", "Pre.seas.", "SPEI", "Fun.rich", "MAT")) +
  xlab("Effect size and 95% CI") +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("gray", "seagreen")) +
  ggtitle("Clima.+Soil (N=483)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(size = 15, hjust = 0.5), text = element_text(size = 18), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))
# for the NEON site
# effect_neon_clima_soil
 # when explanatory variables were grouped with different colors
effect_all_clima_soil=cbind(vabe=rownames(effect_all_clima_soil),effect_all_clima_soil)
effect_all_clima_soil$vabe=factor(effect_all_clima_soil$vabe,levels=c("(Intercept)" , "c"  ,"bio1" , "bio2","bio8" , "bio18" , "bio4" , "bio12" , "bio15" , "spei",  "organicCPercent", "ph" ,"nitrogen" , "sand" , "funrich" ))

ggplot() +
  geom_segment(data = effect_all_clima_soil[3:15, ], color = "gray", size = 0.8, aes(y = vabe, yend = vabe, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_all_clima_soil[3:15, ], pch = 21, aes(x = Estimate, y = vabe, fill = sig), size = 4) +

  scale_y_discrete (breaks = c("bio1" , "bio2","bio8" , "bio18" , "bio4" , "bio12" , "bio15" , "spei",  "organicCPercent", "ph" ,"nitrogen" , "sand" , "funrich" ), label = c("MAT" , "Bio2","Bio8" , "Bio18" , "Tem.seas." , "MAP" , "Pre.seas." , "SPEI",  "SoilC", "pH" ,"SoilN" , "Sand" , "Fun.rich" )) +
  xlab("Effect size and 95% CI") +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("gray", "mediumpurple")) +
  ggtitle("Clima.+Soil (N=483)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(size = 15, hjust = 0.5), text = element_text(size = 18), axis.text.y = element_text(hjust = 0,color=c("blue","blue","blue","blue","blue","blue","blue","blue","tan","tan","tan","tan","seagreen")), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))

  


b <- ggplot() +
  geom_segment(data = effect_neon_clima_soil[c(2, 4:16), ], color = "tan", size = 0.8, aes(y = 1:14, yend = 1:14, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_neon_clima_soil[c(2, 4:16), ], pch = 21, aes(x = Estimate, y = 1:14, fill = sig), size = 4) +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("gray", "seagreen")) +
  scale_y_continuous(breaks = 1:14, label = c("SoilOC", "pH", "SoilN", "Pla.rich.", "Sand", "Bio2", "Bio8", "Bio18", "Tem.seas.", "MAP", "Pre.seas.", "SPEI", "Fun.rich", "MAT")) +
  xlab("Effect size and 95% CI") +
  ggtitle("Clima.+Soil+Plant (N=396)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(size = 15, hjust = 0.5), text = element_text(size = 18), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))

3. ## when explanatory variables were grouped with different colors

effect_neon_clima_soil=cbind(vabe=rownames(effect_neon_clima_soil),effect_neon_clima_soil)
effect_neon_clima_soil$vabe=factor(effect_neon_clima_soil$vabe,levels=c("(Intercept)","c","bio1","bio2","bio8","bio18","bio12","bio15","bio4","spei","organicCPercent","ph","nitrogen","sand","funrich","rich"))
ggplot() +
  geom_segment(data = effect_neon_clima_soil[c(2,4:16), ], color = "gray", size = 0.8, aes(y = vabe, yend = vabe, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_neon_clima_soil[c(2,4:16), ], pch = 21, aes(x = Estimate, y = vabe, fill = sig), size = 4) +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("white", "mediumpurple")) +
  scale_y_discrete(breaks =  c("(Intercept)","c","bio1","bio2","bio8","bio18","bio12","bio15","bio4","spei","organicCPercent","ph","nitrogen","sand","funrich","rich"),label=c("(Intercept1)","c","MAT","Bio2","Bio8","Bio18","MAP","Pre.seas.","Tem.seas.","SPEI","SoilC","pH","SoilN","Sand","Fun.rich","Pla.rich") ) +
  xlab("Effect size and 95% CI") +
  ggtitle("Clima.+Soil+Plant (N=396)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(size = 15, hjust = 0.5), text = element_text(size = 18), axis.text.y = element_text(hjust = 0,color=c("blue","blue","blue","blue","blue","blue","blue","blue","tan","tan","tan","tan","seagreen","seagreen")), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))


4. # when the plotid was included as the only random effect model
# reorder the variables 
effect_neon_clima_soil_plotran=cbind(vabe=rownames(effect_neon_clima_soil_plotran),effect_neon_clima_soil_plotran)
effect_neon_clima_soil_plotran$vabe=factor(effect_neon_clima_soil_plotran$vabe,levels=c("(Intercept)","c","bio1","bio2","bio8","bio18","bio12","bio15","bio4","spei","organicCPercent","ph","nitrogen","sand","funrich","rich"))

ggplot() +
  geom_segment(data = effect_neon_clima_soil_plotran[c(2,4:16), ], color = "gray", size = 0.8, aes(y = vabe, yend = vabe, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_neon_clima_soil_plotran[c(2,4:16), ], pch = 21, aes(x = Estimate, y = vabe, fill = sig), size = 4) +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("white", "mediumpurple")) +
  scale_y_discrete(breaks =  c("(Intercept)","c","bio1","bio2","bio8","bio18","bio12","bio15","bio4","spei","organicCPercent","ph","nitrogen","sand","funrich","rich"),label=c("(Intercept1)","c","MAT","Bio2","Bio8","Bio18","MAP","Pre.seas.","Tem.seas.","SPEI","SoilC","pH","SoilN","Sand","Fun.rich","Pla.rich") ) +
  xlab("Effect size and 95% CI") +
  ggtitle("Clima.+Soil+Plant (N=396)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", plot.title = element_text(size = 15, hjust = 0.5), text = element_text(size = 18), axis.text.y = element_text(hjust = 0,color=c("blue","blue","blue","blue","blue","blue","blue","blue","tan","tan","tan","tan","seagreen","seagreen")), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))



# climate, soil plant and root model

c <- ggplot() +
  geom_segment(data = effect_neon_clima_soil_root[c(2, 4:19), ], color = "tan", size = 0.8, aes(y = 1:17, yend = 1:17, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_neon_clima_soil_root[c(2, 4:19), ], pch = 21, aes(x = Estimate, y = 1:17, fill = sig), size = 4) + # pch=21

  scale_y_continuous(breaks = 1:17, label = c("SoilOC", "pH", "SoilN", "Pla.rich.", "Sand", "Bio2", "Bio8", "Bio18", "MAP", "Pre.seas.", "SPEI", "Fun.rich", "MAT", expression(Root[mass]), expression(Root[d13C]), expression(Root[C]), expression(Root[C:N]))) +
  xlab("Effect size and 95% CI") +
  ggtitle("Clima.+Soil+Plant+Root \n(N=104)") +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("gray", "seagreen")) +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))
# when plotid was included as the only random effect

effect_neon_clima_soil_root_plotran=cbind(vabe=rownames(effect_neon_clima_soil_root_plotran),effect_neon_clima_soil_root_plotran)
effect_neon_clima_soil_root_plotran$vabe=factor(effect_neon_clima_soil_root_plotran$vabe,levels=c("(Intercept)" , "c" , "bio1" ,"bio2", "bio8", "bio18", "bio12" , "bio15", "spei" ,"organicCPercent",  "ph"  , "nitrogen" , "sand" , "funrich" ,"rich" , "fine", "d13C" , "rootc", "rootcn" ))

ggplot() +
  geom_segment(data = effect_neon_clima_soil_root_plotran[c(2,4:19), ], color = "gray", size = 0.8, aes(y = vabe, yend = vabe, x = Estimate - 1.96 * Std..Error, xend = Estimate + 1.96 * Std..Error)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1.2) +
  geom_point(data = effect_neon_clima_soil_root_plotran[c(2,4:19), ], pch = 21, aes(x = Estimate, y = vabe, fill = sig), size = 4) +
  scale_fill_manual("", breaks = c("no", "sig"), values = c("white", "mediumpurple"))+

  scale_y_discrete (breaks =c( "bio1" ,"bio2", "bio8", "bio18", "bio12" , "bio15", "spei" ,"organicCPercent",  "ph"  , "nitrogen" , "sand" ,"funrich" ,"rich" , "fine", "d13C" , "rootc", "rootcn") , label = c("MAT", "Bio2", "Bio8", "Bio18", "MAP", "Pre.seas.", "SPEI", "SoilC", "pH", "SoilN" , "Sand","Fun.rich", "Pla.rich", expression(Root[mass]), expression(Root[d13C]), expression(Root[C]), expression(Root[C:N]))) +
  xlab("Effect size and 95% CI") +
  ggtitle("Clima.+Soil+Plant+Root (N=104)") +
  ylab("") +
  guides(fill = "none") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0,color=c("blue","blue","blue","blue","blue","blue","blue","tan","tan","tan","tan","seagreen","seagreen","orange","orange","orange","orange")), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(face = "italic", size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))

  
