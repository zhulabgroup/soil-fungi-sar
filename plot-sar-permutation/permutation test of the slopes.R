# permutation test for the slope to account for non-independence of the data


full_dob_neon_richness_data_rarefaction%>%filter(!is.na(mean_value)&guild=="all")%>%
  group_by(plotid)%>%summarise(n())%>%filter(`n()`>3)%>%pull(plotid)->plotid_with_four_point

p_value=numeric()
for (i in 1:length(plotid_with_four_point)){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
d=subset(full_dob_neon_richness_data_rarefaction%>%filter(guild=="all"),plotid==plotid_with_four_point[i])

model <- lm(log(mean_value) ~log(area), data = d)

observed_slope <- coef(model)["log(area)"]

num_permutations <- 10000
permuted_slopes <- numeric(num_permutations)
for (j in 1:num_permutations) {
  
  # Permute the independent variable
  perm_data <- d
  perm_data$area <- sample(perm_data$area)
  
  # Fit a regression model and store the slope
  perm_model <- lm(log(mean_value) ~ log(area), data = perm_data)
  permuted_slopes[j] <- coef(perm_model)["log(area)"]
}
p_value[i]=mean(abs(permuted_slopes) >= abs(observed_slope))
}

p.adjust(p_value,"BH")->adjust_p

k=data.frame(p=adjust_p,plot=1:469)

k1=data.frame(p=p_value,plot=1:469)

###

p1=ggplot(data=k1, aes(x =p_value )) +
  geom_histogram(binwidth = 0.0012, fill = "#FFCB05", color = "", alpha = 0.7) +
  ggtitle("") +
  xlab(expression(italic(P)*" value")) +
  ylab("Frequency") +
  geom_vline(xintercept = 0.05,linetype="dashed",color="red")+
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
        panel.border = element_rect(color = "black", size = 1, fill = NA))
  
p2=ggplot(data=k, aes(x =p )) +
  geom_histogram(binwidth = 0.0012, fill = "#FFCB05", color = "black", alpha = 0.7) +
  ggtitle("") +
  xlab(expression(italic(P)*" value")) +
  ylab("Frequency") +
  geom_vline(xintercept = 0.05,linetype="dashed",color="red")+
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
        panel.border = element_rect(color = "black", size = 1, fill = NA))

