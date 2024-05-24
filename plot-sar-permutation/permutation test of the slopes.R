# permutation test for the slope to account for non-indepence of the data

a5=unique(point_number_three_neon$Var1)

p_value=numeric()
for (i in 1:length(a5)){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
d=subset(sar_neon_permutation,plotid==a5[i])

model <- lm(log(richness) ~log(area), data = d)

observed_slope <- coef(model)["log(area)"]

num_permutations <- 1000 
permuted_slopes <- numeric(num_permutations)
for (j in 1:num_permutations) {
  
  # Permute the independent variable
  perm_data <- d
  perm_data$area <- sample(perm_data$area)
  
  # Fit a regression model and store the slope
  perm_model <- lm(log(richness) ~ log(area), data = perm_data)
  permuted_slopes[j] <- coef(perm_model)["log(area)"]
}
p_value[i]=mean(abs(permuted_slopes) >= abs(observed_slope))
}


