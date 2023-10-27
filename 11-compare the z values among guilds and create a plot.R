# step 11 compare the mean of the z values among guilds, which will be based on the data generated in step 10
# look at how many observations for each guild
library(car)
source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")
dim(subset(acm_z_ranall, X1 < 10)) # acm,365 plots,X1<10 means that we excluded plots with <3 soil cores
dim(subset(ecm_z_ranall, X1 < 10)) # acm,493 plots
dim(subset(litsap_z_ranall, X1 < 10)) # litter sap,493 plots
dim(subset(para_z_ranall, X1 < 10)) # parasitic,459 plots
dim(subset(soilsap_z_ranall, X1 < 10)) # parasitic,493 plots
dim(subset(woosap_z_ranall, X1 < 10)) # parasitic,482 plots
dim(subset(epiphy_z_ranall, X1 < 10)) # parasitic,482 plots
dim(subset(plapat_z_ranall, X1 < 10)) # parasitic,482 plots
# we now have eight guilds and need to compare the mean of z among guilds

a <- list()
for (i in 1:dim(acm_z_ranall)[1])
{
  a[[i]] <- data.frame(t(acm_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(acm_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(acm_z_ranall$a1, each = 30)

acm_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

acm_z_ranall_30$X1 <- as.numeric(acm_z_ranall_30$X1) # the data will also be used for model

### for "ecm" guild
a <- list()
for (i in 1:dim(ecm_z_ranall)[1])
{
  a[[i]] <- data.frame(t(ecm_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(ecm_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(ecm_z_ranall$a1, each = 30)

ecm_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "litsap" guild
a <- list()
for (i in 1:dim(litsap_z_ranall)[1])
{
  a[[i]] <- data.frame(t(litsap_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(litsap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(litsap_z_ranall$a1, each = 30)

litsap_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "para" guild
a <- list()
for (i in 1:dim(para_z_ranall)[1])
{
  a[[i]] <- data.frame(t(para_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(para_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(para_z_ranall$a1, each = 30)

para_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "soilsap"
a <- list()
for (i in 1:dim(soilsap_z_ranall)[1])
{
  a[[i]] <- data.frame(t(soilsap_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(soilsap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(soilsap_z_ranall$a1, each = 30)

soilsap_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "woodsap" guild
a <- list()
for (i in 1:dim(woosap_z_ranall)[1])
{
  a[[i]] <- data.frame(t(woosap_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(woosap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(woosap_z_ranall$a1, each = 30)

woosap_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "epiphy" guild
a <- list()
for (i in 1:dim(epiphy_z_ranall)[1])
{
  a[[i]] <- data.frame(t(epiphy_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(epiphy_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(epiphy_z_ranall$a1, each = 30)

epiphy_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# for "plapat" guild
a <- list()
for (i in 1:dim(plapat_z_ranall)[1])
{
  a[[i]] <- data.frame(t(plapat_z_ranall[i, 3:32]))
}

b <- a[[1]]
for (i in 2:dim(plapat_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(plapat_z_ranall$a1, each = 30)

plapat_z_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values
# combine all the data
# guild=rep(c("acm","ecm","litsap","para","soilsap","woosap","epiphy","plapat"),times=c(dim(acm_z_ranall_30)[1],dim(ecm_z_ranall_30)[1],dim(litsap_z_ranall_30)[1],dim(para_z_ranall_30)[1],dim(soilsap_z_ranall_30)[1],dim(woosap_z_ranall_30)[1],dim(epiphy_z_ranall_30)[1],dim(plapat_z_ranall_30)[1]))
names(acm_z_ranall_30)[2] <- "z"
names(ecm_z_ranall_30)[2] <- "z"
names(litsap_z_ranall_30)[2] <- "z"
names(para_z_ranall_30)[2] <- "z"
names(soilsap_z_ranall_30)[2] <- "z"
names(woosap_z_ranall_30)[2] <- "z"
names(epiphy_z_ranall_30)[2] <- "z"
names(plapat_z_ranall_30)[2] <- "z"

# combine all the eight guilds
com_guild <- rbind(acm_z_ranall_30, ecm_z_ranall_30, litsap_z_ranall_30, para_z_ranall_30, soilsap_z_ranall_30, woosap_z_ranall_30, epiphy_z_ranall_30, plapat_z_ranall_30)
guild <- rep(c("acm", "ecm", "litsap", "para", "soilsap", "woosap", "epiphy", "plapat"), times = c(dim(acm_z_ranall_30)[1], dim(ecm_z_ranall_30)[1], dim(litsap_z_ranall_30)[1], dim(para_z_ranall_30)[1], dim(soilsap_z_ranall_30)[1], dim(woosap_z_ranall_30)[1], dim(epiphy_z_ranall_30)[1], dim(plapat_z_ranall_30)[1]))
com_guild <- cbind(com_guild, guild)
com_guild <- subset(com_guild, z < 10) # exclude the plots with <3 core

write.csv(com_guild, "com_guild.csv") # the data

leveneTest(z ~ guild, data = com_guild) # testing variance homogenety, unblanced

oneway.test(z ~ guild, data = com_guild, na.action = na.omit, var.equal = FALSE)

source("http://aoki2.si.gunma-u.ac.jp/R/src/tukey.R", encoding = "euc-jp")

tukey(com_guild$z, com_guild$guild, method = "G") # multiple comparison with 'Games.Howell' Post-hoc Test


boxplot(z ~ guild, data = subset(com_guild, z < 10))


k <- aggregate(z ~ guild, data = com_guild, FUN = mean)
od1 <- k[order(k$z), ] # with the increase trend to display the box plots
com_guild$guild <- factor(com_guild$guild, levels = od1$guild)

b <- ggboxplot(com_guild, x = "guild", y = "z", fill = "guild", outlier.colour = "gray", outlier.shape = NA) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  geom_hline(yintercept = 0.787, linetype = "dashed", color = "red") +
  xlab("") +
  theme(legend.position = c(0.51, 0.85), legend.text = element_text(size = 14), text = element_text(size = 15), axis.text.x = element_blank(), axis.title.y = element_text(face = "italic", size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank()) +
  scale_fill_manual("", breaks = od1$guild, values = c("chocolate1", "gray", "cadetblue1", "lavender", "greenyellow", "forestgreen", "purple", "tan"), labels = c("soil saprotroph(N=493)", "parasite(N=459)", "wood saprotroph(N=482)", "epiphyte(N=482)", "plant pathogen(N=482)", "litter saprotroph(N=493)", "ECM(N=493)", "ACM(N=482)")) +
  annotate(x = 1, y = 1.2, "text", label = "g", size = 6) +
  annotate(x = 2, y = 1.5, "text", label = "f", size = 6) +
  annotate(x = 3, y = 1.56, "text", label = "e", size = 6) +
  annotate(x = 4, y = 1.53, "text", label = "e", size = 6) +
  annotate(x = 5, y = 1.52, "text", label = "de", size = 6) +
  annotate(x = 6, y = 1.4, "text", label = "c", size = 6) +
  annotate(x = 7, y = 1.35, "text", label = "b", size = 6) +
  annotate(x = 8, y = 1.8, "text", label = "a", size = 6) +
  ylim(0, 3)
