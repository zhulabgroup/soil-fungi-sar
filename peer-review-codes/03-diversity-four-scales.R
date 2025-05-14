## use the extrapolation approach to determine the 40 m scale richness for the neon site

## need to confirm that the gx and gy are in the format of numeric
load(file = "rare_all_assign.RData")
load(file = "dob_permutation.RData")

p <- sample_data(rare_all_assign) %>%
  data.frame() %>%
  select(gx, gy)
neon_data <- sample_data(rare_all_assign) %>% data.frame()
neon_data$gx <- NULL
neon_data$gy <- NULL
sample_data(rare_all_assign) <- sample_data(neon_data)
p$gx <- as.numeric(p$gx)
p$gy <- as.numeric(p$gy)
row.names(p) <- row.names(sample_data(rare_all_assign))
p <- sample_data(p)
rare_all_assign <- merge_phyloseq(rare_all_assign, p)
# replace the initial object
#####
ori <- expand.grid(x = seq(0, 10, 0.5), y = seq(0, 10, 0.5))
#

# for different guilds
# to add a new column called guild to the intial data
# these several groups were grouped into one broad group of "para"
taxa_df <- as.data.frame(tax_table(rare_all_assign))
guild <- tax_table(rare_all_assign) %>%
  data.frame() %>%
  select(ta2)
term <- c("protistan_parasite", "lichen_parasite", "algal_parasite", "mycoparasite", "animal_parasite")

guild %>%
  mutate(ta2 = ifelse(ta2 %in% term, "para", ta2)) %>%
  rename(guild = ta2) -> guild
taxa_df$guild <- guild
new_tax_table <- as.matrix(taxa_df)
tax_table(rare_all_assign) <- new_tax_table

guild_select <- c(
  "ectomycorrhizal", "arbuscular_mycorrhizal", "soil_saprotroph", "litter_saprotroph", "plant_pathogen", "wood_saprotroph",
  "para", "epiphyte"
)

## estimating fungal diversity at the 10 m by 10 m scale for the neon plots


data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(rare_all_assign, guild %in% my_list[[m]])


  set.seed(1201)
  times <- 30
  a4 <- sample_data(data_select)
  a4 <- unique(a4$subplotID10)
  richness <- vector("list", length(a4))
  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, subplotID10 == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]
    if (dim1 >= 2) {
      cl <- makeCluster(3)
      registerDoParallel(cl)
      richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
        species <- vector(length = dim1[1]) # create a vector to save diversity
        for (j in 1:1)
        {
          # randomly sample j samples in the plot
          flag <- rep(FALSE, dim1[1])
          flag[sample(1:dim1[1], 1)] <- TRUE
          temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
          species[j] <- sum(otu_table(temp)["TRUE"] > 0)
          return(species[j])
        }
      }
      stopCluster(cl)
    } else {
      richness[[i]] <- table(colSums(otu_table(data_sub)) > 0)["TRUE"] %>% as.numeric()
    }
  }

  # get the diversity for each 10 x 10 m subplot for the NEON sites

  richness_subplot10_neon <- data.frame(nrow = length(a4), ncol = 2) # the mean indicates the mean value for the cores within the same subplot
  for (i in 1:length(a4))
  {
    ak <- dim(richness[[i]])[1]

    if (is.null(ak)) {
      richness_subplot10_neon[i, 1] <- richness[[i]]
      richness_subplot10_neon[i, 2] <- richness[[i]]
    } else {
      richness_subplot10_neon[i, 1] <- mean(richness[[i]])
      richness_subplot10_neon[i, 2] <- sd(richness[[i]])
    }
  }

  richness_subplot10_neon %>%
    mutate(plotid = substr(a4, 1, 8)) %>%
    group_by(plotid) %>%
    summarise(mean_value = mean(nrow, na.rm = TRUE), sd_value = sd(nrow, na.rm = TRUE), area = rep(100, length.out = n)) %>%
    mutate(guild = rep(my_list[[m]], n())) %>%
    select(plotid, mean_value, sd_value, area, guild) -> data[[m]]
}

do.call(rbind, data) -> data_neon_10m_scale

save(data_neon_10m_scale, file = "data_neon_10m_scale.RData")



## diversity estimation at the 20 m by 20 m scale for the NEON sites

data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(rare_all_assign, guild %in% my_list[[m]])

  set.seed(1203)
  times <- 30
  a4 <- sample_data(data_select)
  a4 <- unique(a4$subplotID20)
  richness <- vector("list", length(a4))

  for (i in 1:length(a4)) # note that the i here is large
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, subplotID20 == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]

    tryCatch(
      {
        if (dim1 >= 4) {
          cl <- makeCluster(3)
          registerDoParallel(cl)
          richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
            species <- vector(length = dim1[1]) # create a vector to save diversity
            for (j in 1:4)
            {
              # randomly sample j samples in the plot
              flag <- rep(FALSE, dim1[1])
              flag[sample(1:dim1[1], 4)] <- TRUE # select 3 cores
              temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
              species[j] <- sum(otu_table(temp)["TRUE"] > 0)
              return(species[j])
            }
          }
          stopCluster(cl)
        } else if (dim1 >= 2 & dim1 < 4) {
          sub_samp <- rownames(sample_data(data_sub))
          kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
          data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
          ot <- t(otu_table(data_sub)) %>% data.frame()
          n <- dim(ot)[2] # the number of "sites" for each plot
          ot1 <- rowSums(ot)
          out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 4, length.out = 4)), se = FALSE)
          richness[[i]] <- out3$iNextEst$size_based %>%
            slice_tail() %>%
            pull(qD)
        } else {
          richness[[i]] <- NA
        }
      },
      error = function(e) {
        # Handle the error and continue
        print(paste("Error on element", i, ":", e))
      },
      finally = {
        print(paste("Processed element", i))
      }
    )
  }


  richness_mean_subplot20_neon <- matrix(ncol = 2, nrow = length(a4))
  for (i in 1:length(a4))
  {
    richness_mean_subplot20_neon[i, 1] <- mean(richness[[i]])
    richness_mean_subplot20_neon[i, 2] <- sd(richness[[i]])
  }

  richness_mean_subplot20_neon %>%
    data.frame() %>%
    bind_cols(a4) %>%
    mutate(plotid = substr(a4, 1, 8)) %>%
    group_by(plotid) %>%
    summarise(mean_value = mean(X1, na.rm = TRUE), sd_value = sd(X1, na.rm = TRUE)) %>%
    mutate(area = rep(400, n()), guild = rep(my_list[[m]], n())) -> data[[m]]
}

# bind all the fungal guilds

do.call(rbind, data) %>%
  rename_all(~ paste0(c("plotid", "mean_value", "sd_value", "area", "guild"))) -> data_neon_20m_scale

save(data_neon_20m_scale, file = "data_neon_20m_scale.RData")



# data for the 30 by 30 m scale was computed on the great lakes
# used to search for a 30 m^2 square in the plot
# for individual guilds
ori <- expand.grid(x = seq(0, 10, 0.5), y = seq(0, 10, 0.5))


data <- list()
for (m in 1:9)
{
  data_select <- subset_taxa(rare_all_assign, guild %in% my_list[[m]])

  plot_ID <- sample_data(data_select) %>% data.frame()
  plot_ID <- unique(plot_ID$plotIDM)

  set.seed(520)
  richness_plot_neon <- matrix(ncol = 3, nrow = length(plot_ID))
  for (k in 1:length(plot_ID))
  {
    cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing

    tryCatch(
      {
        core_number <- numeric()
        for (i in 1:dim(ori)[1])
        {
          # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing

          d <- sample_data(data_select) %>%
            data.frame() %>%
            select(gx, gy, plotIDM) %>%
            filter(plotIDM == plot_ID[k])

          core_number[i] <- subset(d, gx >= ori$x[i] & gx <= ori$x[i] + 30 & gy <= ori$y[i] + 30 & gy >= ori$y[i]) %>%
            summarize(row_count = n()) %>%
            as.numeric()
        }
        core_number_8 <- core_number > 8
        core_number_8[!core_number_8] <- NA

        richness <- numeric()
        for (j in 1:length(core_number)) # only select the locations that over 9 cores can be sampled/81 cores
        {
          cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
          tryCatch(
            {
              if (j %in% c(which(!is.na(core_number_8)))) {
                # informs the processing
                sub_samp <- subset(d, gx >= ori$x[j] & gx <= ori$x[j] + 30 & gy <= ori$y[j] + 30 & gy >= ori$y[j])
                sub_samp <- sample(rownames(sub_samp), 9) # select only nine cores
                kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
                richness[j] <- table(colSums(otu_table(kk)) > 0)["TRUE"] %>% as.numeric()
              } else {
                sub_samp <- subset(d, gx >= ori$x[j] & gx <= ori$x[j] + 30 & gy <= ori$y[j] + 30 & gy >= ori$y[j])
                if (dim(sub_samp)[1] < 2) {
                  richness[j] <- NA
                } else {
                  sub_samp <- rownames(sub_samp) # get the sample names
                  kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp) # based on this to extrapolate the richness
                  data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
                  ot <- t(otu_table(data_sub)) %>% data.frame()
                  n <- dim(ot)[2] # the number of "sites" for each plot
                  ot1 <- rowSums(ot)
                  out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 9, length.out = 9)), se = FALSE)
                  richness[j] <- out3$iNextEst$size_based %>%
                    slice_tail() %>%
                    pull(qD)
                }
              }
            },
            error = function(e) {
              # Handle the error and continue
              print(paste("Error on element", j, ":", e))
            },
            finally = {
              print(paste("Processed element", j))
            }
          )
        }

        richness_plot_neon[k, 1] <- plot_ID[k]
        richness_plot_neon[k, 2] <- richness %>% mean(na.rm = TRUE)
        richness_plot_neon[k, 3] <- richness %>% sd(na.rm = TRUE)
      },
      error = function(e) {
        # Handle the error and continue
        print(paste("Error on element", k, ":", e))
      },
      finally = {
        print(paste("Processed element", k))
      }
    )
  }

  cbind(richness_plot_neon, my_list[[m]]) %>%
    data.frame() %>%
    mutate(area = rep(900, n())) %>%
    select(X1, X2, X3, area, X4) %>%
    rename_all(~ paste0(c("plotid", "mean_value", "sd_value", "area", "guild"))) -> data[[m]]
}

do.call(rbind, data) -> data_neon_30m_scale

save(data_neon_30m_scale, file = "data_neon_30m_scale.RData")



## fungal diversity at the 40m by 40m scale for the none plots

data <- list()
for (m in 1:9)
{
  data_select <- subset_taxa(rare_all_assign, guild %in% my_list[[m]])


  set.seed(5679)
  times <- 30
  a4 <- sample_data(data_select)
  a4 <- unique(a4$plotIDM)
  richness <- vector("list", length(a4))

  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, plotIDM == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]

    if (dim1 >= 16) {
      cl <- makeCluster(3)
      registerDoParallel(cl)
      richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
        species <- vector(length = dim1[1]) # create a vector to save diversity
        for (j in 1:16)
        {
          # randomly sample j samples in the plot
          flag <- rep(FALSE, dim1[1])
          flag[sample(1:dim1[1], 16)] <- TRUE
          temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
          species[j] <- sum(otu_table(temp)["TRUE"] > 0)
          return(species[j])
        }
      }
      stopCluster(cl)
    } else if (dim1 >= 2 & dim1 <= 15) {
      sub_samp <- rownames(sample_data(data_sub))
      kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
      data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
      ot <- t(otu_table(data_sub)) %>% data.frame()
      n <- dim(ot)[2] # the number of "sites" for each plot
      ot1 <- rowSums(ot)
      out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 16, length.out = 16)), se = FALSE)
      richness[[i]] <- out3$iNextEst$size_based %>%
        slice_tail() %>%
        pull(qD)
    } else {
      richness[[i]] <- NA
    }
  }

  richness_mean_subplot40_neon <- matrix(ncol = 2, nrow = length(a4))
  for (i in 1:length(a4))
  {
    richness_mean_subplot40_neon[i, 1] <- mean(richness[[i]])
    richness_mean_subplot40_neon[i, 2] <- sd(richness[[i]])
  }


  cbind(
    plotid = a4, mean_value = richness_mean_subplot40_neon[, 1],
    sd_value = richness_mean_subplot40_neon[, 2],
    rep(1600, length.out = n)
  ) %>%
    data.frame() %>%
    mutate(guild = rep(my_list[[m]], n()))
  rename_all(~ paste0(c("plotid", "mean_value", "sd_value", "area"))) -> data[[m]]
}

do.call(rbind, data) -> data_neon_40m_scale

save(data_neon_40m_scale, file = "data_neon_40m_scale.rds")


###################### for the dob plots#######################

# these several groups were grouped into one broad group of "para"
taxa_df <- as.data.frame(tax_table(dob_permutation))
guild <- tax_table(dob_permutation) %>%
  data.frame() %>%
  select(ta2)
term <- c("protistan_parasite", "lichen_parasite", "algal_parasite", "mycoparasite", "animal_parasite")

guild %>%
  mutate(ta2 = ifelse(ta2 %in% term, "para", ta2)) %>%
  rename(guild = ta2) -> guild
taxa_df$guild <- guild
new_tax_table <- as.matrix(taxa_df)
tax_table(dob_permutation) <- new_tax_table


tax_table(dob_permutation) %>%
  data.frame() %>%
  pull(guild) %>%
  unique() -> all

guild_select <- c(
  "ectomycorrhizal", "arbuscular_mycorrhizal", "soil_saprotroph", "litter_saprotroph", "plant_pathogen", "wood_saprotroph",
  "para", "epiphyte"
)

my_list <- c(guild_select[1:8], list(all))

# dob at the 10 m by 10 m plot scale

data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(dob_permutation, guild %in% my_list[[m]])

  set.seed(2202)
  times <- 30
  a4 <- sample_data(data_select)

  a4 <- unique(a4$subplotID10)

  richness <- vector("list", length(a4))

  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(dob, subplotID10 == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]

    if (dim1 > 1) {
      cl <- makeCluster(3)
      registerDoParallel(cl)
      richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
        species <- vector(length = dim1[1]) # create a vector to save diversity
        for (j in 1:1)
        {
          # randomly sample j samples in the plot
          flag <- rep(FALSE, dim1[1])
          flag[sample(1:dim1[1], 1)] <- TRUE # select 3 cores
          temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
          species[j] <- sum(otu_table(temp)["TRUE"] > 0)
          return(species[j])
        }
      }
      stopCluster(cl)
    } else if (dim1 < 1) {
      richness[[i]] <- NA
    } else {
      richness[[i]] <- table(colSums(otu_table(data_sub)) > 0)["TRUE"] %>% as.numeric()
    }
  }

  # get the richness at the 10 by 10 m

  richness_subplot10_dob <- data.frame(nrow = length(a4), ncol = 2) # the mean indicates the mean value for the cores within the same subplot
  for (i in 1:length(a4))
  {
    ak <- dim(richness[[i]])[1]

    if (is.null(ak)) {
      richness_subplot10_dob[i, 1] <- richness[[i]]
      richness_subplot10_dob[i, 2] <- a4[i]
    } else {
      richness_subplot10_dob[i, 1] <- mean(richness[[i]])
      richness_subplot10_dob[i, 2] <- a4[i]
    }
  }

  richness_subplot10_dob %>%
    mutate(plotid = substr(richness_subplot10_dob$ncol, 1, 3)) %>%
    dplyr::rename(richness = nrow) %>%
    group_by(plotid) %>%
    summarise(mean_value = mean(richness, na.rm = TRUE), sd_value = sd(richness, na.rm = TRUE)) %>%
    mutate(area = rep(100, n()), guild = rep(my_list[[m]], n())) -> data[[m]]
}

do.call(rbind, data) %>%
  select(plotid, mean_value, sd_value, area, guild) %>%
  rename_all(~ paste0(c("mean_value", "sd_value", "plotid", "area", "guild"))) -> data_dob_10m_scale

save(data_dob_10m_scale, file = "data_dob_10m_scale.rds")



# dob at the 20m by 20 m scale


data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(dob_permutation, guild %in% my_list[[m]])

  set.seed(1203)
  times <- 30
  a4 <- sample_data(data_select)
  a4 <- unique(a4$subplotID20)
  richness <- vector("list", length(a4))

  for (i in 1:length(a4)) # note that the i here is large
  {
    # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, subplotID20 == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]

    tryCatch(
      {
        if (dim1 >= 4) {
          cl <- makeCluster(3)
          registerDoParallel(cl)
          richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
            species <- vector(length = dim1[1]) # create a vector to save diversity
            for (j in 1:4)
            {
              # randomly sample j samples in the plot
              flag <- rep(FALSE, dim1[1])
              flag[sample(1:dim1[1], 4)] <- TRUE # select 3 cores
              temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
              species[j] <- sum(otu_table(temp)["TRUE"] > 0)
              return(species[j])
            }
          }
          stopCluster(cl)
        } else if (dim1 >= 2 & dim1 < 4) {
          sub_samp <- rownames(sample_data(data_sub))
          kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
          data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
          ot <- t(otu_table(data_sub)) %>% data.frame()
          n <- dim(ot)[2] # the number of "sites" for each plot
          ot1 <- rowSums(ot)
          out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 4, length.out = 4)), se = FALSE)
          richness[[i]] <- out3$iNextEst$size_based %>%
            slice_tail() %>%
            pull(qD)
        } else {
          richness[[i]] <- NA
        }
      },
      error = function(e) {
        # Handle the error and continue
        print(paste("Error on element", i, ":", e))
      },
      finally = {
        print(paste("Processed element", i))
      }
    )
  }

  richness_mean_subplot20_dob <- matrix(ncol = 2, nrow = length(a4))
  for (i in 1:length(a4))
  {
    richness_mean_subplot20_dob[i, 1] <- mean(richness[[i]])
    richness_mean_subplot20_dob[i, 2] <- sd(richness[[i]])
  }

  richness_mean_subplot20_dob %>%
    data.frame() %>%
    bind_cols(a4) %>%
    mutate(plotid = substr(a4, 1, 3)) %>%
    group_by(plotid) %>%
    summarise(mean_value = mean(X1, na.rm = TRUE), sd_value = sd(X1, na.rm = TRUE)) %>%
    mutate(area = rep(400, n()), guild = rep(my_list[[m]], n())) -> data[[m]]
}

do.call(rbind, data) %>%
  select(plotid, mean_value, sd_value, area, guild) %>%
  rename_all(~ paste0(c("mean_value", "sd_value", "plotid", "area", "guild"))) -> data_dob_20m_scale

save(data_dob_20m_scale, file = "data_dob_20m_scale.rds")



# dob plots at the 30m by 30 m scale

location <- sample_data(dob_permutation) %>%
  data.frame() %>%
  select(X1, X2) %>%
  rename_all(~ paste0(c("gx", "gy")))

rownames(location) <- rownames(sample_data(dob_permutation))

location <- sample_data(location)

dob_permutation <- merge_phyloseq(dob_permutation, location)

ori <- expand.grid(x = seq(0, 10, 0.5), y = seq(0, 10, 0.5))

###

term <- c("protistan_parasite", "lichen_parasite", "algal_parasite", "mycoparasite", "animal_parasite")

guild_select <- c(
  "ectomycorrhizal", "arbuscular_mycorrhizal", "soil_saprotroph", "litter_saprotroph", "plant_pathogen", "wood_saprotroph",
  "para", "epiphyte"
)

taxa_df <- as.data.frame(tax_table(dob_permutation))

guild <- tax_table(dob_permutation) %>%
  data.frame() %>%
  select(ta2)

guild %>%
  mutate(ta2 = ifelse(ta2 %in% term, "para", ta2)) %>%
  rename(guild = ta2) -> guild

taxa_df$guild <- guild

new_tax_table <- as.matrix(taxa_df)


## at the 30m by 30 m scale for the dob site
guild_select <- c(
  "ectomycorrhizal", "arbuscular_mycorrhizal", "soil_saprotroph", "litter_saprotroph", "plant_pathogen", "wood_saprotroph",
  "para", "epiphyte"
)
tax_table(dob_permutation) %>%
  data.frame() %>%
  pull(ta2) %>%
  unique() -> all
my_list <- c(guild_select[1:8], list(all))


data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(dob_permutation, ta2 %in% my_list[[m]])


  plot_ID <- sample_data(dob_permutation) %>% data.frame()
  plot_ID <- unique(plot_ID$plotIDM)
  set.seed(520)

  richness_plot_dob <- matrix(ncol = 3, nrow = length(plot_ID))
  for (k in 1:length(plot_ID))
  {
    cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing

    tryCatch(
      {
        core_number <- numeric()
        for (i in 1:dim(ori)[1])
        {
          # cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing

          d <- sample_data(data_select) %>%
            data.frame() %>%
            select(gx, gy, plotIDM) %>%
            filter(plotIDM == plot_ID[k])

          core_number[i] <- subset(d, gx >= ori$x[i] & gx <= ori$x[i] + 30 & gy <= ori$y[i] + 30 & gy >= ori$y[i]) %>%
            summarize(row_count = n()) %>%
            as.numeric()
        }
        core_number_8 <- core_number > 8
        core_number_8[!core_number_8] <- NA


        richness <- numeric()
        for (j in 1:length(core_number)) # only select the locations that over 9 cores can be sampled/81 cores
        {
          cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
          tryCatch(
            {
              if (j %in% c(which(!is.na(core_number_8)))) {
                # informs the processing
                sub_samp <- subset(d, gx >= ori$x[j] & gx <= ori$x[j] + 30 & gy <= ori$y[j] + 30 & gy >= ori$y[j])
                sub_samp <- sample(rownames(sub_samp), 9) # select only nine cores
                kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
                richness[j] <- table(colSums(otu_table(kk)) > 0)["TRUE"] %>% as.numeric()
              } else {
                sub_samp <- subset(d, gx >= ori$x[j] & gx <= ori$x[j] + 30 & gy <= ori$y[j] + 30 & gy >= ori$y[j])
                if (dim(sub_samp)[1] < 2) {
                  richness[j] <- NA
                } else {
                  sub_samp <- rownames(sub_samp) # get the sample names
                  kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp) # based on this to extrapolate the richness
                  data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
                  ot <- t(otu_table(data_sub)) %>% data.frame()
                  n <- dim(ot)[2] # the number of "sites" for each plot
                  ot1 <- rowSums(ot)
                  out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 9, length.out = 9)), se = FALSE)
                  richness[j] <- out3$iNextEst$size_based %>%
                    slice_tail() %>%
                    pull(qD)
                }
              }
            },
            error = function(e) {
              # Handle the error and continue
              print(paste("Error on element", j, ":", e))
            },
            finally = {
              print(paste("Processed element", j))
            }
          )
        }

        richness_plot_dob[k, 1] <- plot_ID[k]
        richness_plot_dob[k, 2] <- richness %>% mean(na.rm = TRUE)
        richness_plot_dob[k, 3] <- richness %>% sd(na.rm = TRUE)
      },
      error = function(e) {
        # Handle the error and continue
        print(paste("Error on element", k, ":", e))
      },
      finally = {
        print(paste("Processed element", k))
      }
    )
  }

  data[[m]] <- cbind(richness_plot_dob, my_list[[m]]) #
}



cbind(richness_plot_dob, my_list[[m]]) %>%
  data.frame() %>%
  mutate(area = rep(900, n())) %>%
  select(X1, X2, X3, area, X4) %>%
  rename_all(~ paste0(c("plotid", "mean_value", "sd_value", "area", "guild"))) -> data_dob_30m_scale

save(data_dob_30m_scale, file = "data_dob_30m_scale.RData")


### for the dob sites at the 40 by 40 m scale

data <- list()
for (m in 1:9)
{
  cat("\r", paste(paste0(rep("*", round(m / 1, 0)), collapse = ""), m, collapse = "")) # informs the processing

  data_select <- subset_taxa(dob_permutation, ta2 %in% my_list[[m]])

  set.seed(5679)
  times <- 30
  a4 <- sample_data(data_select)
  a4 <- unique(a4$plotIDM)
  richness <- vector("list", length(a4))

  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, plotIDM == a4[i])
    dim1 <- dim(sample_data(data_sub))[1]

    tryCatch(
      {
        if (dim1 >= 16) {
          cl <- makeCluster(3)
          registerDoParallel(cl)
          richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "dplyr")) %dopar% {
            species <- vector(length = dim1[1]) # create a vector to save diversity
            for (j in 1:16)
            {
              # randomly sample j samples in the plot
              flag <- rep(FALSE, dim1[1])
              flag[sample(1:dim1[1], 16)] <- TRUE
              temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
              species[j] <- sum(otu_table(temp)["TRUE"] > 0)
              return(species[j])
            }
          }
          stopCluster(cl)
        } else if (dim1 >= 2 & dim1 <= 15) {
          sub_samp <- rownames(sample_data(data_sub))
          kk <- subset_samples(data_select, sample_names(data_select) %in% sub_samp)
          data_sub <- transform_sample_counts(kk, function(x) ifelse(x > 0, 1, 0))
          ot <- t(otu_table(data_sub)) %>% data.frame()
          n <- dim(ot)[2] # the number of "sites" for each plot
          ot1 <- rowSums(ot) # when the sum is o, it incurs errors
          out3 <- iNEXT(c(n, ot1), q = 0, datatype = "incidence_freq", size = round(seq(1, 16, length.out = 16)), se = FALSE)
          richness[[i]] <- out3$iNextEst$size_based %>%
            slice_tail() %>%
            pull(qD)
        } else {
          richness[[i]] <- NA
        }
      },
      error = function(e) {
        # Handle the error and continue
        print(paste("Error on element", i, ":", e))
      },
      finally = {
        print(paste("Processed element", i))
      }
    )
  }


  richness_mean_subplot40_dob <- matrix(ncol = 2, nrow = length(a4))
  for (i in 1:length(a4))
  {
    richness_mean_subplot40_dob[i, 1] <- mean(richness[[i]])
    richness_mean_subplot40_dob[i, 2] <- sd(richness[[i]])
  }

  richness_mean_subplot40_dob %>%
    data.frame() %>%
    mutate(plotid = a4, area = rep(1600, n()), guild = rep(my_list[[m]], n())) %>%
    rename_all(~ paste0(c("mean_value", "sd_value", "plotid", "area", "guild"))) %>%
    select(plotid, mean_value, sd_value, area, guild) -> data[[m]]
}

do.call(rbind, data) -> data_dob_40m_scale

save(data_dob_40m_scale, file = "data_dob_40m_scale.RData")


# combing all the data

bind_rows(
  data_neon_10m_scale,
  data_neon_20m_scale,
  data_neon_30m_scale,
  data_neon_40m_scale,
  data_dob_10m_scale,
  data_dob_20m_scale,
  data_dob_30m_scale,
  data_dob_40m_scale,
) -> full_dob_neon_richness_data
