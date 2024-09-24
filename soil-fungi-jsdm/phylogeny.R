library(tibble)
library(Hmsc)
library(dplyr)
library(readr)
library(raster)
library(ape)
library(phyloseq)
library(pROC)
library(fungaltraits)
library(phytools)
library(stringr)


# the joint species distribution model for the neon data
neon_dob_prevalent_v4.1 <- readRDS("~/soil-sar/data/neon_dob_prevalent_v4.1.Rds")

# to select the guilds of interests

fungal_trait=read_csv("~/soil-sar/plot-sar-permutation/FungalTraits_1.2_ver_16Dec_2020.csv")

fungal_trait%>%select(GENUS,species,primary_lifestyle)

tax_table(neon_dob_prevalent_v4.1)%>%data.frame()%>%dplyr::select(genus,species)%>%rename(GENUS=genus)%>%
  left_join(fungal_trait%>%dplyr::select(GENUS,primary_lifestyle),by="GENUS")%>%as.matrix() ->full_species

# the taxa table must be a matrix
rownames(full_species)=taxa_names(neon_dob_prevalent_v4.1)
new_tax_table <- tax_table(full_species)
tax_table(neon_dob_prevalent_v4.1)=new_tax_table

# the richness for different guilds

new_tax_table%>%data.frame()%>%count(n=primary_lifestyle)
# only 96 AM fungi
# select the species with trait available for the subsequent analysis


#with a thredshold occurence of 10, we have 2893 species
# to see the traits data available for the species

fungal_trait_spore=fungal_traits()

tax_table(data_with_full_name)%>%data.frame()->full_name_list

# to see the coverage of the data
names=colnames(fungal_trait_spore)

trait_data=list()
for (i in 6:101)
{
  fungal_trait_spore%>%left_join(full_name_list%>%dplyr::select(species,primary_lifestyle),by="species")%>%dplyr::filter(!is.na(!!sym(names[i])))%>%
    group_by(species)%>%summarize(mean_value=mean(!!sym(names[i]),rm.na=TRUE))->trait_data[[i]]
}


# see the dimension of each data list

data.cover=numeric()
for (i in 6:101)
{
  data.cover[i]=dim(trait_data[[i]])[1]
}

cbind(data.cover,names)

for (i in 1:101)
{
  print(trait_data[i])
}

# select the spore size and tissure cn as the two traits

#trait_data[[82]]%>%left_join(trait_data[[93]],by="species")%>%filter(!is.na( mean_value.x)&!is.na(mean_value.y))

## if we use spore size and tissue cn content, only 60 species would be selected

## if we use spore size and fruit body size, 690 species would selected

trait_data[[82]]%>%left_join(trait_data[[50]],by="species")%>%filter(!is.na( mean_value.x)&!is.na(mean_value.y))%>%
  filter(!grepl("\\?",species))->body_spore_data

# subset the phyloseq object based on the slected species with two traits
# need to add a column showing the trophic guild of the species

body_spore_data%>%left_join(full_species%>%data.frame(),by="species")%>%filter(!is.na(primary_lifestyle))%>%
  count(species)->species_with_traits

# to see the shared species
# only 88 species would be retained with both traits and 
full_species%>%data.frame()%>%pull()

intersect(body_spore_data$species,full_species%>%data.frame()%>%pull(species))->species_included_model


# guilds             species
#ectomycorrhizal     106
#litter_saprotroph    27
#soil_saprotroph       4
#wood_saprotroph       7


# let's construct the model based on the 88 species

species_subset_phyloseq <- subset_taxa(neon_dob_prevalent_v4.1, species%in%c(species_included_model))

# note that some species are duplicated 

data_otu_table=otu_table(species_subset_phyloseq)%>%data.frame()
data_taxa_table=tax_table(species_subset_phyloseq)%>%data.frame()

colnames(data_otu_table)=data_taxa_table$species

col_names <- colnames(data_otu_table)

grouped_columns <- split(seq_along(col_names), col_names)

df_summed <- sapply(grouped_columns, function(cols) {
  if (length(cols) > 1) {
    # If there are multiple columns, sum them across rows
    rowSums(data_otu_table[, cols], na.rm = TRUE)
  } else {
    # If there's only one column, just return that column
    data_otu_table[, cols]
  }
})

df_summed <- as.data.frame(df_summed)# the modified otu table

# need to modified the taxa_table
# the otu number should be recorded

colnames(df_summed)%>%data.frame()%>%rename_all(~paste0("species"))%>%
  left_join(data_taxa_table%>%distinct(),by="species")%>%
  left_join(data_taxa_table%>%distinct()%>%
              rownames_to_column(var = "OTU")%>%dplyr::select(OTU,species))->taxa_table_with_otu

rownames(taxa_table_with_otu)=taxa_table_with_otu$OTU

# remove the otu column
taxa_table_with_otu=taxa_table_with_otu%>%dplyr::select(-OTU)
# to rename the colnames of the otu table

colnames(df_summed)=rownames(taxa_table_with_otu)

otu_table(df_summed%>%as.matrix(),taxa_are_rows = FALSE)->new_otu_table
# important to set the taxa_are_rows = FALSE
# convert the table as a tax

tax_table(taxa_table_with_otu%>%as.matrix())->taxa_table_with_otu

taxa_table_with_otu=tax_table(taxa_table_with_otu%>%as.matrix())

physeq@tax_table=taxa_table_with_otu
physeq@otu_table=new_otu_table

physeq_88_species=physeq

saveRDS(physeq_88_species,file="physeq_88_species.rds")

#################################

# to get the mean traits of the species

taxa_table_with_otu%>%data.frame()%>%left_join(body_spore_data,by="species")%>%
  rename(spore_size=mean_value.x,body_size= mean_value.y)->species_traits_88

# construct the model
y_data=otu_table(physeq_88_species)%>%data.frame()

x_data=sample_data(physeq_88_species)%>%data.frame()%>%dplyr::select(mat_celsius, map_mm,temp_seasonality)%>%
  scale()%>%data.frame()%>%scale()%>%as.data.frame()

trait_data=species_traits_88

rownames(trait_data)=colnames(y_data)

m=Hmsc(Y=y_data,XData = x_data,TrData=trait_data,
       TrFormula=~primary_lifestyle +spore_size +body_size,
       XFormula=~mat_celsius+ map_mm,distr = "probit",XScale = TRUE)

nChains=2
m = sampleMcmc(m, thin = 10, samples = 1000, transient = 3000,
               nChains = nChains, nParallel = nChains, verbose = TRUE)

mpost = convertToCodaObject(m)
effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate=FALSE)$psrf

par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
###

# what if the R2 is low for most species?
preds = computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds)

#We next evaluate the model’s predictive power through two-fold cross validation.
partition = createPartition(m, nfolds = 2)
preds = computePredictedValues(m, partition = partition, nParallel = nChains)

evaluateModelFit(hM = m, predY = preds)

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95,cex = c(0.3, 0.21, 0.2))

## species trait response to the variables

postgamma = getPostEstimate(m, parName = "Gamma")

plotGamma(m, post = postgamma, param = "Support", supportLevel = 0.85,cex = c(0.7, 0.7, 0.7))








# get the unique table for the initial table
# for the duplicated rows, only the first row would be retained.
data_taxa_table%>%distinct()%>%rownames_to_column(var = "OTU")


# to see other combinations
trait_species_number=numeric()
for (i in 6:101)
{
  trait_data[[82]]%>%left_join(trait_data[[i]],by="species")%>%filter(!is.na( mean_value.x)&!is.na(mean_value.y))->trait_combinations
  trait_species_number[i]=dim(trait_combinations)[1] 
}








# select data with specific guilds

data_ecm_soilsap=subset_taxa(neon_dob_prevalent_v4.1, primary_lifestyle%in%c("ectomycorrhizal","soil_saprotroph"))

# there are 2200 species of the guilds and then we select species with more than 30 sites

otu_table(data_ecm_soilsap)%>%data.frame()%>%colSums()%>%data.frame()%>%rownames_to_column(var = "species")%>%
  rename_all(~paste0(c("species","frequency")))%>%filter(frequency>29)->select_species

#30 occurence will lead to 227 species
#20 occurence will lead to 561 species
#30 occurence will lead to 227 species

data_ecm_soilsap_sub=prune_taxa(select_species$species, data_ecm_soilsap)

# to see how many species are involved
#71 ectomycorrhizal fungai 
#156 soil_saprotroph fungi

data_ecm_soilsap_sub%>%tax_table()%>%data.frame()%>%count(primary_lifestyle)

# start to build the joint species distribution model
# the model without random effect
# adaptNf 

y_data=otu_table(data_ecm_soilsap_sub)%>%data.frame()

x_data=sample_data(data_ecm_soilsap_sub)%>%data.frame()%>%dplyr::select(mat_celsius, map_mm,temp_seasonality)%>%
  scale()%>%data.frame()%>%scale()%>%as.data.frame()

m=Hmsc(Y=y_data,XData = x_data,XFormula=~mat_celsius+ map_mm,distr = "probit",XScale = TRUE)

nChains=2
m = sampleMcmc(m, thin = 10, samples = 100, transient = 3000,
               nChains = nChains, nParallel = nChains, verbose = TRUE)

mpost = convertToCodaObject(m)
effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate=FALSE)$psrf

par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")

# what if the R2 is low for most species?
preds = computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds)

#We next evaluate the model’s predictive power through two-fold cross validation.
partition = createPartition(m, nfolds = 2)
preds = computePredictedValues(m, partition = partition, nParallel = nChains)

evaluateModelFit(hM = m, predY = preds)

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95,cex = c(0.3, 0.21, 0.2))


#If your model used 1000 posterior samples, for instance, it might return predictions for each sample, leading to a large list of predictions.

predictions <- predict(m, XData = extracted_values[1:5,],
                       type = "response",expected=TRUE,
                       posteriorMean = FALSE)# to get the mean of the posterior distribution.

#determine the cutpoint for each species

observed_climate=observed_climate=sample_data(neon_dob_prevalent_v4.1)%>%
  data.frame()%>%dplyr::select(mat_celsius, map_mm)%>%
  scale()%>%data.frame()

# the presence probability for each species

predictions <- predict(m, XData = observed_climate,
                       type = "response",expected=TRUE,
                       posteriorMean = FALSE)# to get the mean of the posterior distribution.

# get the mean probability for the species across the 128 sites
stacked_array <- simplify2array(predictions )
mean_matrix <- apply(stacked_array, c(1, 2), mean)

# to get the cutpoint for each of the 227 species

optimal_cutpoint=numeric()
for (i in 1:dim(mean_matrix)[2]){
  #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  if(i!=63){
    roc_curve <- roc(y_data[,i], mean_matrix[1:dim(mean_matrix)[1],i])
    optimal_cutpoint[i] <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")%>%as.numeric()
    
  }
  else{
    roc_curve <- roc(y_data[,i], mean_matrix[1:dim(mean_matrix)[1],i])
    
    two_threshold=coords(roc_curve, "best", ret = "threshold", best.method = "youden")
    optimal_cutpoint[i] =mean(two_threshold$threshold)
    
  }
}

# make predictions and then based on the cutpoints to dedermine if species would present
# based on 200 plot for example


## use the model to make predictions based on env raster
load("~/soil-sar/SDM/coords_present.RData")
load("~/soil-sar/plot-sar-permutation/r_present_northam.RData")

points <- SpatialPoints(coords_present, proj4string = CRS(projection(r_present_northam)))

extracted_values <- extract(r_present_northam, points)

extracted_values=extracted_values%>%data.frame()%>%dplyr::select(mat_celsius, map_mm)

extracted_values=extracted_values[complete.cases(extracted_values),]


predictions <- predict(m, XData = extracted_values,
                       type = "response",expected=TRUE,
                       posteriorMean = FALSE)# to get the mean of the posterior distribution.

stacked_array <- simplify2array(predictions )
mean_matrix <- apply(stacked_array, c(1, 2), mean)

predicted_presence=matrix(nrow=1000,ncol=227)
for (i in 1:227)
{
  predicted_presence[,i] <- as.integer(mean_matrix [,i] > optimal_cutpoint[i])
  
}







function_cutpoint=function(data)
{
  stacked_array <- simplify2array(data)
  mean_matrix <- apply(stacked_array, c(1, 2), mean)
  roc_curve <- roc(y_data[,1], mean_matrix[1:dim(mean_matrix)[1],1])
  optimal_cutpoint <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")
  return(optimal_cutpoint)
}



## when we set the random

studyDesign = data.frame(sample = as.factor(1:128))

rL = HmscRandomLevel(units = studyDesign$sample)

m_random = Hmsc(Y = y_data, XData = x_data, XFormula = ~mat_celsius+ map_mm,
                studyDesign = studyDesign, ranLevels = list(sample = rL))

m_random = sampleMcmc(m_random, thin = 10, samples = 1000, transient = 3000,
                      nChains = 2, nParallel = 3, verbose = 0)


mpost = convertToCodaObject(m_random)
par(mfrow=c(2,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
hist(effectiveSize(mpost$Omega[[1]]), main="ess(omega)")
hist(gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)")

# the last line code above will cause abortion

# let's take a look at the species association



postBeta = getPostEstimate(m_random, parName="Beta")
plotBeta(m_random, post=postBeta, param="Support", supportLevel = 0.95)


OmegaCor = computeAssociations(m_random)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

corrplot(toPlot, method = "color",type="lower",
         col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("random effect level:", m_random$rLNames[1]), mar=c(0,0,1,0))

# include the phylogeny and traits in the model
# just selct the species with full names to compromise the trait data

neon_dob_prevalent_v4.1%>%tax_table()%>%data.frame()%>%
  filter(!grepl("_sp", species))%>%pull(species)->full_names

gsub("-", "_", full_names)->full_names

neon_dob_prevalent_v4.1%>%tax_table()%>%data.frame()%>%
  filter(!grepl("_sp", species)&!is.na(primary_lifestyle))%>%count(primary_lifestyle)

#2572 species retained

# to look at the occurence of different species

subset_taxa(neon_dob_prevalent_v4.1, species%in%full_names)->data_with_full_name

# there are 2200 species of the guilds and then we select species with more than 30 sites
data_with_full_name%>%tax_table()%>%data.frame()%>%count(primary_lifestyle)
#

otu_table(data_with_full_name)%>%data.frame()%>%colSums()%>%data.frame()%>%rownames_to_column(var = "species")%>%
  rename_all(~paste0(c("species","frequency")))%>%filter(frequency>9)->select_species


# to check the relationship between AMF and ectormycorrhizal fungi associations

data_ecm_amf=subset_taxa(neon_dob_prevalent_v4.1, primary_lifestyle%in%c("ectomycorrhizal","arbuscular_mycorrhizal"))

# there are 2200 species of the guilds and then we select species with more than 30 sites

otu_table(data_ecm_amf)%>%data.frame()%>%colSums()%>%data.frame()%>%rownames_to_column(var = "species")%>%
  rename_all(~paste0(c("species","frequency")))%>%filter(frequency>9)->select_species

#

#30 occurence will lead to 227 species
#20 occurence will lead to 561 species
#30 occurence will lead to 227 species

data_ecm_amf_sub=prune_taxa(select_species$species, data_ecm_amf)

# to see how many species are involved
#71 ectomycorrhizal fungai 
#156 soil_saprotroph fungi

data_ecm_amf_sub%>%tax_table()%>%data.frame()%>%count(primary_lifestyle)

# start to build the joint species distribution model
# the model without random effect
# adaptNf 

y_data=otu_table(data_ecm_amf_sub)%>%data.frame()

x_data=sample_data(data_ecm_amf_sub)%>%data.frame()%>%select(mat_celsius, map_mm,temp_seasonality)%>%
  scale()%>%data.frame()


# build the model with random effect


studyDesign = data.frame(sample = as.factor(1:128))

rL = HmscRandomLevel(units = studyDesign$sample)

m_random = Hmsc(Y = y_data, XData = x_data, XFormula = ~mat_celsius+ map_mm,
                studyDesign = studyDesign, ranLevels = list(sample = rL))

m_random = sampleMcmc(m_random, thin = 10, samples = 1000, transient = 3000,
                      nChains = 2, nParallel = 3, verbose = 0)


mpost = convertToCodaObject(m_random)
par(mfrow=c(2,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
hist(effectiveSize(mpost$Omega[[1]]), main="ess(omega)")
hist(gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)")

# to construct different guild combinations
# to see the number of different species for the sub set data



OmegaCor = computeAssociations(m_random)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         title=paste("random effect level:", m_random$rLNames[1]), mar=c(0,0,1,0))

### for the high-dimensional hmsc Venetian

#------------------------construct a phylogeny for the fungal species data------------------------->



setwd("/Users/luowenqi/soil-sar/plot-sar-permutation")
tree <- read.tree("labelled_supertree_ottnames.tre")
tree_names=tree$tip.label%>%data.frame()%>%rename_all(~paste0("species"))


# select the genus that belongs to fungal species and exclude other groups

patterns=full_species%>%data.frame()%>%dplyr::select(GENUS)%>%pull()%>%unique()%>%sort()
#some genus have a long name, which contains part of a fungal names
patterns=paste0(patterns,"_")

# filter out the species that do not contain the genus names
# this does not ensure that all the species would be included in the tree?

filtered_df<- tree_names[!grepl(paste(patterns, collapse = "|"), tree_names$species), ]
tree_new=drop.tip(tree,filtered_df)


## rename the tip names of the phylogeny and get the binary nomenclature
# to replace the "'" symbol
tree_new$tip.label=gsub("'","",tree_new$tip.label)
name.revised= matrix(nrow=length(tree_new$tip.label),ncol=3)
for(i in 1:length(tree_new$tip.label))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  name.revised[i,]=strsplit(tree_new$tip.label[i], "_")[[1]][1:3]
}


name.revised=name.revised%>%data.frame()%>%
  mutate(new_name=case_when(X2%in% c("cf.", "aff.")~paste(X1,X3,sep="_"),
                            X1%in% c("cf.", "aff.")~paste(X2,X3,sep="_"),
                            TRUE~paste(X1,X2,sep="_")))

tree_new$tip.label=name.revised$new_name

# removed the dumplicated tips

unique_tips <- !duplicated(tree_new$tip.label)
tree_pruned <- drop.tip(tree_new, tree_new$tip.label[!unique_tips])

tip_labels <- tree_pruned$tip.label

# need to find the species that contained in the fungal data
# we have 2512 species in total

# to see how many genus are overlapped between the tree and the fungal data

genus_for_fungal=full_species%>%data.frame()%>%distinct(GENUS)%>%pull()

genus_for_tree=tree_pruned$tip.label%>%data.frame()%>%rename_all(~paste0("species"))%>%
  mutate(genus = word(species, 1, sep = "_"))%>%distinct(genus)%>%pull()

## see the overlap between the two data set
#676 genus shared between the two data sets

shared_genus=intersect(genus_for_fungal,genus_for_tree)

# to see what species shared with the given genus from the fungal data 

congeneric_species=list()
for(i in 1:676)
{
  congeneric_species[[i]] <- tip_labels[grepl(paste0("^", shared_genus[i], "_"), tip_labels)]
  
}
# find the common node number of the the selected species, and the con generic species will be added to the node

node_number=numeric()
for (i in 1:676)
{
  temp=getMRCA(tree_pruned, congeneric_species[[i]])
  if(is.null(temp)){
    node_number[i]=0 # can not find more species to share a node
  }
  else{
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    node_number[i] <- temp
  }
}

# based on the genes, we can see the node number to be add

shared_genus%>%data.frame()%>%bind_cols(node_number)%>%data.frame()%>%
  rename_all(~paste0(c("genus","node")))->node_number_to_add

# bind with the full speces list to determine to which node should the specis be added
full_species%>%data.frame()%>%rename_all(~paste0(c("genus","species","primary_lifestyle")))%>%
  left_join(node_number_to_add,by="genus")->full_species_node

full_species_node_complete=full_species_node%>%filter(!is.na(node)&node!=0)

# add the species to the tree

new_species <- full_species_node_complete$species
node_numbers <- full_species_node_complete$node

try_tree=tree_pruned

for (i in 1:length(new_species)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  try_tree<- bind.tip(try_tree, tip.label = new_species[i], where = node_numbers[i])
}

# for the single species, we can replace the con generic species

full_species_node%>%filter(node==0)%>%pull(genus)->single_genus

tree_pruned$tip.label%>%data.frame()%>%rename_all(~paste0("species"))%>%
  mutate(genus = word(species, 1, sep = "_"))%>%filter(genus%in%single_genus)->tree_single_genus


full_species_node%>%filter(node==0)%>%
  distinct()%>%left_join(tree_single_genus,by="genus")->tree_fungal_singl_genus

try_tree$tip.label <- gsub(tree_fungal_singl_genus$species.y, tree_fungal_singl_genus$species.x, try_tree$tip.label)

#just keep the tips of interest 

final_tree <- keep.tip(try_tree, c(new_species,tree_fungal_singl_genus$species.y))


ggtree(final_tree,ladderize = FALSE,layout = "fan") +
  geom_tiplab(color="blue",size=2)

geom_nodelab(color="purple",size=3)+
  geom_nodelab(aes(label = node), size = 2, hjust = -0.3,color="red")




# this will lead to 1630 species in the tree

full_species_node%>%filter(!is.na(node))%>%distinct(species)



## to construct a tree with all the species available in the fungal data set

subtree <- extract.clade(tree_pruned, "Hygrocybe_ott282216")


ggtree(tree,ladderize = FALSE) +
  geom_tiplab(color="blue",size=4)+
  geom_nodelab(color="purple",size=3)+
  geom_nodelab(aes(label = node), size = 5, hjust = -0.3,color="red")


ggtree(updated_tree ,ladderize = FALSE,layout = "ape") +
  geom_tiplab(color="blue",size=1)+
  geom_nodelab(aes(label = node), size = 3, hjust = -0.3,color="red")


#geom_nodepoint(aes(subset = node == 2312), size = 3, hjust = -0.3,color="purple")
##
tip_labels <- tree_pruned$tip.label
# find the con generic species on the tree for each of the 88 species

genus_to_search=species_traits_88%>%pull(GENUS)

# for each genus to find the con generic species
congeneric_species=list()
for(i in 1:88)
{
  congeneric_species[[i]] <- tip_labels[grepl(paste0("^", genus_to_search[i], "_"), tip_labels)]
  
}
# find the common node number and then add the tip to the tree

mrca_node <- getMRCA(tree_pruned, congeneric_species[[11]])



pruned_tree <- drop.tip(tree, species_to_drop)

keep_tree <- keep.tip(tree_pruned, intersect(tree_pruned$tip.label,species_traits_88$species))

interaction()
# removed the duplicated names



# to see the overlapped in species between phylogeny and with the trait data

tree_species=data.frame(tree_pruned $tip.label)

tree_species%>%rename(species=tree_pruned.tip.label)%>%left_join(body_spore_data,by="species")%>%
  filter(!is.na(mean_value.x)&!is.na(mean_value.y))

# only 29 species shared between the two data sets 
# if all the species genus were retained


body_spore_data%>%mutate(left_side = word(species, 1, sep = "_"))%>%
  distinct(left_side)%>%pull()->unique_genus


tree_species%>%mutate(left_side = word(tree_pruned.tip.label, 1, sep = "_"))%>%
  distinct(left_side)%>%pull()->unique_genus_tree

shared_elements <- intersect(unique_genus, unique_genus_tree)%>%data.frame()%>%mutate(number=1:25)

names(shared_elements)=c("left_side","number")

# to see how many 
body_spore_data%>%mutate(left_side = word(species, 1, sep = "_"))%>%count(left_side)%>%left_join(shared_elements,by="left_side")%>%
  print(n=64)%>%filter(!is.na(n)&!is.na(number))%>%summarize(total = sum(n, na.rm = TRUE))

body_spore_data%>%mutate(left_side = word(species, 1, sep = "_"))%>%count(left_side)%>%left_join(shared_elements,by="left_side")%>%
  print(n=64)->dk

dk%>%filter(!is.na(number))%>%pull(left_side )->genus_of_interest


# we will have 546 species 
# retain the tips with the genus included

tips_to_keep <- unlist(lapply(genus_of_interest, function(genus) {
  grep(genus, tree_pruned$tip.label, value = TRUE)
}))

pruned_tree_with_trait <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, tips_to_keep))

plot(pruned_tree_with_trait ,type="fan")

## to construct a model without phylogeny but with species traits

# trait data #body_spore_data
# to add a column showing the mycorrhizal type of the species

body_spore_dataleft_join(taxa_table%>%data.frame()%>%dplyr::select(species,GENUS,primary_lifestyle),by="species")












#3 to see the mismatch part of the names and the patterns



name.revised= matrix(nrow=length(tree_new$tip.label),ncol=2)

for(i in 1:length(tree_new$tip.label))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  name.revised[i,]=strsplit(tree_new$tip.label[i], "_")[[1]][1:2]
}




tip.names=name.revised$X1


only_in_tipnames <- setdiff(tip.names, patterns)

tree_tip=tree_new$tip.label
# just get the first par

strsplit(tree_tip, "_")

setdiff(tree_tip, patterns)

##to look at the availability of the fungal trait data






## some exercise for the jsdm analysis

library(Hmsc)
library(corrplot)
library(ape)
library(MASS)
library(fields)
library(knitr)

dat <- read_csv("dat.csv")
phy <- ape::read.tree("tree.nex")

Y<-dat[,1:50]#划分物种栖息地情况
XData<-dat[,51:52]

head(XData)

climate <- XData[,1]

habitat <- as.factor(XData[,2])#does not work

habitat$hab=as.factor(habitat$hab)

XData <- data.frame(climate = climate, habitat = habitat)#读

ns = 50#定义采样点数量
n = 200

XFormula = ~habitat + poly(climate,degree = 2,raw = TRUE)

TrFormula = ~habitat.use + thermal.optimum

studyDesign = data.frame(sample = sprintf('sample_%.3d',1:n), stringsAsFactors=TRUE)
#创建一个随机水平对象
rL = HmscRandomLevel(units = studyDesign$sample)
#设置最大允许的自由度
rL$nfMax = 15
#Hmsc()函数进行分层混合模型分析

m = Hmsc(Y = Y, #物种丰度数据         
         XData = XData, #栖息地状况数据         
         XFormula = XFormula,#固定效应         
         TrData = traits, #物种性状数据        
         TrFormula = TrFormula,#随机效应        
         phyloTree = phy,#系统发育树         
         studyDesign = studyDesign, #随机效应水平         
         ranLevels = list(sample = rL))

rownames(traits)=traits$spnames

nChains = 2#设置test.run状态，当为True时，运行速度快，但是参数不可靠test.run = FALSE
#定义参数
test.run = FALSE
if (test.run)
{  
  thin = 1  
  samples = 100 
  transient = 50
}else{
  thin = 10 
  samples = 1000 
  transient = 500
  
}
verbose = 0
#采用sampleMcmc函数解释模型

m = sampleMcmc(m, thin = 10, samples = 1000, 
               transient = 500,              
               nChains = 2,
               nParallel = nChains, verbose = 0)

#评估MCMC收敛性


mpost = convertToCodaObject(m)
#物种生态位

ess.beta = effectiveSize(mpost$Beta)
Betaess.beta = effectiveSize(mpost$Beta)

psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf

hist(ess.beta)
hist(psrf.beta)

#性状对物种生态位的影响Gamma

ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)

sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp))
{  
  tmp[[chain]] = tmp[[chain]][,sppairs]
}

ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

##
print("ess.rho:")
effectiveSize(mpost$Rho)
print("psrf.rho:")
gelman.diag(mpost$Rho)$psrf
#物种特异性解释r2值直方图
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

#50个物种的方差划

head(m$X)

VP = computeVariancePartitioning(m, group = c(1,1,2,2),  groupnames = c("habitat","climate"))

plotVariancePartitioning(m, VP = VP)

#来查看这些特征的解释能力
kable(VP$R2T$Beta)
VP$R2T$Y

postBeta = getPostEstimate(m, parName = "Beta")

plotBeta(m, post = postBeta, param = "Support",         plotTree = TRUE, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))

postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.95)

##
OmegaCor = computeAssociations(m)
supportLevel = 0.65
toPlot = ((OmegaCor[[1]]$support>supportLevel)  + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

corrplot(toPlot, method = "color",col=colorRampPalette(c("blue","white","red"))(200),  tl.cex=.6, tl.col="black",  title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))#

summary(mpost$Rho)

Gradient = constructGradient(m,focalVariable = "climate",non.focalVariables = list("habitat"=list(3,"open")))

Gradient$XDataNew#只在开放栖息地中采样，气候梯度与物种之间的关系

predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,  
                ranLevels=Gradient$rLNew, expected=TRUE)

plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE)#通过设置measure=“Y”，可视化单个物种的相同预测，看index值就代表物种排序
plotGradient(m, Gradient, pred=predY, measure="Y", index = 1, showData = TRUE)
plotGradient(m, Gradient, pred=predY, measure="T", index = 3, showData = TRUE)#这里可以查看气候条件一样，栖息地变化对物种的影响，注意事项跟上方一样

Gradient = constructGradient(m,focalVariable = "habitat",   
                             non.focalVariables = list("climate"=list(1)))

#注意数据是自己构造的，而不是原有的，需要自己替换实际数据
Gradient$XDataNewpredY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew, ranLevels=Gradient$rLNew, expected=TRUE)

#看左侧是物种37,这也说明物种37是最具代表的物种

plotGradient(m, Gradient, pred=predY, measure="Y", 
             index=which.max(m$TrData$habitat.use), 
             showData = TRUE, jigger = 0.2)

#生境对群落加权平均生境利用的影响

plotGradient(m, Gradient, pred=predY, measure="T", index=2, showData = TRUE, jigger = 0.2)



