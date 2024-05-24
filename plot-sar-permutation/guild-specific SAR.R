# determine the SAR based on different fungal trophic guilds
# load the fungal guilds data into the object
ft <- read.csv("FungalTraits_1.2_ver_16Dec_2020.csv", na.strings="")
# the taxa table includes the trophic guild information of the OTUS

tk <- subset_taxa(rare_all, ta2 == "ectomycorrhizal")
# to see the number of different guilds for the data set
tax_table(rare_all)%>%data.frame()%>%dplyr::select(ta2)%>%table()%>%as.data.frame()%>%arrange(Freq)
#
#           mycoparasite   403
#         root_endophyte   431
#        animal_parasite   497
#  arbuscular_mycorrhizal  786
#               epiphyte   930
#          plant_pathogen  1675
#  unspecified_saprotroph  1838
#        wood_saprotroph   2220
#       litter_saprotroph  3792
#    animal_endosymbiont   6062
#  nectar/tap_saprotroph   6287
#         soil_saprotroph  8808
#        ectomycorrhizal   9143

# to select eight fungal guilds

data_EM <- subset_taxa(rare_all_assign, ta2 == "ectomycorrhizal")
data_AM <- subset_taxa(rare_all_assign, ta2 == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(rare_all_assign, ta2 == "soil_saprotroph")
data_littersap <- subset_taxa(rare_all_assign, ta2 == "litter_saprotroph")
data_plapat <- subset_taxa(rare_all_assign, ta2 == "plant_pathogen")
data_woodsap <- subset_taxa(rare_all_assign, ta2 == "wood_saprotroph")
data_para <- subset_taxa(rare_all_assign, ta2%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(rare_all_assign, ta2 == "epiphyte")

# We construct the SAR for each guild


