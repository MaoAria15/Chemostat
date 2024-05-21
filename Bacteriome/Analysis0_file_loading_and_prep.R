

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("2_CSS_phyloseq_bac.R")


############Set data categories - and order



PSB.CSS = prune_samples(sample_sums(PSB) >= 2000, PSB.CSS)
PSB = prune_samples(sample_sums(PSB) >= 2000, PSB)





#
# saveRDS(PSB, file = "data/corr_data/ps_bacterium.rds")

