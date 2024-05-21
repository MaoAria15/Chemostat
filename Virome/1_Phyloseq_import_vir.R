

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")
######################################Construct Phyloseq objects##########################################



##Load raw data 

#OTU table
otu_mat <- read.delim('Data/Virome/tpm_otu_table.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

#04/11-2021 taxonomy from Josue
tax_mat <- read.csv('Data/Virome/Viral_taxonomy.csv' ,row.names = 1)

#04/11-2021 host predictions form Josue
host_mat <- read.csv('Data/Virome/host_taxonomy.csv', row.names = 1)

#Convert to dataframe
tax_mat <- as.matrix.data.frame(tax_mat)

map_mat <- read.delim('Data/Virome/mapping_Chemostat_virome.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

# map_mat <- dplyr::mutate_at(map_mat,1:10,funs(factor)) 
# map_mat <- dplyr::mutate_if(map_mat, is.character, as.numeric)

##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(tax_mat)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSV <- phyloseq(OTU, TAX, samples)

#Make metadata files for easy overview

meta.vir <- as.data.frame(as.matrix(sample_data(PSV)))

#############Make host prediction file


host_mat <- as.matrix.data.frame(host_mat)

HOST <- phyloseq::tax_table(host_mat)

##Combine in phyloseq object

PSV.HOST <- phyloseq(OTU, HOST, samples)

#############
PSV <- subset_samples(PSV, remove_23112022_human %in% "s")
PSV.HOST <- subset_samples(PSV.HOST, remove_23112022_human %in% "s")

#############Remove left over objects

rm(tax_mat, map_mat, otu_mat, samples,OTU, TAX,TAX.no.realm,meta.vir,host,host.clean,host_mat,tax,tax_mat,tax.clean,HOST,order,kingdom,class,family,i,j,phylum,Realm,rownames)



