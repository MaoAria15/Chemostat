
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3_CSS_phyloseq_vir.R")






ps <- subset_samples(PSV.no.realm, Sample_origin %in% "Mice")
############################################Slow#######################################
ps1 <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Chemostat_Slow","Chemostat_Slow_Intermediate"))



taxonomy <- read.csv('Data/Virome/eukaryotic_virome_taxonomy.csv' ,row.names = 1)
votus_to_remove <- as.vector(rownames(taxonomy))
ps1 <- subset_taxa(ps1, !(taxa_names(ps1) %in% votus_to_remove))
vir.phyl <- tax_glom(ps1, "Order", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
#Bacteriome
PSVamp$metadata$Cultivation_phase <- factor(PSVamp$metadata$Cultivation_phase, 
                                            levels = c("Inoculum",
                                                       "Batch", 
                                                       "Chemostat_Vol1.8",
                                                       "Chemostat_Vol3",
                                                       "Chemostat_Vol4",
                                                       "Chemostat_Vol5"))
PSVamp$metadata$Dilution_rate <- factor(PSVamp$metadata$Dilution_rate, levels = c("Inoculum","Slow"))

p_abundance_heatmap_vir_slow <-amp_heatmap(PSVamp,
                                           group_by = "Cultivation_phase",
                                           facet_by = "Dilution_rate",
                                           plot_values = FALSE,
                                           tax_show = 12,
                                           tax_aggregate = "Order" ,
                                           tax_add = "Class",
                                           tax_empty = "best",
                                           showRemainingTaxa = TRUE,
                                           normalise = TRUE,
                                           color_vector = c("white", "red"),
                                           plot_colorscale = "sqrt",
                                           plot_legendbreaks = c(1, 10, 25, 50)
) +
  scale_x_discrete(breaks=c("Inoculum",
                            "Batch", 
                            "Chemostat_Vol1.8",
                            "Chemostat_Vol3",
                            "Chemostat_Vol4",
                            "Chemostat_Vol5"),
                   labels=c("Inoculum",
                            "Batch", 
                            "Chemostat_Vol1.8",
                            "Chemostat_Vol3",
                            "Chemostat_Vol4",
                            "Chemostat_Vol5"))+
  theme_classic() +
  mytheme_abundance_text_vertical
p_abundance_heatmap_vir_slow


############################################Fast#######################################
ps1 <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Fast","Chemostat_Fast","Chemostat_Fast_Intermediate"))



taxonomy <- read.csv('Data/Virome/eukaryotic_virome_taxonomy.csv' ,row.names = 1)
votus_to_remove <- as.vector(rownames(taxonomy))
ps1 <- subset_taxa(ps1, !(taxa_names(ps1) %in% votus_to_remove))
vir.phyl <- tax_glom(ps1, "Order", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
#Bacteriome
PSVamp$metadata$Cultivation_phase <- factor(PSVamp$metadata$Cultivation_phase, 
                                            levels = c("Inoculum",
                                                       "Batch", 
                                                       "Chemostat_Vol1",
                                                       "Chemostat_Vol2",
                                                       "Chemostat_Vol5"))
PSVamp$metadata$Dilution_rate <- factor(PSVamp$metadata$Dilution_rate, levels = c("Inoculum","Fast"))

p_abundance_heatmap_vir_fast <-amp_heatmap(PSVamp,
                                           group_by = "Cultivation_phase",
                                           facet_by = "Dilution_rate",
                                           plot_values = FALSE,
                                           tax_show = 12,
                                           tax_aggregate = "Order" ,
                                           tax_add = "Class",
                                           tax_empty = "best",
                                           showRemainingTaxa = TRUE,
                                           normalise = TRUE,
                                           color_vector = c("white", "red"),
                                           plot_colorscale = "sqrt",
                                           plot_legendbreaks = c(1, 10, 25, 50)
) +
  scale_x_discrete(breaks=c("Inoculum",
                            "Batch", 
                            "Chemostat_Vol1",
                            "Chemostat_Vol2",
                            "Chemostat_Vol5"),
                   labels=c("Inoculum",
                            "Batch", 
                            "Chemostat_Vol1",
                            "Chemostat_Vol2",
                            "Chemostat_Vol5"))+
  theme_classic() +
  mytheme_abundance_text_vertical
p_abundance_heatmap_vir_fast

########################################################

ggarrange(p_abundance_heatmap_vir_slow,
          p_abundance_heatmap_vir_fast,
          nrow = 1,
          ncol = 2,
          labels = c("C","D")
)
