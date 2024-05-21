setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps <- PSB.CSS


#############################Mice###########################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Mice")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


vir.phyl <- tax_glom(ps.rarefied, "Order", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
#Bacteriome
PSVamp$metadata$Groups <- factor(PSVamp$metadata$Groups, levels = groups)
PSVamp$metadata$Dilution_rate <- factor(PSVamp$metadata$Dilution_rate, levels = dil_groups)

p_abundance_heatmap_mice_vir <-amp_heatmap(PSVamp,
                                      group_by = "Groups",
                                      facet_by = "Dilution_rate",
                                      plot_values = FALSE,
                                      tax_show = 12,
                                      tax_aggregate = "Order" ,
                                      tax_add = "Class",
                                      tax_empty = "best",
                                      showRemainingTaxa = TRUE,
                                      normalise = TRUE,
                                      color_vector = c("white", "blue"),
                                      plot_colorscale = "sqrt",
                                      plot_legendbreaks = c(1, 10, 25, 50)
) +
  ggtitle(paste("Mouse Virome")) +
  scale_x_discrete(breaks=groups,
                   labels=groups)+
  theme_classic() +
  mytheme_abundance
p_abundance_heatmap_mice_vir


#############################Human###########################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Human")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))



vir.phyl <- tax_glom(ps.rarefied, "Order", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
#Bacteriome
PSVamp$metadata$Groups <- factor(PSVamp$metadata$Groups, levels = groups)
PSVamp$metadata$Dilution_rate <- factor(PSVamp$metadata$Dilution_rate, levels = dil_groups)

p_abundance_heatmap_human_vir <-amp_heatmap(PSVamp,
                                      group_by = "Groups",
                                      facet_by = "Dilution_rate",
                                      plot_values = FALSE,
                                      tax_show = 12,
                                      tax_aggregate = "Order" ,
                                      tax_add = "Class",
                                      tax_empty = "best",
                                      showRemainingTaxa = TRUE,
                                      normalise = TRUE,
                                      color_vector = c("white", "blue"),
                                      plot_colorscale = "sqrt",
                                      plot_legendbreaks = c(1, 10, 25, 50)
) +
  ggtitle(paste("Human Virome")) +
  scale_x_discrete(breaks=groups,
                   labels=groups)+
  theme_classic() +
  mytheme_abundance
p_abundance_heatmap_human_vir


################################################

ggarrange(p_abundance_heatmap_mice_vir,
          p_abundance_heatmap_human_vir,
          p_abundance_heatmap_mice_bac,    #  ran p_abundance_heatmao_group.R first
          p_abundance_heatmap_human_bac,
          nrow = 2,
          ncol = 2,
          labels = c("A","B","C","D"),
          common.legend = FALSE)
