
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3_CSS_phyloseq_vir.R")



ps <- PSV.no.realm

############################################Mice#######################################
ps1 <- subset_samples(ps, Sample_origin %in% "Mice")
ps1 <- subset_samples(ps1, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

ps2<- ps1

taxonomy <- read.csv('Data/Virome/eukaryotic_virome_taxonomy.csv' ,row.names = 1)
tax_table(ps2) <- as.matrix(taxonomy)


sum_counts1  <- colSums(otu_table(ps1))
sum_counts2  <- colSums(otu_table(ps2))
percentage <- as.data.frame((sum_counts2 / sum_counts1) * 100)

mapping <- as.data.frame(sample_data(ps1))
percentage$Eukaryotic_Percentage<- percentage$`(sum_counts2/sum_counts1) * 100`
percentage$Groups<- mapping$Groups
percentage$Cultivation_phase<- mapping$Cultivation_phase


#stat
stat_cul <- percentage %>%
  wilcox_test(Eukaryotic_Percentage~Cultivation_phase,
              comparisons = cul_compare,
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_cul


stat <- percentage %>%
  wilcox_test(Eukaryotic_Percentage~Groups,
              comparisons = group_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat


sheets <- list("wilcox_2side_cul_compare"=stat_cul,
               "wilcox_fdr_2side_group_compare"=stat)






write_xlsx(sheets, "Virome/stat_result/virome_alpha_eukaryotic_mice.xlsx")


###############################################################################

# Get the sample size 
n1 <-  count(percentage$Groups == "Inoculum")[[2,2]]
n2 <-  count(percentage$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(percentage$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(percentage$Groups == "Chemostat_Slow")[[2,2]]
n5 <-  count(percentage$Groups == "Chemostat_Fast")[[2,2]]





percentage$Groups <- factor(percentage$Groups,levels = groups)
p_mice <- ggplot(percentage, aes(x= Groups, y= Eukaryotic_Percentage)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Groups), coef=1.5) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=groups,
                   labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
  ) +
  scale_y_continuous(limits = c(0, 0.3))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Relative abundance (%)", title="Mice") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha
p_mice




############################################Human#######################################
ps1 <- subset_samples(ps, Sample_origin %in% "Human")
ps1 <- subset_samples(ps1, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

ps2<- ps1

taxonomy <- read.csv('Data/Virome/eukaryotic_virome_taxonomy.csv' ,row.names = 1)
tax_table(ps2) <- as.matrix(taxonomy)


sum_counts1  <- colSums(otu_table(ps1))
sum_counts2  <- colSums(otu_table(ps2))
percentage <- as.data.frame((sum_counts2 / sum_counts1) * 100)

mapping <- as.data.frame(sample_data(ps1))
percentage$Eukaryotic_Percentage<- percentage$`(sum_counts2/sum_counts1) * 100`
percentage$Groups<- mapping$Groups
percentage$Cultivation_phase<- mapping$Cultivation_phase


#stat
stat_cul <- percentage %>%
  wilcox_test(Eukaryotic_Percentage~Cultivation_phase,
              comparisons = cul_compare,
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_cul


stat <- percentage %>%
  wilcox_test(Eukaryotic_Percentage~Groups,
              comparisons = group_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat


sheets <- list("wilcox_2side_cul_compare"=stat_cul,
               "wilcox_fdr_2side_group_compare"=stat)






write_xlsx(sheets, "Virome/stat_result/virome_alpha_eukaryotic_human.xlsx")


###############################################################################

# Get the sample size 
n1 <-  count(percentage$Groups == "Inoculum")[[2,2]]
n2 <-  count(percentage$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(percentage$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(percentage$Groups == "Chemostat_Slow")[[2,2]]
n5 <-  count(percentage$Groups == "Chemostat_Fast")[[2,2]]





percentage$Groups <- factor(percentage$Groups,levels = groups)
p_human <- ggplot(percentage, aes(x= Groups, y= Eukaryotic_Percentage)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Groups), coef=1.5) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=groups,
                   labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
  ) +
  scale_y_continuous(limits = c(0, 0.3))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Relative abundance (%)", title="Human") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha
p_human





########### merge plot#############################

p_euka_vir <-ggarrange(p_mice,
                        p_human,
                        nrow=1, ncol = 2,
                        font.label = list(size = 20),
                        common.legend = TRUE,
                        legend = "right")


p_euka_vir
