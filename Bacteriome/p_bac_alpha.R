setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")


Alpha_obs <- c("Shannon")

set.seed(111) # keep result reproductive

ps = rarefy_even_depth(PSB, rngseed=1, sample.size=2000, replace=F)

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Mice")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

#Normalize to mean read count
#Preparing data sheet for further differential analysis
rich <- estimate_richness(ps.rarefied)
tab <- subset(rich, select = Alpha_obs)
index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
tab$Groups <-sample_data(ps.rarefied)$Groups[index]
tab$Cultivation_phase <- sample_data(ps.rarefied)$Cultivation_phase[index]

tab$Sample <- rownames(rich)

###############################################################################




stat_cul <- tab %>%
  wilcox_test(Shannon~Cultivation_phase,
              comparisons = cul_compare,
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_cul


stat <- tab %>%
  wilcox_test(Shannon~Groups,
              comparisons = group_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat


sheets <- list("wilcox_2side_cul_compare"=stat_cul,
               "wilcox_fdr_2side_group_compare"=stat)




write_xlsx(sheets, "Bacteriome/stat_result/Bacteriome_alpha_Mice.xlsx")


###############################################################################

# Get the sample size 
n1 <-  count(tab$Groups == "Inoculum")[[2,2]]
n2 <-  count(tab$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(tab$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(tab$Groups == "Chemostat_Slow")[[2,2]]
n5 <-  count(tab$Groups == "Chemostat_Fast")[[2,2]]



#no significant exist in kruskal_test

tab$Groups <- factor(tab$Groups,levels = groups)
p_mice <- ggplot(tab, aes(x= Groups, y= Shannon)) +
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
                   )+
  scale_y_continuous(limits = c(2.4, 4.1))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Shannon diversity index", title="α-diversity Mouse") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size = 6, y.position = c(NA,NA,NA))+
  theme_classic() +
  mytheme_alpha

p_mice


######################################### Human ##################################################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Human")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

#Normalize to mean read count
#Preparing data sheet for further differential analysis
rich <- estimate_richness(ps.rarefied)
tab <- subset(rich, select = Alpha_obs)
index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
tab$Groups <-sample_data(ps.rarefied)$Groups[index]
tab$Cultivation_phase <- sample_data(ps.rarefied)$Cultivation_phase[index]

tab$Sample <- rownames(rich)

###############################################################################




stat_cul <- tab %>%
  wilcox_test(Shannon~Cultivation_phase,
              comparisons = cul_compare,
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_cul


stat <- tab %>%
  wilcox_test(Shannon~Groups,
              comparisons = group_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat


sheets <- list("wilcox_2side_cul_compare"=stat_cul,
               "wilcox_fdr_2side_group_compare"=stat)




write_xlsx(sheets, "Bacteriome/stat_result/Bacteriome_alpha_human.xlsx")


###############################################################################

# Get the sample size 
n1 <-  count(tab$Groups == "Inoculum")[[2,2]]
n2 <-  count(tab$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(tab$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(tab$Groups == "Chemostat_Slow")[[2,2]]
n5 <-  count(tab$Groups == "Chemostat_Fast")[[2,2]]



#no significant exist in kruskal_test

tab$Groups <- factor(tab$Groups,levels = groups)
p_human <- ggplot(tab, aes(x= Groups, y= Shannon)) +
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
  )+
  scale_y_continuous(limits = c(1.7, 4))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Shannon diversity index", title="α-diversity Human") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size = 6, y.position = c(NA,NA,NA))+
  theme_classic() +
  mytheme_alpha

p_human



############################################
#width 900 * height 700

ggarrange(p_beta_mice,  #beta diversity ran p_vir_beta_bray.R first
          p_beta_human,  #beta diversity
          p_mice,
          p_human,
          ncol = 2,
          nrow = 2,
          labels = c("A","B","C","D"),
          common.legend = FALSE
)
