setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3_CSS_phyloseq_vir.R")



Alpha_obs <- c("Observed")

set.seed(111) # keep result reproductive

ps = rarefy_even_depth(PSV.no.realm, rngseed=1, sample.size=20000, replace=F)

taxonomy <- read.csv('Data/Virome/eukaryotic_virome_taxonomy.csv' ,row.names = 1)
tax_table(ps) <- as.matrix(taxonomy)

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Mice")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


###############################################################################

# Get the sample size 
n1 <-  count(ps.rarefied@sam_data$Groups == "Inoculum")[[2,2]]
n2 <-  count(ps.rarefied@sam_data$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(ps.rarefied@sam_data$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(ps.rarefied@sam_data$Groups== "Chemostat_Slow")[[2,2]]
n5 <-  count(ps.rarefied@sam_data$Groups == "Chemostat_Fast")[[2,2]]



#no significant exist in kruskal_test

ps.rarefied@sam_data$Groups <- factor(ps.rarefied@sam_data$Groups,levels = groups)
p_mice <- plot_richness(ps.rarefied, x="Groups", measures=c("Observed")) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Groups), coef=1.5) +
  # geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=groups,
                   labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
  )+
  scale_y_continuous(limits = c(0, 6))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Observed eukaryotic viral\n OTUs", title=" Mouse") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size = 6, y.position = c(NA,NA,NA))+
  theme_classic() +
  mytheme_alpha

p_mice


######################################### Human ##################################################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Human")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

###############################################################################

# Get the sample size 
n1 <-  count(ps.rarefied@sam_data$Groups == "Inoculum")[[2,2]]
n2 <-  count(ps.rarefied@sam_data$Groups == "Batch_Slow")[[2,2]]
n3 <-  count(ps.rarefied@sam_data$Groups == "Batch_Fast")[[2,2]]
n4 <-  count(ps.rarefied@sam_data$Groups== "Chemostat_Slow")[[2,2]]
n5 <-  count(ps.rarefied@sam_data$Groups == "Chemostat_Fast")[[2,2]]



#no significant exist in kruskal_test

ps.rarefied@sam_data$Groups <- factor(ps.rarefied@sam_data$Groups,levels = groups)
p_human <- plot_richness(ps.rarefied, x="Groups", measures=c("Observed")) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Groups), coef=1.5) +
  # geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=groups,
                   labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
  )+
  scale_y_continuous(limits = c(0,20))+
  scale_fill_manual(values = cols) +
  labs(x="", y="Observed eukaryotic viral\n OTUs", title=" Human") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size = 6, y.position = c(NA,NA,NA))+
  theme_classic() +
  mytheme_alpha

p_human




############################################
#width 900 * height 700

ggarrange(p_mice,
          p_human,
          ncol = 2,
          nrow = 1,
          labels = c("A","B"),
          common.legend = FALSE
)
