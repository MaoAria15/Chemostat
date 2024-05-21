setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")


set.seed(19950915)




ps <- PSB.CSS

Method="bray"
#############################Mice###########################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Mice")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


GP.ord <- ordinate(ps.rarefied, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps.rarefied, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps.rarefied))
adonis_sub <- adonis2(bray.PSB ~ Groups, method = Method, data = sampledf.PSB, permutations = 999)

adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps.rarefied@sam_data$Groups <- factor(ps.rarefied@sam_data$Groups, levels = groups)
ps.rarefied@sam_data$Dilution_rate <- factor(ps.rarefied@sam_data$Dilution_rate, levels = dil_groups)
p_beta_mice <- phyloseq::plot_ordination(ps.rarefied, GP.ord, axes = 1:2,color="Groups") + 
  geom_point(alpha=2, size=2,aes(color=Groups,shape = Dilution_rate)) +
  scale_color_manual(values= cols)+
  ggtitle(paste(" Mouse Bacteriome")) +
  theme_classic()+
  mytheme_beta 


p_beta_mice


#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps.rarefied))
cbn <- combn(x=unique(metadata$Groups), m = 2)
cbn <- cbn[,c(1,2,5)]

p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.rarefied, Groups %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Groups, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)

p_treatment.adj <- c(p_treatment.adj)
t_mice_group <- data.frame(Group1=c("Inoculum","Inoculum","Chemostat_Slow"),
                        Group2=c("Chemostat_Slow","Chemostat_Fast","Chemostat_Fast"),
                        R2= round(r,digits = 3),
                        F=round(f,digits = 3),
                        p=round(p,digits = 3),
                        p.adj=p_treatment.adj)

#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps.rarefied))
cbn <- combn(x=unique(metadata$Cultivation_phase), m = 2)
cbn <- cbn[,]

p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.rarefied, Cultivation_phase %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Cultivation_phase, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


t_mice_cul <- data.frame(Group1=c("Inoculum"),
                           Group2=c("Batch"),
                           R2= round(r[2],digits = 3),
                           F=round(f[2],digits = 3),
                           p=round(p[2],digits = 3),
                           p.adj="-")

#############################Human###########################################

ps.rarefied <- subset_samples(ps, Sample_origin %in% "Human")
ps.rarefied <- subset_samples(ps.rarefied, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


GP.ord <- ordinate(ps.rarefied, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps.rarefied, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps.rarefied))
adonis_sub <- adonis2(bray.PSB ~ Groups, method = Method, data = sampledf.PSB, permutations = 999)

adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps.rarefied@sam_data$Groups <- factor(ps.rarefied@sam_data$Groups, levels = groups)
ps.rarefied@sam_data$Dilution_rate <- factor(ps.rarefied@sam_data$Dilution_rate, levels = dil_groups)
p_beta_human <- phyloseq::plot_ordination(ps.rarefied, GP.ord, axes = 1:2,color="Groups") + 
  geom_point(alpha=2, size=2,aes(color=Groups,shape = Dilution_rate)) +
  scale_color_manual(values= cols)+
  ggtitle(paste("Human Bacteriome")) +
  theme_classic()+
  mytheme_beta 


p_beta_human


#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps.rarefied))
cbn <- combn(x=unique(metadata$Groups), m = 2)
cbn <- cbn[,c(1,2,5)]

p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.rarefied, Groups %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Groups, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)

p_treatment.adj <- c(p_treatment.adj)
t_human_group <- data.frame(Group1=c("Inoculum","Inoculum","Chemostat_Slow"),
                           Group2=c("Chemostat_Slow","Chemostat_Fast","Chemostat_Fast"),
                           R2= round(r,digits = 3),
                           F=round(f,digits = 3),
                           p=round(p,digits = 3),
                           p.adj=p_treatment.adj)

#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps.rarefied))
cbn <- combn(x=unique(metadata$Cultivation_phase), m = 2)
cbn <- cbn[,]

p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.rarefied, Cultivation_phase %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Cultivation_phase, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


t_human_cul <- data.frame(Group1=c("Inoculum"),
                         Group2=c("Batch"),
                         R2= round(r[2],digits = 3),
                         F=round(f[2],digits = 3),
                         p=round(p[2],digits = 3),
                         p.adj="-")





################################################################################
sheets <- list(
  "adonis_mice_group" =  t_mice_group ,
  "adonis_mice_cul" =  t_mice_cul ,
  "adonis_human_group" =  t_human_group ,
  "adonis_human_cul" =  t_human_cul 
)
# 
write_xlsx(sheets, "Bacteriome/stat_result/bray_permanova.xlsx")
