setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3_CSS_phyloseq_vir.R")

#############################################Mice###################################################################################
ps <- subset_samples(PSV.no.realm, Sample_origin %in% c("Mice"))

ps <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))
ps <- tax_glom(ps, "Family", NArm = FALSE) #select a level to compare
#################################################################################
#################################################################################
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Virome/stat_result/desep2/Mice/Batch_Inoculum.tsv", sep = "\t",header = T, row.names = 1)
tab_2_all <- read.table("Virome/stat_result/desep2/Mice/Chemostat_Slow_Inoculum.tsv", sep = "\t", header = T, row.names = 1)
tab_3_all <- read.table("Virome/stat_result/desep2/Mice/Chemostat_Fast_Inoculum.tsv", sep = "\t",header = T, row.names = 1)
tab_4_all <- read.table("Virome/stat_result/desep2/Mice/Chemostat_Fast_Chemostat_Slow.tsv", sep = "\t", header = T, row.names = 1)


tab_1 <-subset(tab_1_all, padj < 0.05)
# tab_1 <-subset(tab_1, abs(log2FoldChange) > 0.6)
tab_2 <- subset(tab_2_all, padj < 0.05)
# tab_2 <- subset(tab_2, abs(log2FoldChange) >0.6)
tab_3 <- subset(tab_3_all, padj < 0.05)

tab_4 <- subset(tab_4_all, padj < 0.05)


# colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
theme_here <-theme(legend.text = element_text(face = "italic"),
                   panel.background = element_rect(fill = "white", color = "black"),
                   axis.line = element_line(color = "black"),
                   panel.grid.major = element_line(color = "lightgrey"),
                   panel.grid.minor = element_line(color = "lightgrey"),
                   axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5,colour = "black" ,size = 10, face = "bold"),
                   axis.text.y = element_text(size = 9, face = "italic", colour = "black"), title = element_text(face = "bold"),
                   plot.title=element_text(size = 10,face = "bold")) 
###############################################################################

OTU1 <- rownames(tab_1)

ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU1,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]  
a <- intersect(OTU1,species2keep)
tab_1 <- tab_1[a,]
colors <- c("#8dd3c7", "#e41a1c", "#bebada", "#fb8072")
p_mice_1 <-   ggplot(tab_1, aes(y = Family, x = log2FoldChange, color = Phylum)) +
                        geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
                        geom_point(size = 4) +
                        theme_here+
                        ggtitle("Mouse\nLog2FoldChange of Family:\nBatch vs Inoculum ") +
                        geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
                        scale_x_continuous(limits = c(-9, 5), n.breaks = 10) +
                        scale_y_discrete(expand = c(0.00005, 0.8)) +
                        scale_color_manual(values = colors) +
                        ylab("Family")



p_mice_1
###############################################################################
OTU2 <- rownames(tab_2)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU2,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
tab_2 <- tab_2[species2keep,]

colors <- c("#8dd3c7", "#e41a1c", "#fb8072")
p_mice_2 <-   ggplot(tab_2, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Mouse\nLog2FoldChange of Family:\nChemostat-Slow vs Inoculum ") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-17, 14), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_mice_2               

###############################################################################
OTU3 <- rownames(tab_3)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU3,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
a <- intersect(OTU3,species2keep)
tab_3 <- tab_3[a,]

colors <- c("#8dd3c7", "#e41a1c", "#bebada", "#fb8072")
p_mice_3 <-   ggplot(tab_3, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Mouse\nLog2FoldChange of Family:\nChemostat-Fast vs Inoculum") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-13, 16), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_mice_3    

###############################################################################
OTU4 <- rownames(tab_4)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU4,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
a <- intersect(OTU4,species2keep)
tab_4 <- tab_4[a,]

colors <- c("#8dd3c7", "#bebada", "#fb8072")
p_mice_4 <-   ggplot(tab_4, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Mouse\nLog2FoldChange of Family:\nChemostat-Fast vs Chemostat-Slow") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-4, 13), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_mice_4    



#############################################Human###################################################################################
ps <- subset_samples(PSV.no.realm, Sample_origin %in% c("Human"))

ps <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))
ps <- tax_glom(ps, "Family", NArm = FALSE) #select a level to compare
#################################################################################
#################################################################################
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Virome/stat_result/desep2/Human/Batch_Inoculum.tsv", sep = "\t",header = T, row.names = 1)
tab_2_all <- read.table("Virome/stat_result/desep2/Human/Chemostat_Slow_Inoculum.tsv", sep = "\t", header = T, row.names = 1)
tab_3_all <- read.table("Virome/stat_result/desep2/Human/Chemostat_Fast_Inoculum.tsv", sep = "\t",header = T, row.names = 1)
tab_4_all <- read.table("Virome/stat_result/desep2/Human/Chemostat_Fast_Chemostat_Slow.tsv", sep = "\t", header = T, row.names = 1)


tab_1 <-subset(tab_1_all, padj < 0.05)
# tab_1 <-subset(tab_1, abs(log2FoldChange) > 0.6)
tab_2 <- subset(tab_2_all, padj < 0.05)
# tab_2 <- subset(tab_2, abs(log2FoldChange) >0.6)
tab_3 <- subset(tab_3_all, padj < 0.05)

tab_4 <- subset(tab_4_all, padj < 0.05)


# colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
theme_here <-theme(legend.text = element_text(face = "italic"),
                   panel.background = element_rect(fill = "white", color = "black"),
                   axis.line = element_line(color = "black"),
                   panel.grid.major = element_line(color = "lightgrey"),
                   panel.grid.minor = element_line(color = "lightgrey"),
                   axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5,colour = "black" ,size = 10, face = "bold"),
                   axis.text.y = element_text(size = 9, face = "italic", colour = "black"), title = element_text(face = "bold"),
                   plot.title=element_text(size = 10,face = "bold")) 
###############################################################################

OTU1 <- rownames(tab_1)

ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU1,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]  
a <- intersect(OTU1,species2keep)
tab_1 <- tab_1[a,]
colors <- c( "#e41a1c", "#bebada")
p_human_1 <-   ggplot(tab_1, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here+
  ggtitle("Human\nLog2FoldChange of Family:\nBatch vs Inoculum ") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-5, 5), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")



p_human_1
###############################################################################
OTU2 <- rownames(tab_2)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU2,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
tab_2 <- tab_2[species2keep,]

colors <- c( "#e41a1c", "#fb8072")
p_human_2 <-   ggplot(tab_2, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Human\nLog2FoldChange of Family:\nChemostat-Slow vs Inoculum ") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_human_2               

###############################################################################
OTU3 <- rownames(tab_3)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU3,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
a <- intersect(OTU3,species2keep)
tab_3 <- tab_3[a,]

colors <- c("#8dd3c7", "#e41a1c",  "#fb8072")
p_human_3 <-   ggplot(tab_3, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Human\nLog2FoldChange of Family:\nChemostat-Fast vs Inoculum") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 7), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_human_3    

###############################################################################
OTU4 <- rownames(tab_4)
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU4,ps.rel)

# select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]    
a <- intersect(OTU4,species2keep)
tab_4 <- tab_4[a,]

colors <- c("#8dd3c7", "#e41a1c","#fb8072")
p_human_4 <-   ggplot(tab_4, aes(y = Family, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 1, lty = 2) +
  geom_point(size = 4) +
  theme_here +
  ggtitle("Human\nLog2FoldChange of Family:\nChemostat-Fast vs Chemostat-Slow") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-5, 5), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = colors) +
  ylab("Family")


p_human_4    



















ggarrange(p_mice_1,
          p_mice_2,
          p_mice_3,
          p_mice_4,
          p_human_1,
          p_human_2,
          p_human_3,
          p_human_4,
          nrow = 2,
          ncol =4,
          labels = c("A","B","C","D","E","F","G","H")
          )
