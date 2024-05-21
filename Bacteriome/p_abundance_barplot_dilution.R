setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps<- PSB
ps<- subset_samples(ps, Sample_origin %in% "Mice")
##########################Mice#########################################################
ps1<- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Chemostat_Slow","Chemostat_Slow_Intermediate"))
ps1<- subset_samples(ps1, Cultivation_phase %in% c("Inoculum","Batch","Chemostat_Vol1.8","Chemostat_Vol3","Chemostat_Vol4","Chemostat_Vol5"))
# Get the sample size 
n1 <-  count(ps1@sam_data$Cultivation_phase == "Inoculum")[[2,2]]
n2 <-  count(ps1@sam_data$Cultivation_phase  == "Batch")[[2,2]]
n3 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol1.8")[[2,2]]
n4 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol3")[[2,2]]
n5 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol4")[[2,2]]
n6 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol5")[[2,2]]


vir.phyl <- tax_glom(ps1, "Genus", NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Cultivation_phase")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

top<-df %>%
  # filter(level=="Phylum")%>%
  group_by(tax)%>%
  dplyr::summarize(mean_abund=mean(Abundance), .groups = "drop")%>%
  arrange(desc(mean_abund))

#Find top 30 tax
top

top20 <- top$tax[1:20] #for species


df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = c("Inoculum",
                                            "Batch", 
                                            "Chemostat_Vol1.8",
                                            "Chemostat_Vol3",
                                            "Chemostat_Vol4",
                                            "Chemostat_Vol5"))


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_slow <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name="Genus")+
  scale_x_discrete(labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\n(N={n2})"),
                            glue("Chemostat\nVol1.8\n(N={n3})"),
                            glue("Chemostat\nVol3\n(N={n4})"),
                            glue("Chemostat\nVol4\n(N={n5})"),
                            glue("Chemostat\nVol5\n(N={n6})"))
  )+
  coord_cartesian(ylim = c(0, 1)) +
  mytheme_alpha+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("Groups")+
  ggtitle("Mouse Bacteriome Slow")

p_slow

#####################################FAST#####################################


ps1<- subset_samples(ps, Groups %in% c("Inoculum","Batch_Fast","Chemostat_Fast","Chemostat_Fast_Intermediate"))
ps1<- subset_samples(ps1, Cultivation_phase %in% c("Inoculum","Batch","Chemostat_Vol1","Chemostat_Vol2","Chemostat_Vol5"))
# Get the sample size 
n1 <-  count(ps1@sam_data$Cultivation_phase == "Inoculum")[[2,2]]
n2 <-  count(ps1@sam_data$Cultivation_phase  == "Batch")[[2,2]]
n3 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol1")[[2,2]]
n4 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol2")[[2,2]]
n5 <-  count(ps1@sam_data$Cultivation_phase  == "Chemostat_Vol5")[[2,2]]


vir.phyl <- tax_glom(ps1, "Genus", NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Cultivation_phase")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps3)



# Check the first few rows of the merged dataframe
head(df)



#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

top<-df %>%
  # filter(level=="Phylum")%>%
  group_by(tax)%>%
  dplyr::summarize(mean_abund=mean(Abundance), .groups = "drop")%>%
  arrange(desc(mean_abund))

#Find top 30 tax
top

top20 <- top$tax[1:20] #for species


df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = c("Inoculum",
                                            "Batch", 
                                            "Chemostat_Vol1",
                                            "Chemostat_Vol2",
                                            "Chemostat_Vol5"))


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_fast <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name="Genus")+
  scale_x_discrete(labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\n(N={n2})"),
                            glue("Chemostat\nVol1\n(N={n3})"),
                            glue("Chemostat\nVol2\n(N={n4})"),
                            glue("Chemostat\nVol5\n(N={n5})"))
  )+
  coord_cartesian(ylim = c(0, 1)) +
  mytheme_alpha+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("Groups")+
  ggtitle("Mouse Bacteriome Fast")

p_fast



#############################################################################
p_barplot <-ggarrange(p_slow,
                      p_fast,
                      labels = c("A","B"),
                      nrow=1, ncol = 2)


p_barplot

