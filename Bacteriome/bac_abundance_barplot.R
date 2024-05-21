setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps<- PSB
##########################################################################

ps1<- subset_samples(ps, Sample_origin %in% "Mice")
ps1<- subset_samples(ps1, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

# Get the sample size 
n1 <-  count(ps1@sam_data$Groups == "Inoculum")[[2,2]]
n2 <-  count(ps1@sam_data$Groups  == "Batch_Slow")[[2,2]]
n3 <-  count(ps1@sam_data$Groups  == "Batch_Fast")[[2,2]]
n4 <-  count(ps1@sam_data$Groups  == "Chemostat_Slow")[[2,2]]
n5 <-  count(ps1@sam_data$Groups  == "Chemostat_Fast")[[2,2]]



vir.phyl <- tax_glom(ps1, "Genus", NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Groups")

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
df0$Sample <- factor(df0$Sample, levels = groups)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_mice <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name="Genus")+
  scale_x_discrete(labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
                            )+
  coord_cartesian(ylim = c(0, 1)) +
  mytheme_alpha+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("Mouse")

p_mice


###################################Human#######################################

ps1<- subset_samples(ps, Sample_origin %in% "Human")
ps1<- subset_samples(ps1, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))

# Get the sample size 
n1 <-  count(ps1@sam_data$Groups == "Inoculum")[[2,2]]
n2 <-  count(ps1@sam_data$Groups  == "Batch_Slow")[[2,2]]
n3 <-  count(ps1@sam_data$Groups  == "Batch_Fast")[[2,2]]
n4 <-  count(ps1@sam_data$Groups  == "Chemostat_Slow")[[2,2]]
n5 <-  count(ps1@sam_data$Groups  == "Chemostat_Fast")[[2,2]]



vir.phyl <- tax_glom(ps1, "Genus", NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Groups")

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
df0$Sample <- factor(df0$Sample, levels = groups)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_human <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name="Genus")+
  scale_x_discrete(labels=c(glue("Inoculum\n(N={n1})"),
                            glue("Batch\nSlow\n(N={n2})"),
                            glue("Batch\nFast\n(N={n3})"),
                            glue("Chemostat\nSlow\n(N={n4})"),
                            glue("Chemostat\nFast\n(N={n5})"))
  )+
  coord_cartesian(ylim = c(0, 1)) +
  mytheme_alpha+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("Human")

p_human






#############################################################################
p_barplot <-ggarrange(p_mice,
                      p_human,
                      nrow=1, ncol = 2)

p_barplot

