
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3_CSS_phyloseq_vir.R")



ps <- subset_samples(PSV.no.realm, Sample_origin %in% c("Mice"))

ps <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


ps <- tax_glom(ps, "Family", NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Groups)

Groups
Groups <- Groups[Groups!="Inoculum"]
path_table <- "Virome/stat_result/desep2/Mice/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Groups %in% c(group, "Inoculum"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Groups)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Inoculum.tsv"), sep="\t", col.names = NA)
}


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Groups)

Groups
Groups <- Groups[Groups!="Chemostat_Slow"]
path_table <- "Virome/stat_result/desep2/Mice/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Groups %in% c(group, "Chemostat_Slow"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Groups)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Chemostat_Slow.tsv"), sep="\t", col.names = NA)
}


### loop for all the grouping compared with Saline
Cultivation_phase <- unique(sample_data(ps)$Cultivation_phase)

Cultivation_phase
Cultivation_phase <- Cultivation_phase[Cultivation_phase!="Inoculum"]
path_table <- "Virome/stat_result/desep2/Mice/"
dir.create(path_table)

for (group in Cultivation_phase){
  ps.sub <- prune_samples(sample_data(ps)$Cultivation_phase %in% c(group, "Inoculum"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Cultivation_phase)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Inoculum.tsv"), sep="\t", col.names = NA)
}





#################################################################################

ps <- subset_samples(PSV.no.realm, Sample_origin %in% c("Human"))

ps <- subset_samples(ps, Groups %in% c("Inoculum","Batch_Slow","Batch_Fast","Chemostat_Slow","Chemostat_Fast"))


ps <- tax_glom(ps, "Family", NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Groups)

Groups
Groups <- Groups[Groups!="Inoculum"]
path_table <- "Virome/stat_result/desep2/Human/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Groups %in% c(group, "Inoculum"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Groups)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Inoculum.tsv"), sep="\t", col.names = NA)
}






### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Groups)

Groups
Groups <- Groups[Groups!="Chemostat_Slow"]
path_table <- "Virome/stat_result/desep2/Human/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Groups %in% c(group, "Chemostat_Slow"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Groups)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Chemostat_Slow.tsv"), sep="\t", col.names = NA)
}





### loop for all the grouping compared with Saline
Cultivation_phase <- unique(sample_data(ps)$Cultivation_phase)

Cultivation_phase
Cultivation_phase <- Cultivation_phase[Cultivation_phase!="Inoculum"]
path_table <- "Virome/stat_result/desep2/Human/"
dir.create(path_table)

for (group in Cultivation_phase){
  ps.sub <- prune_samples(sample_data(ps)$Cultivation_phase %in% c(group, "Inoculum"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Cultivation_phase)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_Inoculum.tsv"), sep="\t", col.names = NA)
}


