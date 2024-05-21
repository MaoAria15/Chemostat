library(ggplot2)
library(MESS)
library(forcats)
library(tidyverse)
library(psycho)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ampvis2)
library(DESeq2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(microbiome)
library(forestmangr)
library(writexl)
library(viridis)
library(readxl)
library(phyloseq)      # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(writexl)
library(directlabels)
library(glue)
library(ggpmisc)
library(ComplexHeatmap)
library(magick)
library(colorRamp2)

library(circlize)
library(microDecon)
library(psych)
library(pheatmap)
library(igraph)
library(ggraph)
library(influential)
library(showtext)
library(gridExtra)
# library(MicrobiotaProcess)
cols <- c("#0073C2FF" ,"#EFC000FF" ,"#868686FF" ,"#CD534CFF" ,"#7AA6DCFF")


group_col <- c(Inoculum="#0073C2FF", Batch_Slow="#EFC000FF", Batch_Fast="#868686FF",Chemostat_Slow="#CD534CFF",Chemostat_Fast="#7AA6DCFF")

group_compare <-list(c("Inoculum","Chemostat_Slow"),
                    c("Inoculum","Chemostat_Fast"),
                    c("Chemostat_Slow","Chemostat_Fast")
                    )

cul_compare <- list(c("Inoculum","Batch"))

groups <- c("Inoculum", "Batch_Slow", "Batch_Fast","Chemostat_Slow", "Chemostat_Fast")


dil_groups <- c("Inoculum","Slow","Fast")






mytheme_with_x <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 axis.ticks.x = element_blank(),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_blank(),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme_alpha <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 axis.ticks.x = element_blank(),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 # axis.text.x=element_blank(),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)

mytheme_alpha_noy <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                       axis.line=element_line(size=0.5),
                       #panel.border = element_blank(),
                       axis.text=element_text(size = 8, colour = "Black"),
                       axis.ticks=element_line(size=1, colour = "Black"),
                       axis.ticks.x = element_blank(),
                       axis.text.y=element_blank(), 
                       axis.ticks.y=element_blank(), 
                       axis.title.y=element_blank(),
                       axis.line.y=element_blank(),
                       strip.background = element_rect(colour = "white", fill = "white"),
                       # axis.text.x=element_blank(),
                       axis.title = element_text(size = 8, face = "bold"),
                       strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(8, "pt"),
                       legend.title = element_text(size = 8,face = "bold"),
                       title = element_text(size =8, face = "bold")
)
mytheme_beta <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                       axis.line=element_line(size=0.5),
                       #panel.border = element_blank(),
                       axis.text=element_text(size = 8, colour = "Black"),
                       axis.ticks=element_line(size=1, colour = "Black"),
                       strip.background = element_rect(colour = "white", fill = "white"),
                       axis.title = element_text(size = 8, face = "bold"),
                       strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(8, "pt"),
                       legend.title = element_text(size = 8,face = "bold"),
                       title = element_text(size =8, face = "bold")
)
mytheme_abundance <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.text.x=element_text(size = 8,colour = "Black"),
                      axis.text.y=element_text(size = 8,colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.title.x  = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                      legend.text = element_text(size=8,face = "italic"),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold")
)
mytheme_abundance_text_vertical <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                           axis.line=element_line(size=0.5),
                           #panel.border = element_blank(),
                           axis.text=element_text(size = 8, colour = "Black"),
                           axis.text.x=element_text(size = 8,angle = -45,hjust=0.1,vjust=0.2,colour = "Black"),
                           axis.text.y=element_text(size = 8,colour = "Black",face = "italic"),
                           axis.ticks=element_line(size=1, colour = "Black"),
                           strip.background = element_rect(colour = "white", fill = "white"),
                           axis.title.x  = element_text(size = 8, face = "bold"),
                           strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                           legend.text = element_text(size=8,face = "italic"),
                           legend.key.size = unit(8, "pt"),
                           legend.title = element_text(size = 8,face = "bold"),
                           title = element_text(size =8, face = "bold")
)
mytheme_abundance_noy <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                           axis.line=element_line(size=0.5),
                           axis.text=element_text(size = 8, colour = "Black"),
                           axis.text.x=element_text(size = 8,angle = -45,hjust=0.1,vjust=0.2,colour = "Black"),
                           axis.text.y=element_blank(), 
                           axis.ticks.y=element_blank(), 
                           axis.title.y=element_blank(),
                           axis.line.y=element_blank(),
                           axis.ticks=element_line(size=1, colour = "Black"),
                           strip.background = element_rect(colour = "white", fill = "white"),
                           axis.title = element_text(size = 8, face = "bold"),
                           strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                           legend.text = element_text(size=8,face = "italic"),
                           legend.key.size = unit(8, "pt"),
                           legend.title = element_text(size = 8,face = "bold"),
                           title = element_text(size =8, face = "bold"))
mytheme_his <-  theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                      axis.title = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                      legend.text = element_text(size=8),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold"),
                      legend.background = element_rect(color = "black",linewidth = 1))


filter_and_replace <- function(data, threshold = 0.05) {
  data %>%
    filter(p < threshold) %>%
    mutate(p_signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )) %>%
    select(-p)  # Remove the original p-value column
}


