# Load libraries ----
  BiocManager::install('DOSE')
  BiocManager::install('clusterProfiler')
  
  if (!require("DOSE")) install.packages("DOSE")
  library(DOSE)
  if(!require("clusterProfiler")) install.packages("clusterProfiler")
  library(clusterProfiler)
  if (!require("enrichplot")) install.packages("enrichplot")
  library(enrichplot)
  if (!require("forcats")) install.packages('forcats')
  library('forcats')
  if (!require("ggplot2")) install.packages("ggplot2")
  library('ggplot2')
  if (!require("ggstance")) install.packages("ggstance")
  library('ggstance')
  if (!require("plyr")) install.packages("plyr")
  library('plyr')
  if (!require("dplyr")) install.packages("dplyr")
  library('dplyr')
  if (!require("viridis")) install.packages("viridis")
  library("viridis")
    
  # Split plot for PTM enriched pathways ----
      NES<-read.csv("PTM_pathway_enrichment/Input_data/split_plot_upload.csv", header = TRUE, stringsAsFactors = FALSE)
      split_NES <- arrange(NES, abs(nes)) %>% 
        group_by(sign(nes)) %>% 
        slice(1:30)
      
      split <- ggplot(split_NES, aes(nes, fct_reorder(Description, nes), fill=pvalue), showCategory=10) + 
        geom_barh(stat='identity') + 
        scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
        theme_minimal() + ylab(NULL) + ggtitle("PTM-signature enrichment", subtitle = "PI3Ki resistant tumors") + xlab("NES")
      
      ggsave(split, units = "in", height = 10, width = 6, dpi = 300, compression="lzw", filename = "PTM_pathway_enrichment/Output_plots/PTM_enrich_splitplot.tiff")
      
      
# lolliplot ---
      colnames(net10)
      net10<-read.csv("Kinase_prediction/Output_data/networkin_input.csv", header = TRUE, stringsAsFactors = FALSE)
      lolli <- ggplot(net10, showCategory = 50, 
             aes(richFactor, fct_reorder(Description, richFactor))) + 
        geom_segment(aes(xend=0, yend = Description)) +
        geom_point(aes(color=log2FC, size = Substrate_count)) +
        scale_color_viridis_c(guide=guide_colorbar(reverse=FALSE), direction = -1) +
        scale_size_continuous(range=c(2, 10)) +
        theme_minimal() + 
        xlab("Kinase Prediction Score") +
        ylab(NULL) + 
        labs(caption = "NetworKIN 3.0; Nat Methods 11, 603–604 (2014)") +
        theme(axis.text = element_text(size = 12, colour = 'black')) +
        theme(axis.title = element_text(size = 14)) +
        theme(title = element_text(size = 18)) +
        theme(plot.subtitle = element_text(size = 10)) +
        ggtitle("Top predicted kinases", subtitle = "Phosphopeptides uniquely upregulated in PI3Ki resistant tumors (n=1786)") +
        theme(plot.caption = element_text(size = 10))
      ggsave(lolli, units = "in", height = 10, width = 6, dpi = 300, compression="lzw", filename = "Kinase_prediction/Output_data/kinase_lilliplot_UP.tiff")
      
      
      