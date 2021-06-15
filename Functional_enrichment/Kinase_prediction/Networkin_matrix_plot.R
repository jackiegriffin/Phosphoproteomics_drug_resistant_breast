
# Load libraries ----
  if (!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)
  if (!require("viridis")) install.packages("viridis")
  library(viridis)           
  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)

  display.brewer.pal(n=8, name = "YlGnBu")
  brewer.pal(n = 8, name = "YlGnBu")


# Heat matrices of substrate (raw peptide abundance) and predicted kinase (Networkin score>10) ----
    # Up-regulated (in Pi3Ki resistant tumors) only ----
      net_10 <- read.csv("Kinase_prediction/Input_data/up_net_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      net_10_positive <- net_10[net_10$log2fc > 0,]
      
      net_10_positive$uniprot <- gsub('HUMAN','', net_10_positive$uniprot)
      net_10_positive$uniprot <- gsub('_','', net_10_positive$uniprot)
      net_10_positive$prot_mod <- paste(net_10_positive$uniprot,net_10_positive$modsite, sep = "_")
      
      u <- ggplot(net_10_positive, aes(x=prot_mod, y = kinase, color = log2fc)) + 
        geom_point(shape=15, size=1.5) +
        scale_colour_gradient(low = "#7FCDBB", high = "#0C2C84") +
        coord_equal() +
        theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size = 3.5, color = "black")) +
        theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size = 3.5, color="black")) +
        
        theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
        theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
        theme(plot.title = element_text(size = 8)) +
        xlab("Phosphorylation Site") +
        ylab("Predicted Kinase \nNetworkin score > 10") +
        ggtitle("Phosphorylation events driving resistance, based on upregulated phosphopeptides \nin PI3Ki resistant tumors")
     
       ggsave(u, units = "in", height = 10, width = 6, dpi = 300, compression="lzw", filename = "Kinase_prediction/Output_data/Networkin_matrix_up.tiff")

    # Down-regulated (in PI3Ki resistant tumors) only ----
      net_down_final <- read.csv("Kinase_prediction/Input_data/down_net_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      #net_n <- net_n_p[net_n_p$log2fc < 0,]
      display.brewer.pal(n=8, name = "YlOrRd")
      brewer.pal(n = 8, name = "YlOrRd")
      
      net_down_final$uniprot <- gsub('HUMAN','', net_down_final$uniprot)
      net_down_final$uniprot <- gsub('_','', net_down_final$uniprot)
      net_down_final$prot_mod <- paste(net_down_final$uniprot,net_down_final$modsite, sep = "_")
      
      d <-ggplot(net_down_final, aes(x=prot_mod, y = kinase, color = log2fc)) + 
        geom_point(shape=15, size = 1) +
        scale_colour_gradient(low = "#B10026", high = "#FEB24C") +
        coord_equal() +
        theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size=3, color="black")) +
        theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size=3, color="black")) +
        
        theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
        theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
        theme(plot.title = element_text(size = 8)) +
        theme(plot.subtitle = element_text(size = 8)) +
        
        xlab("Phosphorylation Site") +
        ylab("Predicted Kinase \nNetworkin score > 10") +
        ggtitle("Phosphorylation events driving resistance, based on downregulated phosphopeptides \nin PI3Ki resistant tumors")

        ggsave(d, units = "in", height = 10, width = 6, dpi = 300, compression="lzw", filename = "Kinase_prediction/Output_data/Networkin_matrix_down.tiff")

    # Bi-directional plot ----
      display.brewer.pal(n=8, name = "YlOrRd")
      brewer.pal(n = 8, name = "YlOrRd")
    
        net_n <- read.csv("Kinase_prediction/Input_data/net_threshold_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
        #net_n <- net_n_p[net_n_p$log2fc < 0,]
        
        net_n$uniprot <- gsub('HUMAN','', net_n$uniprot)
        net_n$uniprot <- gsub('_','', net_n$uniprot)
        net_n$prot_mod <- paste(net_n$uniprot,net_n$modsite, sep = "_")
        
       combo<- ggplot(net_n, aes(x=prot_mod, y = kinase, color = log2fc)) + 
          geom_point(shape=15, size = 1.5) +
          scale_colour_gradient(low = "blue", high = "red") +
          coord_equal() +
          theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size=3.5, color="black")) +
          theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size=3.5, color = "black")) +
          
          theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
          theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 6)) +
          theme(plot.title = element_text(size = 6)) +
          theme(plot.subtitle = element_text(size = 6)) +
          
          xlab("Phosphopeptide") +
          ylab("Predicted Kinase") +
          ggtitle("Predicted Kinases in PI3Ki Resistant ZR75-1 Tumors", subtitle = "Phosphopeptides uniquely Up & DOWN regulated in PI3Ki resistant tumors")
    
        ggsave(combo, units = "in", height = 10, width = 6, dpi = 300, compression="lzw", filename = "Kinase_prediction/Output_data/Networkin_matrix.tiff")
        



