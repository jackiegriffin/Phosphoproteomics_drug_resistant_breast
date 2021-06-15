
# Load libraries ----
  if(!require("openxlsx")) {install.packages("openxlsx")}
  library(openxlsx)
  if(!require("pheatmap")) {install.packages("pheatmap")}
  library(pheatmap)
  if(!require("ggplot2")) {install.packages("ggplot2")}
  library(ggplot2)
  if(!require("calibrate")) {install.packages("calibrate")}
  library(calibrate)
  if(!require("tidyr")) {install.packages("tidyr")}
  library(tidyr)
  if(!require("tidyverse")) {install.packages("tidyverse")}
  library(tidyverse)
  if(!require("dplR")) {install.packages("dplR")}
  library(dplR)
  if(!require("plyr")) {install.packages("plyr")}
  library(plyr)
  if(!require("dplyr")) {install.packages("dplyr")}
  library(dplyr)
  if(!require("reshape2")) {install.packages("reshape2")}
  library(reshape2)
  if(!require("data.table")) {install.packages("data.table")}
  library(data.table)
  if(!require("RColorBrewer")) {install.packages("RcolorBrewer")}
  library(RColorBrewer)
  if(!require("gplots")) {install.packages("gplots")}
  library(gplots)
  if(!require("Lahman")) {install.packages("Lahman")}
  library(Lahman)
  if(!require("VennDiagram")) {install.packages("VennDiagram")}
  library(VennDiagram)
  if(!require("plotly")) {install.packages("plotly")}
  library(plotly)
  if(!require("manhattanly")) {install.packages("manhattanly")}
  library(manhattanly)
  if(!require("gridExtra")) {install.packages("gridExtra")}
  library(gridExtra)
  if(!require("tiff")) {install.packages("tiff")}
  library(tiff)


#######################  OBJECTIVE  ############################################
#                                                                              #
#    PLATFORM : LC-MS/MS phosphoproteomics                                     #
#    OBJECTIVE : 1. Clean data                                                 #
#                2. Calculate differential peptide abundance (DPA) = log2 FC   #
#                   Control = treatment naive, proliferating                   #
#                   3d = short term PI3Ki treated; cytotoxic, cell death       #
#                   56 = long term PI3Ki treated; resistant, proliferating     #
#                3. T-test, FDR adjusted P-value                               #
#                4. Reformat feature IDs for downstream software compatibility #
#                   Protein_modsite = NetworKin kinase prediction              #
#                   Gene_modsite-p  = PTM-SEA enrichment                       #
#                                                                              #
################################################################################

# Load data ----
  raw <- read.xlsx(xlsxFile = "Input_data/with cetnred phosphosites Copy of ZR75-1FR tumor phosphoproteomics (Dec 2017)(9578) (1).xlsx", sheet = 1) 
  raw_trim <- raw[-c(5:12,19,23:29)]
  raw_trim <- setNames(raw_trim,c("raw_sequence", "centered_sequence", "modsites", "referece_id", "ctrl_a", "ctrl_b", "ctrl_c",
                                    "short_a", "short_b", "short_3", "long_a", "long_b","long_c")) 


        
  
# Calculate treatment group means means ----
  ctrl_means <- rep(NA, nrow(raw_trim)) 
  for (i in 1:length(ctrl_means)) { 
    ctrl_means[i] <- mean(as.numeric(raw_trim[i,5:7]), na.rm = TRUE)
  }
  sum(is.na(ctrl_means))
  raw_trim$ctrl_mean <- ctrl_means
  
  short_means <- rep(NA, nrow(raw_trim))
  for (i in 1:length(short_means)) {
    short_means[i] <- mean(as.numeric(raw_trim[i,8:10]), na.rm = TRUE)
  }
  sum(is.na(short_means))
  raw_trim$short_mean <- short_means
  
  long_means <- rep(NA, nrow(raw_trim))
  for (i in 1:length(long_means)) {
    long_means[i] <- mean(as.numeric(raw_trim[i,11:13]), na.rm = TRUE)
  }
  sum(is.na(long_means))
  raw_trim$long_mean <- long_means

# Remove incomplete data ----
  raw_trim_complete <- raw_trim[complete.cases(raw_trim[ , 5:16]),] # keep complete data


# Calculate log2 fold change ----
  fc_56v3<- raw_trim_complete$long_mean - raw_trim_complete$short_mean    # 56d - 3d
  raw_trim_complete$fc_56v3<-fc_56v3
  fc_56vctrl<- raw_trim_complete$long_mean - raw_trim_complete$ctrl_mean  # 56d - control
  raw_trim_complete$fc_56vctrl<-fc_56vctrl
  fc_3vctrl<- raw_trim_complete$short_mean - raw_trim_complete$ctrl_mean  # 3d  - control
  raw_trim_complete$fc_3vctrl<-fc_3vctrl


# Independent t-test & BH adjustment ----

  # 56d vs 3d 
    t_test_56v3 <- rep(NA, nrow(raw_trim_complete))
    t_p_value <- function(...) {
      obj<-try(t.test(...), silent=TRUE) 
      if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
    }
    for (i in  1:length(t_test_56v3)) { 
      t_test_56v3[i] <- t_p_value(raw_trim_complete[i, 8:10], raw_trim_complete[i, 11:13])
    }
    raw_trim_complete$pvalue_56v3 <- t_test_56v3
    padjust_BH_56v3<-p.adjust(t_test_56v3, method = 'BH')
    raw_trim_complete$padjust_BH_56v3<-padjust_BH_56v3

  # 56d vs control
    t_test_56vctrl <- rep(NA, nrow(raw_trim_complete))
    t_p_value <- function(...) {
      obj<-try(t.test(...), silent=TRUE) 
      if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
    }
    for (i in  1:length(t_test_56vctrl)) { 
      t_test_56vctrl[i] <- t_p_value(raw_trim_complete[i, 5:7], raw_trim_complete[i, 11:13])
    }
    raw_trim_complete$pvalue_56vctrl <- t_test_56vctrl
    padjust_BH_56vctrl<-p.adjust(t_test_56vctrl, method = 'BH')
    raw_trim_complete$padjust_BH_56vctrl<-padjust_BH_56vctrl

  # 3d vs control
    t_test_3vctrl <- rep(NA, nrow(raw_trim_complete))
    t_p_value <- function(...) {
      obj<-try(t.test(...), silent=TRUE) 
      if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
    }
    for (i in  1:length(t_test_3vctrl)) { 
      t_test_3vctrl[i] <- t_p_value(raw_trim_complete[i, 5:7], raw_trim_complete[i, 8:10])
    }
    raw_trim_complete$pvalue_3vctrl <- t_test_3vctrl
    padjust_BH_3vctrl<-p.adjust(t_test_3vctrl, method = 'BH')
    raw_trim_complete$padjust_BH_3vctrl<-padjust_BH_3vctrl

    
    # Heat map ----
    
        # Cluster raw phosphopeptide abundance data
          colnames(raw_trim_complete)
          
          heat <- raw_trim_complete[, 5:13]
          metadata <- read.csv("Input_data/metadata.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
          heat_palette = colorRampPalette(c("navy", "white", "firebrick3"))(200)
          out <- pheatmap(heat,show_rownames=F, cluster_cols=T, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", 
                          cex=1, clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
                          annotation_col=metadata, legend_labels = NA, show_colnames = F, color = heat_palette,
                          units = "in", height = 10, width = 6, res = 300, compression = "lzw", filename = "Output_data/heatmap.tiff")   

    
# Reformat feature ids ----

    # GENPEPT_ACCESSION 
      head(raw_trim_complete)
      ref_protein_id_3 <- raw_trim_complete$referece_id
      names_i_want_3 <- rep(NA, length(ref_protein_id_3))
      for (i in 1:length(names_i_want_3)) {                                                             # Isolate GENPEPT_ACCESSION from reference_ids
        split_line <- strsplit(ref_protein_id_3[i], '\\|')
        names_i_want_3[i] <- split_line[[1]][2]
      }
      raw_trim_complete$gene_id <- names_i_want_3                                                       # add GENPEPT_ACCESSION col to df
      head(raw_trim_complete)
    
    # UNIPROT_ID 
      ref_protein_id_3 <-raw_trim_complete$referece_id
      names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
      for (i in 1:length(names_i_want_prot_3)){                                                          # Isolate UNIPROT_ID from reference_id 
        split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
        names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
      }
      raw_trim_complete$uniprot_id<- names_i_want_prot_3                                                 # add UNIPROT_ID col to df
      raw_trim_complete$uniprot_id <- gsub(',.*','', raw_trim_complete$uniprot_id)
      head(raw_trim_complete)
    
    # UNIPROT_ID (w/o _HUMAN)
      ref_protein_id_3 <-raw_trim_complete$referece_id
      names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
      for (i in 1:length(names_i_want_prot_3)){                                                          # Isolate UNIPROT_ID from reference_id 
        split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
        names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
      }
      raw_trim_complete$protein_id <- names_i_want_prot_3                                                # add UNIPROT_ID col to df
      raw_trim_complete$protein_id <- gsub('HUMAN','', raw_trim_complete$protein_id)                     # Remove _HUMAN & replace with nothing
      raw_trim_complete$protein_id <- gsub('_','', raw_trim_complete$protein_id)     
      raw_trim_complete$protein_id <- gsub(',.*','', raw_trim_complete$protein_id)   
      head(raw_trim_complete)
      
    # split multiple modsites by :
      raw_trim_complete_modsplit <-separate_rows(raw_trim_complete, modsites, sep = ":")
      raw_trim_complete_modsplit$modsites <- gsub(",.*","", raw_trim_complete_modsplit$modsites)         # remove everything after comma
      
    # Paste protein_id and modsite with _ (NETWORKIN FORMAT)
      raw_trim_complete_modsplit$protein_modsites <- paste(raw_trim_complete_modsplit$protein_id,
                                                          raw_trim_complete_modsplit$modsites, sep = "_")
      
    # Split individual modsites by residue (chr) and position (num)
      raw_trim_complete_modsplit$position <- gsub('[^0-9]','', raw_trim_complete_modsplit$modsites)      # remove all non-numeric string
      raw_trim_complete_modsplit$residue <- gsub('[0-9]','', raw_trim_complete_modsplit$modsites)        # remove numeric strings
      
    # Paste gene_id and modsite with ; (PTMSEA FORMAT)
      raw_trim_complete_modsplit$ptmsea_format <- paste(raw_trim_complete_modsplit$gene_id,
                                                       raw_trim_complete_modsplit$modsites, sep = ";")
      raw_trim_complete_modsplit$ptmsea_format <- paste0(raw_trim_complete_modsplit$ptmsea_format, "-p") # add '-p' to the end of ptmsea_format col
      

      
      
    
    # Volcano plot (plotly)----

          voly_plot <- raw_trim_complete_modsplit
          cols_i_want <-voly_plot[c("protein_modsites", "fc_56vctrl", "pvalue_56vctrl")] 
          cols_i_want <- setNames(cols_i_want, c("external_gene_name", "Fold", "FDR"))
        
          cols_i_want["group"] <- "Not Significant" # add a grouping column; default value is "not significant"
          cols_i_want[which(cols_i_want['FDR'] < 0.05 & abs(cols_i_want['Fold']) < 1.5 ),"group"] <- "Adj. P-value <0.05"
          cols_i_want[which(cols_i_want['FDR'] > 0.05 & abs(cols_i_want['Fold']) > 1.5 ),"group"] <- "Log2 Fold Change >1.5"
          cols_i_want[which(cols_i_want['FDR'] < 0.05 & abs(cols_i_want['Fold']) > 1.5 ),"group"] <- "Significant"
            
          top_peaks <- cols_i_want[with(cols_i_want, order(Fold, FDR)),][1:5,] # Find and label the top peaks
          top_peaks <- rbind(top_peaks, cols_i_want[with(cols_i_want, order(-Fold, FDR)),][1:5,])
            
          a <- list()                           # create empty list
          for (i in seq_len(nrow(top_peaks))) { # fill list with entries for each row in the df
            m <- top_peaks[i, ]                 
            a[[i]] <- list(
                          x = m[['Fold']],
                          y = -log10(m[['FDR']]),
                          text = m[['external_gene_name']],
                          xref = "x",
                          yref = "y",
                          showarrow = FALSE,
                          arrowhead = 0.0000001,
                          ax = 20,
                          ay = -40
              )
            }
            
          # Make plot
            p <- plot_ly(data = cols_i_want, type = "scatter", x = cols_i_want$Fold, y = -log10(cols_i_want$FDR), 
                         text = cols_i_want$external_gene_name, mode = "markers", color = cols_i_want$group, size=1) %>% 
                 layout(title = "Phosphopeptide Dysregulation in PI3K Inhibitor Resistant ZR75-1 Tumors", 
                        xaxis=list(title='Log2FC (PI3Ki Resistant - Control)'), 
                        yaxis = list(title='-log10(p-value)')) %>% 
                 layout(annotations = a)
            p 
            
          # Save plot as HTML file
            htmlwidgets::saveWidget(as_widget(p), "Output_data/C:volcano_ploty.html")
            

      
      
# Subset significant data ----
      
      # 56d vs control
        raw_trim_complete_modsplit_56vctrl <- raw_trim_complete_modsplit[raw_trim_complete_modsplit$padjust_BH_56vctrl < 0.05 & abs(raw_trim_complete_modsplit$fc_56vctrl > 1.5), ]
        colnames(raw_trim_complete_modsplit_56vctrl)
        raw_trim_complete_modsplit_56vctrl <- raw_trim_complete_modsplit_56vctrl[-c(8:10,15,17,19:21,24:25)]
        colnames(raw_trim_complete_modsplit_56vctrl)
      
      # 56d vs 3d 
        raw_trim_complete_modsplit_56v3 <- raw_trim_complete_modsplit[raw_trim_complete_modsplit$padjust_BH_56v3 < 0.05 & abs(raw_trim_complete_modsplit$fc_56v3 > 1.5), ]
        colnames(raw_trim_complete_modsplit_56v3)
        raw_trim_complete_modsplit_56v3 <- raw_trim_complete_modsplit_56v3[-c(5:7,14,18:19,22:25)]
        colnames(raw_trim_complete_modsplit_56v3)
        
      # 3d vs control: NS
        # raw_trim_complete_modsplit_3vctrl <- raw_trim_complete_modsplit[raw_trim_complete_modsplit$padjust_BH_3vctrl < 0.05 & abs(raw_trim_complete_modsplit$fc_3vctrl3 > 1.5), ]
        # colnames(raw_trim_complete_modsplit_3vctrl)
        # raw_trim_complete_modsplit_3vctrl <- raw_trim_complete_modsplit_3vctrl[-c(5:7,14,18:19,22:25)]
        # colnames(raw_trim_complete_modsplit_3vctrl)
        
      # Venn Diagram 
        set1<-raw_trim_complete_modsplit_56vctrl$ptmsea_format
        set2<-raw_trim_complete_modsplit_56v3$ptmsea_format
          
        venn.diagram(
          x = list(set1, set2),
          category.names = c("56 v control", "56 v 3"),
          filename = 'venn_diagram/Venn.tiff',
          imagetype="tiff" ,
          resolution = 110,
          compression = "lzw",
          cex = 8,
          cat.cex=8,
          fill = c("#F0027F", "#386CB0"),
          cat.pos = c(-50,60), 
          cat.default.pos = "outer",
          cat.dist = c(0.000005, 0.008)
        )
        
      
# Subset significant data uniquely found in 56d vs control----
      unique_DPAs_56vctrl <- anti_join (raw_trim_complete_modsplit_56vctrl, raw_trim_complete_modsplit_56v3, by='ptmsea_format')
      write.csv(unique_DPAs_56vctrl, file = "Output_data/unique_DPAs_56vctrl.csv")
      
      
      
      
      
      






