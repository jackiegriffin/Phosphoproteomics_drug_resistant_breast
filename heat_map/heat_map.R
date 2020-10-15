
## Heat map ----
library(pheatmap)
library(plotly)
library(RColorBrewer)

heat <- read.csv("ZR751_phosphoproteomics_complete_stats_1.csv", stringsAsFactors = FALSE, row.names = 1, header = TRUE)
heat <- heat[, 5:13]
colnames(heat)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
metadata

heat_palette = colorRampPalette(c("navy", "white", "firebrick3"))(200)

out <- pheatmap(heat, 
                show_rownames=F, cluster_cols=T, cluster_rows=T,
                cex=1, clustering_distance_rows="euclidean", cex=1,
                clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
                annotation_col=metadata, legend_labels = NA, show_colnames = F, color = heat_palette)


# Correct publication formatting below ---- 

tiff ('filename.tiff', 
         units = "in", 
           height = 10,
           width = 6,
         res = 300, compression = "lzw")
dev.off ()
         
## done ----