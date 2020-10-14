# Install  ----
library(plotly)
library(manhattanly)
library(gridExtra)
library(tiff)

# Upload data ----
voly_plot <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
cols_i_want <-voly_plot[c("protein_modsites", "fc_56vctrl", "pvalue_56vctrl")] # gene names, FC, P-value
cols_i_want <- setNames(cols_i_want, c("external_gene_name", "Fold", "FDR"))

# add a grouping column; default value is "not significant"
cols_i_want["group"] <- "Not Significant"
cols_i_want[which(cols_i_want['FDR'] < 0.05 & abs(cols_i_want['Fold']) < 1.5 ),"group"] <- "Adj. P-value <0.05"
cols_i_want[which(cols_i_want['FDR'] > 0.05 & abs(cols_i_want['Fold']) > 1.5 ),"group"] <- "Log2 Fold Change >1.5"
cols_i_want[which(cols_i_want['FDR'] < 0.05 & abs(cols_i_want['Fold']) > 1.5 ),"group"] <- "Significant"

# Find and label the top peaks..
top_peaks <- cols_i_want[with(cols_i_want, order(Fold, FDR)),][1:5,]
top_peaks <- rbind(top_peaks, cols_i_want[with(cols_i_want, order(-Fold, FDR)),][1:5,])

a <- list() # creast empty list
for (i in seq_len(nrow(top_peaks))) { # fill list with entries for each row in the df
  m <- top_peaks[i, ]                 # eady list entry is another list with named items that iwll be used by Plot_ly
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

p <- plot_ly(data = cols_i_want, type = "scatter", x = cols_i_want$Fold, y = -log10(cols_i_want$FDR), text = cols_i_want$external_gene_name, mode = "markers", 
             color = cols_i_want$group, size=1) %>% 
  layout(title = "Phosphopeptide Dysregulation in PI3K Inhibitor Resistant ZR75-1 Tumors", 
         xaxis=list(title='Log2FC (PI3Ki Resistant - Control)'), 
         yaxis = list(title='-log10(p-value)')) %>% 
           layout(annotations = a)
p 


# to save plot to a HTML file:
htmlwidgets::saveWidget(as_widget(p), "C:volcano_ploty_thick_bubbles.html")

## END ----
