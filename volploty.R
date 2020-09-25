##  plotly vol plot ----

library(plotly)
library(manhattanly) #HapMap = reference data
install.packages("plotly")
library(plotly)
library(gridExtra)
voly <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
head(voly)

# Upload data ----
diff_df <-read.delim(file = "ZR751_phosphoproteomics_complete_stats_reformatted_2.csv", header = TRUE, sep = ",")

#     keep only the fields needed for this plot:
#     Gene names
#     Fold change value 
#     P-value
colnames(diff_df)

head(diff_df)
diff_df <- diff_df[c("protein_modsites", "fc_56vctrl", "pvalue_56vctrl")]

colnames(diff_df)
diff_df <- setNames(diff_df, c("external_gene_name", "Fold", "FDR"))
head(diff_df)


# add a grouping column; default value is "not significant"
diff_df["group"] <- "NotSignificant"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['Fold']) < 1.5 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df[which(diff_df['FDR'] > 0.05 & abs(diff_df['Fold']) > 1.5 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['Fold']) > 1.5 ),"group"] <- "Significant&FoldChange"

# Find and label the top peaks..
top_peaks <- diff_df[with(diff_df, order(Fold, FDR)),][1:5,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-Fold, FDR)),][1:5,])

# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks))) {
  m <- top_peaks[i, ]
  a[[i]] <- list(
    x = m[['Fold']],
    y = -log10(m[['FDR']]),
    text = m[['external_gene_name']],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

# make the Plot.ly plot
p <- plot_ly(data = diff_df, x = diff_df$Fold, y = -log10(diff_df$FDR), text = diff_df$external_gene_name, mode = "markers", 
             color = diff_df$group) %>% 
  layout(title = "Phosphopeptide Dysregulation in PI3K Inhibitor Resistant ZR75-1 Tumors", xaxis=list(title='Log2FC (PI3KiR - Control)'), 
         yaxis = list(title='-log10(p-value)')) %>%
  layout(annotations = a)
p 

# to save plot to a HTML file:
htmlwidgets::saveWidget(as_widget(p), "volploty_zr751.html")
