
## Install packages ----
library(ggplot2)
library(dplR)
library(tidyverse)
library(tidyr)
library(data.table)
library(reshape2)
library(calibrate)
library(Lahman)
library(plyr)
library(RColorBrewer)
library(VennDiagram)

## Volcano plot ----
# Upload _1.csv  
vol_plot_test <-read.csv('ZR751_phosphoproteomics_complete_stats_1.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1) # 16933

# 56 v ctrl
with(vol_plot_test, plot(fc_56vctrl, -log10(pvalue_56vctrl), ylab = NA, xlab = NA, pch=20, 
                         xlim=c(-7,7), ylim = c(0,6))) #plot
#mtext('Phosphopeptide Dysregulation in PI3K Inhibitor Resistant ZR75-1 Tumors', side = 1, line = -26, cex = 1.3) #title
mtext('-log10 pvalue', side = 2, line = 2.5, cex = 1.1) #y-axis
mtext('Log2FC (PI3Ki Resistant - Control)', side = 1, line = 2.5, cex = 1.1) #x-axis
mtext('Phosphopeptide Dysregulation in PI3K Inhibitor Resistant ZR75-1 Tumors', side = 1, line = -30, cex = 1.4) #title

#mtext('N=16933', side = 3, line = 0, cex = 1) #total peptides / sig (p<0.05 & fc>2): n=3076

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl< -1.5), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#E7298A")) #1515
text('n=16933', x = 5.5, y = 0.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1.1, font = NULL)

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl>1.5), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#1B9E77")) #1561
text(' ', x = 6, y = 5.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 2, col="#1B9E77", font = NULL)
