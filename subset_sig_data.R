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


## Upload data ----
sig <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(sig)

## Subset significant data ----

# 56 v control
sig_56vctrl <- sig[sig$padjust_BH_56vctrl < 0.05 & abs(sig$fc_56vctrl > 1.5), ]
colnames(sig_56vctrl)
sig_56vctrl <- sig_56vctrl[-c(8:10,15,17,19:21,24:25)]
colnames(sig_56vctrl)

# 56 v 3 day 
sig_56v3 <- sig[sig$padjust_BH_56v3 < 0.05 & abs(sig$fc_56v3 > 1.5), ]
colnames(sig_56v3)
sig_56v3 <- sig_56v3[-c(5:7,14,18:19,22:25)]
colnames(sig_56v3)

# 3 day v control: NS
# sig_3vctrl <- sig[sig$padjust_BH_3vctrl < 0.05 & abs(sig$fc_3vctrl3 > 1.5), ]
# colnames(sig_3vctrl)
# sig_3vctrl <- sig_3vctrl[-c(5:7,14,18:19,22:25)]
# colnames(sig_3vctrl)

## Venn Diagram ----

set1<-sig_56vctrl$ptmsea_format
set2<-sig_56v3$ptmsea_format

venn.diagram(
  x = list(set1, set2),
  category.names = c("LT-Control", "LT-ST"),
  filename = 'venn_diagram/Venn.tiff',
  imagetype="tiff" ,
  resolution = 100,
  compression = "lzw",
  cex = 8,
  cat.cex=8,
  fill = c("#F0027F", "#386CB0"),
  cat.pos = c(-50,60), 
  cat.default.pos = "outer",
  cat.dist = c(0.005, 0.005)
)


## subset 798 unique DPAs in 56 v ctrl group ----
# anti_join(x, y) return all rows from x with no matching values in y, keeping x
unique_DPAs_56vctrl <- anti_join (sig_56vctrl, sig_56v3, by='ptmsea_format')
write.csv(unique_DPAs_56vctrl, file = "unique_DPAs_56vctrl.csv")




