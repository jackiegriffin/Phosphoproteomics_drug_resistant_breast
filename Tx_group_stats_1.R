<<<<<<< HEAD:Tx_group_stats_1.R
## install packages ----
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
if(!require("tidyverse")) {install.packages("tidyverse")}
library(tidyverse)


## Upload data ----
raw <- read.xlsx(xlsxFile = "with cetnred phosphosites Copy of ZR75-1FR tumor phosphoproteomics (Dec 2017)(9578) (1).xlsx",
                 sheet = 1) # tab 1 containinig summary phosphoproteomic data
## Remove & re-name columns  ----
vars_eliminate <- c("charge",
                    "match",
                    "tp",
                    "X8",
                    "xcorr_max",
                    "dcn_max",
                    "ppm_min",
                    "localization_max",
                    "longterm1A",
                    "pvalue.3D_control",
                    "average.log2.ratio.3D/control", 
                    "pvalue.longterm_control",
                    "average.log2.ratio.longterm/control",
                    "pvalue.longterm_3D",
                    "average.log2.ratio.longterm/3D",
                    "X29") # naming columns to remove
raw_filter <- raw[,!colnames(raw) %in% vars_eliminate] # remove named columns
raw_filter <- setNames(raw_filter,c("raw_sequence", "centered_sequence", "modsites", 
                                    "referece_id", "ctrl_a", "ctrl_b", "ctrl_c",
                                    "short_a", "short_b", "short_3", "long_a", "long_b",
                                    "long_c")) # change column names
colnames(raw_filter)
## Treatment group means ----
ctrl_means <- rep(NA, nrow(raw_filter)) # create vector filled with NA
for (i in 1:length(ctrl_means)) { # i in 1 - start at row 1 and loop every row for specified columns
  ctrl_means[i] <- mean(as.numeric(raw_filter[i,5:7]), na.rm = TRUE)
}
sum(is.na(ctrl_means))
raw_filter$ctrl_mean <- ctrl_means

short_means <- rep(NA, nrow(raw_filter))
for (i in 1:length(short_means)) {
  short_means[i] <- mean(as.numeric(raw_filter[i,8:10]), na.rm = TRUE)
}
sum(is.na(short_means))
raw_filter$short_mean <- short_means

long_means <- rep(NA, nrow(raw_filter))
for (i in 1:length(long_means)) {
  long_means[i] <- mean(as.numeric(raw_filter[i,11:13]), na.rm = TRUE)
}
sum(is.na(long_means))
raw_filter$long_mean <- long_means

## Remove incomplete data ----
raw_filter_complete <- raw_filter[complete.cases(raw_filter[ , 5:16]),] # keep complete data


## Fold change log2(B)-log2(A) ----
fc_56v3<- raw_filter_complete$long_mean - raw_filter_complete$short_mean # 56 - 3
raw_filter_complete$fc_56v3<-fc_56v3
fc_56vctrl<- raw_filter_complete$long_mean - raw_filter_complete$ctrl_mean # 56 - control
raw_filter_complete$fc_56vctrl<-fc_56vctrl
fc_3vctrl<- raw_filter_complete$short_mean - raw_filter_complete$ctrl_mean # 3 - control
raw_filter_complete$fc_3vctrl<-fc_3vctrl


## Independent t-test & BH adjustment ----
colnames(raw_filter_complete)
# 56 v 3
t_test_56v3 <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_56v3)) { 
  t_test_56v3[i] <- t_p_value(raw_filter_complete[i, 8:10], raw_filter_complete[i, 11:13])
}
raw_filter_complete$pvalue_56v3 <- t_test_56v3
padjust_BH_56v3<-p.adjust(t_test_56v3, method = 'BH')
raw_filter_complete$padjust_BH_56v3<-padjust_BH_56v3

# 56 v ctrl
t_test_56vctrl <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_56vctrl)) { 
  t_test_56vctrl[i] <- t_p_value(raw_filter_complete[i, 5:7], raw_filter_complete[i, 11:13])
}
raw_filter_complete$pvalue_56vctrl <- t_test_56vctrl
padjust_BH_56vctrl<-p.adjust(t_test_56vctrl, method = 'BH')
raw_filter_complete$padjust_BH_56vctrl<-padjust_BH_56vctrl

# 3 v ctrl
t_test_3vctrl <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_3vctrl)) { 
  t_test_3vctrl[i] <- t_p_value(raw_filter_complete[i, 5:7], raw_filter_complete[i, 8:10])
}
raw_filter_complete$pvalue_3vctrl <- t_test_3vctrl
padjust_BH_3vctrl<-p.adjust(t_test_3vctrl, method = 'BH')
raw_filter_complete$padjust_BH_3vctrl<-padjust_BH_3vctrl


## Write .csv file  ----
write.csv(raw_filter_complete, file="ZR751_phosphoproteomics_complete_stats_1.csv")
## Volcano plot ----
# Upload _1.csv  
vol_plot_test <-read.csv('ZR751_phosphoproteomics_complete_stats_1.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1) # 16933

# 56 v ctrl
with(vol_plot_test, plot(fc_56vctrl, -log10(pvalue_56vctrl), ylab = NA, xlab = NA, pch=20, 
                    xlim=c(-7,7), ylim = c(0,6))) #plot
mtext('Differential Peptide Abundance in PI3Ki Resistant Tumors', side = 1, line = -30, cex = 1.3) #title
mtext('-log10 pvalue', side = 2, line = 2.5, cex = 1.1) #y-axis
mtext('Log2 fold change (PI3KiR - Control)', side = 1, line = 2.5, cex = 1.1) #x-axis
#mtext('N=16933', side = 3, line = 0, cex = 1) #total peptides / sig (p<0.05 & fc>2): n=3076

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl< -1), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#E7298A")) #1515
text('n=16933', x = 5.5, y = 0.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1.1, font = NULL)

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl>1), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#1B9E77")) #1561
text(' ', x = 6, y = 5.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 2, col="#1B9E77", font = NULL)





# System Information
sessionInfo()










## Heat map ----





pal <- wes_palette("Zissou1", 100, type = "continuous")
install.packages("viridis")
library(viridis)
heat_palette = colorRampPalette(c("navy", "white", "firebrick3"))(200)
heat_palette = wes_palette
head(vol_plot)
colnames(vol_plot, color=heat_palette("Zissou1", 100, type="continuous"))

pheatmap(heat_plot_matrix)
heat_plot <- vol_plot[, 5:13]
heat_plot_matrix <- data.matrix(heat_plot)
pheatmap(heat_plot_matrix, main = 'Phosphopeptide abundance (Log2)
Independent t-test with pvalue_56vctrl-adjusted (BH) p-values < 0.05', 
         show_rownames = FALSE,
         color = colorRampPalette(c("blue", "white", "firebrick3"))(200), labels_col = c('','PI3ki; untreated 
(proliferating; control)',
                                                                                         '','PI3Ki; 3 days 
(cytotoxic; drug adaptation)',
                                                                                         '','','PI3Ki; 8 weeks 
(proliferating; resistant)','',''), angle_col = 0)
=======
## install packages ----
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
if(!require("tidyverse")) {install.packages("tidyverse")}
library(tidyverse)


## Upload data ----
raw <- read.xlsx(xlsxFile = "with cetnred phosphosites Copy of ZR75-1FR tumor phosphoproteomics (Dec 2017)(9578) (1).xlsx",
                 sheet = 1) # tab 1 containinig summary phosphoproteomic data
## Remove & re-name columns  ----
vars_eliminate <- c("charge",
                    "match",
                    "tp",
                    "X8",
                    "xcorr_max",
                    "dcn_max",
                    "ppm_min",
                    "localization_max",
                    "longterm1A",
                    "pvalue.3D_control",
                    "average.log2.ratio.3D/control", 
                    "pvalue.longterm_control",
                    "average.log2.ratio.longterm/control",
                    "pvalue.longterm_3D",
                    "average.log2.ratio.longterm/3D",
                    "X29") # naming columns to remove
raw_filter <- raw[,!colnames(raw) %in% vars_eliminate] # remove named columns
raw_filter <- setNames(raw_filter,c("raw_sequence", "centered_sequence", "modsites", 
                                    "referece_id", "ctrl_a", "ctrl_b", "ctrl_c",
                                    "short_a", "short_b", "short_3", "long_a", "long_b",
                                    "long_c")) # change column names
colnames(raw_filter)
## Treatment group means ----
ctrl_means <- rep(NA, nrow(raw_filter)) # create vector filled with NA
for (i in 1:length(ctrl_means)) { # i in 1 - start at row 1 and loop every row for specified columns
  ctrl_means[i] <- mean(as.numeric(raw_filter[i,5:7]), na.rm = TRUE)
}
sum(is.na(ctrl_means))
raw_filter$ctrl_mean <- ctrl_means

short_means <- rep(NA, nrow(raw_filter))
for (i in 1:length(short_means)) {
  short_means[i] <- mean(as.numeric(raw_filter[i,8:10]), na.rm = TRUE)
}
sum(is.na(short_means))
raw_filter$short_mean <- short_means

long_means <- rep(NA, nrow(raw_filter))
for (i in 1:length(long_means)) {
  long_means[i] <- mean(as.numeric(raw_filter[i,11:13]), na.rm = TRUE)
}
sum(is.na(long_means))
raw_filter$long_mean <- long_means

## Remove incomplete data ----
raw_filter_complete <- raw_filter[complete.cases(raw_filter[ , 5:16]),] # keep complete data


## Fold change log2(B)-log2(A) ----
fc_56v3<- raw_filter_complete$long_mean - raw_filter_complete$short_mean # 56 - 3
raw_filter_complete$fc_56v3<-fc_56v3
fc_56vctrl<- raw_filter_complete$long_mean - raw_filter_complete$ctrl_mean # 56 - control
raw_filter_complete$fc_56vctrl<-fc_56vctrl
fc_3vctrl<- raw_filter_complete$short_mean - raw_filter_complete$ctrl_mean # 3 - control
raw_filter_complete$fc_3vctrl<-fc_3vctrl


## Independent t-test & BH adjustment ----
colnames(raw_filter_complete)
# 56 v 3
t_test_56v3 <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_56v3)) { 
  t_test_56v3[i] <- t_p_value(raw_filter_complete[i, 8:10], raw_filter_complete[i, 11:13])
}
raw_filter_complete$pvalue_56v3 <- t_test_56v3
padjust_BH_56v3<-p.adjust(t_test_56v3, method = 'BH')
raw_filter_complete$padjust_BH_56v3<-padjust_BH_56v3

# 56 v ctrl
t_test_56vctrl <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_56vctrl)) { 
  t_test_56vctrl[i] <- t_p_value(raw_filter_complete[i, 5:7], raw_filter_complete[i, 11:13])
}
raw_filter_complete$pvalue_56vctrl <- t_test_56vctrl
padjust_BH_56vctrl<-p.adjust(t_test_56vctrl, method = 'BH')
raw_filter_complete$padjust_BH_56vctrl<-padjust_BH_56vctrl

# 3 v ctrl
t_test_3vctrl <- rep(NA, nrow(raw_filter_complete))
t_p_value <- function(...) {
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
}
for (i in  1:length(t_test_3vctrl)) { 
  t_test_3vctrl[i] <- t_p_value(raw_filter_complete[i, 5:7], raw_filter_complete[i, 8:10])
}
raw_filter_complete$pvalue_3vctrl <- t_test_3vctrl
padjust_BH_3vctrl<-p.adjust(t_test_3vctrl, method = 'BH')
raw_filter_complete$padjust_BH_3vctrl<-padjust_BH_3vctrl


## Write .csv file  ----
write.csv(raw_filter_complete, file="ZR751_phosphoproteomics_complete_stats_1.csv")
## Volcano plot ----
# Upload _1.csv  
vol_plot_test <-read.csv('ZR751_phosphoproteomics_complete_stats_1.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1) # 16933

# 56 v ctrl
with(vol_plot_test, plot(fc_56vctrl, -log10(pvalue_56vctrl), ylab = NA, xlab = NA, pch=20, 
                    xlim=c(-7,7), ylim = c(0,6))) #plot
mtext('Differential Peptide Abundance in PI3Ki Resistant Tumors', side = 1, line = -30, cex = 1.3) #title
mtext('-log10 pvalue', side = 2, line = 2.5, cex = 1.1) #y-axis
mtext('Log2 fold change (PI3KiR - Control)', side = 1, line = 2.5, cex = 1.1) #x-axis
#mtext('N=16933', side = 3, line = 0, cex = 1) #total peptides / sig (p<0.05 & fc>2): n=3076

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl< -1), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#E7298A")) #1515
text('n=16933', x = 5.5, y = 0.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1.1, font = NULL)

with(subset(vol_plot_test, padjust_BH_56vctrl<.05 & fc_56vctrl>1), 
     points(fc_56vctrl, -log10(pvalue_56vctrl), pch=20, col="#1B9E77")) #1561
text(' ', x = 6, y = 5.5, adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 2, col="#1B9E77", font = NULL)





# System Information
sessionInfo()










## Heat map ----





pal <- wes_palette("Zissou1", 100, type = "continuous")
install.packages("viridis")
library(viridis)
heat_palette = colorRampPalette(c("navy", "white", "firebrick3"))(200)
heat_palette = wes_palette
head(vol_plot)
colnames(vol_plot, color=heat_palette("Zissou1", 100, type="continuous"))

pheatmap(heat_plot_matrix)
heat_plot <- vol_plot[, 5:13]
heat_plot_matrix <- data.matrix(heat_plot)
pheatmap(heat_plot_matrix, main = 'Phosphopeptide abundance (Log2)
Independent t-test with pvalue_56vctrl-adjusted (BH) p-values < 0.05', 
         show_rownames = FALSE,
         color = colorRampPalette(c("blue", "white", "firebrick3"))(200), labels_col = c('','PI3ki; untreated 
(proliferating; control)',
                                                                                         '','PI3Ki; 3 days 
(cytotoxic; drug adaptation)',
                                                                                         '','','PI3Ki; 8 weeks 
(proliferating; resistant)','',''), angle_col = 0)
>>>>>>> ad152c3c1fbf3a5c7d63fe25d911afcc078811bb:data_processing_1.R
