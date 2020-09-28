
## install packages ----
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("viridis")) install.packages("viridis")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(ggplot2)
library(viridis)           
library(RColorBrewer)

display.brewer.pal(n=8, name = "YlGnBu")
brewer.pal(n = 8, name = "YlGnBu")


## Heatplot UP > 10 ----
# Phosphopeptides upregulated in PI3K inhibitor resistant tumors were input into 
# kinase prediction software (Networkin) kinase - substrate associations with 
# Networkin score > 10 are plotted

net_10 <- read.csv("keep/up_net_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
net_10_positive <- net_10[net_10$log2fc > 0,]

net_10_positive$uniprot <- gsub('HUMAN','', net_10_positive$uniprot)
net_10_positive$uniprot <- gsub('_','', net_10_positive$uniprot)
net_10_positive$prot_mod <- paste(net_10_positive$uniprot,net_10_positive$modsite, sep = "_")

ggplot(net_10_positive, aes(x=prot_mod, y = kinase, color = log2fc)) + 
  geom_point(shape=15, size=3.1) +
  scale_colour_gradient(low = "#7FCDBB", high = "#0C2C84") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size = 9)) +
  theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size = 9)) +
  
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 12)) +
  theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 12)) +
  theme(plot.title = element_text(size = 14)) +
  xlab("Phosphopeptide") +
  ylab("Predicted Kinase") +
  ggtitle("Predicted Kinases Driving PI3Ki Resistance in ZR75-1 Tumors")


## Heatplot DOWN > 10 ----
# Phosphopeptides downreglated in PI3K inhibitor resistant tumors were input into 
# kinase prediction software (Networkin) kinase - substrate associations with 
# Networkin score > 10 are plotted

net_down_final <- read.csv("keep/down_net_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#net_n <- net_n_p[net_n_p$log2fc < 0,]
display.brewer.pal(n=8, name = "YlOrRd")
brewer.pal(n = 8, name = "YlOrRd")

net_down_final$uniprot <- gsub('HUMAN','', net_down_final$uniprot)
net_down_final$uniprot <- gsub('_','', net_down_final$uniprot)
net_down_final$prot_mod <- paste(net_down_final$uniprot,net_down_final$modsite, sep = "_")

ggplot(net_down_final, aes(x=prot_mod, y = kinase, color = log2fc)) + 
  geom_point(shape=15, size = 3) +
  scale_colour_gradient(low = "#B10026", high = "#FEB24C") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size=8)) +
  theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size=8)) +
  
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 12)) +
  theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 12)) +
  theme(plot.title = element_text(size = 14)) +
  theme(plot.subtitle = element_text(size = 12)) +
  
  xlab("Phosphopeptides") +
  ylab("Predicted Kinase") +
  ggtitle("Predicted Kinases Downregulated in PI3Ki Resistant ZR75-1 Tumors")



## Heatplot down & up > 10 ----

display.brewer.pal(n=8, name = "YlOrRd")
brewer.pal(n = 8, name = "YlOrRd")

net_n <- read.csv("keep/net_threshold_10.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#net_n <- net_n_p[net_n_p$log2fc < 0,]

net_n$uniprot <- gsub('HUMAN','', net_n$uniprot)
net_n$uniprot <- gsub('_','', net_n$uniprot)
net_n$prot_mod <- paste(net_n$uniprot,net_n$modsite, sep = "_")

ggplot(net_n, aes(x=prot_mod, y = kinase, color = log2fc)) + 
  geom_point(shape=15, size = 2.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5, size=8)) +
  theme(axis.text.y = element_text(hjust =1, vjust = 0.5, size=8)) +
  
  theme(axis.title.y = element_text(hjust = 0.5, vjust = 1.5, size = 10)) +
  theme(axis.title.x = element_text(hjust = 0.5, vjust = 1.5, size = 10)) +
  theme(plot.title = element_text(size = 12)) +
  theme(plot.subtitle = element_text(size = 10)) +
  
  xlab("Phosphopeptide") +
  ylab("Predicted Kinase") +
  ggtitle("Predicted Kinases in PI3Ki Resistant ZR75-1 Tumors", subtitle = "Phosphopeptides uniquely Up & DOWN regulated in PI3Ki resistant tumors")





