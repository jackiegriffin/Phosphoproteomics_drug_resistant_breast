

## morpheous ----
morpheous <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(complete_reform)
## write file to save and upload to morpheus to convert to gct ----
complete_morph_idv <- complete_reform[c(5:13, 28, 32)] # subset columns i want
complete_morph_idv <- complete_morph_idv[c(11, 10, 1:9)] # reorder columns for gct format uploat in morpheous

write.csv(complete_morph_idv, file="complete_morph_idv.csv")



## three day analyses ----

threed <- read.csv('complete_reform.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
threed_sig <- threed[threed$pvalue_3vctrl <= 0.05, ]
colnames(threed)

write.csv(threed_sig, file = "ctrl_3_sig.csv")












