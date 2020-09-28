<<<<<<< HEAD


#************************    REFORMAT FEATURE IDENTIFIERS FOR DOWNSTREAM SOFTWARE COMPATIBILITY    **********************
#
#   
#    PROTEIN_MODSITE format needed for NETWORKIN kinase prediction software 
#
#    GENE_MODSITE-p format need for PTM-SEA pathway enrichment analysis 
#
# ***********************************************************************************************************************


## Upload data ----
complete_stats_1 <- read.csv('ZR751_phosphoproteomics_complete_stats_1.csv',
                           header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(complete_stats_1)

## REFERENCE_ID ----
# GENPEPT_ACCESSION 
ref_protein_id_3 <- complete_stats_1$referece_id
names_i_want_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_3)) {                  # Isolate GENPEPT_ACCESSION from reference_ids
  split_line <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_3[i] <- split_line[[1]][2]
}
complete_stats_1$gene_id <- names_i_want_3             # add GENPEPT_ACCESSION col to df
head(complete_stats_1)

# UNIPROT_ID 
ref_protein_id_3 <-complete_stats_1$referece_id
names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_prot_3)){               # Isolate UNIPROT_ID from reference_id 
  split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
}
complete_stats_1$uniprot_id<- names_i_want_prot_3       # add UNIPROT_ID col to df
complete_stats_1$uniprot_id <- gsub(',.*','', complete_stats_1$uniprot_id)
head(complete_stats_1)

# UNIPROT_ID (w/o _HUMAN)
ref_protein_id_3 <-complete_stats_1$referece_id
names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_prot_3)){               # Isolate UNIPROT_ID from reference_id 
  split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
}
complete_stats_1$protein_id <- names_i_want_prot_3       # add UNIPROT_ID col to df
complete_stats_1$protein_id <- gsub('HUMAN','', complete_stats_1$protein_id) # Remove _HUMAN & replace with nothing
complete_stats_1$protein_id <- gsub('_','', complete_stats_1$protein_id)     
complete_stats_1$protein_id <- gsub(',.*','', complete_stats_1$protein_id)   
head(complete_stats_1)

## MODSITES ----

# split multiple modsites by :
complete_stats_1_modsplit <-separate_rows(complete_stats_1, modsites, sep = ":")
complete_stats_1_modsplit$modsites <- gsub(",.*","", complete_stats_1_modsplit$modsites) # remove everything after comma

# Paste protein_id and modsite with _ (NETWORKIN FORMAT)
complete_stats_1_modsplit$protein_modsites <- paste(complete_stats_1_modsplit$protein_id,
                                                    complete_stats_1_modsplit$modsites, sep = "_")

# Split individual modsites by residue (chr) and position (num)
complete_stats_1_modsplit$position <- gsub('[^0-9]','', complete_stats_1_modsplit$modsites) # remove all non-numeric string
complete_stats_1_modsplit$residue <- gsub('[0-9]','', complete_stats_1_modsplit$modsites) # remove numeric strings

# Paste gene_id and modsite with ; (PTMSEA FORMAT)
complete_stats_1_modsplit$ptmsea_format <- paste(complete_stats_1_modsplit$gene_id,
                                                complete_stats_1_modsplit$modsites, sep = ";")
complete_stats_1_modsplit$ptmsea_format <- paste0(complete_stats_1_modsplit$ptmsea_format, "-p") # add '-p' to the end of ptmsea_format col

## Write .csv file ----
write.csv(complete_stats_1_modsplit, file="ZR751_phosphoproteomics_complete_stats_reformatted_2.csv")




complete_reform <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
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











=======


#************************    REFORMAT FEATURE IDENTIFIERS FOR DOWNSTREAM SOFTWARE COMPATIBILITY    **********************
#
#   
#    PROTEIN_MODSITE format needed for NETWORKIN kinase prediction software 
#
#    GENE_MODSITE-p format need for PTM-SEA pathway enrichment analysis 
#
# ***********************************************************************************************************************


## Upload data ----
complete_stats_1 <- read.csv('ZR751_phosphoproteomics_complete_stats_1.csv',
                           header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(complete_stats_1)

## REFERENCE_ID ----
# GENPEPT_ACCESSION 
ref_protein_id_3 <- complete_stats_1$referece_id
names_i_want_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_3)) {                  # Isolate GENPEPT_ACCESSION from reference_ids
  split_line <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_3[i] <- split_line[[1]][2]
}
complete_stats_1$gene_id <- names_i_want_3             # add GENPEPT_ACCESSION col to df
head(complete_stats_1)

# UNIPROT_ID 
ref_protein_id_3 <-complete_stats_1$referece_id
names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_prot_3)){               # Isolate UNIPROT_ID from reference_id 
  split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
}
complete_stats_1$uniprot_id<- names_i_want_prot_3       # add UNIPROT_ID col to df
complete_stats_1$uniprot_id <- gsub(',.*','', complete_stats_1$uniprot_id)
head(complete_stats_1)

# UNIPROT_ID (w/o _HUMAN)
ref_protein_id_3 <-complete_stats_1$referece_id
names_i_want_prot_3 <- rep(NA, length(ref_protein_id_3))
for (i in 1:length(names_i_want_prot_3)){               # Isolate UNIPROT_ID from reference_id 
  split_line_prot_3 <- strsplit(ref_protein_id_3[i], '\\|')
  names_i_want_prot_3[i] <- split_line_prot_3[[1]][3]
}
complete_stats_1$protein_id <- names_i_want_prot_3       # add UNIPROT_ID col to df
complete_stats_1$protein_id <- gsub('HUMAN','', complete_stats_1$protein_id) # Remove _HUMAN & replace with nothing
complete_stats_1$protein_id <- gsub('_','', complete_stats_1$protein_id)     
complete_stats_1$protein_id <- gsub(',.*','', complete_stats_1$protein_id)   
head(complete_stats_1)

## MODSITES ----

# split multiple modsites by :
complete_stats_1_modsplit <-separate_rows(complete_stats_1, modsites, sep = ":")
complete_stats_1_modsplit$modsites <- gsub(",.*","", complete_stats_1_modsplit$modsites) # remove everything after comma

# Paste protein_id and modsite with _ (NETWORKIN FORMAT)
complete_stats_1_modsplit$protein_modsites <- paste(complete_stats_1_modsplit$protein_id,
                                                    complete_stats_1_modsplit$modsites, sep = "_")

# Split individual modsites by residue (chr) and position (num)
complete_stats_1_modsplit$position <- gsub('[^0-9]','', complete_stats_1_modsplit$modsites) # remove all non-numeric string
complete_stats_1_modsplit$residue <- gsub('[0-9]','', complete_stats_1_modsplit$modsites) # remove numeric strings

# Paste gene_id and modsite with ; (PTMSEA FORMAT)
complete_stats_1_modsplit$ptmsea_format <- paste(complete_stats_1_modsplit$gene_id,
                                                complete_stats_1_modsplit$modsites, sep = ";")
complete_stats_1_modsplit$ptmsea_format <- paste0(complete_stats_1_modsplit$ptmsea_format, "-p") # add '-p' to the end of ptmsea_format col

## Write .csv file ----
write.csv(complete_stats_1_modsplit, file="ZR751_phosphoproteomics_complete_stats_reformatted_2.csv")




complete_reform <- read.csv('ZR751_phosphoproteomics_complete_stats_reformatted_2.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
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











>>>>>>> ad152c3c1fbf3a5c7d63fe25d911afcc078811bb
