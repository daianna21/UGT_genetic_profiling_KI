

library(here)
library(readr)
library(rlang)
library(ggplot2)
library(cowplot)

####################################################################################################
##                                1. Data Processing and Correction
####################################################################################################

## UGT genes
UGT_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10',
               'UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28',
               'UGT3A1', 'UGT3A2', 
               'UGT8')


## Process and build variables for each gene 
for (gene in UGT_genes){
  
  ## Read data
  UGT_gene_info <- read.csv(paste0("raw-data/", gene, "_gnomAD_v2.1.1.csv"))
  
 
  ## Change variable names
  colnames(UGT_gene_info) <- gsub('_$+', '', gsub('\\.+', '_', colnames(UGT_gene_info)))
  ## Specify 'or' in Latino/Admixed American and African/African American population names
  colnames(UGT_gene_info) <- gsub('ino_', 'ino_or_', colnames(UGT_gene_info))
  colnames(UGT_gene_info) <- gsub('can_Af', 'can_or_Af', colnames(UGT_gene_info))
  
  ## Create unique variant IDs with Chr-Position-Reference-Alternate info
  UGT_gene_info$Variant_ID <- paste(UGT_gene_info$Chromosome, UGT_gene_info$Position, UGT_gene_info$Reference, UGT_gene_info$Alternate, sep='-')
  
  assign(paste0(gene, '_data'), UGT_gene_info)
  #save(UGT_gene_info, file=paste0('processed-data/', gene, '_processed_data.Rdata'))

}



################################################################################
##           1.1  Correct variant counts within each gene family
################################################################################

# _______________________________________________________________________________
#  1.1.1 Detect overlapping variants between different genes of the same family 
# _______________________________________________________________________________

############################
####      UGT1 genes 
############################
UGT1_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10')

## Confirm all variants in UGT1 genes are in the same chr
sapply(UGT1_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT1A1  UGT1A3  UGT1A4  UGT1A5  UGT1A6  UGT1A7  UGT1A8  UGT1A9 UGT1A10 
#      2       2       2       2       2       2       2       2       2 


## For each variant, determine in what UGT1 genes' txs they are present 

########## *** UGT1A1 and UGT1A8 txs are swapped; take only the canonical tx of these genes for all of their variants *** ##########

UGT1A1_data$Transcript <- rep('ENST00000305208.5', dim(UGT1A1_data)[1])
UGT1A8_data$Transcript <- rep('ENST00000373450.4', dim(UGT1A8_data)[1])

## Unique variants in UGT1 genes
unique_UGT1_variants <- unique(c(UGT1A1_data$Variant_ID, UGT1A3_data$Variant_ID, UGT1A4_data$Variant_ID, UGT1A5_data$Variant_ID, UGT1A6_data$Variant_ID, UGT1A7_data$Variant_ID, UGT1A8_data$Variant_ID, UGT1A9_data$Variant_ID, UGT1A10_data$Variant_ID))

## Compute in which transcript each variant appears for each gene 
## (if present in a gene, variants are expected to appear only in one transcript since they are all unique in each gene dataset)
UGT1_variants <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(UGT1_variants) <- c(UGT1_genes, 'VEP_Annotation', 'Position')

i=1
for (variant in unique_UGT1_variants){
  
  ## Search variant in each gene; if present, annotate tx
  variant_txs <- sapply(UGT1_genes, function(gene){
    tx <- eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'Transcript']
    if (length(tx)==0) {NA}
    else {tx}
    })
  
  UGT1_variants[i, 1:9] <- variant_txs
  
  ## VEP annotation for variant in each gene that contains it
  variant_anno <- unlist(sapply(UGT1_genes, function(gene){
                           eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
                         }))
  ## Evaluate if the anno of a shared variable is the same across all genes
  if (length(unique(variant_anno))==1){
    UGT1_variants$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT1_variants$VEP_Annotation[i] <- toString(variant_anno)
  }
  
  ## Add variant position
  UGT1_variants$Position[i] <- as.numeric(strsplit(variant, '-')[[1]][2])
  
  i=i+1
}

rownames(UGT1_variants) <- unique_UGT1_variants


############################
####    UGT2 genes 
############################
UGT2_genes <- c('UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28')

## Confirm all UGT2 variants are in the same chr
sapply(UGT2_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT2A1  UGT2A2  UGT2A3  UGT2B4  UGT2B7 UGT2B10 UGT2B11 UGT2B15 UGT2B17 UGT2B28 
#      4       4       4       4       4       4       4       4       4       4 

unique_UGT2_variants <- unique(c(UGT2A1_data$Variant_ID, UGT2A2_data$Variant_ID, UGT2A3_data$Variant_ID, 
                                          UGT2B4_data$Variant_ID, UGT2B7_data$Variant_ID, UGT2B10_data$Variant_ID, UGT2B11_data$Variant_ID, 
                                          UGT2B15_data$Variant_ID, UGT2B17_data$Variant_ID, UGT2B28_data$Variant_ID))
UGT2_variants <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(UGT2_variants) <- c(UGT2_genes, 'VEP_Annotation', 'Position')

i=1
for (variant in unique_UGT2_variants){
  
  variant_txs <- sapply(UGT2_genes, function(gene){
    tx <- eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'Transcript']
    if (length(tx)==0) {NA}
    else {tx}
  })
  
  UGT2_variants[i, 1:10] <- variant_txs
  
  variant_anno <- unlist(sapply(UGT2_genes, function(gene){
    eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
  }))
  
  if (length(unique(variant_anno))==1){
    UGT2_variants$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT2_variants$VEP_Annotation[i] <- toString(variant_anno)
  }
  
  UGT2_variants$Position[i] <- as.numeric(strsplit(variant, '-')[[1]][2])
  
  i=i+1
}

rownames(UGT2_variants) <- unique_UGT2_variants


############################
####    UGT3 genes 
############################
UGT3_genes <- c('UGT3A1', 'UGT3A2')

## Confirm all UGT3 variants are in the same chr
sapply(UGT3_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT3A1 UGT3A2 
#      5      5 

unique_UGT3_variants <- unique(c(UGT3A1_data$Variant_ID, UGT3A2_data$Variant_ID))
UGT3_variants <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(UGT3_variants) <- c(UGT3_genes, 'VEP_Annotation', 'Position')

i=1
for (variant in unique_UGT3_variants){
  
  variant_txs <- sapply(UGT3_genes, function(gene){
    tx <- eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'Transcript']
    if (length(tx)==0) {NA}
    else {tx}
  })
  
  UGT3_variants[i, 1:2] <- variant_txs
  
  variant_anno <- unlist(sapply(UGT3_genes, function(gene){
    eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
  }))

  if (length(unique(variant_anno))==1){
    UGT3_variants$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT3_variants$VEP_Annotation[i] <- toString(variant_anno)
  }
  
  UGT3_variants$Position[i] <- as.numeric(strsplit(variant, '-')[[1]][2])
  
  i=i+1
}

rownames(UGT3_variants) <- unique_UGT3_variants


############################
####    UGT8 genes 
############################
unique(UGT8_data$Chromosome)
#    4 

## Variant_ID's are unique so they all appear once in UGT8
which(duplicated(UGT8_data$Variant_ID))
#  integer(0)


# ________________________________________________________________________________________
#  1.1.1.1 Verify there are no overlapping variants between genes from different families

length(intersect(unique_UGT1_variants, unique_UGT2_variants))
#  0
length(intersect(unique_UGT1_variants, unique_UGT3_variants))
#  0
length(intersect(unique_UGT1_variants, unique_UGT8_variants))
#  0
length(intersect(unique_UGT2_variants, unique_UGT3_variants))
#  0
length(intersect(unique_UGT2_variants, unique_UGT8_variants))
#  0
length(intersect(unique_UGT3_variants, unique_UGT8_variants))
#  0




# _______________________________________________________________________________________________
#  1.1.2 Quantify the proportion of variants present in the canonical transcripts* of UGT genes
# _______________________________________________________________________________________________

##  *Canonical transcript in GRCh38.p14 or the most frequent/unique transcript for the gene

## Define canonical/most common tx for each gene
canonical_UGT1_txs <- list('UGT1A1'= 'ENST00000305208.5', 'UGT1A3'='ENST00000482026.1', 'UGT1A4'='ENST00000373409.3', 
                           'UGT1A5'='ENST00000373414.3', 'UGT1A6'='ENST00000305139.6', 'UGT1A7'='ENST00000373426.3', 
                           'UGT1A8'= 'ENST00000373450.4', 'UGT1A9'= 'ENST00000354728.4', 'UGT1A10'='ENST00000344644.5')

canonical_UGT2_txs <- list('UGT2A1'= 'ENST00000503640.1', 'UGT2A2'='ENST00000457664.2', 'UGT2A3'='ENST00000251566.4', 
                           'UGT2B4'='ENST00000305107.6', 'UGT2B7'='ENST00000305231.7', 'UGT2B10'='ENST00000265403.7', 
                           'UGT2B11'= 'ENST00000446444.1', 'UGT2B15'= 'ENST00000338206.5', 'UGT2B17'='ENST00000317746.2', 
                           'UGT2B28'='ENST00000335568.5')

canonical_UGT3_txs <- list('UGT3A1'= 'ENST00000274278.3', 'UGT3A2'='ENST00000282507.3')

canonical_UGT8_txs <- list('UGT8'= 'ENST00000310836.6')


## % of variants in a gene that are present in the canonical tx of the gene

###################
####  UGT1 genes
###################

## For each gene
sapply(UGT1_genes, function(x){
  length(which(UGT1_variants[[x]] == canonical_UGT1_txs[[x]]))/length(which(!is.na(UGT1_variants[[x]]))) * 100}
  )
# UGT1A1    UGT1A3    UGT1A4    UGT1A5    UGT1A6    UGT1A7    UGT1A8    UGT1A9   UGT1A10 
# 100.00000 100.00000 100.00000 100.00000  96.03524  95.40682 100.00000 100.00000  97.13262 


## In the whole family
## (For shared variants, if they are in canonical txs of all genes in which they are)
for (i in 1:dim(UGT1_variants)[1]){
  var_txs <- UGT1_variants[i, which(!is.na(UGT1_variants[i, 1:9]))] 
  ## If all txs of a variant are canonical
  if (length(which(var_txs %in% unlist(canonical_UGT1_txs))) == length(var_txs)){
    UGT1_variants$Canonical_txs[i] <- TRUE
  }
  else{
    UGT1_variants$Canonical_txs[i] <- FALSE
  }
}
table(UGT1_variants$Canonical_txs)['TRUE'] / dim(UGT1_variants)[1] * 100
## % of total UGT1 variants that appear in canonical txs
# 98.31384 


###################
####  UGT2 genes
###################

sapply(UGT2_genes, function(x){
  length(which(UGT2_variants[[x]] == canonical_UGT2_txs[[x]]))/length(which(!is.na(UGT2_variants[[x]]))) * 100}
)
# UGT2A1    UGT2A2    UGT2A3    UGT2B4    UGT2B7   UGT2B10   UGT2B11   UGT2B15   UGT2B17   UGT2B28 
# 97.41268  99.86339  98.91041  94.64752  99.47507  99.88208 100.00000 100.00000 100.00000 100.00000  

for (i in 1:dim(UGT2_variants)[1]){
  var_txs <- UGT2_variants[i, which(!is.na(UGT2_variants[i, 1:10]))] 
  if (length(which(var_txs %in% unlist(canonical_UGT2_txs))) == length(var_txs)){
    UGT2_variants$Canonical_txs[i] <- TRUE
  }
  else{
    UGT2_variants$Canonical_txs[i] <- FALSE
  }
}
table(UGT2_variants$Canonical_txs)['TRUE'] / dim(UGT2_variants)[1] * 100
## % of total UGT2 variants in canonical txs
# 98.95475


###################
####  UGT3 genes
###################

sapply(UGT3_genes, function(x){
  length(which(UGT3_variants[[x]] == canonical_UGT3_txs[[x]]))/length(which(!is.na(UGT3_variants[[x]]))) * 100}
)
# UGT3A1   UGT3A2 
# 93.36735 99.56012 

for (i in 1:dim(UGT3_variants)[1]){
  var_txs <- UGT3_variants[i, which(!is.na(UGT3_variants[i, 1:2]))] 
  if (length(which(var_txs %in% unlist(canonical_UGT3_txs))) == length(var_txs)){
    UGT3_variants$Canonical_txs[i] <- TRUE
  }
  else{
    UGT3_variants$Canonical_txs[i] <- FALSE
  }
}
table(UGT3_variants$Canonical_txs)['TRUE'] / dim(UGT3_variants)[1] * 100
## % of total UGT3 variants in canonical txs
# 96.24829 


###################
####  UGT8 gene
###################

sapply(UGT8_genes, function(x){
  length(which(UGT8_data$Transcript == canonical_UGT8_txs[[x]]))/length(which(!is.na(UGT8_data$Transcript))) * 100}
)
# UGT8 
# 97.22222 


######################
####  All UGT genes
######################

## % of total UGT variants in canonical txs
vars_in_canonical_txs <- sum(table(UGT1_variants$Canonical_txs)['TRUE'], table(UGT2_variants$Canonical_txs)['TRUE'], table(UGT3_variants$Canonical_txs)['TRUE'], 
    length(which(UGT8_data$Transcript == canonical_UGT8_txs)))
total_vars <- sum(dim(UGT1_variants)[1], dim(UGT2_variants)[1], dim(UGT3_variants)[1], dim(UGT8_data)[1])

vars_in_canonical_txs / total_vars * 100
# [1] 98.39432




# _____________________________________________________________________________
#  1.1.3 Manual annotation of variants based on boundaries of canonical txs
# _____________________________________________________________________________

## Conserve variants in canonical txs only

UGT1_variants_canonical <- UGT1_variants[UGT1_variants$Canonical_txs==TRUE, ]
UGT2_variants_canonical <- UGT2_variants[UGT2_variants$Canonical_txs==TRUE, ]
UGT3_variants_canonical <- UGT3_variants[UGT3_variants$Canonical_txs==TRUE, ]
UGT8_variants_canonical <- UGT8_data[UGT8_data$Canonical_txs==TRUE, ]
rownames(UGT8_variants_canonical) <- UGT8_variants_canonical$Variant_ID

## For a given variant, examine where it's located with respect to txs boundaries

## Function to evaluate the location of a variant according to the boundaries of a gene transcript
location_determination <- function(variant_pos, tx, feature){

  ## Read tx data
  tx_seq_data <- as.data.frame(read.csv(paste0("raw-data/Tx_seq_data/", tx, "_seq_data.csv")))
  tx_seq_data <- tx_seq_data[-c(1, dim(tx_seq_data)[1]),]
  ## Column names
  colnames(tx_seq_data)[-1] <- gsub('\\.+', '_', colnames(tx_seq_data)[-1])
  ## Char to integer for End, Start and Length
  tx_seq_data$Start <- as.numeric(gsub(',', '', tx_seq_data$Start))
  tx_seq_data$End <- as.numeric(gsub(',', '', tx_seq_data$End))
  tx_seq_data$Length <- as.numeric(gsub(',', '', tx_seq_data$Length))
  
  
  ## Overall transcript composition:  
  ## (5' upstream seq) ... [5'UTR + Exon A] + Intron A-B + Exon B + Intron B-C + Exon C + Intron C-D + Exon D + ... + [Exon Z + 3'UTR] ... (3' downstream seq)
  rownames(tx_seq_data) <- tx_seq_data$Exon_Intron
  
  ## Add exon/intron info
  tx_seq_data$Exon_Intron <- sapply(tx_seq_data$Exon_Intron, function(x){
    if(length(grep('Intron', x))==1){x}  else{paste0('Exon ', tx_seq_data[x, 'No.'])}})
  
  ## Specify we have 1st Exon + 5'UTR and last Exon + 3'UTR
  tx_seq_data$Exon_Intron[1] <- '5-UTR + First Exon'
  tx_seq_data$Exon_Intron[dim(tx_seq_data)[1]] <- 'Last Exon + 3-UTR'
  rownames(tx_seq_data) <- tx_seq_data$Exon_Intron
  
  ## Evaluate if the variant is within a tx exon/intron/UTR or outside
  location_within_tx <- unlist(sapply(rownames(tx_seq_data), function(x){
                                  if(variant_pos %in% tx_seq_data[x, 'Start']:tx_seq_data[x, 'End']){x}}))
  
  ## If the variant is not within tx region then is either in the 5' upstream or 3' downstream sequence
  if (is.null(location_within_tx))  {
    if (variant_pos < tx_seq_data["5-UTR + First Exon", 'Start']){
      location <- '5\' upstream'
    }
    
    if (variant_pos > tx_seq_data["Last Exon + 3-UTR", 'End']){
      location <- '3\' downstream'
    }
  }
  
  else {
    location <- location_within_tx
  }
  
  ## Extract boundaries of a tx feature
  if (! is.null(feature)){
    range <- tx_seq_data[feature, c('Start', 'End')]
  }
  else{
    range <- NULL
  }
  
  return(list(location, range))
}


####################################
####  Annotation of UGT1 variants
####################################

## For each variant, obtain in which txs it is present and its location for each one

for (i in 1:dim(UGT1_variants_canonical)[1]){
 pos <- UGT1_variants_canonical$Position[i]
 txs <- UGT1_variants_canonical[i,which(!is.na(UGT1_variants_canonical[i, 1:9]))]
 locations <- unlist(sapply(txs, function(x){location_determination(pos, x, NULL)[[1]][1]}))
 if (length(unique(locations))==1){
   UGT1_variants_canonical$Location_in_txs[i] <- unique(locations)
 }
 else
   UGT1_variants_canonical$Location_in_txs[i] <- toString(locations)
}

table(UGT1_variants_canonical[, c('VEP_Annotation', 'Location_in_txs')])

#########################  All variants  #########################
## Compare VEP annotation vs manual annotation
#                                                                                              Location_in_txs
# VEP_Annotation                        5-UTR + First Exon      5-UTR + First Exon, Intron 1-2       5' upstream        5' upstream, 5-UTR + First Exon
# 3_prime_UTR_variant                                    0                                   0                 0                                     0
# 5_prime_UTR_variant                                   60                                   0                28                                     7
# frameshift_variant                                    56                                   8                19                                    47
# inframe_deletion                                       7                                   2                 1                                     5
# inframe_insertion                                      1                                   0                 0                                     0
# intron_variant                                         0                                   0                64                                     0
# intron_variant, 5_prime_UTR_variant                    0                                  11                 0                                    12
# missense_variant                                     821                                 164               208                                   609
# splice_acceptor_variant                                0                                   0                 0                                     0
# splice_donor_variant                                   1                                   0                 1                                     0
# splice_region_variant                                  1                                   1                 1                                     0
# start_lost                                             4                                   0                 2                                     6
# stop_gained                                           51                                   5                18                                    27
# synonymous_variant                                   349                                  62               101                                   269

#                                                                                              Location_in_txs
# VEP_Annotation                        5' upstream, Intron 1-2      Exon 2  Exon 3   Exon 4   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5   Last Exon + 3-UTR
#   3_prime_UTR_variant                                       0         0         0        0            0            0            0            0                  10
#   5_prime_UTR_variant                                       0         0         0        0            0            0            0            0                   0
#   frameshift_variant                                        0         0         2        3            0            0            0            0                   3
#   inframe_deletion                                          0         0         1        0            0            0            0            0                   1
#   inframe_insertion                                         0         0         0        0            0            0            0            0                   0
#   intron_variant                                           54         0         0        0          112           30           25           50                   0
#   intron_variant, 5_prime_UTR_variant                       1         0         0        0            0            0            0            0                   0
#   missense_variant                                          0        22        24       40            0            0            0            0                  79
#   splice_acceptor_variant                                   0         0         0        0            1            0            0            1                   0
#   splice_donor_variant                                      2         0         0        0            6            0            1            0                   0
#   splice_region_variant                                    11         0         1        0           16            2            3            3                   0
#   start_lost                                                0         0         0        0            0            0            0            0                   0
#   stop_gained                                               0         1         4     .  0            0            0            0            0                   1
#   synonymous_variant                                        0         2        13       23            0            0            0            0                  38
  

#########################  Tx anno of variants common in all UGT1 genes  #########################
UGT1_shared_variants_allGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==9)]
table(UGT1_variants_canonical[UGT1_shared_variants_allGenes, 'Location_in_txs'])
# Exon 2            Exon 3            Exon 4        Intron 1-2        Intron 2-3        Intron 3-4        Intron 4-5       Last Exon + 3-UTR 
#     25                45                66                12                32                29                33                     132 

## Position of shared variants in Intron 1-2 region
pos <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs=='Intron 1-2' & rownames(UGT1_variants_canonical) %in% UGT1_shared_variants_allGenes), 'Position']
## Confirm this overlapping intronic region starts after exon 1 of UGT1A1 and finishes before exon 2 of all UGT1 genes
min(pos) > location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], '5-UTR + First Exon')[[2]]['End']
# TRUE
max(pos) < location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], 'Exon 2')[[2]]['Start']
# TRUE


#########################  Tx anno of variants common in 4 UGT1 genes  #########################
UGT1_shared_variants_fourGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){ length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==4)]
## Which 4 genes share those variants?
table(apply(sapply(UGT1_shared_variants_fourGenes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,1:9]))]}), 2, toString))
# UGT1A1, UGT1A4, UGT1A6, UGT1A10 
#                              21
table(UGT1_variants_canonical[UGT1_shared_variants_fourGenes, c('VEP_Annotation', 'Location_in_txs')])
#                Location_in_txs
# VEP_Annotation   Intron 4-5
# intron_variant         21


#########################  Tx anno of variants common in 2 UGT1 genes  #########################
UGT1_shared_variants_twoGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT1_shared_variants_twoGenes, function(x){colnames(UGT1_variants_canonical)[which(!is.na(UGT1_variants_canonical[x,1:9]))]}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A5     UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#            392                317                279                340 
table(UGT1_variants_canonical[UGT1_shared_variants_twoGenes, 'Location_in_txs'])
# 5-UTR + First Exon, Intron 1-2          5' upstream, 5-UTR + First Exon         5' upstream, Intron 1-2           Intron 1-2 
#                            253                                      982                              68                    25 
 
# ------------------------------------------------------
## Genes with 5-UTR + First Exon / Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5-UTR + First Exon, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A8 
# 253 

## Verify these variants are in the 1st Exon of UGT1A1 = Intron 1-2 of UGT1A8
pos <- UGT1_variants_canonical[vars, 'Position']
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
# 5-UTR + First Exon 
#               253
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#        253 

## Check 1st Exon of UGT1A1 is within Intron 1-2 of UGT1A8
exon_boundaries <- location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], '5-UTR + First Exon')[[2]]
exon_range <- exon_boundaries[['Start']]:exon_boundaries[['End']]
intron_boundaries <- location_determination(pos[1], canonical_UGT1_txs[['UGT1A8']], 'Intron 1-2')[[2]]
intron_range <- intron_boundaries[['Start']]:intron_boundaries[['End']]

length(which(exon_range %in% intron_range)) == length(exon_range)

# ------------------------------------------------------
## Genes with 5' upstream / 5-UTR + First Exon variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, 5-UTR + First Exon')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A3       UGT1A1, UGT1A5      UGT1A1, UGT1A9 
#            368                 292                  322 

## Verify these variants are in the 1st Exon of UGT1A[3,5,9] = 5' upstream seq of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
# 5' upstream 
#         982 

sapply(c('UGT1A3', 'UGT1A5', 'UGT1A9'), function(y){
  vars_pos <- UGT1_variants_canonical[which(!is.na(UGT1_variants_canonical[,y]) & UGT1_variants_canonical$Location_in_txs=='5\' upstream, 5-UTR + First Exon'), 'Position']
  table(sapply(vars_pos, function(x){location_determination(x, canonical_UGT1_txs[[y]], NULL)[[1]]}))
})
# UGT1A3.5-UTR + First Exon      UGT1A5.5-UTR + First Exon      UGT1A9.5-UTR + First Exon 
#                       368                            292                            322 

# ------------------------------------------------------
## Genes with 5' upstream / Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A3      UGT1A1, UGT1A5      UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#             24                  25                   1                 18

## Confirm these variants are in Intron 1-2 of UGT1A[3,5,8,9] = 5' upstream seq of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
# 5' upstream 
#          68

sapply(c('UGT1A3', 'UGT1A5', 'UGT1A8', 'UGT1A9'), function(y){
  vars_pos <- UGT1_variants_canonical[which(!is.na(UGT1_variants_canonical[,y]) & UGT1_variants_canonical$Location_in_txs=='5\' upstream, Intron 1-2'), 'Position']
  table(sapply(vars_pos, function(x){location_determination(x, canonical_UGT1_txs[[y]], NULL)[[1]]}))
})
# UGT1A3.Intron 1-2    UGT1A5.Intron 1-2    UGT1A8.Intron 1-2    UGT1A9.Intron 1-2 
#                24                   25                    1                   18

## Check the intronic region of UGT1A[3,5,8,9] goes before exon 1 of UGT1A1
max(pos) < location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], '5-UTR + First Exon')[[2]][['Start']]
# [1] TRUE

# ------------------------------------------------------
## Genes with Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == 'Intron 1-2' & rownames(UGT1_variants_canonical) %in% UGT1_shared_variants_twoGenes)
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A8 
#             25 

## Corroborate the common intronic region goes after 1st exon and before 2nd exon of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
min(pos) > location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], '5-UTR + First Exon')[[2]][['End']]
# [1] TRUE
max(pos) < location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], 'Exon 2')[[2]][['Start']]
# [1] TRUE



##############################################
####  Tx location of variants in UGT2 genes
##############################################

for (i in 1:dim(UGT2_variants_canonical)[1]){
  pos <- UGT2_variants_canonical$Position[i]
  txs <- UGT2_variants_canonical[i,which(!is.na(UGT2_variants_canonical[i, 1:10]))]
  locations <- unlist(sapply(txs, function(tx){location_determination(pos, tx, NULL)}))
  if (length(unique(locations))==1){
    UGT2_variants_canonical$Location_in_txs[i] <- unique(locations)
  }
  else
    UGT2_variants_canonical$Location_in_txs[i] <- toString(locations)
}


#########################  All variants  #########################
## Compare VEP annotation and tx location of variants in UGT2 genes
table(UGT2_variants_canonical[, c('VEP_Annotation', 'Location_in_txs')])
#                                                                                    Location_in_txs
# VEP_Annotation            5-UTR + First Exon Exon 2 Exon 3 Exon 4 Exon 5 Intron 1-2 Intron 2-3 Intron 3-4 Intron 4-5 Intron 5-6 Last Exon + 3-UTR
# 3_prime_UTR_variant                        0      0      0      0      0          0          0          0          0          0               155
# 5_prime_UTR_variant                       78      0      0      0      0          0          0          0          0          0                 0
# frameshift_variant                       111     16     16      8     29          0          0          0          0          0                44
# inframe_deletion                           7      0      3      1      2          0          0          1          0          0                 5
# inframe_insertion                          2      0      1      0      0          0          0          0          0          0                 0
# intron_variant                             0      0      0      0      0        326        314        236        270        291                 0
# missense_variant                        1709    284    253    200    521          0          0          0          0          0               600
# splice_acceptor_variant                    0      0      0      0      0          5          4          7          4          5                 0
# splice_donor_variant                       1      0      0      0      0          8          6          7          3         11                 0
# splice_region_variant                      1      8      2      4      8         47         39         46         37         32                 4
# start_lost                                 8      0      0      0      0          0          0          0          0          0                 0
# stop_gained                               80     24     11     18     12          0          0          0          0          0                41
# stop_lost                                  0      0      0      0      0          0          0          0          0          0                13
# stop_retained_variant                      0      0      0      0      0          0          0          0          0          0                 2
# synonymous_variant                       584    102     93     62    178          0          0          0          0          0               195


#########################  Tx anno of variants common in 2 UGT2 genes  #########################
UGT2_shared_variants_twoGenes <- rownames(UGT2_variants_canonical)[which(sapply(1:dim(UGT2_variants_canonical)[1], function(x){ length(which(!is.na(UGT2_variants_canonical[x,1:10]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT2_shared_variants_twoGenes, function(x){colnames(UGT2_variants_canonical)[which(!is.na(UGT2_variants_canonical[x,1:10]))]}), 2, toString))
# UGT2A1, UGT2A2 
#            482
table(UGT2_variants_canonical[UGT2_shared_variants_twoGenes, 'Location_in_txs'])
# Exon 2            Exon 3            Exon 4            Exon 5        Intron 1-2        Intron 2-3        Intron 3-4        Intron 4-5        Intron 5-6     Last Exon + 3-UTR 
#     45                39                28                99                14                34                28                35                38                   122

## Check Intron 1-2 variants span a region between exon 2 of UGT2A1 and UGT2A2 and before exon 1 of UGT2A2
pos <- UGT2_variants_canonical[which(UGT2_variants_canonical$Location_in_txs=='Intron 1-2' & rownames(UGT2_variants_canonical) %in% UGT2_shared_variants_twoGenes), 'Position']
min(pos) > location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], 'Exon 2')[[2]][['End']]
# [1] TRUE
max(pos) < location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], '5-UTR + First Exon')[[2]][['Start']]
# [1] TRUE



##############################################
####  Tx location of variants in UGT3 genes
##############################################

## No shared variants between UGT3A1 and UGT3A2
which(sapply(1:dim(UGT3_variants_canonical)[1], function(x){ length(which(!is.na(UGT3_variants_canonical[x,1:2]))) })==2)
#  integer(0)

for (i in 1:dim(UGT3_variants_canonical)[1]){
  pos <- UGT3_variants_canonical$Position[i]
  txs <- UGT3_variants_canonical[i,which(!is.na(UGT3_variants_canonical[i, 1:2]))]
  locations <- unlist(sapply(txs, function(tx){location_determination(pos, tx, NULL)}))
  if (length(unique(locations))==1){
    UGT3_variants_canonical$Location_in_txs[i] <- unique(locations)
  }
  else
    UGT3_variants_canonical$Location_in_txs[i] <- toString(unlist(sapply(txs, function(tx){location_determination(pos, tx)})))
}

#########################  All variants  #########################
## Compare VEP annotation and tx location of variants in UGT3 genes
table(UGT3_variants_canonical[, c('VEP_Annotation', 'Location_in_txs')])

#                                                                                            Location_in_txs
# VEP_Annotation            5-UTR + First Exon Exon 2 Exon 3 Exon 4 Exon 5 Exon 6 Intron 1-2 Intron 2-3 Intron 3-4 Intron 4-5 Intron 5-6 Intron 6-7 Last Exon + 3-UTR
# 3_prime_UTR_variant                        0      0      0      0      0      0          0          0          0          0          0          0                31
# 5_prime_UTR_variant                       49      0      0      0      0      0          0          0          0          0          0          0                 0
# frameshift_variant                         1      1      4     16      5      5          0          0          0          0          0          0                 6
# inframe_deletion                           0      0      0      3      0      1          0          0          0          0          0          0                 0
# inframe_insertion                          0      0      0      1      0      0          0          0          0          0          0          0                 0
# intron_variant                             0      0      0      0      0      0         96         47         34         66         55         49                 0
# missense_variant                          36     39     35    175     97     94          0          0          0          0          0          0               122
# splice_acceptor_variant                    0      0      0      0      0      0          0          0          1          1          2          1                 0
# splice_donor_variant                       0      0      0      0      0      0          1          5          0          0          2          1                 0
# splice_region_variant                      1      2      1      0      0      0          8          2          8          6          8         11                 0
# start_lost                                 4      0      0      0      0      0          0          0          0          0          0          0                 0
# stop_gained                                0      1      1     14      2      3          0          0          0          0          0          0                 8
# synonymous_variant                        29      9      7     86     34     29          0          0          0          0          0          0                55



##############################################
####  Tx location of variants in UGT8 gene 
##############################################

for (i in 1:dim(UGT8_variants_canonical)[1]){
  pos <- UGT8_variants_canonical$Position[i]
  tx <- UGT8_variants_canonical$Transcript[i]
  location <- location_determination(pos, tx, NULL)[[1]]
  UGT8_variants_canonical[i, 'Location_in_txs'] <- location
}

#########################  All variants  #########################
## Compare VEP annotation and tx location of variants in UGT8 genes
table(UGT8_variants_canonical[, c('VEP_Annotation', 'Location_in_txs')])
#                                                                         Location_in_txs
# VEP_Annotation            Exon 2 Exon 3 Exon 4 Exon 5 Intron 1-2 Intron 2-3 Intron 3-4 Intron 4-5 Intron 5-6 Last Exon + 3-UTR
# 3_prime_UTR_variant            0      0      0      0          0          0          0          0          0                18
# frameshift_variant             1      0      0      3          0          0          0          0          0                 0
# inframe_deletion               1      1      0      0          0          0          0          0          0                 0
# inframe_insertion              1      0      0      0          0          0          0          0          0                 0
# intron_variant                 0      0      0      0          0         28         29         28         28                 0
# missense_variant             120     15      4     24          0          0          0          0          0                57
# splice_acceptor_variant        0      0      0      0          0          0          0          0          1                 0
# splice_donor_variant           0      0      0      0          0          0          1          0          1                 0
# splice_region_variant          1      1      1      2          2          3          2          3          8                 0
# start_lost                     1      0      0      0          0          0          0          0          0                 0
# stop_gained                    2      0      0      0          0          0          0          0          0                 1
# synonymous_variant            77     10      3     15          0          0          0          0          0                32




# __________________________________________________________________________
#  1.1.4  Variant quantification and comparison between genes and families
# __________________________________________________________________________

## Plot total number of variants in each gene/family

# ___________________________________________________________________________
#  1.1.4.1  All variants per gene

## Define colors for variant categories
var_colors <- list('Exon'= 'coral1',
                   'Intron'= 'lightsteelblue3',
                   '5-UTR + First Exon'='plum',
                   'Last Exon + 3-UTR'= 'slateblue4',
                   '5\' upstream'= 'mediumseagreen', 
                   '3\' downstream'='khaki1')

## Function to create barplot for all variants in genes of a certain family
barplot_gene_fam<- function(gene_family, which_variants){
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants_canonical')))
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    ## Variants in canonical tx only
    gene_data <- gene_data[which(gene_data$Variant_ID %in% rownames(UGT_variants)), ]
    gene_data$Location_in_txs <- sapply(gene_data$Position, function(x){
      location_determination(x, unique(gene_data$Transcript), NULL)[[1]]})
    ## Categorize in exonic and intronic 
    gene_data$Location_in_txs <- sapply(gene_data$Location_in_txs, function(x){if(length(grep('Intron', x))==1){'Intron'} else if (length(grep('^Exon', x))==1){'Exon'} else {x}})
    assign(paste0(gene, '_data'), gene_data)
  }
  
  annotations <- c('5-UTR + First Exon', 'Exon', 'Intron', 'Last Exon + 3-UTR',  '5\' upstream', '3\' downstream') 
  var_data <- data.frame(matrix(ncol = 3))
  colnames(var_data) <- c('gene', 'annotation', 'frequency')
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    data <- data.frame(matrix(ncol = 3))
    colnames(data) <- c('gene', 'annotation', 'frequency')
    ylim = 2170
    ylab='Frequency of variants in canonical tx'
    
    if(which_variants=='exonic'){
      annotations <- c('5-UTR + First Exon', 'Exon', 'Last Exon + 3-UTR') 
      ## Exonic variants within each gene
      gene_data <- gene_data[grep('Exon', gene_data$Location_in_txs), ]
      ylim = 750
      ylab='Frequency of exonic variants in canonical tx'
    }
    
    for (i in 1:length(annotations)){
      data[i,'gene'] <- gene
      data[i,'annotation'] <- annotations[i]
      data[i,'frequency'] <- table(gene_data$Location_in_txs)[annotations[i]]
    }
    var_data <- rbind(var_data, data)
  }
  
  var_data$frequency <- replace(var_data$frequency, which(is.na(var_data$frequency)), 0)
  var_data$gene <- factor(var_data$gene, levels = unique(var_data$gene))
  var_data <- var_data[-1,]
  
  p <- ggplot(var_data, aes(fill=annotation, y=frequency, x=gene)) + 
          geom_bar(position="stack", stat="identity") + 
          labs(x=paste(gene_family, 'genes', sep=' '), y=ylab, fill='Predicted effect') +
          theme_bw() + 
          ylim(0, ylim) +
          scale_fill_manual(values = var_colors) +
          theme(axis.text = element_text(size = 8),
                legend.text = element_text(size=9))
  
  return(p)
}


#####################
####   UGT1 genes
#####################
# Stacked barplots
p1 <- barplot_gene_fam('UGT1', 'all') + theme(legend.position="none")

#####################
####   UGT2 genes
#####################
p2 <- barplot_gene_fam('UGT2', 'all') + theme(legend.position="none")

#####################
####   UGT3 genes
#####################
p3 <- barplot_gene_fam('UGT3', 'all') + theme(legend.position="none")

#####################
####  UGT8 gene
#####################
p4 <- barplot_gene_fam('UGT8', 'all')

plot_grid(p1, p2, p3, p4, nrow=1, rel_widths = c(1,1.02,0.28, 0.45))
ggsave(filename=paste0('plots/01_Data_Processing/All_variants_genes.pdf'), width = 19, height = 7)


# ___________________________________________________________________________
#  1.1.4.2  Exonic variants per gene 

#####################
####   UGT1 genes
#####################

p1 <- barplot_gene_fam('UGT1', 'exonic') + theme(legend.position="none")

#####################
####   UGT2 genes
#####################
p2 <- barplot_gene_fam('UGT2', 'exonic') + theme(legend.position="none")

#####################
####   UGT3 genes
#####################
p3 <- barplot_gene_fam('UGT3', 'exonic') + theme(legend.position="none")

#####################
####  UGT8 gene
#####################
p4 <- barplot_gene_fam('UGT8', 'exonic')

plot_grid(p1, p2, p3, p4, nrow=1, rel_widths = c(1,1.02,0.28, 0.45))
ggsave(filename=paste0('plots/01_Data_Processing/Exonic_variants_genes.pdf'), width = 19, height = 7)


# _____________________________________________________________________________
#  1.1.4.3  Exonic variants per gene family

## Map exonic variants to individual genes


## Function to create barplot for exonic variants for each gene family
exonic_vars_gene_fam<- function(gene_family){
  
  UGT_exonic_variants <- eval(parse_expr(paste0(gene_family, '_exonic_variants')))
  
  var_data <- data.frame(matrix(ncol = 3))
  colnames(var_data) <- c('gene_family', 'annotation', 'frequency')
  
  for (gene in genes){
    data <- data.frame(matrix(ncol = 3))
    colnames(data) <- c('gene', 'annotation', 'frequency')
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    gene_data <- gene_data[which(gene_data$Variant_ID %in% rownames(UGT_variants_canonical)), ]
    for (i in 1:length(annotations)){
      data[i,'gene'] <- gene
      data[i,'annotation'] <- annotations[i]
      data[i,'frequency'] <- table(gene_data$VEP_Annotation)[annotations[i]]
    }
    var_data <- rbind(var_data, data)
  }
  var_data$frequency <- replace(var_data$frequency, which(is.na(var_data$frequency)), 0)
  var_data$gene <- factor(var_data$gene, levels = unique(var_data$gene))
  var_data <- var_data[-1,]
  
  p <- ggplot(var_data, aes(fill=annotation, y=frequency, x=gene)) + 
    geom_bar(position="stack", stat="identity") + 
    labs(x=paste(gene_family, 'genes', sep=' '), y='Frequency of variants in canonical tx', fill='Predicted effect') +
    theme_bw() + 
    scale_fill_manual(values = var_colors) +
    theme(axis.text = element_text(size = 8),
          legend.text = element_text(size=9))
  
  return(p)
}




# _________________________________________________________________________
#  1.1.4.1 Shared variants per gene family














#################################################################################
#################################################################################
#################################################################################

## TODO'S:

# - Compare my anno with VEP
# - Identify weird cases and plot


# - Shared variants anno: compare my anno with VEP
# - What about shared variants with different VEP anno or manual locations? -> plot



#################################################################################
#################################################################################
#################################################################################



## FUTURE

# 
# for (gene in UGT_genes){
#   
#   ## Read data
#   UGT_gene_info <- read.csv(paste0("Processed-data/", gene, "_processed.csv"))
#   
#   ## Add alternate Allele Frequency within each population
#   populations <- c('Other', 'Latino_or_Admixed_American', 'European_Finnish', 'Amish', 'East_Asian',
#                    'Middle_Eastern', 'African_or_African_American', 'South_Asian', 'Ashkenazi_Jewish',
#                    'European_non_Finnish')
#   
#   for(population in populations){
#     ## Frequency given by the variant counts in the population over the number of alleles in the same population 
#     UGT_gene_info[, paste0('Allele_Frequency_', population)] <- eval(parse_expr(paste0('UGT_gene_info$Allele_Count_', population)))/eval(parse_expr(paste0('UGT_gene_info$Allele_Number_', population)))
#     
#   }
#   
#   ## Add overall Allele Frequency for each population
#   for(population in populations){
#     ## Frequency given by the variant counts over the total number of alleles (considering all populations)
#     UGT_gene_info[, paste0('Overall_Allele_Frequency_', population)] <- eval(parse_expr(paste0('UGT_gene_info$Allele_Count_', population)))/UGT_gene_info$Allele_Number
#     
#   }
#   
# }
