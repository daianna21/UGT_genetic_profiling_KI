

library(here)
library(readr)
library(rlang)

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

##### *** UGT1A1 and UGT1A8 txs are swapped; take only the main tx of these genes for all of their variants ***

UGT1A1_data$Transcript <- rep('ENST00000305208.5', dim(UGT1A1_data)[1])
UGT1A8_data$Transcript <- rep('ENST00000373450.4', dim(UGT1A8_data)[1])


## Unique variants in UGT1 genes
unique_UGT1_variants <- unique(c(UGT1A1_data$Variant_ID, UGT1A3_data$Variant_ID, UGT1A4_data$Variant_ID, UGT1A5_data$Variant_ID, UGT1A6_data$Variant_ID, UGT1A7_data$Variant_ID, UGT1A8_data$Variant_ID, UGT1A9_data$Variant_ID, UGT1A10_data$Variant_ID))

## Compute in what transcript each variant appears for each gene 
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




# __________________________________________________________________________________________
#  1.1.2 Verify there are no overlapping variants between genes from different subfamilies
# __________________________________________________________________________________________

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




# _______________________________________________________
#  1.1.3 Examine annotation and genes of shared variants
# _______________________________________________________

##############################
####  UGT1 shared variants 
##############################

###########  Annotation for shared variants between all UGT1 genes  ###########
table(UGT1_variants[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,-c(10,11)]))) })==9), 'VEP_Annotation'])
# 3_prime_UTR_variant      frameshift_variant        inframe_deletion          intron_variant        missense_variant splice_acceptor_variant    splice_donor_variant 
#                  10                       8                       2                      94                     165                       2                       1 
# splice_region_variant             stop_gained      synonymous_variant 
#                    10                       6                      76 


###########  Annotation for variants shared between 4 UGT1 genes  ###########
table(UGT1_variants[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,-c(10,11)]))) })==4), 'VEP_Annotation'])
# 3_prime_UTR_variant, intron_variant, 3_prime_UTR_variant, 3_prime_UTR_variant 
#                                                                            16 
# intron_variant 
#             21 
# missense_variant, intron_variant, missense_variant, missense_variant 
#                                                                    2 
# splice_acceptor_variant, intron_variant, splice_acceptor_variant, splice_acceptor_variant 
#                                                                                         2 
# stop_retained_variant, intron_variant, stop_retained_variant, stop_retained_variant 
#                                                                                   1 
# synonymous_variant, intron_variant, synonymous_variant, synonymous_variant 
#                                                                         3 

## Which 4 genes share those variants?
UGT1_shared_variants_fourGenes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,-c(10,11)]))) })==4)]
table(apply(sapply(UGT1_shared_variants_fourGenes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,-c(10,11)]))]}), 2, toString))
# UGT1A1, UGT1A4, UGT1A6, UGT1A10 
#                              45 
## Txs
table(apply(sapply(UGT1_shared_variants_fourGenes, function(x){UGT1_variants[x,which(!is.na(UGT1_variants[x, -c(10,11)]))]}), 2, toString))
# ENST00000305208.5, ENST00000373409.3, ENST00000305139.6, ENST00000344644.5        ENST00000305208.5, ENST00000373409.3, ENST00000406651.1, ENST00000373445.1 
#                                                                         21                                                                                24
 

###########  Annotation for variants shared between 2 UGT1 genes  ###########
table(UGT1_variants[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,-c(10,11)]))) })==2), 'VEP_Annotation'])
# 5_prime_UTR_variant                  frameshift_variant                    inframe_deletion                      intron_variant 
#                   7                                  55                                   7                                  75 
# intron_variant, 5_prime_UTR_variant                    missense_variant                splice_donor_variant               splice_region_variant 
#                                 24                                 773                                   4                                  14 
# start_lost                         stop_gained                  synonymous_variant 
#          6                                  32                                 331 

## Which 2 genes share those variants?
UGT1_shared_variants_twoGenes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,-c(10,11)]))) })==2)]
table(apply(sapply(UGT1_shared_variants_twoGenes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,-c(10,11)]))]}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A5     UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#            392                317                279                340 

## Txs
table(apply(sapply(UGT1_shared_variants_twoGenes, function(x){UGT1_variants[x,which(!is.na(UGT1_variants[x, -c(10,11)]))]}), 2, toString))
# ENST00000305208.5, ENST00000354728.4      ENST00000305208.5, ENST00000373414.3      ENST00000305208.5, ENST00000373450.4      ENST00000305208.5, ENST00000482026.1 
#                                  340                                       317                                       279                                       392 



##############################
####  UGT2 shared variants 
##############################

###########  Annotation for variants shared between 2 UGT2 genes  ###########
table(UGT2_variants[which(sapply(1:dim(UGT2_variants)[1], function(x){ length(which(!is.na(UGT2_variants[x,-c(11,12)]))) })==2), 'VEP_Annotation'])
# 3_prime_UTR_variant      frameshift_variant        inframe_deletion          intron_variant        missense_variant splice_acceptor_variant    splice_donor_variant 
#                  10                      18                       1                     122                     217                       3                       3 
# splice_region_variant             stop_gained   stop_retained_variant      synonymous_variant 
#                    23                      15                       1                      69 

## Which 2 genes share those variants?
UGT2_shared_variants_twoGenes <- rownames(UGT2_variants)[which(sapply(1:dim(UGT2_variants)[1], function(x){ length(which(!is.na(UGT2_variants[x,-c(11,12)]))) })==2)]
table(apply(sapply(UGT2_shared_variants_twoGenes, function(x){colnames(UGT2_variants)[which(!is.na(UGT2_variants[x,-c(11,12)]))]}), 2, toString))
# UGT2A1, UGT2A2 
#            482



##############################
####  UGT3 shared variants 
##############################

## No shared variants between UGT3A1 and UGT3A2
which(sapply(1:dim(UGT3_variants)[1], function(x){ length(which(!is.na(UGT3_variants[x,-c(3,4)]))) })==2)
#  integer(0)



##############################
####  UGT8 shared variants 
##############################

## Only one gene for UGT8




# ________________________________________________________________________________
#  1.1.4 Re-evaluate annotation of shared variants based on genes' txs boundaries
# ________________________________________________________________________________

## For a given shared variant, examine where it's located with respect to genes' txs boundaries

########################
####  UGT1 genes' txs
########################

## Function to evaluate the location of a variant according to the boundaries of a gene transcript
location_determination <- function(variant_pos, gene_tx){
  
  ## Read tx data
  tx_seq_data <- as.data.frame(read.csv(paste0("raw-data/Tx_seq_data/", tx, "_seq_data.csv")))
  tx_seq_data <- tx_seq_data[-c(1, dim(tx_seq_data)[1]),]
  ## Column names
  colnames(tx_seq_data)[-1] <- gsub('\\.+', '_', colnames(tx_seq_data)[-1])
  
  ## Transcript composition: 
  ## (5' upstream seq) ... [5'UTR + Exon 1] + Intron 1-2 + Exon 2 + Intron 2-3 + Exon 3 + Intron 3-4 + Exon 4 + Intron 4-5 + [Exon 5 + 3'UTR] ... (3' downstream seq)
  
  ## Add feature variable
  rownames(tx_seq_data)<- sapply(1:dim(tx_seq_data)[1], function(x){
    if(length(grep('Intron', tx_seq_data$Exon_Intron[x]))==1){tx_seq_data$Exon_Intron[x]}
    else{paste0('Exon-', tx_seq_data$No.[x])}
  })
  ## Specify Exon-1 + 5'UTR and Exon-5 + 3'UTR
  rownames(tx_seq_data)[1] <- '5-UTR + Exon-1'
  rownames(tx_seq_data)[dim(tx_seq_data)[1]] <- 'Exon-5 + 3-UTR'
  
  ## Char to integer for End, Start and Length
  tx_seq_data$Start <- as.numeric(gsub(',', '', tx_seq_data$Start))
  tx_seq_data$End <- as.numeric(gsub(',', '', tx_seq_data$End))
  tx_seq_data$Length <- as.numeric(gsub(',', '', tx_seq_data$Length))
  
  
  ## Evaluate if the variant is within a tx exon/intron/UTR region
  location_within_tx <- unlist(sapply(rownames(tx_seq_data), function(x){
                                  if(variant_pos %in% tx_seq_data[x, 'Start']:tx_seq_data[x, 'End']){x}}))
  
  ## If the variant is not within tx region then is either in 5' or 3' upstream sequences
  if (is.null(location_within_tx))  {
    if (variant_pos < tx_seq_data["5-UTR + Exon-1", 'Start']){
      location <- '5\' upstream'
    }
    
    if (variant_pos > tx_seq_data["Exon-5 + 3-UTR", 'End']){
      location <- '3\' downstream'
    }
  }
  
  else {
    location <- location_within_tx
  }
  
  return(location)
}



## Shared variants by 4 genes
UGT1_shared_variants_fourGenes <- UGT1_variants_counts[rownames(UGT1_variants_counts)[which(apply(UGT1_variants_counts[,-10], 1, sum)==4)],]
## Genes that share the variants
apply(sapply(UGT1_shared_variants_fourGenes, function(x){colnames(UGT1_variants_counts)[which(UGT1_variants_counts[x,]==1)]}), 2, toString)


## Add location of each variant for the genes in which it's present
for (i in 1:dim(UGT1_variants_counts)[1]){
  ## Genes with the variant
  genes_tx <- colnames(UGT1_variants_counts)[which(UGT1_variants_counts[i,]==1)]
  ## Variant location within each gene
  sapply(genes, function(x){location_determination(UGT1_variants_counts[i,'Position'], x)})
}





























for (gene in UGT_genes){
  
  ## Read data
  UGT_gene_info <- read.csv(paste0("Processed-data/", gene, "_processed.csv"))
  
  ## Add alternate Allele Frequency within each population
  populations <- c('Other', 'Latino_or_Admixed_American', 'European_Finnish', 'Amish', 'East_Asian',
                   'Middle_Eastern', 'African_or_African_American', 'South_Asian', 'Ashkenazi_Jewish',
                   'European_non_Finnish')
  
  for(population in populations){
    ## Frequency given by the variant counts in the population over the number of alleles in the same population 
    UGT_gene_info[, paste0('Allele_Frequency_', population)] <- eval(parse_expr(paste0('UGT_gene_info$Allele_Count_', population)))/eval(parse_expr(paste0('UGT_gene_info$Allele_Number_', population)))
    
  }
  
  ## Add overall Allele Frequency for each population
  for(population in populations){
    ## Frequency given by the variant counts over the total number of alleles (considering all populations)
    UGT_gene_info[, paste0('Overall_Allele_Frequency_', population)] <- eval(parse_expr(paste0('UGT_gene_info$Allele_Count_', population)))/UGT_gene_info$Allele_Number
    
  }
  
}
