

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

## For each variant, determine in what UGT genes' txs they are present 

create_gene_fam_table <- function(gene_family){
  
  ## Genes of gene family
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  
  ## Unique variants in gene family
  unique_UGT_variants <- unique(unlist(sapply(paste0(genes, '_data$Variant_ID'), function(x){eval(parse_expr(x))})))
  
  ## Compute in which transcript each variant appears for each gene 
  ## (if present in a gene, variants are expected to appear only in one transcript since they are all unique in each gene dataset)
  UGT_variants <- data.frame(matrix(ncol = length(genes)+2, nrow = 0))
  colnames(UGT_variants) <- c(genes, 'VEP_Annotation', 'Position')
  
  i=1
  for (variant in unique_UGT_variants){
    
    ## Search variant in each gene; if present, annotate tx
    variant_txs <- sapply(genes, function(gene){
      tx <- eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'Transcript']
      if (length(tx)==0) {NA}
      else {tx}
    })
    
    UGT_variants[i, 1:length(genes)] <- variant_txs
    
    ## VEP annotation for variant in each gene that contains it
    variant_anno <- unlist(sapply(genes, function(gene){
      eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
    }))
    ## Evaluate if the anno of a shared variable is the same across all genes
    if (length(unique(variant_anno))==1){
      UGT_variants$VEP_Annotation[i] <- unique(variant_anno)
    }
    else{
      UGT_variants$VEP_Annotation[i] <- toString(variant_anno)
    }
    
    ## Add variant position
    UGT_variants$Position[i] <- as.numeric(strsplit(variant, '-')[[1]][2])
    
    i=i+1
  }
  
  rownames(UGT_variants) <- unique_UGT_variants
  
  assign(paste0('unique_', gene_family, '_variants'), unique_UGT_variants, envir = parent.frame())
  assign(paste0(gene_family, '_variants'), UGT_variants, envir = parent.frame())
}


############################
####      UGT1 genes 
############################

UGT1_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10')

## Confirm all variants in UGT1 genes are in the same chr
sapply(UGT1_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT1A1  UGT1A3  UGT1A4  UGT1A5  UGT1A6  UGT1A7  UGT1A8  UGT1A9 UGT1A10 
#      2       2       2       2       2       2       2       2       2 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## !! UGT1A1 and UGT1A8 txs are swapped; take only the unique / canonical tx of these 
##                     genes for all of their variants !!

UGT1A1_data$Transcript <- rep('ENST00000305208.5', dim(UGT1A1_data)[1])
UGT1A8_data$Transcript <- rep('ENST00000373450.4', dim(UGT1A8_data)[1])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

create_gene_fam_table('UGT1')


############################
####    UGT2 genes 
############################
UGT2_genes <- c('UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28')

## Confirm all UGT2 variants are in the same chr
sapply(UGT2_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT2A1  UGT2A2  UGT2A3  UGT2B4  UGT2B7 UGT2B10 UGT2B11 UGT2B15 UGT2B17 UGT2B28 
#      4       4       4       4       4       4       4       4       4       4 

create_gene_fam_table('UGT2')


############################
####    UGT3 genes 
############################
UGT3_genes <- c('UGT3A1', 'UGT3A2')

## Confirm all UGT3 variants are in the same chr
sapply(UGT3_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT3A1 UGT3A2 
#      5      5 

create_gene_fam_table('UGT3')


############################
####    UGT8 genes 
############################

UGT8_genes <- c('UGT8')
unique(UGT8_data$Chromosome)
#    4 

## Variant_ID's are unique so they all appear once in UGT8
which(duplicated(UGT8_data$Variant_ID))
#  integer(0)
create_gene_fam_table('UGT8')
unique_UGT8_variants <- UGT8_data$Variant_ID

  
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

for (gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')) {
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants')))
  canonical_UGT_txs <- eval(parse_expr(paste0('canonical_', gene_family, '_txs')))
  
  for(gene in gene_family){
    
    ## Genes of gene family
    genes <- eval(parse_expr(paste0(gene_family, '_genes')))
    
    ## For each gene
    print(sapply(genes, function(x){
      length(which(UGT_variants[[x]] == canonical_UGT_txs[[x]]))/length(which(!is.na(UGT_variants[[x]]))) * 100}
    ))
    
    ## In the whole family
    ## (For shared variants, if they are in canonical txs of all genes in which they are)
    for (i in 1:dim(UGT_variants)[1]){
      var_txs <- UGT_variants[i, which(!is.na(UGT_variants[i, 1:length(genes)]))] 
      ## If all txs of a variant are canonical
      if (length(which(var_txs %in% unlist(canonical_UGT_txs))) == length(var_txs)){
        UGT_variants$Canonical_txs[i] <- TRUE
      }
      else{
        UGT_variants$Canonical_txs[i] <- FALSE
      }
    }
    
    assign(paste0(gene, '_variants'), UGT_variants)
    print(table(UGT_variants$Canonical_txs)['TRUE'] / dim(UGT_variants)[1] * 100)
  }
}

###################
####  UGT1 genes
###################

## % of variants in each gene that are in the canonical tx

# UGT1A1    UGT1A3    UGT1A4    UGT1A5    UGT1A6    UGT1A7    UGT1A8    UGT1A9   UGT1A10 
# 100.00000 100.00000 100.00000 100.00000  96.03524  95.40682 100.00000 100.00000  97.13262 

## % of total UGT1 variants that appear in canonical txs
# 98.31384 


###################
####  UGT2 genes
###################

# UGT2A1    UGT2A2    UGT2A3    UGT2B4    UGT2B7   UGT2B10   UGT2B11   UGT2B15   UGT2B17   UGT2B28 
# 97.41268  99.86339  98.91041  94.64752  99.47507  99.88208 100.00000 100.00000 100.00000 100.00000  

## % of total UGT2 variants in canonical txs
# 98.95475


###################
####  UGT3 genes
###################

# UGT3A1   UGT3A2 
# 93.36735 99.56012 

## % of total UGT3 variants in canonical txs
# 96.24829 


###################
####  UGT8 gene
###################

# UGT8 
# 97.22222 


######################
####  All UGT genes
######################

gene_families <- c('UGT1', 'UGT2', 'UGT3', 'UGT8')

## Total number of UGT variants in canonical txs
vars_in_canonical_txs <- sum(unlist(sapply(paste0('table(', gene_families, '_variants$Canonical_txs)[\'TRUE\']'), function(x){eval(parse_expr(x))})))
## Total variants in all UGT families
total_vars <- sum(unlist(sapply(paste0('dim(', gene_families, '_variants)[1]'), function(x){eval(parse_expr(x))})))

## % of total UGT variants in canonical txs
vars_in_canonical_txs / total_vars * 100
# [1] 98.39432




# _____________________________________________________________________________
#  1.1.3 Manual annotation of variants based on boundaries of canonical txs
# _____________________________________________________________________________

## Conserve variants in canonical txs only
UGT1_variants_canonical <- UGT1_variants[UGT1_variants$Canonical_txs==TRUE, ]
UGT2_variants_canonical <- UGT2_variants[UGT2_variants$Canonical_txs==TRUE, ]
UGT3_variants_canonical <- UGT3_variants[UGT3_variants$Canonical_txs==TRUE, ]
UGT8_variants_canonical <- UGT8_variants[UGT8_variants$Canonical_txs==TRUE, ]

## Define end of 5'-UTR of each tx 
fiveUTR_end <- list('ENST00000305208.5'= 234668933,
                    'ENST00000482026.1'= 234637772,
                    'ENST00000373409.3'= 234627466, 
                    'ENST00000373414.3'= 234621638,
                    'ENST00000305139.6'= 234601650,
                    'ENST00000373426.3'= 234590584,
                    'ENST00000373450.4'= 234526353, 
                    'ENST00000354728.4'= 234580580,
                    'ENST00000344644.5'= 234545168,
                    'ENST00000503640.1'= 70513363, 
                    'ENST00000457664.2'= 70505358, 
                    'ENST00000251566.4'= 69817479, 
                    'ENST00000305107.6'= 70361580, 
                    'ENST00000305231.7'= 69962238, 
                    'ENST00000265403.7'= 69681737, 
                    'ENST00000446444.1'= 70080441, 
                    'ENST00000338206.5'= 69536337, 
                    'ENST00000317746.2'= 69434203, 
                    'ENST00000335568.5'= 70146218, 
                    'ENST00000274278.3'= 35991343, 
                    'ENST00000282507.3'= 36066892, 
                    'ENST00000310836.6'= 115520130)

## Define start of 3'-UTR of each tx (same for all UGT1 genes)
threeUTR_start <- list('ENST00000305208.5'= 234681206,
                       'ENST00000482026.1'= 234681206,
                       'ENST00000373409.3'= 234681206,
                       'ENST00000373414.3'= 234681206,
                       'ENST00000305139.6'= 234681206,
                       'ENST00000373426.3'= 234681206,
                       'ENST00000373450.4'= 234681206, 
                       'ENST00000354728.4'= 234681206,
                       'ENST00000344644.5'= 234681206,
                       'ENST00000503640.1' = 70455089, 
                       'ENST00000457664.2'= 70455089, 
                       'ENST00000251566.4'= 69795530, 
                       'ENST00000305107.6'= 70346351, 
                       'ENST00000305231.7'= 69978455, 
                       'ENST00000265403.7'= 69696598, 
                       'ENST00000446444.1'= 70066157, 
                       'ENST00000338206.5'= 69512821, 
                       'ENST00000317746.2'= 69403342, 
                       'ENST00000335568.5'= 70160528, 
                       'ENST00000274278.3'= 35954303, 
                       'ENST00000282507.3'= 36035799,
                       'ENST00000310836.6'= 115597445)

## For a given variant, examine where it's located with respect to txs boundaries

## Function to evaluate the location of a variant according to the boundaries of a gene transcript
location_determination <- function(variant_pos, tx, feature){

  ## Read tx data
  tx_seq_data <- as.data.frame(read.csv(paste0("raw-data/Tx_seq_data/", tx, "_seq_data.csv")))
  tx_seq_data <- tx_seq_data[-c(1, dim(tx_seq_data)[1]),]
  ## Column names
  colnames(tx_seq_data)[-1] <- gsub('\\.+', '_', colnames(tx_seq_data)[-1])
  ## Char to integer for End and Start 
  tx_seq_data$Start <- as.numeric(gsub(',', '', tx_seq_data$Start))
  tx_seq_data$End <- as.numeric(gsub(',', '', tx_seq_data$End))
  tx_seq_data <- tx_seq_data[, c('No.', 'Exon_Intron', 'Start', 'End')]
  
  
  ## Overall transcript composition:  
  ## (5' upstream seq) ... (5'UTR) + Exon A + Intron A-B + Exon B + Intron B-C + Exon C + Intron C-D + Exon D + ... + Exon Z + (3'UTR) ... (3' downstream seq)
  colnames(tx_seq_data)[2] <- 'Location'
  
  ## Delimit UTR regions
  ## 5'-UTR = first position of tx to end of 5'UTR 
  fiveUTR <- tx_seq_data[1,'Start']:fiveUTR_end[[tx]]
  ## If UTR is not defined (length 1 from start to end)
  if (length(fiveUTR)==1){fiveUTR=0}
  ## 3'-UTR = start of 3'-UTR to end of tx
  threeUTR <- threeUTR_start[[tx]]:tx_seq_data[nrow(tx_seq_data), 'End']
  if (length(threeUTR)==1){threeUTR=0}
  
  ## Delimit real boundaries of first exon
  first_exon <- tx_seq_data[1, 'Start']:tx_seq_data[1, 'End']
  ## Remove 5'-UTR region
  if(!length(which(first_exon %in% fiveUTR))==0){
    first_exon <- first_exon[-which(first_exon %in% fiveUTR)]
  }
  ## Delimit real boundaries of last exon
  last_exon <- tx_seq_data[nrow(tx_seq_data), 'Start']:tx_seq_data[nrow(tx_seq_data), 'End']
  ## Remove 3'-UTR region
  last_exon <- last_exon[-which(last_exon %in% threeUTR)]
  
  ## Replace exon 1 limits
  if (length(first_exon)==0){first_exon <- c(0,0)}
  tx_seq_data[1,'Start'] <- first_exon[1]
  tx_seq_data[1,'End'] <- first_exon[length(first_exon)]
  
  ## Replace last exon limits
  if (length(last_exon)==0){last_exon <- c(0,0)}
  tx_seq_data[nrow(tx_seq_data),'Start'] <- last_exon[1]
  tx_seq_data[nrow(tx_seq_data),'End'] <- last_exon[length(last_exon)]
  
  ## Add UTRs
  tx_seq_data <- rbind(c(NA, '5\'UTR', fiveUTR[1], fiveUTR[length(fiveUTR)]), tx_seq_data)
  tx_seq_data <- rbind(tx_seq_data, c(NA, '3\'UTR', threeUTR[1], threeUTR[length(threeUTR)]))
  
  ## Location categories
  tx_seq_data$Location <- sapply(1:nrow(tx_seq_data), 
                                 function(i){if(length(grep('ENSE', tx_seq_data$Location[i]))==1){paste0('Exon ', tx_seq_data$No.[i])}
                                             else{tx_seq_data$Location[i]}})
  rownames(tx_seq_data) <- tx_seq_data$Location
  
  ## Evaluate if the variant is within a tx exon/intron/UTR or outside
  location_within_tx <- unlist(sapply(rownames(tx_seq_data), function(x){
                                  if(variant_pos %in% tx_seq_data[x, 'Start']:tx_seq_data[x, 'End']){x}}))
  
  ## If the variant is not within tx region then is either in the 5' upstream or 3' downstream sequence
  if (is.null(location_within_tx))  {
    
    ## If no 5'-UTR present, take start of first exon
    if (! as.numeric(tx_seq_data["5\'UTR", 'Start'])==0) {start <- as.numeric(tx_seq_data["5\'UTR", 'Start'])}
    else {start <- as.numeric(tx_seq_data["Exon 1", 'Start'])}
    
    ## If tx goes from 5'-3'
    if (tx_seq_data$Start[3]<tx_seq_data$End[3]){
      
      if (variant_pos < start){
        location <- '5\' upstream'
      }
      
      if (variant_pos > as.numeric(tx_seq_data["3\'UTR", 'End'])){
        location <- '3\' downstream'
      }
    }
    
  ## If tx goes from 3'-5'
   else{
     if (variant_pos > start){
       location <- '5\' upstream'
     }
     
     if (variant_pos < as.numeric(tx_seq_data["3\'UTR", 'End'])){
       location <- '3\' downstream'
     }
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


## Function to add variant location data in all genes of a family
add_tx_location <- function(gene_family){
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants_canonical')))
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    ## Variants in canonical tx only
    gene_data <- gene_data[which(gene_data$Variant_ID %in% rownames(UGT_variants)), ]
    ## Location
    gene_data$Location_in_txs <- sapply(gene_data$Position, function(x){
      location_determination(x, unique(gene_data$Transcript), NULL)[[1]]})
    assign(paste0(gene, '_data'), gene_data, envir = parent.frame())
  }
}

## Add tx location for variants in each gene dataset
add_tx_location('UGT1')
add_tx_location('UGT2')
add_tx_location('UGT3')
add_tx_location('UGT8')


## Add tx location info for variants in each gene family
for (gene_family in gene_families){
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants_canonical')))
  
  for (i in 1:dim(UGT_variants)[1]){
    pos <- UGT_variants$Position[i]
    txs <- UGT_variants[i,which(!is.na(UGT_variants[i, 1:length(genes)]))]
    locations <- unlist(sapply(txs, function(tx){location_determination(pos, tx, NULL)}))
    if (length(unique(locations))==1){
      UGT_variants$Location_in_txs[i] <- unique(locations)
    }
    else
      UGT_variants$Location_in_txs[i] <- toString(locations)
  }
  assign(paste0(gene_family, '_variants_canonical'), UGT_variants)
}


# __________________________________________________________________________________
#  1.1.3.1  Check manual annotation of all and shared variants in each gene family

######################
####  UGT1 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT1_variants_canonical[, 'Location_in_txs'])
# 3'UTR             5' upstream      5' upstream, 5'UTR     5' upstream, Exon 1       5' upstream, Intron 1-2           5'UTR 
#    10                     443                      19                     963                            68              60 
# 5'UTR, Intron 1-2           Exon 1      Exon 1, Intron 1-2             Exon 2             Exon 3              Exon 4 
#                11             1291                     242                 25                 45                  66 
# Exon 5              Intron 1-2              Intron 2-3              Intron 3-4              Intron 4-5 
#    122                     135                      32                      29                      54 


################  Tx anno of variants common in all UGT1 genes  ################ 
UGT1_shared_variants_allGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==9)]

table(UGT1_variants_canonical[UGT1_shared_variants_allGenes, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5 
#     10         25         45         66        122           12           32           29           33 

# ------------------------------------------------------
## Position of shared variants in Intron 1-2 region
pos <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs=='Intron 1-2' & rownames(UGT1_variants_canonical) %in% UGT1_shared_variants_allGenes), 'Position']
## Confirm this overlapping intronic region starts after Exon 1 of UGT1A1 and finishes before exon 2 of all UGT1 genes
min(pos) > location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], 'Exon 1')[[2]]['End']
# TRUE
max(pos) < location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], 'Exon 2')[[2]]['Start']
# TRUE

# ------------------------------------------------------
## Check all variants in Exon 2 - Exon 5 region are shared by all UGT1 genes
exon2to5_vars <- rownames(UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5')), ])
which(sapply(exon2to5_vars, function(x){length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })!=9)
# named integer(0)


#################  Tx anno of variants common in 4 UGT1 genes  #################
UGT1_shared_variants_fourGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){ length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==4)]
## Which 4 genes share those variants?
table(apply(sapply(UGT1_shared_variants_fourGenes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,1:9]))]}), 2, toString))
# UGT1A1, UGT1A4, UGT1A6, UGT1A10 
#                              21
table(UGT1_variants_canonical[UGT1_shared_variants_fourGenes, 'Location_in_txs'])
# Intron 4-5 
#         21 


#################  Tx anno of variants common in 2 UGT1 genes  #################
UGT1_shared_variants_twoGenes <- rownames(UGT1_variants_canonical)[which(sapply(1:dim(UGT1_variants_canonical)[1], function(x){length(which(!is.na(UGT1_variants_canonical[x,1:9]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT1_shared_variants_twoGenes, function(x){colnames(UGT1_variants_canonical)[which(!is.na(UGT1_variants_canonical[x,1:9]))]}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A5     UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#            392                317                279                340 
table(UGT1_variants_canonical[UGT1_shared_variants_twoGenes, 'Location_in_txs'])
# 5' upstream, 5'UTR     5' upstream, Exon 1      5' upstream, Intron 1-2       5'UTR, Intron 1-2      Exon 1, Intron 1-2         Intron 1-2 
#                 19                     963                           68                      11                     242                 25 

# ------------------------------------------------------------------------------
##  1.  Genes with 5' upstream / 5'-UTR variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, 5\'UTR')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A9 
#              8                 11 

## Verify these variants are in the 5'-UTR of UGT1A[3,9] = 5' upstream seq of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
pos_A3 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, 5\'UTR' & !is.na(UGT1_variants_canonical[['UGT1A3']])), 'Position']
pos_A9 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, 5\'UTR' & !is.na(UGT1_variants_canonical[['UGT1A9']])), 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
# 5' upstream 
#          19
table(sapply(pos_A3, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A3']], NULL)[[1]]}))
#  5'UTR 
#      8 
table(sapply(pos_A9, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A9']], NULL)[[1]]}))
#  5'UTR 
#     11 

# ------------------------------------------------------------------------------
## 2.  Genes with 5' upstream / Exon 1 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Exon 1')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A3    UGT1A1, UGT1A5     UGT1A1, UGT1A9 
#           360                292                311 

## Corroborate these variants are in Exon 1 of UGT1A[3,5,9] = 5' upstream seq of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
pos_A3 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Exon 1' & !is.na(UGT1_variants_canonical[['UGT1A3']])), 'Position']
pos_A5 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Exon 1' & !is.na(UGT1_variants_canonical[['UGT1A5']])), 'Position']
pos_A9 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Exon 1' & !is.na(UGT1_variants_canonical[['UGT1A9']])), 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
#  5' upstream 
#          963
table(sapply(pos_A3, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A3']], NULL)[[1]]}))
#  Exon 1 
#     360 
table(sapply(pos_A5, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A5']], NULL)[[1]]}))
#  Exon 1 
#     292
table(sapply(pos_A9, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A9']], NULL)[[1]]}))
#  Exon 1
#     311 

# ------------------------------------------------------------------------------
## 3.  Genes with 5' upstream / Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A5     UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#             24                 25                  1                 18 

## Check these variants are in Intron 1-2 of UGT1A[3,5,8,9] = 5' upstream seq of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
pos_A3 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A3']])), 'Position']
pos_A5 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A5']])), 'Position']
pos_A8 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A8']])), 'Position']
pos_A9 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\' upstream, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A9']])), 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
#  5' upstream 
#           68
table(sapply(pos_A3, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A3']], NULL)[[1]]}))
# Intron 1-2 
#         24 
table(sapply(pos_A5, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A5']], NULL)[[1]]}))
# Intron 1-2 
#         25
table(sapply(pos_A8, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#          1
table(sapply(pos_A9, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A9']], NULL)[[1]]}))
# Intron 1-2 
#         18 

# -------------------------------------------------------------------------------
## 4.  Genes with 5'-UTR / Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == '5\'UTR, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A8 
#             11

## Check these variants are in 5'-UTR of UGT1A1 = Intron 1-2 of UGT1A8
pos <- UGT1_variants_canonical[vars, 'Position']
pos_A8 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == '5\'UTR, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A8']])), 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
# 5'UTR 
#    11
table(sapply(pos_A8, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#         11

# ------------------------------------------------------------------------------
## 5.  Genes with Exon 1 / Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == 'Exon 1, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A8 
#            242

## Verify these variants are in Exon 1 of UGT1A1 = Intron 1-2 of UGT1A8
pos <- UGT1_variants_canonical[vars, 'Position']
pos_A8 <- UGT1_variants_canonical[which(UGT1_variants_canonical$Location_in_txs == 'Exon 1, Intron 1-2' & !is.na(UGT1_variants_canonical[['UGT1A8']])), 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]}))
#  Exon 1 
#     242 
table(sapply(pos_A8, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#        242

# ------------------------------------------------------
## 6.  Genes with Intron 1-2 variants
vars <- which(UGT1_variants_canonical$Location_in_txs == 'Intron 1-2' & rownames(UGT1_variants_canonical) %in% UGT1_shared_variants_twoGenes)
table(apply(sapply(vars, function(x){colnames(UGT1_variants_canonical[x, which(!is.na(UGT1_variants_canonical[x, 1:9]))])}), 2, toString))
# UGT1A1, UGT1A8 
#             25 

## Corroborate the common intronic region goes after 1st exon and before 2nd exon of UGT1A1
pos <- UGT1_variants_canonical[vars, 'Position']
min(pos) > location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], 'Exon 1')[[2]][['End']]
# [1] TRUE
max(pos) < location_determination(pos[1], canonical_UGT1_txs[['UGT1A1']], 'Exon 2')[[2]][['Start']]
# [1] TRUE



######################
####  UGT2 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT2_variants_canonical[, 'Location_in_txs'])
# 3'UTR      5'UTR     Exon 1     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6 Intron 1-2 Intron 2-3 Intron 3-4 Intron 4-5 Intron 5-6 
#   155         78       2503        434        379        293        750        904        386        363        297        314        339 


#################  Tx anno of variants common in 2 UGT2 genes  #################
UGT2_shared_variants_twoGenes <- rownames(UGT2_variants_canonical)[which(sapply(1:dim(UGT2_variants_canonical)[1], function(x){ length(which(!is.na(UGT2_variants_canonical[x,1:10]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT2_shared_variants_twoGenes, function(x){colnames(UGT2_variants_canonical)[which(!is.na(UGT2_variants_canonical[x,1:10]))]}), 2, toString))
# UGT2A1, UGT2A2 
#            482
table(UGT2_variants_canonical[UGT2_shared_variants_twoGenes, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5   Intron 5-6 
#     10         45         39         28         99        112           14           34           28           35           38 

# --------------------------------------------------------------------------------------------------------
## Check Intron 1-2 variants span a region between exon 2 of UGT2A1 and UGT2A2 and before exon 1 of UGT2A2
pos <- UGT2_variants_canonical[which(UGT2_variants_canonical$Location_in_txs=='Intron 1-2' & rownames(UGT2_variants_canonical) %in% UGT2_shared_variants_twoGenes), 'Position']
min(pos) > location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], 'Exon 2')[[2]][['End']]
# [1] TRUE
max(pos) < location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], 'Exon 1')[[2]][['Start']]
# [1] TRUE

# --------------------------------------------------------------------------------------------------------
## Check all variants in Exon 2 - Exon 6 of UGT2A1 and UGT2A2 are always shared 
exon2to6_varsA1 <- rownames(UGT2A1_data[which(UGT2A1_data$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6')),])
exon2to6_varsA2 <- rownames(UGT2A2_data[which(UGT2A2_data$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6')),])
identical(exon2to6_varsA1, exon2to6_varsA2)
# [1] TRUE



######################
####  UGT3 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT3_variants_canonical[, 'Location_in_txs'])
# 3'UTR      5'UTR     Exon 1     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6     Exon 7   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5   Intron 5-6   Intron 6-7 
#    31         49         71         52         48        295        138        132        191          105           54           43           73           67           62 

## No shared variants between UGT3A1 and UGT3A2
which(sapply(1:dim(UGT3_variants_canonical)[1], function(x){length(which(!is.na(UGT3_variants_canonical[x,1:2])))})==2)
#  integer(0)



######################
####  UGT8 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT8_variants_canonical[, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6    Intron 1-2    Intron 2-3    Intron 3-4    Intron 4-5    Intron 5-6 
#     18        204         27          8         44         90             2            31            32            31            38 


# ------------------------------------------------------------------------------
## Plot number of variants from each category in each gene

## Define colors for variant categories
var_colors <- list('Exon 1'= 'mistyrose2',
                   'Exon 2'= 'lightpink1',
                   'Exon 3'= 'rosybrown3',
                   'Exon 4'= 'palevioletred2',
                   'Exon 5'= 'palevioletred3',
                   'Exon 6'= 'plum3',
                   'Exon 7'= 'orchid3',
                   'Exon'= 'salmon',
                   'Exon 1, Intron 1-2'= 'tan1',
                   '5\'UTR, Intron 1-2'= 'firebrick4',
                   '5\' upstream, Intron 1-2'= 'olivedrab1',
                   'Intron 1-2'= 'lightblue2',
                   'Intron 2-3'= 'lightsteelblue1', 
                   'Intron 3-4'= 'lightsteelblue2',
                   'Intron 4-5'= 'lightsteelblue',
                   'Intron 5-6'= 'lightsteelblue3',
                   'Intron 6-7'= 'lightsteelblue4',
                   '5\' upstream, Exon 1'= 'cornsilk3',
                   '5\' upstream, 5\'UTR'= 'mediumspringgreen',
                   '5\'UTR'='royalblue2',
                   '3\'UTR'= 'slateblue4',
                   '5\' upstream'= 'thistle3', 
                   '3\' downstream'='khaki1',
                   'Exonic variants'= 'salmon',
                   'All variants'= 'darkslategray3')


## Function to create barplot for all variants in genes of a certain family
barplot_gene_fam<- function(gene_family){
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  var_data <- data.frame(matrix(ncol = 3))
  colnames(var_data) <- c('gene', 'location', 'number')
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    data <- data.frame(matrix(ncol = 3))
    colnames(data) <- c('gene', 'location', 'number')
    
    if (gene_family=='UGT1'){
      ordered_locations <- c('5\' upstream', '5\'UTR', 'Intron 4-5', 'Intron 3-4', 'Intron 2-3', 'Intron 1-2',
                             'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', '3\'UTR')
    }
    
    else if(gene_family=='UGT2' | gene_family=='UGT8') {
      ordered_locations <- c('5\'UTR', 'Intron 5-6', 'Intron 4-5', 'Intron 3-4', 'Intron 2-3', 'Intron 1-2',
                             'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', '3\'UTR')
    }
    
    else {
      ordered_locations <- c('5\'UTR', 'Intron 6-7', 'Intron 5-6', 'Intron 4-5', 'Intron 3-4', 'Intron 2-3', 'Intron 1-2',
                             'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', 'Exon 7', '3\'UTR')
    }
    
    locations <- unique(gene_data$Location_in_txs)
    
    for (i in 1:length(ordered_locations)){
      data[i,'gene'] <- gene
      data[i,'location'] <- ordered_locations[i]
      data[i,'number'] <- table(gene_data$Location_in_txs)[ordered_locations[i]]
    }
    var_data <- rbind(var_data, data)
  }
  
  var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
  var_data$gene <- factor(var_data$gene, levels = unique(var_data$gene))
  var_data <- var_data[-1,]
  
  p <- ggplot(var_data, aes(fill=factor(location, levels=ordered_locations), y=number, x=gene)) + 
    geom_bar(position="stack", stat="identity") + 
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of variants in canonical tx', fill='Location') +
    theme_bw() + 
    ylim(0, 2170) +
    scale_fill_manual(values = var_colors) +
    theme(axis.text = element_text(size = 8),
          legend.text = element_text(size=9))
  
  return(p)
}

p1 <- barplot_gene_fam('UGT1')
p2 <- barplot_gene_fam('UGT2')
p3 <- barplot_gene_fam('UGT3')
p4 <- barplot_gene_fam('UGT8')

plot_grid(p1, p2, p3, p4, nrow=1, rel_widths = c(1,1.1, 0.48, 0.40))
ggsave(filename=paste0('plots/01_Data_Processing/All_variants_genes.pdf'), width = 20, height = 6)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##  ** Note there are no Exon 1 variants in UGT1A8 because there's no real tx info. 
##     For the canonical tx most variants are taken as Intron 1-2.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# ______________________________________________________________________________
#  1.1.3.2 Evaluate VEP annotation of exonic variants in each gene

## Extract exonic variants per gene
exonic_variants_per_gene <- function(gene_family){
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    exonic_vars <- gene_data[grep('Exon', gene_data$Location_in_txs), ]
    assign(paste0(gene, '_exonic_data'), exonic_vars, envir = parent.frame())
  }
}  

exonic_variants_per_gene('UGT1')
exonic_variants_per_gene('UGT2')
exonic_variants_per_gene('UGT3')
exonic_variants_per_gene('UGT8')


## Check annotation of exonic variants is consistent with exonic location

############################
####  UGT1 exonic variants
############################

## Unique annotations of exonic variants in all genes of a family
unique(unlist(sapply(UGT1_genes, function(gene){names(table(eval(parse_expr(paste0(gene, '_exonic_data')))[,'VEP_Annotation']))})))

# [1] "frameshift_variant"    "inframe_deletion"      "missense_variant"      "splice_region_variant" "stop_gained"           "synonymous_variant"   
# [7] "start_lost"            "inframe_insertion"     "splice_donor_variant" 

############################
####  UGT2 exonic variants
############################

unique(unlist(sapply(UGT2_genes, function(gene){names(table(eval(parse_expr(paste0(gene, '_exonic_data')))[,'VEP_Annotation']))})))

# [1] "frameshift_variant"    "inframe_deletion"      "missense_variant"      "splice_region_variant" "start_lost"            "stop_gained"          
# [7] "stop_retained_variant" "synonymous_variant"    "stop_lost"             "splice_donor_variant"  "inframe_insertion"    

############################
####  UGT3 exonic variants
############################

unique(as.vector(unlist(sapply(UGT3_genes, function(gene){names(table(eval(parse_expr(paste0(gene, '_exonic_data')))[,'VEP_Annotation']))}))))

# [1] "frameshift_variant"    "inframe_deletion"      "missense_variant"      "splice_region_variant" "start_lost"            "stop_gained"          
# [7] "synonymous_variant"    "inframe_insertion"   

############################
####  UGT8 exonic variants
############################

as.vector(unique(unlist(sapply(UGT8_genes, function(gene){names(table(eval(parse_expr(paste0(gene, '_exonic_data')))[,'VEP_Annotation']))}))))

# [1] "frameshift_variant"    "inframe_deletion"      "inframe_insertion"     "missense_variant"      "splice_region_variant" "start_lost"           
# [7] "stop_gained"           "synonymous_variant" 


## Plot number of exonic variants for each gene
## Function to create barplot for exonic variants in genes of a certain family
barplot_exonic_variants <- function(gene_family){
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  var_data <- data.frame(matrix(ncol = 2))
  colnames(var_data) <- c('gene', 'number')
  
  ## Number of exonic variants per gene
  i=1
  for (gene in genes){
    exonic_gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
    number <- dim(exonic_gene_data)[1]
    var_data[i,] <- c(gene, number)
    i=i+1
  }
  var_data$gene <- factor(var_data$gene, levels = unique(var_data$gene))
  var_data$number <- as.numeric(var_data$number)
  
  
  p <- ggplot(var_data, aes(y=number, x=gene)) + 
    geom_bar(stat="identity", fill=var_colors[['Exon']]) + 
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of exonic variants in canonical tx') +
    theme_bw() + 
    ylim(0, 750) +
    theme(axis.text = element_text(size = 10),
          legend.position="none")
  
  ## Plot with exon categories

  return(p)
}

p1 <- barplot_exonic_variants('UGT1') + theme(plot.margin = margin(30, 5, 30, 5))
p2 <- barplot_exonic_variants('UGT2') + theme(plot.margin = margin(30, 5, 30, 5))
p3 <- barplot_exonic_variants('UGT3') + theme(plot.margin = margin(30, 5, 30, 5))
p4 <- barplot_exonic_variants('UGT8') + theme(plot.margin = margin(30, 5, 30, 5))

plot_grid(p1, p2, p3, p4, nrow=1, rel_widths = c(1,1.1, 0.3, 0.20), align = "hv")
ggsave(filename=paste0('plots/01_Data_Processing/Exonic_variants_genes.pdf'), width = 20, height = 7)


## Plot number of exonic variants from each anno for each gene

exonic_vars_anno_colors <- list(
  "frameshift_variant" ='lightsalmon',
  "inframe_deletion" ='chartreuse2',  
  "inframe_insertion" = 'red1',
  "missense_variant" ='lightskyblue',
  "synonymous_variant" = 'orchid1',
  "splice_region_variant" ='mediumblue',
  "splice_donor_variant" = 'yellow',
  "start_lost" ='deepskyblue',   
  "stop_lost" = 'magenta3',
  "stop_gained" ='wheat2',          
  "stop_retained_variant" ='darkslategray'
)

barplot_exonic_variants_anno <- function(gene_family){
  
  ordered_annotations <- c("frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant",
                           "splice_donor_variant", "splice_region_variant", "synonymous_variant", "start_lost",   
                            "stop_lost", "stop_gained", "stop_retained_variant")
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  exonic_var_data <- data.frame(matrix(ncol = 3))
  colnames(exonic_var_data) <- c('gene', 'annotation', 'number')
  
  ## Number of exonic variants from each annotation per gene
  i=1
  for (gene in genes){
    exonic_gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
    for (annotation in ordered_annotations){
      number <- table(exonic_gene_data$VEP_Annotation)[annotation]
      annotation <- annotation
      exonic_var_data[i,] <- c(gene, annotation, number)
      i=i+1
    }
  }
  
  exonic_var_data$gene <- factor(exonic_var_data$gene, levels = unique(exonic_var_data$gene))
  exonic_var_data$number <- as.numeric(exonic_var_data$number)
  exonic_var_data$number <- replace(exonic_var_data$number, which(is.na(exonic_var_data$number)), 0)
  
  p <- ggplot(exonic_var_data, aes(y=number, x=gene, fill=factor(annotation, levels = ordered_annotations))) + 
    geom_bar(position="stack", stat="identity") + 
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of exonic variants in canonical tx',
         fill='Predicted effect') +
    theme_bw() + 
    scale_fill_manual(values=exonic_vars_anno_colors) +
    ylim(0, 750) +
    theme(axis.text = element_text(size = 11),
          legend.text = element_text(size = 11))
  
  return(p)
}

p1 <- barplot_exonic_variants_anno('UGT1') + theme(legend.position="none")
p2 <- barplot_exonic_variants_anno('UGT2') + theme(legend.position="none")
p3 <- barplot_exonic_variants_anno('UGT3') + theme(legend.position="none")
p4 <- barplot_exonic_variants_anno('UGT8')

plot_grid(p1, p2, p3, p4, ncol=4, rel_widths = c(1,1.1, 0.3, 0.48))
ggsave(filename=paste0('plots/01_Data_Processing/Exonic_variants_genes_anno.pdf'), width = 22, height = 7)




# _____________________________________________________________________________
#  1.1.4 Variant quantification in each gene family
# _____________________________________________________________________________

################# Number of variants in each location category #################

## Create stacked barplot for each gene family
var_data <- data.frame(matrix(ncol = 3))
colnames(var_data) <- c('gene_family', 'location', 'number')

locations <- c('5\' upstream', '5\'UTR', '5\' upstream, 5\'UTR',  'Intron 6-7', 'Intron 5-6',  'Intron 4-5', 'Intron 3-4', 'Intron 2-3', 'Intron 1-2', 
                       'Exon 1, Intron 1-2', '5\'UTR, Intron 1-2', '5\' upstream, Intron 1-2', '5\' upstream, Exon 1', 
                       'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', 'Exon 7', '3\'UTR')

for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  
    UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants_canonical')))
    data <- data.frame(matrix(ncol = 3))
    colnames(data) <- c('gene_family', 'location', 'number')
    
    for (i in 1:length(locations)){
      data[i,'gene_family'] <- gene_family
      data[i,'location'] <- locations[i]
      data[i,'number'] <- table(UGT_variants$Location_in_txs)[locations[i]]
    }
  var_data <- rbind(var_data, data)
}
  
var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
var_data$gene_family <- factor(var_data$gene_family, levels = unique(var_data$gene_family))
var_data <- var_data[-1,]
  
p1 <- ggplot(var_data, aes(fill=factor(location, levels=locations), y=number, x=gene_family)) + 
  geom_bar(position="stack", stat="identity") + 
  labs(x='UGT gene family', y='Total number of variants in canonical tx of genes', fill='Location') +
  theme_bw() + 
  ylim(0, 7200) +
  scale_fill_manual(values = var_colors) +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size=9),
        plot.margin = margin(10, 2, 10, 2))
  

####################### Number of total/exonic variants #######################

## Plot total number of variants and number of exonic variants per gene family
data <- data.frame(matrix(ncol = 3))
colnames(data) <- c('gene_family',  'category' ,'number')

i=1
for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants_canonical')))
  ## Exonic variants only: shared variants with one exonic location are exonic by definition (in at least one gene)
  UGT_exonic_variants <- UGT_variants[grep('Exon', UGT_variants$Location_in_txs),]
  
  data[i, ] <- c(gene_family, 'All variants', dim(UGT_variants)[1])
  data[i+1, ] <- c(gene_family, 'Exonic variants', dim(UGT_exonic_variants)[1])
  i=i+2
}
data$number <- as.numeric(data$number)

p2 <- ggplot(data, aes(y=number, x=gene_family, fill=category)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  labs(x='UGT gene family', y='Number of variants in canonical tx of genes', fill='Type') +
  theme_bw() + 
  scale_fill_manual(values = var_colors) +
  ylim(0, 7200) +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size=9),
        plot.margin = margin(10, 2, 10, 2))


################# Number of exonic variants of each annotation #################

## Function to create table of exonic variants and their location in all genes of each family
create_exonic_gene_fam_table <- function(gene_family){
  
  ## Genes of gene family
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  
  ## Unique exonic variants in gene family
  unique_UGT_exonic_variants <- unique(unlist(sapply(paste0(genes, '_exonic_data$Variant_ID'), function(x){eval(parse_expr(x))})))

  UGT_exonic_variants <- data.frame(matrix(ncol = length(genes)+3, nrow = 0))
  colnames(UGT_exonic_variants) <- c(genes,  'Position', 'VEP_Annotation', 'Location_in_txs')
  
  i=1
  for (variant in unique_UGT_exonic_variants){
    variant_txs <- sapply(genes, function(gene){
      tx <- eval(parse_expr(paste0(gene, '_exonic_data')))[which(eval(parse_expr(paste0(gene, '_exonic_data$Variant_ID'))) == variant), 'Transcript']
      if (length(tx)==0) {NA}
      else {tx}
    })
    
    UGT_exonic_variants[i, 1:length(genes)] <- variant_txs
    
    ## Add variant position
    UGT_exonic_variants$Position[i] <- as.numeric(strsplit(variant, '-')[[1]][2])
    
    ## VEP annotation for variant in each gene 
    variant_anno <- unlist(sapply(genes, function(gene){
      eval(parse_expr(paste0(gene, '_exonic_data')))[which(eval(parse_expr(paste0(gene, '_exonic_data$Variant_ID'))) == variant), 'VEP_Annotation']
    }))
    if (length(unique(variant_anno))==1){
      UGT_exonic_variants$VEP_Annotation[i] <- unique(variant_anno)
    }
    else{
      UGT_exonic_variants$VEP_Annotation[i] <- toString(variant_anno)
    }
    
    ## Add tx location 
    location <- unlist(sapply(genes, function(gene){
      eval(parse_expr(paste0(gene, '_exonic_data')))[which(eval(parse_expr(paste0(gene, '_exonic_data$Variant_ID'))) == variant),
                                                     'Location_in_txs']}))
    if (length(unique(variant_anno))==1){
      UGT_exonic_variants$Location_in_txs[i] <- unique(location)
    }
    else{
      UGT_exonic_variants$Location_in_txs[i] <- toString(location)
    }
    
    i=i+1
  }
  
  rownames(UGT_exonic_variants) <- unique_UGT_exonic_variants
  
  assign(paste0('unique_', gene_family, '_exonic_variants'), unique_UGT_exonic_variants, envir = parent.frame())
  assign(paste0(gene_family, '_exonic_variants'), UGT_exonic_variants, envir = parent.frame())
}

create_exonic_gene_fam_table('UGT1')
create_exonic_gene_fam_table('UGT2')
create_exonic_gene_fam_table('UGT3')
create_exonic_gene_fam_table('UGT8')


## Stacked barplot
annotations <- c("frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant",
                 "splice_donor_variant", "splice_region_variant", "synonymous_variant", "start_lost",   
                 "stop_lost", "stop_gained", "stop_retained_variant")

var_data <- data.frame(matrix(ncol = 3))
colnames(var_data) <- c('gene_family', 'annotation', 'number')

for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  
  UGT_exonic_variants <- eval(parse_expr(paste0(gene_family, '_exonic_variants')))
  data <- data.frame(matrix(ncol = 3))
  colnames(data) <- c('gene_family', 'annotation', 'number')
  
  for (i in 1:length(annotations)){
    data[i,'gene_family'] <- gene_family
    data[i,'annotation'] <- annotations[i]
    data[i,'number'] <- table(UGT_exonic_variants$VEP_Annotation)[annotations[i]]
  }
  var_data <- rbind(var_data, data)
}

var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
var_data$gene_family <- factor(var_data$gene_family, levels = unique(var_data$gene_family))
var_data <- var_data[-1,]

p3 <- ggplot(var_data, aes(fill=factor(annotation, levels=annotations), y=number, x=gene_family)) + 
  geom_bar(position="stack", stat="identity") + 
  labs(x='UGT gene family', y='Number of exonic variants in canonical tx of genes', fill='Predicted effect') +
  theme_bw() + 
  ylim(0, 7200) +
  scale_fill_manual(values = exonic_vars_anno_colors) +
  theme(axis.text = element_text(size = 8),
        legend.text = element_text(size=9),
        plot.margin = margin(10, 2, 10, 2))


## Plots for gene families
plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.9, 1, 0.68))
ggsave(filename=paste0('plots/01_Data_Processing/Num_variants_per_gene_fam.pdf'), width = 20, height = 7)














## FUTURE

## Define functional categories of predicted effects of exonic variants

LoF <- c('splice_donor', 'splice_acceptor', 'stop_lost', 'stop_gained', 'frameshift')
neutral <- c('synonymous')
## missense functions with ANOVA
# deletereous?
  
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
