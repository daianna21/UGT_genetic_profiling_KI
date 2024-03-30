
library(here)
library(readr)
library(rlang)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
library(sessioninfo)


####################################################################################################
##                                1. Data Processing and Filtering
####################################################################################################

## UGT genes
UGT_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10',
               'UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28',
               'UGT3A1', 'UGT3A2', 
               'UGT8')


## Process and valid variables for each gene 
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

}


################################################################################
##             1.1  Obtain variants within each gene feature
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

############################
####    UGT2 genes 
############################
UGT2_genes <- c('UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28')

## Confirm all UGT2 variants are in the same chr
sapply(UGT2_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT2A1  UGT2A2  UGT2A3  UGT2B4  UGT2B7 UGT2B10 UGT2B11 UGT2B15 UGT2B17 UGT2B28 
#      4       4       4       4       4       4       4       4       4       4 

############################
####    UGT3 genes 
############################
UGT3_genes <- c('UGT3A1', 'UGT3A2')

## Confirm all UGT3 variants are in the same chr
sapply(UGT3_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT3A1 UGT3A2 
#      5      5 

############################
####    UGT8 genes 
############################

UGT8_genes <- c('UGT8')
unique(UGT8_data$Chromosome)
#    4 


# _______________________________________________________________________________
#  1.1.1.1 Subset to variants present in the canonical transcripts* of UGT genes

##  * Canonical gene transcript in Ensembl (under GRCh37 / hg19 genome build; under  GRCh38 for UGT1A8)


## Define canonical tx for each gene
canonical_UGT1_txs <- list('UGT1A1'= 'ENST00000305208.5', 'UGT1A3'='ENST00000482026.1', 'UGT1A4'='ENST00000373409.3', 
                           'UGT1A5'='ENST00000373414.3', 'UGT1A6'='ENST00000305139.6', 'UGT1A7'='ENST00000373426.3', 
                           'UGT1A8'= 'ENST00000373450.4','UGT1A9'= 'ENST00000354728.4', 'UGT1A10'='ENST00000344644.5')

canonical_UGT2_txs <- list('UGT2A1'= 'ENST00000503640.1', 'UGT2A2'='ENST00000457664.2', 'UGT2A3'='ENST00000251566.4', 
                           'UGT2B4'='ENST00000305107.6', 'UGT2B7'='ENST00000305231.7', 'UGT2B10'='ENST00000265403.7', 
                           'UGT2B11'= 'ENST00000446444.1', 'UGT2B15'= 'ENST00000338206.5', 'UGT2B17'='ENST00000317746.2', 
                           'UGT2B28'='ENST00000335568.5')

canonical_UGT3_txs <- list('UGT3A1'= 'ENST00000274278.3', 'UGT3A2'='ENST00000282507.3')

canonical_UGT8_txs <- list('UGT8'= 'ENST00000310836.6')

save(canonical_UGT1_txs, file = paste0('processed-data/01_Data_Processing/canonical_UGT1_txs.Rdata'))
save(canonical_UGT2_txs, file = paste0('processed-data/01_Data_Processing/canonical_UGT2_txs.Rdata'))
save(canonical_UGT3_txs, file = paste0('processed-data/01_Data_Processing/canonical_UGT3_txs.Rdata'))
save(canonical_UGT8_txs, file = paste0('processed-data/01_Data_Processing/canonical_UGT8_txs.Rdata'))


## % of variants in a gene that are present in the canonical tx of the gene

for (gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')) {
  canonical_UGT_txs <- eval(parse_expr(paste0('canonical_', gene_family, '_txs')))
  
  ## Genes of gene family
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  
  ## For each gene
  for(gene in genes){
    
    gene_data <- eval(parse_expr(paste0(gene, '_data')))
    gene_data$Canonical_txs <- sapply(gene_data$Transcript, function(x){if (x==canonical_UGT_txs[[gene]]){TRUE} else {FALSE}})
    gene_canonical_data <- gene_data[gene_data$Canonical_txs==TRUE,]
    
    print(paste0(gene, ': ', length(which(gene_data$Canonical_txs==TRUE))/dim(gene_data)[1]*100, '%,  ', 
                 length(which(gene_data$Canonical_txs==TRUE)), ' variants'))
    
    rownames(gene_data) <- gene_data$Variant_ID
    rownames(gene_canonical_data) <- gene_canonical_data$Variant_ID
    
    assign(paste0(gene, '_data'), gene_data)
    assign(paste0(gene, '_canonical_data'), gene_canonical_data)
    save(gene_data, file = paste0('processed-data/01_Data_Processing/', gene, '_data.Rdata'))
  }
  
}

## % and number of variants in each gene that are in the canonical tx

###################
####  UGT1 genes
###################

# [1] "UGT1A1: 100%  653 variants"
# [1] "UGT1A3: 100%  766 variants"
# [1] "UGT1A4: 100%  838 variants"
# [1] "UGT1A5: 100%  691 variants"
# [1] "UGT1A6: 96.0352422907489%  654 variants"
# [1] "UGT1A7: 95.4068241469816%  727 variants"
# [1] "UGT1A8: 21.7351598173516%,  476 variants"
# [1] "UGT1A9: 100%  714 variants"
# [1] "UGT1A10: 97.1326164874552%  813 variants"

###################
####  UGT2 genes
###################

# [1] "UGT2A1: 97.4126778783959%,  753 variants"
# [1] "UGT2A2: 99.8633879781421%,  731 variants"
# [1] "UGT2A3: 98.910411622276%,  817 variants"
# [1] "UGT2B4: 94.6475195822454%,  725 variants"
# [1] "UGT2B7: 99.4750656167979%,  758 variants"
# [1] "UGT2B10: 99.8820754716981%,  847 variants"
# [1] "UGT2B11: 100%,  951 variants"
# [1] "UGT2B15: 100%,  715 variants"
# [1] "UGT2B17: 100%,  562 variants"
# [1] "UGT2B28: 100%,  818 variants"

###################
####  UGT3 genes
###################

# [1] "UGT3A1: 93.3673469387755%,  732 variants"
# [1] "UGT3A2: 99.5601173020528%,  679 variants"

###################
####  UGT8 gene
###################

# [1] "UGT8: 97.2222222222222%,  525 variants"


# _______________________________________________________________________________
#  1.1.1.2 Find variants shared in canonical txs of the genes in a family

## Determine in which canonical gene txs each variant is present

create_gene_fam_table <- function(gene_family){
  
  ## Genes of gene family
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  
  ## Unique canonical variants in gene family
  unique_UGT_canonical_variants <- unique(unlist(sapply(paste0(genes, '_canonical_data$Variant_ID'), function(x){eval(parse_expr(x))})))
  ## Total unique variants in gene family
  unique_UGT_variants <- unique(unlist(sapply(paste0(genes, '_data$Variant_ID'), function(x){eval(parse_expr(x))})))
  
  ## Annotate canonical transcript in which each variant appears for each gene 
  ## (if present in a gene, variants are expected to appear only in one transcript since they are all unique in each gene dataset)
  UGT_variants <- data.frame(matrix(ncol = length(genes)+2, nrow = 0))
  colnames(UGT_variants) <- c(genes, 'VEP_Annotation', 'Position')
  
  i=1
  for (variant in unique_UGT_canonical_variants){
    
    ## Search variant in each gene; if present, annotate tx
    variant_txs <- sapply(genes, function(gene){
      tx <- eval(parse_expr(paste0(gene, '_canonical_data')))[which(eval(parse_expr(paste0(gene, '_canonical_data$Variant_ID'))) == variant), 'Transcript']
      if (length(tx)==0) {NA}
      else {tx}
    })
    
    UGT_variants[i, 1:length(genes)] <- variant_txs
    
    ## VEP annotation for variant in each gene that contains it
    variant_anno <- unlist(sapply(genes, function(gene){
      eval(parse_expr(paste0(gene, '_canonical_data')))[which(eval(parse_expr(paste0(gene, '_canonical_data$Variant_ID'))) == variant), 'VEP_Annotation']
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
  
  rownames(UGT_variants) <- unique_UGT_canonical_variants
  
  assign(paste0('unique_', gene_family, '_canonical_variants'), unique_UGT_canonical_variants, envir = parent.frame())
  assign(paste0('unique_', gene_family, '_total_variants'), unique_UGT_variants, envir = parent.frame())
  assign(paste0(gene_family, '_variants'), UGT_variants, envir = parent.frame())
}


############################
####      UGT1 genes 
############################

create_gene_fam_table('UGT1')

############################
####      UGT2 genes 
############################

create_gene_fam_table('UGT2')

############################
####      UGT3 genes 
############################

create_gene_fam_table('UGT3')

############################
####      UGT8 genes 
############################

create_gene_fam_table('UGT8')


## Verify there are no overlapping variants between genes from different families

length(intersect(unique_UGT1_total_variants, unique_UGT2_total_variants))
#  0
length(intersect(unique_UGT1_total_variants, unique_UGT3_total_variants))
#  0
length(intersect(unique_UGT1_total_variants, unique_UGT8_total_variants))
#  0
length(intersect(unique_UGT2_total_variants, unique_UGT3_total_variants))
#  0
length(intersect(unique_UGT2_total_variants, unique_UGT8_total_variants))
#  0
length(intersect(unique_UGT3_total_variants, unique_UGT8_total_variants))
#  0


## Percentage and number of variants in each gene family that are in canonical txs
## (For shared variants, if they are in at least one canonical tx)

for (gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')) {
  
  unique_canonical_variants <- eval(parse_expr(paste0('unique_', gene_family, '_canonical_variants')))
  unique_total_variants <- eval(parse_expr(paste0('unique_', gene_family, '_total_variants')))
  
  print(paste0(gene_family, ': ', length(unique_canonical_variants)/length(unique_total_variants)*100, ', ', 
               length(unique_canonical_variants), ' variants in canonical txs'))

}

# [1] "UGT1: 98.6129997280392, 3626 variants in canonical txs"
# [1] "UGT2: 98.9547517535415, 7195 variants in canonical txs"
# [1] "UGT3: 96.2482946793997, 1411 variants in canonical txs"
# [1] "UGT8: 97.2222222222222, 525 variants in canonical txs"



## Total number of UGT variants in canonical txs

gene_families <- c('UGT1', 'UGT2', 'UGT3', 'UGT8')
vars_in_canonical_txs <- sum(unlist(sapply(paste0('length(unique_', gene_families, '_canonical_variants)'), function(x){eval(parse_expr(x))})))
## Total variants in all UGT families
total_vars <- sum(unlist(sapply(paste0('length(unique_', gene_families, '_total_variants)'), function(x){eval(parse_expr(x))})))

## % of total UGT variants in canonical txs
vars_in_canonical_txs / total_vars * 100
# [1] 98.47923




# _____________________________________________________________________________
#  1.1.2 Manual annotation of variants based on boundaries of canonical txs
# _____________________________________________________________________________

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
                       'ENST00000503640.1'= 70455089, 
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
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_canonical_data')))
    ## Location
    gene_data$Location_in_txs <- sapply(1:dim(gene_data)[1], function(i){
      location_determination(gene_data$Position[i], gene_data$Transcript[i], NULL)[[1]]})
    assign(paste0(gene, '_canonical_data'), gene_data, envir = parent.frame())
    save(gene_data, file = paste0('processed-data/01_Data_Processing/', gene, '_canonical_data.Rdata'))
  }
}

## Add tx location of canonical variants in each gene dataset
add_tx_location('UGT1')
add_tx_location('UGT2')
add_tx_location('UGT3')
add_tx_location('UGT8')


## Add tx location of canonical variants in each gene family dataset
for (gene_family in gene_families){
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants')))
  
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
  assign(paste0(gene_family, '_variants'), UGT_variants)
  save(UGT_variants, file = paste0('processed-data/01_Data_Processing/', gene_family, '_variants.Rdata'))
}


# __________________________________________________________________________________
#  1.1.2.1  Check manual annotation of all and shared variants in each gene family

######################
####  UGT1 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT1_variants[, 'Location_in_txs'])
#   3'UT       5' upstream              5'UTR       5'UTR, Intron 1-2                  Exon 1                  Exon 2 
#      10                1                106                       8                    2845                      25 
#   Exon 3                     Exon 4                  Exon 5              Intron 1-2         Intron 1-2, 5'UTR     Intron 2-3              Intron 3-4 
#       45                         66                     122                     255                         4             32                      29 
#   Intron 4-5 
#           78 


################  Tx anno of variants common in 8 UGT1 genes  ################ 
UGT1_shared_variants_8Genes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })==8)]

table(UGT1_variants[UGT1_shared_variants_8Genes, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5 
#     10         25         45         66        122           12           32           29           33 

## Which 8 genes share those variants
table(apply(sapply(UGT1_shared_variants_8Genes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,1:9]))]}), 2, toString))
# UGT1A1, UGT1A3, UGT1A4, UGT1A5, UGT1A6, UGT1A7, UGT1A9, UGT1A10 
#                                                             374 

# ------------------------------------------------------
## 1.  Position of shared variants in Intron 1-2 region
pos <- UGT1_variants[which(UGT1_variants$Location_in_txs=='Intron 1-2' & rownames(UGT1_variants) %in% UGT1_shared_variants_8Genes), 'Position']
## Confirm this overlapping intronic region starts after Exon 1 of UGT1A1 and finishes before exon 2 of all UGT1 genes
min(pos) > location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], 'Exon 1')[[2]]['End']
# TRUE
max(pos) < location_determination(max(pos), canonical_UGT1_txs[['UGT1A1']], 'Exon 2')[[2]]['Start']
# TRUE

# ------------------------------------------------------
## 2.  Check all variants in Exon 2 - Exon 5 + 3'UTR region are shared by all these 8 UGT1 genes
exon2to5_vars <- rownames(UGT1_variants[which(UGT1_variants$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', '3\'UTR')), ])
which(sapply(exon2to5_vars, function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })!=8)
# named integer(0)

## Check if variants in shared introns (Intron 4-5, Intron 3-4 and Intron 2-3) are shared by all 8 UGT1 genes
intron4to5_vars <- rownames(UGT1_variants[which(UGT1_variants$Location_in_txs == 'Intron 4-5'), ])
length(which(sapply(intron4to5_vars, function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })!=8))
# 45

intron3to4_vars <- rownames(UGT1_variants[which(UGT1_variants$Location_in_txs == 'Intron 3-4'), ])
length(which(sapply(intron3to4_vars, function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })!=8))
# 0

intron2to3_vars <- rownames(UGT1_variants[which(UGT1_variants$Location_in_txs == 'Intron 2-3'), ])
length(which(sapply(intron2to3_vars, function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })!=8))
# 0

## Add the shared variants from exon 2 to 3' UTR to UGT1A8 
UGT1A8_canonical_data <- rbind(UGT1A8_canonical_data, UGT1A1_canonical_data[UGT1_shared_variants_8Genes,])
UGT1A8_canonical_data$Transcript <- canonical_UGT1_txs[['UGT1A8']]
save(UGT1A8_canonical_data, file = 'processed-data/01_Data_Processing/UGT1A8_canonical_data.Rdata')   
UGT1_variants[UGT1_shared_variants_8Genes, 'UGT1A8'] <- canonical_UGT1_txs[['UGT1A8']]


## Verify
################  Tx anno of variants common in all UGT1 genes  ################ 
UGT1_shared_variants_9Genes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })==9)]
table(UGT1_variants[UGT1_shared_variants_9Genes, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5 
#     10         25         45         66        122           12           32           29           33 


#################  Tx anno of variants common in 3 UGT1 genes  #################
UGT1_shared_variants_3Genes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){ length(which(!is.na(UGT1_variants[x,1:9]))) })==3)]
## Which 4 genes share those variants?
table(apply(sapply(UGT1_shared_variants_3Genes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,1:9]))]}), 2, toString))
# UGT1A4, UGT1A6, UGT1A10 
#                      21
table(UGT1_variants[UGT1_shared_variants_3Genes, 'Location_in_txs'])
# Intron 4-5 
#         21 


#################  Tx anno of variants common in 2 UGT1 genes  #################
UGT1_shared_variants_2Genes <- rownames(UGT1_variants)[which(sapply(1:dim(UGT1_variants)[1], function(x){length(which(!is.na(UGT1_variants[x,1:9]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT1_shared_variants_2Genes, function(x){colnames(UGT1_variants)[which(!is.na(UGT1_variants[x,1:9]))]}), 2, toString))
# UGT1A3, UGT1A8      UGT1A5, UGT1A8      UGT1A8, UGT1A9 
#              8                  21                  17 
table(UGT1_variants[UGT1_shared_variants_2Genes, 'Location_in_txs'])
# 5'UTR, Intron 1-2        Intron 1-2         Intron 1-2, 5'UTR 
#                 8                34                         4 

# ------------------------------------------------------------------------------
## 1.  Genes with Intron 1-2 / 5' UTR variants
vars <- which(UGT1_variants$Location_in_txs == 'Intron 1-2, 5\'UTR')
table(apply(sapply(vars, function(x){colnames(UGT1_variants[x, which(!is.na(UGT1_variants[x, 1:9]))])}), 2, toString))
# UGT1A8, UGT1A9 
#              4 

## Check this variant is in 5'-UTR of UGT1A9 = Intron 1-2 of UGT1A8
pos <- UGT1_variants[vars, 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A9']], NULL)[[1]]}))
#  5'UTR 
#      4
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#          4 

# -------------------------------------------------------------------------------
## 2.  Genes with 5'-UTR / Intron 1-2 variants
vars <- which(UGT1_variants$Location_in_txs == '5\'UTR, Intron 1-2')
table(apply(sapply(vars, function(x){colnames(UGT1_variants[x, which(!is.na(UGT1_variants[x, 1:9]))])}), 2, toString))
# UGT1A3, UGT1A8 
#              8

## Check these variants are in 5'-UTR of UGT1A3 = Intron 1-2 of UGT1A8
pos <- UGT1_variants[vars, 'Position']

table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A3']], NULL)[[1]]}))
# 5'UTR 
#     8
table(sapply(pos, function(x){location_determination(x, canonical_UGT1_txs[['UGT1A8']], NULL)[[1]]}))
# Intron 1-2 
#          8

# ------------------------------------------------------
## 3.  Genes with Intron 1-2 variants
vars <- which(UGT1_variants$Location_in_txs == 'Intron 1-2' & rownames(UGT1_variants) %in% UGT1_shared_variants_2Genes)
table(apply(sapply(vars, function(x){colnames(UGT1_variants[x, which(!is.na(UGT1_variants[x, 1:9]))])}), 2, toString))
# UGT1A5, UGT1A8      UGT1A8, UGT1A9 
#             21                  13 

## Corroborate the common intronic region of UGT1A5-UGT1A8 and UGT1A9-UGT1A8 goes after 1st exon and before 2nd exon of UGT1A5 and UGT1A9, respectively
pos <- UGT1_variants[vars, 'Position']
posA5 <- UGT1_variants[vars,][which(!is.na(UGT1_variants[vars,'UGT1A5'])), 'Position']
min(posA5) > location_determination(0, canonical_UGT1_txs[['UGT1A5']], 'Exon 1')[[2]][['End']]
# [1] TRUE
max(posA5) < location_determination(0, canonical_UGT1_txs[['UGT1A5']], 'Exon 2')[[2]][['Start']]
# [1] TRUE

posA9 <- UGT1_variants[vars,][which(!is.na(UGT1_variants[vars,'UGT1A9'])), 'Position']
min(posA9) > location_determination(0, canonical_UGT1_txs[['UGT1A9']], 'Exon 1')[[2]][['End']]
# [1] TRUE
max(posA9) < location_determination(0, canonical_UGT1_txs[['UGT1A9']], 'Exon 2')[[2]][['Start']]
# [1] TRUE



######################
####  UGT2 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT2_variants[, 'Location_in_txs'])
# 3'UTR      5'UTR     Exon 1     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6 Intron 1-2 Intron 2-3 Intron 3-4 Intron 4-5 Intron 5-6 
#   155         78       2503        434        379        293        750        904        386        363        297        314        339 


#################  Tx anno of variants common in 2 UGT2 genes  #################
UGT2_shared_variants_2Genes <- rownames(UGT2_variants)[which(sapply(1:dim(UGT2_variants)[1], function(x){ length(which(!is.na(UGT2_variants[x,1:10]))) })==2)]
## Which 2 genes share those variants?
table(apply(sapply(UGT2_shared_variants_2Genes, function(x){colnames(UGT2_variants)[which(!is.na(UGT2_variants[x,1:10]))]}), 2, toString))
# UGT2A1, UGT2A2 
#            482
table(UGT2_variants[UGT2_shared_variants_2Genes, 'Location_in_txs'])
#  3'UTR     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5   Intron 5-6 
#     10         45         39         28         99        112           14           34           28           35           38 

# --------------------------------------------------------------------------------------------------------
## Check Intron 1-2 variants span a region between exon 2 of UGT2A1 and UGT2A2 and before exon 1 of UGT2A2
pos <- UGT2_variants[which(UGT2_variants$Location_in_txs=='Intron 1-2' & rownames(UGT2_variants) %in% UGT2_shared_variants_2Genes), 'Position']
min(pos) > location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], 'Exon 2')[[2]][['End']]
# [1] TRUE
max(pos) < location_determination(pos[1], canonical_UGT2_txs[['UGT2A1']], 'Exon 1')[[2]][['Start']]
# [1] TRUE

# --------------------------------------------------------------------------------------------------------
## Check all variants in Exon 2 - Exon 6 of UGT2A1 and UGT2A2 are always shared 
exon2to6_varsA1 <- rownames(UGT2A1_canonical_data[which(UGT2A1_canonical_data$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6')),])
exon2to6_varsA2 <- rownames(UGT2A2_canonical_data[which(UGT2A2_canonical_data$Location_in_txs %in% c('Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6')),])
identical(exon2to6_varsA1, exon2to6_varsA2)
# [1] TRUE



######################
####  UGT3 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT3_variants[, 'Location_in_txs'])
# 3'UTR      5'UTR     Exon 1     Exon 2     Exon 3     Exon 4     Exon 5     Exon 6     Exon 7   Intron 1-2   Intron 2-3   Intron 3-4   Intron 4-5   Intron 5-6   Intron 6-7 
#    31         49         71         52         48        295        138        132        191          105           54           43           73           67           62 

## No shared variants between UGT3A1 and UGT3A2
which(sapply(1:dim(UGT3_variants)[1], function(x){length(which(!is.na(UGT3_variants[x,1:2])))})==2)
#  integer(0)



######################
####  UGT8 variants
######################

#####################  Manual annotation of all variants  ######################
table(UGT8_variants[, 'Location_in_txs'])
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
                   '5\'UTR, Intron 1-2'= 'firebrick4',
                   'Intron 1-2, 5\'UTR'= 'olivedrab1',
                   'Intron 1-2'= 'lightblue2',
                   'Intron 2-3'= 'lightsteelblue1', 
                   'Intron 3-4'= 'lightsteelblue2',
                   'Intron 4-5'= 'lightsteelblue',
                   'Intron 5-6'= 'lightsteelblue3',
                   'Intron 6-7'= 'lightsteelblue4',
                   '5\'UTR'='royalblue2',
                   '3\'UTR'= 'slateblue4',
                   '5\' upstream'= 'thistle3', 
                   '3\' downstream'='khaki1',
                   'Exonic variants'= 'salmon',
                   'Non-exonic variants'= 'darkseagreen3')


## Function to create barplot for all variants in genes of a certain family
barplot_gene_fam<- function(gene_family){
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  var_data <- data.frame(matrix(ncol = 3))
  colnames(var_data) <- c('gene', 'location', 'number')
  total_num <- list()
  
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_canonical_data')))
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
    total_num[[gene]] <- sum(data[!is.na(data$number), 'number'])
  }
  
  var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
  var_data$gene <- factor(var_data$gene, levels = unique(var_data$gene))
  var_data <- var_data[-1,]
  total_num <- melt(total_num)
  colnames(total_num) <- c('n', 'gene')
  
  p <- ggplot(var_data, aes(fill=factor(location, levels=ordered_locations), y=number, x=gene)) + 
    geom_bar(position="stack", stat="identity", width = 0.8) + 
    geom_text(data=total_num, aes(label=n, y=n+2, x=gene, fill=NULL), vjust=-0.25, size=3.2) +
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of reported variants in canonical transcript', fill='Location') +
    theme_classic() +
    scale_fill_manual(values = var_colors) +
    scale_y_continuous(limits = c(0,990), expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, face='italic'),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size=11), 
          legend.title = element_text(size =12, face='bold'),
          axis.title = element_text(size = (11.5), face='bold'))
  
  return(p)
}

p1 <- barplot_gene_fam('UGT1')
p2 <- barplot_gene_fam('UGT2')
p3 <- barplot_gene_fam('UGT3')
p4 <- barplot_gene_fam('UGT8')

plot_grid(p1, p2, p3, p4, nrow=1, rel_widths = c(1,1.08, 0.52, 0.46))
ggsave(filename=paste0('plots/01_Data_Processing/All_variants_genes.pdf'), width = 16, height = 5.1)



# ______________________________________________________________________________
#  1.1.2.2 Evaluate VEP annotation of exonic variants in each gene

## Extract exonic variants per gene
exonic_variants_per_gene <- function(gene_family){
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  for (gene in genes){
    gene_data <- eval(parse_expr(paste0(gene, '_canonical_data')))
    exonic_vars <- gene_data[grep('Exon', gene_data$Location_in_txs), ]
    assign(paste0(gene, '_exonic_data'), exonic_vars, envir = parent.frame())
    save(exonic_vars, file = paste0('processed-data/01_Data_Processing/', gene, '_exonic_data.Rdata'))
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
    geom_bar(stat="identity", fill=var_colors[['Exon']], width = 0.8) + 
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of exonic variants in canonical transcript') +
    geom_text(aes(label=number), vjust=-0.25, size=3.4) +
    theme_classic() + 
    scale_y_continuous(limits = c(0,740), expand = c(0,0)) +
    theme(plot.margin = margin(30, 20, 30, 20),
          axis.text = element_text(size = 10),
          legend.position="none", 
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 10, face='italic'),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = (11.5), face='bold'))
  
  ## Plot with exon categories

  return(p)
}

p1 <- barplot_exonic_variants('UGT1') 
p2 <- barplot_exonic_variants('UGT2')
p3 <- barplot_exonic_variants('UGT3') 
p4 <- barplot_exonic_variants('UGT8') 

plot_grid(p1, p2, p3, p4, NULL, nrow=1, labels = c("A", "B", "C", "D"), rel_widths = c(1,1.08, 0.41, 0.32, 0.01))
ggsave(filename=paste0('plots/01_Data_Processing/Exonic_variants_genes.pdf'), width = 13, height = 6)


## Plot number of exonic variants from each anno for each gene

exonic_vars_anno_colors <- list(
  "frameshift_variant" ='lightsalmon2',
  "missense_variant" ='deeppink2',
  "synonymous_variant" = 'darkslategray3',
  "stop_gained" ='wheat2',
  "other"='gray60'
)

barplot_exonic_variants_anno <- function(gene_family){
  
  ordered_annotations <- c("other", "stop_gained", "frameshift_variant", "synonymous_variant", "missense_variant")
  
  ## Location of variants in each gene of the family
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  exonic_var_data <- data.frame(matrix(ncol = 3))
  colnames(exonic_var_data) <- c('gene', 'annotation', 'number')
  
  ## Number of exonic variants from each annotation per gene
  total_num <- list()
  i=1
  for (gene in genes){
    exonic_gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
    ## Define 5 categories
    exonic_gene_data$anno <- sapply(exonic_gene_data$VEP_Annotation, function(x){if(x %in% ordered_annotations){x}else{'other'}})
    for (annotation in ordered_annotations){
      number <- table(exonic_gene_data$anno)[annotation]
      annotation <- annotation
      exonic_var_data[i,] <- c(gene, annotation, number)
      i=i+1
    }
    total_num[[gene]] <- dim(exonic_gene_data)[1]
  }
  total_num <- melt(total_num)
  colnames(total_num) <- c('n', 'gene')
  
  exonic_var_data$gene <- factor(exonic_var_data$gene, levels = unique(exonic_var_data$gene))
  exonic_var_data$annotation <- factor(exonic_var_data$annotation, levels=ordered_annotations)
  exonic_var_data$number <- as.numeric(exonic_var_data$number)
  exonic_var_data$number <- replace(exonic_var_data$number, which(is.na(exonic_var_data$number)), 0)
  
  p <- ggplot(exonic_var_data, aes(y=number, x=gene, fill=factor(annotation, levels = ordered_annotations))) + 
    geom_bar(position="stack", stat="identity", width=0.8) +
    geom_text(data=total_num, aes(label=n, y=n, x=gene, fill=NULL), vjust=-0.25, size=3.5) +
    labs(x=paste(gene_family, 'genes', sep=' '), y='Number of exonic variants in canonical transcript',
         fill='Variant annotation') +
    theme_classic() + 
    scale_fill_manual(values=exonic_vars_anno_colors, labels=c('Other', 'Stop gained', 'Frameshift', 'Synonymous', 'Missense')) +
    scale_y_continuous(limits = c(0,740), expand = c(0,0)) +
    theme(legend.text = element_text(size = 11),
          legend.title = element_text(size =12, face='bold'),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 10, face='italic'),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = (11.5), face='bold'))
  
 
  
  return(p)
}

p1 <- barplot_exonic_variants_anno('UGT1') + theme(legend.position="none")
p2 <- barplot_exonic_variants_anno('UGT2') + theme(legend.position="none")
p3 <- barplot_exonic_variants_anno('UGT3') + theme(legend.position="none")
p4 <- barplot_exonic_variants_anno('UGT8')

plot_grid(p1, p2, p3, p4, ncol=4, rel_widths = c(1,1.1, 0.34, 0.6))
ggsave(filename=paste0('plots/01_Data_Processing/Exonic_variants_genes_anno.pdf'), width = 15, height = 5.5)




# _____________________________________________________________________________
#  1.1.3 Variant quantification in each gene family
# _____________________________________________________________________________

################# Number of variants in each location category #################

## Create stacked barplot for each gene family
var_data <- data.frame(matrix(ncol = 3))
colnames(var_data) <- c('gene_family', 'location', 'number')

locations <- c('5\' upstream', '5\'UTR',  'Intron 6-7', 'Intron 5-6',  'Intron 4-5', 'Intron 3-4', 'Intron 2-3', 'Intron 1-2', 
                       '5\'UTR, Intron 1-2', 'Intron 1-2, 5\'UTR', 
                       'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', 'Exon 7', '3\'UTR')

total_num <- list()
for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  
    UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants')))
    data <- data.frame(matrix(ncol = 3))
    colnames(data) <- c('gene_family', 'location', 'number')
    
    for (i in 1:length(locations)){
      data[i,'gene_family'] <- gene_family
      data[i,'location'] <- locations[i]
      data[i,'number'] <- table(UGT_variants$Location_in_txs)[locations[i]]
    }
  var_data <- rbind(var_data, data)
  total_num[[gene_family]] <- dim(UGT_variants)[1]
}
  
var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
var_data$gene_family <- factor(var_data$gene_family, levels = unique(var_data$gene_family))
var_data <- var_data[-1,]
total_num <- as.data.frame(melt(total_num))
colnames(total_num) <- c('n', 'gene_family')

p1 <- ggplot(var_data, aes(fill=factor(location, levels=locations), y=number, x=gene_family)) + 
  geom_bar(position="stack", stat="identity", width=0.8) + 
  labs(x='UGT gene family', y='Number of variants in canonical transcripts of genes', fill='Location') +
  theme_classic() + 
  geom_text(data=total_num, aes(label=n, y=n+30, x=gene_family, fill=NULL), vjust=-0.25, size=3.4) +
  scale_y_continuous(limits = c(0,7600), expand = c(0,0)) +
  scale_fill_manual(values = var_colors) +
  theme(plot.margin = margin(10, 2, 10, 2),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=11), 
        legend.title = element_text(size =12, face='bold'),
        axis.title = element_text(size = (11.5), face='bold'))


####################### Number of total/exonic variants #######################

## Plot total number of variants and number of exonic variants per gene family
data <- data.frame(matrix(ncol = 3))
colnames(data) <- c('gene_family',  'category' ,'number')

i=1
for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  UGT_variants <- eval(parse_expr(paste0(gene_family, '_variants')))
  ## Exonic variants only
  UGT_exonic_variants <- UGT_variants[grep('Exon', UGT_variants$Location_in_txs),]
  
  data[i, ] <- c(gene_family, 'Non-exonic variants', dim(UGT_variants)[1]-dim(UGT_exonic_variants)[1])
  data[i+1, ] <- c(gene_family, 'Exonic variants', dim(UGT_exonic_variants)[1])
  i=i+2
}
data$number <- as.numeric(data$number)
data$category <- factor(data$category, levels = c('Exonic variants', 'Non-exonic variants'))

p2 <- ggplot(data, aes(y=number, x=gene_family, fill=category)) + 
  geom_bar(stat="identity", position = position_stack(reverse = TRUE), width = 0.8) + 
  labs(x='UGT gene family', y='Number of variants in canonical transcripts of genes', fill='Category') +
  theme_classic() + 
  geom_text(data=total_num, aes(label=n, x=gene_family, y=n+150, fill=NULL), size=3.5) +
  geom_text(data=subset(data, category=='Exonic variants'), aes(label=number, x=gene_family, y=number/2, fill=NULL), size=3.5) +
  scale_y_continuous(limits = c(0,7600), expand = c(0,0)) +
  scale_fill_manual(values = var_colors) +
  theme(plot.margin = margin(10, 2, 10, 2),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=11), 
        legend.title = element_text(size =12, face='bold'),
        axis.title = element_text(size = (11.5), face='bold'))


################# Number of exonic variants of each annotation #################

## Function to create table of exonic variants and their location in all genes of each family
create_exonic_gene_fam_table <- function(gene_family){
  
  ## Genes of gene family
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  
  ## Unique exonic variants in gene family
  unique_UGT_exonic_variants <- unique(unlist(sapply(paste0(genes, '_exonic_data$Variant_ID'), function(x){eval(parse_expr(x))})))

  UGT_exonic_variants <- data.frame(matrix(ncol = length(genes)+4, nrow = 0))
  colnames(UGT_exonic_variants) <- c(genes,  'Position', 'VEP_Annotation', 'Location_in_txs', 'MAF')
  
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
    
    ## GMAF for variant in each gene 
    variant_maf <- unlist(sapply(genes, function(gene){
      eval(parse_expr(paste0(gene, '_exonic_data')))[which(eval(parse_expr(paste0(gene, '_exonic_data$Variant_ID'))) == variant), 'Allele_Frequency']
    }))
    if (length(unique(variant_maf))==1){
      UGT_exonic_variants$MAF[i] <- unique(variant_maf)
    }
    else{
      UGT_exonic_variants$MAF[i] <- toString(variant_maf)
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
  save(UGT_exonic_variants, file = paste0('processed-data/01_Data_Processing/', gene_family, '_exonic_variants.Rdata'))
}

create_exonic_gene_fam_table('UGT1')
create_exonic_gene_fam_table('UGT2')
create_exonic_gene_fam_table('UGT3')
create_exonic_gene_fam_table('UGT8')


## Stacked barplot
total_num <- list()
var_data <- data.frame(matrix(ncol = 3))
colnames(var_data) <- c('gene_family', 'annotation', 'number')

for(gene_family in c('UGT1', 'UGT2', 'UGT3', 'UGT8')){
  
  ordered_annotations <- c("other", "stop_gained", "frameshift_variant", "synonymous_variant", "missense_variant")
  UGT_exonic_variants <- eval(parse_expr(paste0(gene_family, '_exonic_variants')))
  UGT_exonic_variants$anno <- sapply(UGT_exonic_variants$VEP_Annotation, function(x){if(x %in% ordered_annotations){x}else{'other'}})
  
  data <- data.frame(matrix(ncol = 3))
  colnames(data) <- c('gene_family', 'annotation', 'number')
  
  for (i in 1:length(ordered_annotations)){
    data[i,'gene_family'] <- gene_family
    data[i,'annotation'] <- ordered_annotations[i]
    data[i,'number'] <- table(UGT_exonic_variants$anno)[ordered_annotations[i]]
  }
  var_data <- rbind(var_data, data)
  total_num[[gene_family]] <- dim(UGT_exonic_variants)[1]
}

var_data$number <- replace(var_data$number, which(is.na(var_data$number)), 0)
var_data$gene_family <- factor(var_data$gene_family, levels = unique(var_data$gene_family))
var_data$annotation <- factor(var_data$annotation, levels = ordered_annotations)
var_data <- var_data[-1,]
total_num <- as.data.frame(melt(total_num))
colnames(total_num) <- c('n', 'gene')

p3 <- ggplot(var_data, aes(fill=factor(annotation, levels=ordered_annotations), y=number, x=gene_family)) + 
  geom_bar(position="stack", stat="identity", width=0.8) + 
  labs(x='UGT gene family', y='Number of exonic variants in canonical transcripts of genes', fill='Variant annotation') +
  theme_classic() + 
  geom_text(data=total_num, aes(label=n, y=n+35, x=gene, fill=NULL), vjust=-0.25, size=3.5) +
  scale_fill_manual(values=exonic_vars_anno_colors, labels=c('Other', 'Stop gained', 'Frameshift', 'Synonymous', 'Missense')) +
  scale_y_continuous(limits = c(0,7600), expand = c(0,0)) +
  theme(plot.margin = margin(10, 2, 10, 2),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=11), 
        legend.title = element_text(size =12, face='bold'),
        axis.title = element_text(size = (11.5), face='bold'))

## Plots for gene families
plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.9, 1, 0.89))
ggsave(filename=paste0('plots/01_Data_Processing/Num_variants_per_gene_fam.pdf'), width = 17, height = 5.5)




# _____________________________________________________________________________
#  1.1.4 Examination of the minor allele frequency of all UGT exonic variants
# _____________________________________________________________________________

## GMAF of all exonic variants
all_exonic_vars <- rbind(UGT1_exonic_variants[c("VEP_Annotation", "MAF")], UGT2_exonic_variants[c("VEP_Annotation", "MAF")], UGT3_exonic_variants[c("VEP_Annotation", "MAF")], UGT8_exonic_variants[c("VEP_Annotation", "MAF")])
all_UGT_vars_MAF <- all_exonic_vars$MAF

## Number of total exonic variants
dim(all_exonic_vars)[1]
# [1] 9666

## Number of exonic variants that are rare (MAF≤0.01)
rare_exonic_vars <- subset(all_exonic_vars, MAF<=0.01)
dim(rare_exonic_vars)[1]
# [1] 9557

## % of rare exonic variants that are missense 
table(rare_exonic_vars$VEP_Annotation)['missense_variant'] / dim(rare_exonic_vars)[1] *100
#  missense_variant 
#          65.85749 

## Number of exonic variants that are singletons (MAF≤0.00001)
dim(subset(all_exonic_vars, MAF<=0.00001))[1] 
# [1] 5777

## Quantiles
quantiles <- as.data.frame(quantile(all_UGT_vars_MAF, probs=seq(from=0, to=1, length=7000)))
quantiles$Quantile <- rownames(quantiles)
colnames(quantiles) <- c('MAF', 'Quantile')
quantiles$Quantile <- as.numeric(gsub('%', '', quantiles$Quantile))

## Quantile for MAF≤0.00001 (Singletons)
singletons_q <- max(subset(quantiles, MAF<=0.00001)$Quantile)

## Quantile for 0.001<MAF≤0.01 (Rare)
rare_q <- max(subset(quantiles, MAF<=0.001)$Quantile)

## Quantile for MAF>0.01 (Common)
common_q <- min(subset(quantiles, MAF>0.01)$Quantile)

## Labels for such points
quantiles$label <- sapply(quantiles$Quantile, function(x){if(x==singletons_q | x==rare_q | x==common_q){signif(x, digits=3)} else{NA}})

## Categorize based on log10MAF
quantiles$maf_cat <- sapply(log10(quantiles$MAF), function(x){if(x<=(-5)){'Singleton'} 
                                         else if(x>(-5) & x<=(-3)){'Very rare'}
                                         else if(x>(-3) & x<=(-2)){'Rare'}
                                         else {'Common'}})
quantiles$maf_cat <- factor(quantiles$maf_cat, levels=c('Singleton', 'Very rare', 'Rare', 'Common'))

ggplot(data = quantiles, aes(x = Quantile, ymin = 0, ymax = MAF, fill = maf_cat)) +
  geom_ribbon(alpha=0.75) +
  scale_y_continuous(trans='log10', breaks=c(1e+00, 1e-01, 1e-02,  1e-03, 1e-04, 1e-05, 1e-06),
                     labels=c('0', '-1', '-2', '-3', '-4', '-5', '-6'), expand = c(0, 0)) +
  scale_x_continuous(breaks=c(0,25,50,75,100),
                     labels=c('0%', '25%', '50%', '75%', '100%')) +
  scale_fill_manual(values=c('Singleton'='honeydew2', 'Very rare'='lightcyan3', 'Rare'='lightblue4', 'Common'='midnightblue')) +
  labs(x='Percentage of variants', y='log10(MAF) of exonic variants', fill='Variant frequency') +
  theme_classic() +  
  geom_line(aes(y = MAF)) +
  ## Singletons
  geom_point(aes(x=singletons_q, y=1e-05, fill=NULL), color='tomato4', size=1, show.legend = FALSE) + 
  geom_hline(yintercept = 1e-05, color = 'indianred3', linetype='dashed', linewidth=0.6) +
  ## Rare variants
  geom_point(aes(x=rare_q, y=1e-03, fill=NULL), color='tomato4', size=1, show.legend = FALSE) + 
  geom_hline(yintercept = 1e-03, color = 'indianred3', linetype='dashed', linewidth=0.6) +
  ## Common variants
  geom_point(aes(x=common_q, y=1e-02, fill=NULL), color='tomato4', size=1, show.legend = FALSE) + 
  geom_hline(yintercept = 1e-02, color = 'indianred3', linetype='dashed', linewidth=0.6) +
  geom_text_repel(aes(x=label, y=MAF, fill=NULL, label=paste0(label ,'%')), 
                  size=3, color='grey40', min.segment.length = unit(0, 'lines'), hjust=1.5, 
                  box.padding = 0.3, lineheight=unit(2, 'lines'), nudge_x=0, nudge_y=0.2, 
                  fontface='bold') + 
  geom_text(aes(label='Singleton', x=20, y=(1e-05)/3*2), size=2.9, color='gray30') +
  geom_text(aes(label='Very rare', x=20, y=1e-04), size=2.9, color='gray30') +
  geom_text(aes(label='Rare', x=20, y=0.003), size=2.9, color='gray30') +
  geom_text(aes(label='Common', x=20, y=1e-01), size=2.9, color='gray30') +
  theme(legend.key.size = unit(0.35, units = 'cm'),
        axis.text = element_text(size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size =10, face='bold'),
        axis.title = element_text(size = (11), face='bold'))

ggsave(filename=paste0('plots/01_Data_Processing/MAF_all_vars.pdf'), width = 6, height = 4.2)






################################################################################
##    1.2  Integration of punctual variants of pharmacogenomic relevance
################################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
## Add missing deleterious promotor variants rs34983651 for UGT1A1, including the regulatory variant '2-234668879-C-CAT' (UGT1A1*28)                     # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
#                                                                                                                                            # |      
## Identify variants in datasets                                                                                                             # |
unlist(sapply(UGT1_genes, function(gene){if('rs34983651' %in% eval(parse_expr(paste0(gene,'_data$rsIDs')))){gene}                            # |
                                         else {NULL}}))                                                                                      # |
#  UGT1A8                                                                                                                                    # |
# "UGT1A8"                                                                                                                                   # |
                                                                                                                                             # | 
## Extract variant info from such gene dataset                                                                                               # | 
var_data <- UGT1A8_data[UGT1A8_data$rsIDs=='rs34983651',]                                                                                    # | 
#                                                                                                                                            # | 
#   #---------------------------------------------------------------------------------#                                                      # |
#   #                             !!! Warning !!!                                     #                                                      # |
#   #  Consider HGVS, transcript and protein consequence, as well as VEP annotation   #                                                      # |
#   #          can differ between genes based on their tx boundaries.                 #                                                      # |
#   #                                                                                 #                                                      # |
#   #---------------------------------------------------------------------------------#                                                      # |
#                                                                                                                                            # | 
## Change Transcript and VEP Annotation data from UGT1A8 to A1                                                                               # |
var_data$Transcript <- canonical_UGT1_txs[['UGT1A1']]                                                                                        # |
var_data$Canonical_txs <- TRUE                                                                                                               # |
var_data$VEP_Annotation <- '5\' upstream'                                                                                                    # |
## Add location in UGT1A1 tx of the variants (5' upstream)                                                                                   # |
var_data$Location_in_txs <- sapply(var_data$Position, function(x){location_determination(x,canonical_UGT1_txs[['UGT1A1']], NULL)[[1]]})      # |
                                                                                                                                             # | 
## Add variant data to UGT1A1 data                                                                                                           # |
UGT1A1_data <- rbind(UGT1A1_data, var_data[,-56])                                                                                            # |
## Remove 2 variants without functional evidence in ClinVar and a 3rd that is benign                                                         # |
UGT1A1_data <- UGT1A1_data[-which(UGT1A1_data$Variant_ID=='2-234668879-C-CATATAT' | UGT1A1_data$Variant_ID=='2-234668879-C-CATATATAT' |      # |    
                                  UGT1A1_data$Variant_ID=='2-234668879-CAT-C'), ]                                                            # |
save(UGT1A1_data, file = 'processed-data/01_Data_Processing/UGT1A1_data.Rdata')                                                              # |
                                                                                                                                             # | 
## Add variant data to UGT1A1 canonical and exonic dataset (though these are not exonic variants)                                            # |
UGT1A1_canonical_data <- rbind(UGT1A1_canonical_data, var_data)                                                                              # |
UGT1A1_canonical_data <- UGT1A1_canonical_data[-which(UGT1A1_canonical_data$Variant_ID=='2-234668879-C-CATATAT' |                            # |
                                                        UGT1A1_canonical_data$Variant_ID=='2-234668879-C-CATATATAT' |                        # |
                                                        UGT1A1_canonical_data$Variant_ID=='2-234668879-CAT-C'), ]                            # |
save(UGT1A1_canonical_data , file = 'processed-data/01_Data_Processing/UGT1A1_canonical_data.Rdata')                                         # |
                                                                                                                                             # |  
UGT1A1_exonic_data <- rbind(UGT1A1_exonic_data, var_data)                                                                                    # |
UGT1A1_exonic_data <- UGT1A1_exonic_data[-which(UGT1A1_exonic_data$Variant_ID=='2-234668879-C-CATATAT' |                                     # |
                                                        UGT1A1_exonic_data$Variant_ID=='2-234668879-C-CATATATAT' |                           # |
                                                        UGT1A1_exonic_data$Variant_ID=='2-234668879-CAT-C'), ]                               # |
save(UGT1A1_exonic_data , file = 'processed-data/01_Data_Processing/UGT1A1_exonic_data.Rdata')                                               # |
#                                                                                                                                            # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234669144-G-A' (UGT1A1*6) is present in exonic UGT1A1 data                                               # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234669144-G-A' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
                                         else {NULL}}))                                                                               # |
#  UGT1A1                                                                                                                             # |
# "UGT1A1"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A1_exonic_data[UGT1A1_exonic_data$Variant_ID=='2-234669144-G-A', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234669144   rs4148323            G           A   ENST00000305208.5            p.Gly71Arg                 c.211G>A      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Examine missing regulatory variant '2-234637707-T-C' for UGT1A3                                                              # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                               # |                 
## Identify variant                                                                                                             # |
unlist(sapply(UGT1_genes, function(gene){if('2-234637707-T-C' %in% eval(parse_expr(paste0(gene,'_data$Variant_ID')))){gene}     # |
  else {NULL}}))                                                                                                                # |
#                                                                                                                               # | 
#  UGT1A8                                                                                                                       # |
# "UGT1A8"                                                                                                                      # |
var_data <- UGT1A8_data[UGT1A8_data$Variant_ID=='2-234637707-T-C',]                                                             # |
                                                                                                                                # |
## Variant not added as it has no functional evidence in ClinVar                                                                # |
var_data[,c('ClinVar_Clinical_Significance', 'ClinVar_Variation_ID')]                                                           # |
                                                                                                                                # | 
# ClinVar_Clinical_Significance ClinVar_Variation_ID                                                                            # |
# 2-234637707-T-C                                                 NA                                                            # |
                                                                                                                                # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234602191-A-G' (UGT1A6*2) is present in exonic UGT1A6 data                                               # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234602191-A-G' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
                                         else {NULL}}))                                                                               # |
#  UGT1A6                                                                                                                             # |
# "UGT1A6"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A6_exonic_data[UGT1A6_exonic_data$Variant_ID=='2-234602191-A-G', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234602191   rs2070959            A           G   ENST00000305139.6           p.Thr181Ala                 c.541A>G      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234602202-A-C' (UGT1A6*2) is present in exonic UGT1A6 data                                               # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234602202-A-C' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
  else {NULL}}))                                                                               # |
#  UGT1A6                                                                                                                             # |
# "UGT1A6"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A6_exonic_data[UGT1A6_exonic_data$Variant_ID=='2-234602202-A-C', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234602202   rs1105879            A           C   ENST00000305139.6           p.Arg184Ser                 c.552A>C      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234601669-T-G' is present in exonic UGT1A6 data                                                          # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234601669-T-G' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
                                         else {NULL}}))                                                                               # |
#  UGT1A6                                                                                                                             # |
# "UGT1A6"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A6_exonic_data[UGT1A6_exonic_data$Variant_ID=='2-234601669-T-G', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234601669   rs6759892            T           G   ENST00000305139.6             p.Ser7Ala                  c.19T>G      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234591205-T-C' is present in exonic UGT1A7 data                                                          # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234591205-T-C' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
                                         else {NULL}}))                                                                               # |
#  UGT1A7                                                                                                                             # |
# "UGT1A7"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A7_exonic_data[UGT1A7_exonic_data$Variant_ID=='2-234591205-T-C', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234591205   rs11692021           T           C   ENST00000373426.3           p.Trp208Arg                 c.622T>C      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '2-234526871-C-G' is present in exonic UGT1A8 data                                                          # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT1_genes, function(gene){if('2-234526871-C-G' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}    # |
                                         else {NULL}}))                                                                               # |
#  UGT1A8                                                                                                                             # |
# "UGT1A8"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT1A8_exonic_data[UGT1A8_exonic_data$Variant_ID=='2-234526871-C-G', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',   # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position       rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence      # |
#          2   234526871   rs1042597            C           G   ENST00000373450.4           p.Ala173Gly                 c.518C>G      # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '4-70346565-A-T' is present in exonic UGT2B4 data                                                           # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT2_genes, function(gene){if('4-70346565-A-T' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}     # |
                                        else {NULL}}))                                                                                # |
#  UGT2B4                                                                                                                             # |
# "UGT2B4"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT2B4_exonic_data[UGT2B4_exonic_data$Variant_ID=='4-70346565-A-T', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',    # |
                                                                       'Transcript', 'Protein_Consequence', 'Transcript_Consequence', # |
                                                                       'VEP_Annotation', 'Location_in_txs')]                          # | 
#                                                                                                                                     # |
# Chromosome    Position        rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence     # |
#          4    70346565   rs13119049            A           T   ENST00000305107.6           p.Asp458Glu                c.1374T>A     # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 6                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '4-69964338-T-C' is present in exonic UGT2B7 data                                                           # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT2_genes, function(gene){if('4-69964338-T-C' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}     # |
                                         else {NULL}}))                                                                               # |
#  UGT2B7                                                                                                                             # |
# "UGT2B7"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT2B7_exonic_data[UGT2B7_exonic_data$Variant_ID=='4-69964338-T-C', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',    # |
                                                                      'Transcript', 'Protein_Consequence', 'Transcript_Consequence',  # |
                                                                      'VEP_Annotation', 'Location_in_txs')]                           # | 
#                                                                                                                                     # |
# Chromosome    Position        rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence     # |
#          4    69964338    rs7439366            T           C   ENST00000305231.7           p.Tyr268His                 c.802T>C     # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 2                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check LoF variant '4-69962449-G-T' is present in exonic UGT2B7 data                                                                # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT2_genes, function(gene){if('4-69962449-G-T' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}     # |
                                         else {NULL}}))                                                                               # |
#  UGT2B7                                                                                                                             # |
# "UGT2B7"                                                                                                                            # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT2B7_exonic_data[UGT2B7_exonic_data$Variant_ID=='4-69962449-G-T', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',    # |
                                                                      'Transcript', 'Protein_Consequence', 'Transcript_Consequence',  # |
                                                                      'VEP_Annotation', 'Location_in_txs')]                           # | 
#                                                                                                                                     # |
# Chromosome    Position        rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence     # |
#          4    69962449    rs12233719           G           T   ENST00000305231.7            p.Ala71Ser                 c.211G>T     # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## Check missense variant '4-69536084-A-C' is present in exonic UGT2B15 data                                                          # |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#                                                                                                                                     # |               
## Identify variant in datasets                                                                                                       # |
unlist(sapply(UGT2_genes, function(gene){if('4-69536084-A-C' %in% eval(parse_expr(paste0(gene,'_exonic_data$Variant_ID')))){gene}     # |
                                         else {NULL}}))                                                                               # |
#  UGT2B15                                                                                                                            # |
# "UGT2B15"                                                                                                                           # |
#                                                                                                                                     # | 
## Explore data                                                                                                                       # | 
UGT2B15_exonic_data[UGT2B15_exonic_data$Variant_ID=='4-69536084-A-C', c('Chromosome', 'Position', 'rsIDs', 'Reference', 'Alternate',  # |
                                                                      'Transcript', 'Protein_Consequence', 'Transcript_Consequence',  # |
                                                                      'VEP_Annotation', 'Location_in_txs')]                           # | 
#                                                                                                                                     # |
# Chromosome    Position        rsIDs    Reference   Alternate          Transcript   Protein_Consequence   Transcript_Consequence     # |
#          4    69536084    rs1902023            A           C   ENST00000338206.5            p.Tyr85Asp                 c.253T>G     # |
#    VEP_Annotation    Location_in_txs                                                                                                # |
#  missense_variant             Exon 1                                                                                                # |
#                                                                                                                                     # | 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/Stockholm
# date     2023-10-12
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# bit           4.0.5   2022-11-15 [1] CRAN (R 4.3.0)
# bit64         4.0.5   2020-08-30 [1] CRAN (R 4.3.0)
# cli           3.6.1   2023-03-23 [1] CRAN (R 4.3.0)
# colorspace    2.1-0   2023-01-23 [1] CRAN (R 4.3.0)
# corrplot    * 0.92    2021-11-18 [1] CRAN (R 4.3.0)
# cowplot     * 1.1.1   2020-12-30 [1] CRAN (R 4.3.0)
# crayon        1.5.2   2022-09-29 [1] CRAN (R 4.3.0)
# dplyr         1.1.2   2023-04-20 [1] CRAN (R 4.3.0)
# fansi         1.0.4   2023-01-22 [1] CRAN (R 4.3.0)
# farver        2.1.1   2022-07-06 [1] CRAN (R 4.3.0)
# generics      0.1.3   2022-07-05 [1] CRAN (R 4.3.0)
# ggplot2     * 3.4.2   2023-04-03 [1] CRAN (R 4.3.0)
# ggrepel     * 0.9.3   2023-02-03 [1] CRAN (R 4.3.0)
# glue          1.6.2   2022-02-24 [1] CRAN (R 4.3.0)
# gtable        0.3.3   2023-03-21 [1] CRAN (R 4.3.0)
# here        * 1.0.1   2020-12-13 [1] CRAN (R 4.3.0)
# hms           1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
# labeling      0.4.2   2020-10-20 [1] CRAN (R 4.3.0)
# lattice       0.21-8  2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.3.0)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
# Matrix        1.6-0   2023-07-08 [1] CRAN (R 4.3.0)
# mgcv          1.9-0   2023-07-11 [1] CRAN (R 4.3.0)
# munsell       0.5.0   2018-06-12 [1] CRAN (R 4.3.0)
# nlme          3.1-162 2023-01-31 [1] CRAN (R 4.3.0)
# pillar        1.9.0   2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.3.0)
# plyr          1.8.8   2022-11-11 [1] CRAN (R 4.3.0)
# pROC        * 1.18.4  2023-07-06 [1] CRAN (R 4.3.0)
# R6            2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
# ragg          1.2.5   2023-01-12 [1] CRAN (R 4.3.0)
# Rcpp          1.0.11  2023-07-06 [1] CRAN (R 4.3.0)
# readr       * 2.1.4   2023-02-10 [1] CRAN (R 4.3.0)
# reshape2    * 1.4.4   2020-04-09 [1] CRAN (R 4.3.0)
# rlang       * 1.1.1   2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.3.0)
# rstudioapi    0.15.0  2023-07-07 [1] CRAN (R 4.3.0)
# scales        1.2.1   2022-08-20 [1] CRAN (R 4.3.0)
# sessioninfo * 1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
# stringi       1.7.12  2023-01-11 [1] CRAN (R 4.3.0)
# stringr       1.5.0   2022-12-02 [1] CRAN (R 4.3.0)
# systemfonts   1.0.4   2022-02-11 [1] CRAN (R 4.3.0)
# textshaping   0.3.6   2021-10-13 [1] CRAN (R 4.3.0)
# tibble        3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.3.0)
# tzdb          0.4.0   2023-05-12 [1] CRAN (R 4.3.0)
# utf8          1.2.3   2023-01-31 [1] CRAN (R 4.3.0)
# vctrs         0.6.3   2023-06-14 [1] CRAN (R 4.3.0)
# vroom         1.6.3   2023-04-28 [1] CRAN (R 4.3.0)
# withr         2.5.0   2022-03-03 [1] CRAN (R 4.3.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────