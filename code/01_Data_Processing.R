

library(here)
library(readr)
library(rlang)

################################################################################
##                    01. Data Processing and Correction
################################################################################

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



##############################################
##  1.1 Correct allele number for all genes???????
##############################################

## Detect overlapping variants between different genes of the same family

############################
####      UGT1 genes 
############################
UGT1_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10')

## Confirm all variants in UGT1 genes are in the same chr
sapply(UGT1_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT1A1  UGT1A3  UGT1A4  UGT1A5  UGT1A6  UGT1A7  UGT1A8  UGT1A9 UGT1A10 
#      2       2       2       2       2       2       2       2       2 


## For each variant, determine in what UGT1 genes they are present 

## Unique variants in UGT1 genes
unique_UGT1_variants <- unique(c(UGT1A1_data$Variant_ID, UGT1A3_data$Variant_ID, UGT1A4_data$Variant_ID, UGT1A5_data$Variant_ID, UGT1A6_data$Variant_ID, UGT1A7_data$Variant_ID, UGT1A8_data$Variant_ID, UGT1A9_data$Variant_ID, UGT1A10_data$Variant_ID))

## Compute the number of times each variant appears in each gene
UGT1_variants_counts <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(UGT1_variants_counts) <- UGT1_genes

i=1
for (variant in unique_UGT1_variants){
  variant_counts <- sapply(UGT1_genes, function(gene){
                      length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
                    })
  
  UGT1_variants_counts[i,] <- variant_counts
  i=i+1
}

rownames(UGT1_variants_counts) <- unique_UGT1_variants




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
UGT2_variants_counts <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(UGT2_variants_counts) <- UGT2_genes

i=1
for (variant in unique_UGT2_variants){
  variant_counts <- sapply(UGT2_genes, function(gene){
    length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
  })
  
  UGT2_variants_counts[i,] <- variant_counts
  i=i+1
}

rownames(UGT2_variants_counts) <- unique_UGT2_variants



############################
####    UGT3 genes 
############################
UGT3_genes <- c('UGT3A1', 'UGT3A2')

## Confirm all UGT3 variants are in the same chr
sapply(UGT3_genes, function(x){unique(eval(parse_expr(paste0(x, '_data$Chromosome'))))})
# UGT3A1 UGT3A2 
#      5      5 

unique_UGT3_variants <- unique(c(UGT3A1_data$Variant_ID, UGT3A2_data$Variant_ID))
UGT3_variants_counts <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(UGT3_variants_counts) <- UGT3_genes

i=1
for (variant in unique_UGT3_variants){
  variant_counts <- sapply(UGT3_genes, function(gene){
    length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
  })
  
  UGT3_variants_counts[i,] <- variant_counts
  i=i+1
}

rownames(UGT3_variants_counts) <- unique_UGT3_variants




############################
####    UGT8 genes 
############################
unique(UGT8_data$Chromosome)
#    4 

## Variant_ID's are unique
unique_UGT8_variants <- UGT8_data$Variant_ID
## Number of times each unique variant appears in UGT8 gene
UGT8_variants_counts <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(UGT8_variants_counts) <- 'UGT8'

i=1
for (variant in unique_UGT8_variants){
  variant_counts <- length(which(UGT8_data$Variant_ID %in% variant))
  
  UGT8_variants_counts[i,] <- variant_counts
  i=i+1
}

rownames(UGT8_variants_counts) <- unique_UGT8_variants





## Verify there are not overlappging variants between genes from different subfamilies

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
