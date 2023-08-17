

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


## For each variant, determine in what UGT1 genes they are present 

## Unique variants in UGT1 genes
unique_UGT1_variants <- unique(c(UGT1A1_data$Variant_ID, UGT1A3_data$Variant_ID, UGT1A4_data$Variant_ID, UGT1A5_data$Variant_ID, UGT1A6_data$Variant_ID, UGT1A7_data$Variant_ID, UGT1A8_data$Variant_ID, UGT1A9_data$Variant_ID, UGT1A10_data$Variant_ID))

## Compute the number of times each variant appears in each gene 
## (if present in a gene, variants are expected to appear only once since they are all unique in each gene)
UGT1_variants_counts <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(UGT1_variants_counts) <- c(UGT1_genes, 'VEP_Annotation')

i=1
for (variant in unique_UGT1_variants){
  ## Search variant in each gene
  variant_counts <- sapply(UGT1_genes, function(gene){
                      length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
                    })
  ## VEP annotation for variant in each gene
  variant_anno <- unlist(sapply(UGT1_genes, function(gene){
                           eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
                         }))
  ## Evaluate if the anno of a shared variable is the same across all genes
  UGT1_variants_counts[i, 1:9] <- variant_counts
  if (length(unique(variant_anno))==1){
    UGT1_variants_counts$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT1_variants_counts$VEP_Annotation[i] <- toString(variant_anno)
  }
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
UGT2_variants_counts <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(UGT2_variants_counts) <- c(UGT2_genes, 'VEP_Annotation')

i=1
for (variant in unique_UGT2_variants){
  variant_counts <- sapply(UGT2_genes, function(gene){
    length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
  })
  
  variant_anno <- unlist(sapply(UGT2_genes, function(gene){
    eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
  }))
  
  UGT2_variants_counts[i, 1:10] <- variant_counts
  if (length(unique(variant_anno))==1){
    UGT2_variants_counts$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT2_variants_counts$VEP_Annotation[i] <- toString(variant_anno)
  }
  
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
UGT3_variants_counts <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(UGT3_variants_counts) <- c(UGT3_genes, 'VEP_Annotation')

i=1
for (variant in unique_UGT3_variants){
  variant_counts <- sapply(UGT3_genes, function(gene){
    length(which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) %in% variant))
  })
  
  variant_anno <- unlist(sapply(UGT3_genes, function(gene){
    eval(parse_expr(paste0(gene, '_data')))[which(eval(parse_expr(paste0(gene, '_data$Variant_ID'))) == variant), 'VEP_Annotation']
  }))
  
  UGT3_variants_counts[i, 1:2] <- variant_counts
  if (length(unique(variant_anno))==1){
    UGT3_variants_counts$VEP_Annotation[i] <- unique(variant_anno)
  }
  else{
    UGT3_variants_counts$VEP_Annotation[i] <- toString(variant_anno)
  }
  
  i=i+1
}

rownames(UGT3_variants_counts) <- unique_UGT3_variants


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
table(UGT1_variants_counts[which(apply(UGT1_variants_counts[,-10], 1, sum)==9), 'VEP_Annotation'])
# 3_prime_UTR_variant      frameshift_variant        inframe_deletion          intron_variant        missense_variant splice_acceptor_variant    splice_donor_variant 
#                  10                       8                       2                      94                     165                       2                       1 
# splice_region_variant             stop_gained      synonymous_variant 
#                    10                       6                      76 


###########  Annotation for variants shared between 4 UGT1 genes  ###########
table(UGT1_variants_counts[which(apply(UGT1_variants_counts[,-10], 1, sum)==4), 'VEP_Annotation'])
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
shared_variants <- rownames(UGT1_variants_counts)[which(apply(UGT1_variants_counts[,-10], 1, sum)==4)]
table(apply(sapply(shared_variants, function(x){colnames(UGT1_variants_counts)[which(UGT1_variants_counts[x,]==1)]}), 2, toString))
# UGT1A1, UGT1A4, UGT1A6, UGT1A10 
#                              45 


###########  Annotation for variants shared between 2 UGT1 genes  ###########
table(UGT1_variants_counts[which(apply(UGT1_variants_counts[,-10], 1, sum)==2), 'VEP_Annotation'])
# 5_prime_UTR_variant                  frameshift_variant                    inframe_deletion                      intron_variant 
#                   7                                  55                                   7                                  75 
# intron_variant, 5_prime_UTR_variant                    missense_variant                splice_donor_variant               splice_region_variant 
#                                 24                                 773                                   4                                  14 
# start_lost                         stop_gained                  synonymous_variant 
#          6                                  32                                 331 

## Which 2 genes share those variants?
shared_variants <- rownames(UGT1_variants_counts)[which(apply(UGT1_variants_counts[,-10], 1, sum)==2)]
table(apply(sapply(shared_variants, function(x){colnames(UGT1_variants_counts)[which(UGT1_variants_counts[x,]==1)]}), 2, toString))
# UGT1A1, UGT1A3     UGT1A1, UGT1A5     UGT1A1, UGT1A8     UGT1A1, UGT1A9 
#            392                317                279                340 
 


##############################
####  UGT2 shared variants 
##############################

###########  Annotation for variants shared between 2 UGT2 genes  ###########
table(UGT2_variants_counts[which(apply(UGT2_variants_counts[,-11], 1, sum)==2), 'VEP_Annotation'])
# 3_prime_UTR_variant      frameshift_variant        inframe_deletion          intron_variant        missense_variant splice_acceptor_variant    splice_donor_variant 
#                  10                      18                       1                     122                     217                       3                       3 
# splice_region_variant             stop_gained   stop_retained_variant      synonymous_variant 
#                    23                      15                       1                      69 

## Which 2 genes share those variants?
shared_variants <- rownames(UGT2_variants_counts)[which(apply(UGT2_variants_counts[,-11], 1, sum)==2)]
table(apply(sapply(shared_variants, function(x){colnames(UGT2_variants_counts)[which(UGT2_variants_counts[x,]==1)]}), 2, toString))
# UGT2A1, UGT2A2 
#            482



##############################
####  UGT3 shared variants 
##############################

## No shared variants between UGT3A1 and UGT3A2
which(apply(UGT3_variants_counts[,-3], 1, sum)==2)
#  named integer(0)



##############################
####  UGT8 shared variants 
##############################

## Only one gene for UGT8























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
