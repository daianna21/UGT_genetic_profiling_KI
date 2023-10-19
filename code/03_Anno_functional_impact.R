
library(here)
library(rlang)
library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(corrplot)
library(pROC)
library(ggrepel)
library(sessioninfo)


####################################################################################################
##                         3. Annotation of the functional impact of variants
####################################################################################################

gene_families <- c('UGT1', 'UGT2', 'UGT3', 'UGT8')
UGT_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10',
               'UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28',
               'UGT3A1', 'UGT3A2', 
               'UGT8')
UGT1_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10')
UGT2_genes <- c('UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28')
UGT3_genes <- c('UGT3A1', 'UGT3A2')
UGT8_genes <- c('UGT8')

## Canonical txs of genes 
canonical_txs <- list('UGT1A1'= 'ENST00000305208.5', 'UGT1A3'='ENST00000482026.1', 'UGT1A4'='ENST00000373409.3', 
                           'UGT1A5'='ENST00000373414.3', 'UGT1A6'='ENST00000305139.6', 'UGT1A7'='ENST00000373426.3', 
                           'UGT1A8'= 'ENST00000373450.4','UGT1A9'= 'ENST00000354728.4', 'UGT1A10'='ENST00000344644.5',
                           'UGT2A1'= 'ENST00000503640.1', 'UGT2A2'='ENST00000457664.2', 'UGT2A3'='ENST00000251566.4', 
                           'UGT2B4'='ENST00000305107.6', 'UGT2B7'='ENST00000305231.7', 'UGT2B10'='ENST00000265403.7', 
                           'UGT2B11'= 'ENST00000446444.1', 'UGT2B15'= 'ENST00000338206.5', 'UGT2B17'='ENST00000317746.2', 
                           'UGT2B28'='ENST00000335568.5', 'UGT3A1'= 'ENST00000274278.3', 'UGT3A2'='ENST00000282507.3',
                           'UGT8'= 'ENST00000310836.6')


## Load exonic data for each gene 
non_minor_alleles_allGenes <- vector()
for (gene in UGT_genes){
  exonic_vars <- eval(parse_expr(load(here(paste0('~/Desktop/UGT_genetic_profiling_KI/processed-data/01_Data_Processing/', gene, '_exonic_data.Rdata')),
       verbose=TRUE)))
  ## Variants with global allele freqs > 0.5 (not minor allele)
  non_minor_alleles <- exonic_vars[which(exonic_vars$Allele_Frequency>0.5), ]
  if (dim(non_minor_alleles)[1]!=0){
    non_minor_alleles <- cbind(gene, exonic_vars[which(exonic_vars$Allele_Frequency>0.5), ])
    non_minor_alleles_allGenes <- rbind(non_minor_alleles_allGenes, non_minor_alleles)
  }
  assign( paste0(gene, '_exonic_data'), exonic_vars)
}



################################################################################
##              3.1  Functional prediction of missense variants
################################################################################

## Subset to missense variants for all genes

for (gene in UGT_genes){
  UGT_exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  missense_vars <- UGT_exonic_data[which(UGT_exonic_data$VEP_Annotation=='missense_variant'),]
  assign(paste0(gene, '_missense_vars'), missense_vars)
}

# _______________________________________________________________________________
#  3.1.1  ANNOVAR input file preparation
# _______________________________________________________________________________

## Necessary columns per variant: Chromosome   Start    End    Reference    Alternate
## Start and End are the same for missense variants (single nt mutations)

## Annotate and predict missense variants

for (gene in UGT_genes){
  UGT_missense_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))
  
  ## Evaluate if reference and alternate alleles are single valid characters (no indels nor missing/null) for all missense variants in each gene
  print(paste0(gene, ': ', names(table(
    sapply(1:dim(UGT_missense_vars)[1], function(i){
    if (UGT_missense_vars$Reference[i] %in% c('A', 'T', 'C', 'G')){TRUE}
    else {FALSE}
  }))), ' ',
  
  names(table(
    sapply(1:dim(UGT_missense_vars)[1], function(i){
      if (UGT_missense_vars$Alternate[i] %in% c('A', 'T', 'C', 'G')){TRUE}
      else {FALSE}
    })))
  ))
  
}

# [1] "UGT1A1: TRUE TRUE"
# [1] "UGT1A3: TRUE TRUE"
# [1] "UGT1A4: TRUE TRUE"
# [1] "UGT1A5: TRUE TRUE"
# [1] "UGT1A6: TRUE TRUE"
# [1] "UGT1A7: TRUE TRUE"
# [1] "UGT1A8: TRUE TRUE"
# [1] "UGT1A9: TRUE TRUE"
# [1] "UGT1A10: TRUE TRUE"
# [1] "UGT2A1: TRUE TRUE"
# [1] "UGT2A2: TRUE TRUE"
# [1] "UGT2A3: TRUE TRUE"
# [1] "UGT2B4: TRUE TRUE"
# [1] "UGT2B7: TRUE TRUE"
# [1] "UGT2B10: TRUE TRUE"
# [1] "UGT2B11: TRUE TRUE"
# [1] "UGT2B15: TRUE TRUE"
# [1] "UGT2B17: TRUE TRUE"
# [1] "UGT2B28: TRUE TRUE"
# [1] "UGT3A1: TRUE TRUE"
# [1] "UGT3A2: TRUE TRUE"
# [1] "UGT8: TRUE TRUE"


## Sub datasets with specified columns for ANNOVAR

for (gene in UGT_genes){
  UGT_missense_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))
  missense_vars_ANNOVAR_format <- data.frame(matrix(ncol=5, nrow=nrow(UGT_missense_vars)))
  colnames(missense_vars_ANNOVAR_format) <- c('Chromosome', 'Start', 'End', 'Ref', 'Obs')
  
  missense_vars_ANNOVAR_format$Chromosome <- UGT_missense_vars$Chromosome
  missense_vars_ANNOVAR_format$Start <- missense_vars_ANNOVAR_format$End <- UGT_missense_vars$Position
  missense_vars_ANNOVAR_format$Ref <- UGT_missense_vars$Reference
  missense_vars_ANNOVAR_format$Obs <- UGT_missense_vars$Alternate
  
  assign(paste0(gene, '_missense_vars_ANNOVAR_format'), missense_vars_ANNOVAR_format)
  save(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/', 
                                                                                      gene, '_missense_vars_ANNOVAR_format.Rdata'))
  write.table(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/', 
                                                   gene, '_missense_vars_ANNOVAR_format.txt'), row.names = FALSE, col.names = FALSE, sep = '\t')
  write.table(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/', 
                                                          gene, '_missense_vars_ANNOVAR_format.csv'), row.names = FALSE, col.names = FALSE, sep = '\t')
}
 


# _______________________________________________________________________________
#  3.1.2  Examination of ANNOVAR gene-based annotation output 
# _______________________________________________________________________________

## Download ANNOVAR output for each gene
for (gene in UGT_genes){
  myanno <- read.csv(paste0('processed-data/03_Anno_functional_impact/ANNOVAR_output/myanno', gene, '.hg19_multianno.csv'))
  assign(paste0('myanno_', gene), myanno)
  
}

############################ 1. Confirm all missense variants are annotated as exonic  ############################# 

sapply(UGT_genes, function(gene){names(table(eval(parse_expr(paste0('myanno_', gene, '$Func.refGene')))))})

# $UGT1A1
# [1] "exonic"
# 
# $UGT1A3
# [1] "exonic"
# 
# $UGT1A4
# [1] "exonic"
# 
# $UGT1A5
# [1] "exonic"
# 
# $UGT1A6
# [1] "exonic"
# 
# $UGT1A7
# [1] "exonic"
# 
# $UGT1A8
# [1] "exonic"
# 
# $UGT1A9
# [1] "exonic"
# 
# $UGT1A10
# [1] "exonic"
# 
# $UGT2A1
# [1] "exonic"
# 
# $UGT2A2
# [1] "exonic"
# 
# $UGT2A3
# [1] "exonic"
# 
# $UGT2B4
# [1] "exonic"
# 
# $UGT2B7
# [1] "exonic"
# 
# $UGT2B10
# [1] "exonic"          "exonic;splicing"
# 
# $UGT2B11
# [1] "exonic"          "exonic;splicing"
# 
# $UGT2B15
# [1] "exonic"
# 
# $UGT2B17
# [1] "exonic"
# 
# $UGT2B28
# [1] "exonic"
# 
# $UGT3A1
# [1] "exonic"
# 
# $UGT3A2
# [1] "exonic"
# 
# $UGT8
# [1] "exonic"


# -----------------------------------------------------------  ?  -----------------------------------------------------------
#                                       Evaluate variants found as exonic and splicing                                      |
# -----------------------------------------------------------  ?  -----------------------------------------------------------

## In UGT2B10: var 4-69681891-G-A is within exon 1 of NM_001075 and NM_001144767 txs, not in the exon boundaries -> exonic
## This variant is in 5'-UTR of another UGT2B10 tx (NM_001290091)
myanno_UGT2B10[which(myanno_UGT2B10$Func.refGene=='exonic;splicing'), 1:10]
# Chr      Start        End   Ref   Alt      Func.refGene      Gene.refGene        GeneDetail.refGene   ExonicFunc.refGene
#   4   69681891   69681891     G     A   exonic;splicing   UGT2B10;UGT2B10   NM_001290091:exon1:UTR5    nonsynonymous SNV
#                                                                 AAChange.refGene
# UGT2B10:NM_001075:exon1:c.G154A:p.V52M,UGT2B10:NM_001144767:exon1:c.G154A:p.V52M

## Corroborate location of variant in NM_001075, the UGT2B10 canonical tx (ENST00000265403.7) and check Exon 1 boundaries 
## (Use function location_determination() in ../01_Data_Processing.R)
location_determination(69681891, canonical_UGT2_txs[['UGT2B10']], 'Exon 1')
# "Exon 1"
#
#           Start       End
# Exon 1 69681738  69682455

# - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## In UGT2B11: variants 4-70066297-G-C and 4-70066297-G-A are both within exon 6 of UGT2B11 NM_001073 tx
## but in exon 3 of LOC105377267 tx
myanno_UGT2B11[which(myanno_UGT2B11$Func.refGene=='exonic;splicing'), 1:10]
# Chr     Start       End  Ref   Alt     Func.refGene           Gene.refGene          GeneDetail.refGene   ExonicFunc.refGene                           AAChange.refGene
#   4  70066297  70066297    G    C   exonic;splicing   UGT2B11;LOC105377267  NR_136191:exon3:c.484+1G>C    nonsynonymous SNV   UGT2B11:NM_001073:exon6:c.C1451G:p.T484S
#   4  70066297  70066297    G    A   exonic;splicing   UGT2B11;LOC105377267  NR_136191:exon3:c.484+1G>A    nonsynonymous SNV   UGT2B11:NM_001073:exon6:c.C1451T:p.T484I

location_determination(70066297, canonical_UGT2_txs[['UGT2B11']], 'Exon 6')
# "Exon 6" 
#
#           Start      End
# Exon 6 70066437 70066158
# ---------------------------------------------------------------------------------------------------------------------------


################################ 2.  Check all missense variants are non-synonymous  ################################

sapply(UGT_genes, function(gene){names(table(eval(parse_expr(paste0('myanno_', gene, '$ExonicFunc.refGene')))))})

# $UGT1A1
# [1] "nonsynonymous SNV"
# 
# $UGT1A3
# [1] "nonsynonymous SNV"
# 
# $UGT1A4
# [1] "nonsynonymous SNV"
# 
# $UGT1A5
# [1] "nonsynonymous SNV"
# 
# $UGT1A6
# [1] "nonsynonymous SNV"
# 
# $UGT1A7
# [1] "nonsynonymous SNV"
# 
# $UGT1A8
# [1] "nonsynonymous SNV"
# 
# $UGT1A9
# [1] "nonsynonymous SNV"
# 
# $UGT1A10
# [1] "nonsynonymous SNV"
# 
# $UGT2A1
# [1] "nonsynonymous SNV"
# 
# $UGT2A2
# [1] "nonsynonymous SNV"
# 
# $UGT2A3
# [1] "nonsynonymous SNV"
# 
# $UGT2B4
# [1] "nonsynonymous SNV"
# 
# $UGT2B7
# [1] "nonsynonymous SNV"
# 
# $UGT2B10
# [1] "nonsynonymous SNV" "startloss"        
# 
# $UGT2B11
# [1] "nonsynonymous SNV"
# 
# $UGT2B15
# [1] "nonsynonymous SNV"
# 
# $UGT2B17
# [1] "nonsynonymous SNV"
# 
# $UGT2B28
# [1] "nonsynonymous SNV"
# 
# $UGT3A1
# [1] "nonsynonymous SNV"
# 
# $UGT3A2
# [1] "nonsynonymous SNV"
# 
# $UGT8
# [1] "nonsynonymous SNV"


# -----------------------------------------------------------  ?  -----------------------------------------------------------
#                                               Evaluate startloss variants:                                                |
# -----------------------------------------------------------  ?  -----------------------------------------------------------

## Variant 4-69683774-T-G is non-synonymous (aa change in AAChange.refGene) for NM_001075 and NM_001144767 txs of UGT2B10
## but is start lost for NM_001290091 tx (p.Met1? as protein consequence)
myanno_UGT2B10[which(myanno_UGT2B10$ExonicFunc.refGene=='startloss'), 1:10]
# Chr      Start        End   Ref  Alt   Func.refGene   Gene.refGene   GeneDetail.refGene ExonicFunc.refGene
#   4   69683774   69683774     T    G         exonic        UGT2B10                    .          startloss
#                                                                                                         AAChange.refGene
# UGT2B10:NM_001075:exon2:c.T746G:p.M249R,UGT2B10:NM_001144767:exon2:c.T494G:p.M165R,UGT2B10:NM_001290091:exon2:c.T2G:p.M1?

## Corroborate this variant is within exon 2 of canonical tx (NM_001075), not at the beginning of exon 1 (Met)
location_determination(69683774, canonical_UGT2_txs[['UGT2B10']], 'Exon 2')
# "Exon 2"
#
#           Start      End
# Exon 2 69683747 69683895
# ---------------------------------------------------------------------------------------------------------------------------


########################### 3.  Check in which genes the variants are present  #############################
sapply(UGT_genes, function(gene){names(table(eval(parse_expr(paste0('myanno_', gene, '$Gene.refGene')))))})
# $UGT1A1
# [1] "UGT1A1"   "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"
# 
# $UGT1A3
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A3"                                                         
# 
# $UGT1A4
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A4"                                                         
# 
# $UGT1A5
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A5"                                                         
# 
# $UGT1A6
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A6"                                                         
# 
# $UGT1A7
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A7"                                                         
# 
# $UGT1A8
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A8"                                                         
# 
# $UGT1A9
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A9"                                                         
# 
# $UGT1A10
# [1] "UGT1A1;UGT1A10;UGT1A3;UGT1A4;UGT1A5;UGT1A6;UGT1A7;UGT1A8;UGT1A9"   "UGT1A10"                                                        
# 
# $UGT2A1
# [1] "UGT2A1"    "UGT2A1;UGT2A2"
# 
# $UGT2A2
# [1] "UGT2A1;UGT2A2"   "UGT2A2"       
# 
# $UGT2A3
# [1] "UGT2A3"
# 
# $UGT2B4
# [1] "UGT2B4"
# 
# $UGT2B7
# [1] "UGT2B7"
# 
# $UGT2B10
# [1] "UGT2B10"    "UGT2B10;UGT2B10"
# 
# $UGT2B11
# [1] "UGT2B11"    "UGT2B11;LOC105377267"
# 
# $UGT2B15
# [1] "UGT2B15"
# 
# $UGT2B17
# [1] "UGT2B17"
# 
# $UGT2B28
# [1] "UGT2B28"
# 
# $UGT3A1
# [1] "UGT3A1"
# 
# $UGT3A2
# [1] "UGT3A2"
# 
# $UGT8
# [1] "UGT8"



# _______________________________________________________________________________
#  3.1.3  Apply ADME-optimized functionality prediction framework 
# _______________________________________________________________________________

## Scores of interest
scores_algorithms <- c('SIFT_score', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'LRT_score','MutationAssessor_score', 'FATHMM_score', 
                       'fathmm.MKL_coding_score', 'PROVEAN_score', 'VEST3_score', 'VEST4_score', 'CADD_phred', 'DANN_score', 'MetaSVM_score', 'MetaLR_score', 
                       'REVEL_score', 'PrimateAI_score', 'M.CAP_score', 'ClinPred_score', 'Eigen.PC.raw_coding', 'MutPred_score', 'MVP_score')
## Categorical predictions 
cat_pred_algorithms <- c('SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'MutationAssessor_pred', 'FATHMM_pred', 
                        'fathmm.MKL_coding_pred', 'PROVEAN_pred', 'MetaSVM_pred', 'MetaLR_pred', 'M.CAP_pred', 'ClinPred_pred', 'LRT_pred', 'PrimateAI_pred')

## Create single dataset for all missense UGT variants and their predicted effects and scores
## Unique variant IDs
for (gene in UGT_genes){
  myanno <- eval(parse_expr(paste0('myanno_', gene)))
  myanno$Variant_ID <- paste(myanno$Chr, myanno$Start, myanno$Ref, myanno$Alt, sep='-')
  assign(paste0('myanno_', gene), myanno)
}

## Total unique UGT variants
unique_UGT_vars <- unique(unlist(sapply(paste0('myanno_', UGT_genes, '$Variant_ID'), function(x){eval(parse_expr(x))})))

variants_predictions <- data.frame(matrix(nrow=0, ncol=37))
colnames(variants_predictions) <- c('Variant_ID', 'Gene.refGene', scores_algorithms, cat_pred_algorithms)

for (UGT_variant in unique_UGT_vars){
  ## Search variant in each gene dataset
  for(gene in UGT_genes){
    myanno <- eval(parse_expr(paste0('myanno_', gene)))
    if (UGT_variant %in% myanno$Variant_ID){
      ## For shared variants, assume the predictions are the same across all genes and take the ones reported in the first one
      variants_predictions  <- rbind(variants_predictions, myanno[which(myanno$Variant_ID==UGT_variant), c('Variant_ID', 'Gene.refGene', scores_algorithms, cat_pred_algorithms)])
      break
    }
  }
}

## Correct category names for FATHMM_pred and PrimateAI_pred (TRUE -> 'T')
variants_predictions$FATHMM_pred <- replace(variants_predictions$FATHMM_pred, 
                                            which(variants_predictions$FATHMM_pred==TRUE), 'T')
variants_predictions$PrimateAI_pred <- replace(variants_predictions$PrimateAI_pred, 
                                               which(variants_predictions$PrimateAI_pred==TRUE), 'T')
## Standardize variable names for raw scores
colnames(variants_predictions)[c(9, 13, 21, 29)] <- c('fathmm.MKL_score', 'CADD_phred_score', 'Eigen.PC_score', 'fathmm.MKL_pred')


## Add ADME-optimized model scores

## ADME-optimized algorithm thresholds to categorize variants 
ADME_thresholds <- list('LRT'='<0.0025',
                        'MutationAssessor'='>2.0566',
                        'PROVEAN'='< -3.286',
                        'VEST3'='>0.4534',
                        'CADD_phred'='>19.19')

## Categorize variants in D/N by each algorithm using these thresholds
ADME_categorical_predictions <- data.frame(matrix(nrow=dim(variants_predictions)[1], ncol=length(ADME_thresholds)+1))
colnames(ADME_categorical_predictions) <- c('Variant_ID', paste0(names(ADME_thresholds), '_pred'))
ADME_categorical_predictions$Variant_ID <- variants_predictions$Variant_ID

for(algorithm in names(ADME_thresholds)){
  ## Evaluate if the algorithm score of each variant passes cutoff (1) or not (0)
  ADME_categorical_predictions[paste0(algorithm, '_pred')] <- apply(variants_predictions, 1, 
                                                               function(x){if (x[paste0(algorithm, '_score')]=='.'){'.'}
                                                                 else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', ADME_thresholds[[algorithm]]))) ){1}
                                                                 else{0} })
}

## Add global ADME scores for each variant
variants_predictions$ADME_score <- signif(apply(ADME_categorical_predictions[,-1], 1, function(x){mean(as.numeric(x[which(x!='.')]))}), digits=3)
## Put '.' of no algorithm used in ADME model returned scores for a variant
variants_predictions$ADME_score[which(is.nan(variants_predictions$ADME_score))] <- '.'



# _______________________________________________________________________________
#  3.1.4  Extract AlphaMissense (AM) predictions for all UGT variants
# _______________________________________________________________________________

## Retrieve AM scores and predictions for all variants of a gene
for (gene in UGT_genes){
  tx <- canonical_txs[[gene]]
  data <- read_tsv(paste0('~/Desktop/UGT_genetic_profiling_KI/raw-data/AlphaMissense_data/', tx, '_AlphaMissense_data'), show_col_types = FALSE)
  ## Add variant IDs
  colnames(data) <- c('chr', 'pos', 'ref', 'alt', 'genome', 'uniprot_id', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class')
  data$chr <- sapply(data$chr, function(x){strsplit(x, 'chr')[[1]][2]})
  data$Variant_ID <- paste(data$chr, data$pos, data$ref, data$alt, sep='-')
  assign(paste0(gene, '_AlphaMissense_data'), data)
}


## Search UGT missense variants in AM datasets
AM_score_pred <- function(variant_id){
  ## Data for each variant
  var_data <- vector()
  ## Search variant within each gene 
  for(gene in UGT_genes){
    AM_data <- eval(parse_expr(paste0(gene, '_AlphaMissense_data')))
    if (variant_id %in% AM_data$Variant_ID){
      var_data  <- rbind(var_data, AM_data[which(AM_data$Variant_ID==variant_id), c('Variant_ID', 'am_pathogenicity', 'am_class')])
    }
    
    if (!is.null(dim(var_data))){
      ## For shared variants with differing scores across genes, take the most pathogenic one (the biggest) and its corresponding class
      score <- var_data[order(var_data$am_pathogenicity, decreasing = TRUE),][1,'am_pathogenicity']
      class <- var_data[order(var_data$am_pathogenicity, decreasing = TRUE),][1,'am_class']
    }
    ## If variant was not found in any gene dataset (no AM score/pred)
    else{
      score <- class <- '.'
    }
  }
  return(c(score, class))
}

## Search AM score/pred for all UGT variants
for (i in 1:dim(variants_predictions)[1]){
  AM_output <- AM_score_pred(variants_predictions$Variant_ID[i])
  variants_predictions$AlphaMissense_score[i] <- AM_output[[1]]
  variants_predictions$AlphaMissense_pred[i] <- AM_output[[2]]
}

variants_predictions$AlphaMissense_pred <- replace(replace(replace(variants_predictions$AlphaMissense_pred, 
                                                   which(variants_predictions$AlphaMissense_pred=='benign'), 'N'), 
                                                   which(variants_predictions$AlphaMissense_pred=='pathogenic'), 'D'),
                                                   which(variants_predictions$AlphaMissense_pred=='ambiguous'), 'U')

save(variants_predictions, file='processed-data/03_Anno_functional_impact/variants_scores_and_predictions.Rdata')


# ____________________________________________________________________________________________
#  3.1.5  Comparison and evaluation of predictive algorithms 
# ____________________________________________________________________________________________

####################  3.1.5.1 Compare predictions of different algorithms  ####################

## Function to create density plot of raw scores for variants in the different functional categories

score_density_plot <- function(algorithm, predicted_cat_type){
  
  algorithm_score <- paste0(algorithm, '_score')
  algorithm_pred <- paste0(algorithm, '_pred')
  
  ## Define colors for categories of predicted effect
  colors <- list('D'='tomato2', 'T'='skyblue1', 'N'='skyblue1', 'B'='skyblue1', 'P'='lightsalmon', 
                 'H'= 'red4', 'M'='red3', 'L'='dodgerblue3', 'U'='grey90')
  
  if (predicted_cat_type=='new'){
    
    threshold <- algorithms_thresholds[[algorithm]]
    numeric_threshold <- as.numeric(gsub('[<, =, >]', '', threshold))
    
    ## Subset to variants with algorithm scores (and predictions)
    data <- new_variants_predictions[which(eval(parse_expr(paste0('new_variants_predictions$',algorithm_score)))!='.'),]

    ## Percentage of variants with missing scores/predictions from each algorithm (percentage of missingness)
    missingness <- apply(new_variants_predictions, 2, function(x){100*length(which(x=='.'))/dim(new_variants_predictions)[1]})[[algorithm_score]]
    ## Number of variants with valid output from each algorithm
    num_vars <- apply(new_variants_predictions, 2, function(x){length(which(x!='.'))})[[algorithm_score]]
    
    ## Density function with score limits 
    density <- density(as.numeric(data[,algorithm_score]), 
                       from = min(as.numeric(data[,algorithm_score])), 
                       to = max(as.numeric(data[,algorithm_score]))) 
    ## Categorize variant scores by defined threshold 
    df <- data.frame(x = density$x, y = density$y, pred = replace(replace(density$x, which(eval(parse_expr(paste0('density$x', threshold)))), 'D'), 
                                                                  which(!eval(parse_expr(paste0('density$x', threshold)))), 'N')) 
    if (algorithm=='FATHMM'){
      hjust = 0.3
    }
    else if(algorithm=='VEST4' | algorithm=='DANN' | algorithm=='PrimateAI' | algorithm=='CADD'){
      hjust = 1.5
    }
    else {
      hjust = -0.3
    }
    
    if (algorithm=='CADD_phred'){score_type <- ' scores'}
    else {score_type <- ' raw scores'}
    
    p1 <- ggplot(data = df, aes(x = x, ymin = 0, ymax = y, fill = pred)) +
      geom_ribbon(alpha=0.7) +
      theme_bw() +
      scale_fill_manual(values=colors[names(table(df$pred))]) +
      labs(x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm)), score_type), y= 'Density', fill='Predicted effect',
           subtitle=paste0('Missingness: ', signif(as.numeric(missingness), digits=3), '%', '\n', 
                           num_vars, ' variants')) +
      geom_line(aes(y = y)) +
      geom_vline(xintercept = numeric_threshold, color = 'indianred3', linetype='dashed', linewidth=0.6) +
      geom_label(aes(x = numeric_threshold, y = max(df$y), color = 'indianred3', label = numeric_threshold), 
               hjust = hjust, vjust = 3, fontface = 2, fill = "white", show.legend = FALSE) +
      theme(legend.key = element_rect(fill = "white", colour = "black"),
            plot.subtitle = element_text(size = 10, color = "gray30"))
    
    return(p1)
  }
  
  else{
    
    data <- variants_predictions[which(eval(parse_expr(paste0('variants_predictions$',algorithm_score)))!='.'),]
    
    p2 <- ggplot(data = data, aes(x = as.numeric(eval(parse_expr(algorithm_score))))) +
      geom_density(alpha=0.6, fill='grey')+
      theme_bw() +
      labs(x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm)), ' raw scores'), y= 'Density') 
    
    p3 <- ggplot(data = data, aes(x = as.numeric(eval(parse_expr(algorithm_score))), 
                                  fill=eval(parse_expr(algorithm_pred))))+
      geom_density(alpha=0.6)+
      scale_fill_manual(values=colors[names(table(eval(parse_expr(paste0('data$', algorithm_pred)))))]) +
      theme_bw() +
      labs(x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm)), ' raw scores'), y= 'Density', fill='Predicted effect') 
    
    return(list(p2, p3))
  }

    
}
  

##################  a) Raw scores of variants already predicted as D and N/B/T  ##################

## Algorithms already returning categorical predictions
algorithms <- c('SIFT', 'Polyphen2_HDIV', 'Polyphen2_HVAR', 'MutationAssessor', 'FATHMM', 
                'fathmm.MKL', 'PROVEAN', 'MetaSVM', 'MetaLR', 'M.CAP', 'ClinPred', 'LRT', 'PrimateAI', 'AlphaMissense')

plots <- list()
j=1
for (i in 1:length(algorithms)){
  plots[[j]] <- score_density_plot(algorithms[i], 'returned')[[1]]
  plots[[j+1]] <- score_density_plot(algorithms[i], 'returned')[[2]]
  j=j+2
}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]], plots[[22]], 
          plots[[23]], plots[[24]], plots[[25]], plots[[26]], plots[[27]], plots[[28]],
          ncol=6, rel_widths = c(rep(c(0.74,1), 14)))

ggsave(filename='plots/03_Anno_functional_impact/Returned_RawScores_density_plots.pdf', width = 30, height = 18)


##################  b) Raw scores of variants categorized by defined score thresholds  ##################

## Reported/conventional threshold to categorize variants as deleterious (D) (or neutral (N) otherwise) in each algorithm 
algorithms_thresholds <- list('SIFT'='<=0.05',               
                              'Polyphen2_HDIV'='>0.452',     
                              'Polyphen2_HVAR'='>0.446',     
                              'MutationAssessor'='>1.9',     
                              'FATHMM'= '<= -1.5',           
                              'fathmm.MKL'='>0.5',   
                              'PROVEAN'= '< -2.282',         
                              'MetaSVM'='>=0',               
                              'MetaLR'='>=0.5',             
                              'M.CAP'='>=0.025',            
                              'ClinPred'='>=0.5',               
                              'CADD_phred'='>15',              
                              'DANN'='>0.96', 
                              'REVEL'='>0.5',               
                              'Eigen.PC'='>=0', 
                              'MVP'= '>0.75',
                              'LRT'='<0.001',
                              'MutPred'='>0.5',
                              'PrimateAI'='>=0.803',
                              'VEST4'='>0.5',
                              'ADME'='>0.5',
                              'AlphaMissense'='>=0.564')

## Categorize variants with these thresholds
categorical_predictions <- data.frame(matrix(nrow=dim(variants_predictions)[1], ncol=length(names(algorithms_thresholds))+1))
colnames(categorical_predictions) <- c('Variant_ID', paste0(names(algorithms_thresholds), '_pred'))
categorical_predictions$Variant_ID <- variants_predictions$Variant_ID

for(algorithm in names(algorithms_thresholds)){
  ## Check if the algorithm score of each variant is above/below/equal to threshold
  categorical_predictions[paste0(algorithm, '_pred')] <- apply(variants_predictions, 1, 
                                                               function(x){if (x[paste0(algorithm, '_score')]=='.' | x[paste0(algorithm, '_score')]=='-'){'.'}
                                                                          else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', algorithms_thresholds[[algorithm]]))) ){'D'}
                                                                          else{'N'} })
}
## Bind scores and new binary predictions per algorithm 
new_variants_predictions <- cbind(categorical_predictions, variants_predictions[,paste0(names(algorithms_thresholds), '_score')])
new_variants_predictions$MutPred_score <- replace(new_variants_predictions$MutPred_score, which(new_variants_predictions$MutPred_score=='-'), '.')


## Distribution of raw scores in D and N variants defined by cutoffs
plots <- list()
for (i in 1:length(names(algorithms_thresholds))){
  plots[[i]] <- score_density_plot(names(algorithms_thresholds)[i], 'new')
}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]], plots[[22]],
          ncol=5)

ggsave(filename='plots/03_Anno_functional_impact/New_RawScores_density_plots.pdf', width = 20, height = 12)
save(new_variants_predictions, file='processed-data/03_Anno_functional_impact/new_variants_predictions.Rdata')


# ------------------------------------------------------------------------------
## Correlation between raw scores from each pair of methods
raw_scores <- as.data.frame(apply(new_variants_predictions[,paste0(names(algorithms_thresholds), '_score')], 2, as.numeric))
corr <- matrix(nrow=22, ncol = 22)
colnames(corr) <- rownames(corr) <- colnames(raw_scores)

for (i in 1:length(colnames(raw_scores))){
  for (j in 1:length(colnames(raw_scores))){
    ## Subset to variants with valid scores in both algorithms
    raw_scores_subset <- raw_scores[which(!is.na(raw_scores[,colnames(raw_scores)[i]]) & !is.na(raw_scores[,colnames(raw_scores)[j]])),]
    corr[i, j] <- signif(cor(raw_scores_subset[,colnames(raw_scores)[i]], raw_scores_subset[,colnames(raw_scores)[j]], method = 'pearson'), digits=2)
  }
}

whole_corr <- corr
colnames(corr) <- rownames(corr) <- gsub('phred', '', gsub('\\.','-', gsub('_', ' ', gsub('_score', '', colnames(corr)))))
## Half matrix
corr[lower.tri(corr)] <- NA
## Take absolute corr
corr <- abs(corr)
half_corr_data <- melt(corr, na.rm = TRUE)
half_corr_data$value <- signif(as.numeric(half_corr_data$value), digits = 3)

## Mean corr coeff
unique_half_corr_data <- half_corr_data[which(half_corr_data$value!=1), ]
mean(unique_half_corr_data$value)
# [1] 0.5680519

## Highest corr coeffs
unique_half_corr_data[order(unique_half_corr_data$value, decreasing = TRUE), ][1:4,]
#               Var1             Var2    value
#     Polyphen2 HDIV   Polyphen2 HVAR     0.97
#            MetaSVM           MetaLR     0.93
#         CADD phred         Eigen-PC     0.93
#            MetaSVM            REVEL     0.88


## Percentage of high coeffs (>0.5)
length(which(unique_half_corr_data$value>0.5))/dim(unique_half_corr_data)[1]*100
# [1] 65.36797

## Percentage of medium coeffs (0.3=<r=<0.5)
length(which(unique_half_corr_data$value>=0.3 & unique_half_corr_data$value<=0.5))/dim(unique_half_corr_data)[1]*100
# [1] 25.97403

## Percentage of low coeffs (<0.3)
length(which(unique_half_corr_data$value<0.3))/dim(unique_half_corr_data)[1]*100
# [1] 8.658009


# ------------------------------------------------------------------------------
## Agreement proportion between predictions from each pair of methods
predictions <- as.data.frame(new_variants_predictions[,paste0(names(algorithms_thresholds), '_pred')])
agreement_prop <- matrix(ncol=22, nrow=22)
colnames(agreement_prop) <- rownames(agreement_prop) <- colnames(predictions)

for (i in 1:length(colnames(predictions))){
  for(j in 1:length(colnames(predictions))){
    ## Subset to variants with predictions from both algorithms
    predictions_subset <- predictions[which(!is.na(predictions[,colnames(predictions)[i]]) & !is.na(predictions[,colnames(predictions)[j]])),]
    ## Evaluate if predictions from algorithm 1 and 2 are different or the same for each variant
    comparisons <- apply(predictions_subset[,c(colnames(predictions)[i], colnames(predictions)[j])], 1, function(x){if(x[1]==x[2]){'equal'}else{'diff'}})
    agreement_prop[i, j]<- table(comparisons)['equal']/dim(predictions_subset)[1]
  }
}

whole_agreement_prop <- agreement_prop
## Half matrix
colnames(agreement_prop) <- rownames(agreement_prop) <- gsub('phred', '', gsub('\\.','-', gsub('_', ' ', gsub('_pred', '', colnames(agreement_prop)))))
agreement_prop[lower.tri(agreement_prop)] <- NA 
half_agreement_prop <- melt(agreement_prop, na.rm = TRUE)
half_agreement_prop$value <- signif(as.numeric(half_agreement_prop$value), digits = 3)


pdf(file = "plots/03_Anno_functional_impact/Compare_algorithms_outcomes.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
corrplot(corr, method="color",
         diag=FALSE,
         type="upper",
         addCoef.col = "black",
         tl.cex = 0.4,
         tl.col = 'black',
         number.cex = 0.3,
         cl.cex = 0.4,
         col.lim = c(0,1),
         col=colorRampPalette(c("white","white","white", 'mistyrose2',"tomato2", 'firebrick4'))(100)
)         
corrplot(agreement_prop, method="color", 
         diag=FALSE, 
         type="upper",
         addCoef.col = "black",
         tl.cex = 0.4,
         tl.col = 'black',
         number.cex = 0.3,
         cl.cex = 0.4,
         col.lim = c(0,1),
         col=colorRampPalette(c("white","white", 'darkseagreen2', "darkgreen"))(100)
         
)

dev.off()

## Mean prop
unique_half_agreement_prop <- half_agreement_prop[which(half_agreement_prop$value!=1), ]
mean(unique_half_agreement_prop$value)
# [1] 0.6372771

## Highest corr coeffs
unique_half_agreement_prop[order(unique_half_agreement_prop$value, decreasing = TRUE), ][1:4,]
#             Var1             Var2   value
#   Polyphen2 HDIV   Polyphen2 HVAR   0.949
#           MetaSVM          MetaLR   0.930
#           MetaSVM           REVEL   0.926
#            MetaLR           REVEL   0.905

## Percentage of high prop (>0.5)
length(which(unique_half_agreement_prop$value>0.5))/dim(unique_half_agreement_prop)[1]*100
# [1] 82.25108


# ------------------------------------------------------------------------------
## Scatterplot of raw scores and categorical predictions of 2 methods

scatterplot_compare_2methods <- function(algorithm1, algorithm2){
  
  algorithm1_score <- paste0(algorithm1, '_score')
  algorithm2_score <- paste0(algorithm2, '_score')
  algorithm1_pred <- paste0(algorithm1, '_pred')
  algorithm2_pred <- paste0(algorithm2, '_pred')
  
  if (algorithm1=='CADD_phred' | algorithm2=='CADD_phred'){score_type <- ' score'}
  else {score_type <- ' raw score'}
  
  ## Corr
  correlation <- whole_corr[algorithm1_score, algorithm2_score]
  ## Agreement
  agreement<- signif(whole_agreement_prop[algorithm1_pred, algorithm2_pred], digits=2)
  
  ## Scores
  data <- new_variants_predictions[which(new_variants_predictions[,algorithm1_score]!='.' & new_variants_predictions[,algorithm2_score]!='.'),
                                        c(algorithm1_score, algorithm2_score)]
  colnames(data) <- gsub('_score', '', colnames(data))
  data <- as.data.frame(apply(data, 2, as.numeric))
  
  ## Categories
  data_pred <- new_variants_predictions[which(new_variants_predictions[,algorithm1_pred]!='.' & new_variants_predictions[,algorithm2_pred]!='.'),
                           c(algorithm1_pred, algorithm2_pred)]
  colnames(data_pred) <- colnames(data)
  data$categories <- apply(data_pred, 1, function(x){if(x[1]=='D' & x[2]=='D'){'D in both'}
                            else if(x[1]=='D' & x[2]=='N'){paste0('D in ', gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[1]))), '; ', 'N in ', 
                                                                  gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[2]))))}
                            else if(x[1]=='N' & x[2]=='D'){paste0('D in ', gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[2]))), '; ', 'N in ', 
                                                                  gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[1]))))}
                            else {'N in both'}
                      })
  colors <- list()
  colors[['D in both']]='indianred'
  colors[[paste0('D in ', gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[1]))), '; ', 'N in ', 
                 gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[2]))))]]='thistle2'
  colors[[paste0('D in ', gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[2]))), '; ', 'N in ', 
                 gsub(' phred', '', gsub('_',' ',gsub('\\.', '-', colnames(data)[1]))))]]='lightpink1'
  colors[['N in both']]='lightblue3'
  
  ggplot(data, aes(x=eval(parse_expr(algorithm1)), y=eval(parse_expr(algorithm2)), color = categories)) +
    geom_point(size=2) +
    stat_smooth(geom = "line", alpha = 0.9, size = 1, method = lm, color = "orangered3", fullrange = FALSE) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      subtitle = paste0("Corr: ", correlation, '; Agreement: ', 100*agreement, '%'),
      x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm1)), score_type),
      y = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm2)), score_type),
      color='Predicted effect'
    ) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.title = element_text(size = (11)),
      axis.text = element_text(size = (10)),
      plot.subtitle = element_text(size = 9.6, color = "gray30"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size =10))
  ggsave(paste0('plots/03_Anno_functional_impact/Corr_', algorithm1, '_vs_', algorithm2, '.pdf'), height = 4.5, width = 6.5)
}

## Compare methods of interest
## Top; similar algorithms
scatterplot_compare_2methods('Polyphen2_HDIV', 'Polyphen2_HVAR')
scatterplot_compare_2methods('MetaSVM', 'MetaLR')
scatterplot_compare_2methods('FATHMM', 'fathmm.MKL')
## Top on agreement
scatterplot_compare_2methods('MetaSVM', 'REVEL')
## Low corr, high agreement
scatterplot_compare_2methods('FATHMM', 'VEST4') 
scatterplot_compare_2methods('FATHMM', 'PrimateAI')
scatterplot_compare_2methods('Eigen.PC', 'M.CAP') 
## High corr, high agreement
scatterplot_compare_2methods('Eigen.PC', 'CADD_phred')
## ADME
scatterplot_compare_2methods('ADME', 'CADD_phred')
scatterplot_compare_2methods('ADME', 'VEST4')
scatterplot_compare_2methods('ADME', 'LRT')
scatterplot_compare_2methods('ADME', 'MutationAssessor')
scatterplot_compare_2methods('ADME', 'PROVEAN')
scatterplot_compare_2methods('ADME', 'ClinPred')
scatterplot_compare_2methods('ADME', 'REVEL')
scatterplot_compare_2methods('ADME', 'MetaLR')
## AlphaMissense
scatterplot_compare_2methods('AlphaMissense', 'MetaSVM')
scatterplot_compare_2methods('AlphaMissense', 'MetaLR')
scatterplot_compare_2methods('AlphaMissense', 'REVEL')
scatterplot_compare_2methods('AlphaMissense', 'VEST4')
## FATHMM
scatterplot_compare_2methods('FATHMM', 'MutPred')
scatterplot_compare_2methods('FATHMM', 'CADD_phred')
scatterplot_compare_2methods('FATHMM', 'DANN')



####################  3.1.5.2 Evaluate predictions of different algorithms  ####################

## Plot the number of predicted D, N and missing variants per algorithm

vars_per_method <- melt(sapply(colnames(new_variants_predictions)[2:23], function(x){table(new_variants_predictions[, x])}))
colnames(vars_per_method) <- c('Prediction', 'Method', 'Numbers')
vars_per_method$Method <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred', '', vars_per_method$Method))))
# Order methods by number of D variants
numD<- sapply(colnames(new_variants_predictions)[2:23], function(x){table(new_variants_predictions[, x])['D']})
numD <- numD[order(numD, decreasing = TRUE)]
names(numD) <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred.D', '', names(numD)))))
vars_per_method$Method <- factor(vars_per_method$Method, levels=names(numD))
## Order to have D variants first in each bar
cat_order <- c('.', 'N', 'D')
vars_per_method$Prediction <- factor(vars_per_method$Prediction, levels=cat_order)

ggplot(vars_per_method, aes(x=Method, y=Numbers, fill=Prediction)) + 
  geom_bar(position="stack", stat="identity", colour = 'black', width=.7) + 
  theme_classic() +
  labs(
    y = 'Number of predicted missense variants',
    x= '',
    subtitle = paste0(dim(new_variants_predictions)[1], ' total missense variants across all UGT genes')
  ) +
  scale_fill_manual(values = c("grey80", "skyblue2", "tomato"), labels = c("Missing", "Neutral", "Deleterious")) + 
  #scale_y_discrete(limits= c(0, 1500, 3000, 4500, dim(new_variants_predictions)[1])) +
  coord_cartesian(ylim = c(0, dim(new_variants_predictions)[1]),clip = 'off') +
  geom_text(data=subset(vars_per_method, Prediction=='.'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+730, 
                                                               fill=NULL), hjust = 0.5, size = 2.2) +
  geom_text(data=subset(vars_per_method, Prediction=='N'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+460, 
                                                               fill=NULL), hjust = 0.5, size = 2.2) +
  geom_text(data=subset(vars_per_method, Prediction=='D'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+190, 
                                                               fill=NULL), hjust = 0.5, size = 2.2) +
  theme(plot.subtitle = element_text(size = (10), vjust = 9.8, hjust=0, color="gray50", face='bold'), 
        legend.direction = "vertical",
        legend.position = c(1.07, 1.05),
        legend.key.size = unit(0.3, units = 'cm'),
        plot.margin = unit(c(3, 3, 0.5, 0.5), "cm"),
        axis.title = element_text(size = (11), face='bold'),
        axis.text = element_text(size = (10)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text = element_text(size = 9),
        legend.title = element_text(size =10, face='bold'))

ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method.pdf', width = 8, height = 6)

# -------------------------------------------------------------------------------------------------

## Plot the number of D, N and missing variants per algorithm per gene

DNM_vars_per_method_per_gene <- function(gene){
  gene_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))$Variant_ID
  gene_vars <- new_variants_predictions[which(new_variants_predictions$Variant_ID %in% gene_vars),]
  
  DNM_numbers <- sapply(colnames(gene_vars)[2:23], function(x){ c('D'=length(which(gene_vars[,x]=='D')),
                                                                           'N'=length(which(gene_vars[,x]=='N')),
                                                                           '.'=length(which(gene_vars[,x]=='.'))) })
  vars_per_method <- melt(DNM_numbers)
  colnames(vars_per_method) <- c('Prediction',  'Method', 'Numbers')
  vars_per_method$Method <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred', '', vars_per_method$Method))))
  
  # Order methods by number of D variants first, and N variants second
  DNM_numbers <- as.data.frame(t(DNM_numbers))
  DNM_numbers <- DNM_numbers[order(DNM_numbers$D, DNM_numbers$N,  decreasing = TRUE),]
  numD <- rownames(DNM_numbers)
  numD <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred', '', numD))))
  vars_per_method$Method <- factor(vars_per_method$Method, levels=numD)
  
  ## Order to have D variants first in each bar
  cat_order <- c('.', 'N', 'D')
  vars_per_method$Prediction <- factor(vars_per_method$Prediction, levels=cat_order)
  vars_per_method$Numbers <- as.numeric(vars_per_method$Numbers)
  
  if (dim(gene_vars)[1]>450){
    s=21
  }
  else if(dim(gene_vars)[1]>400){
    s=18
  }
  else if (dim(gene_vars)[1]>=344){
    s=15
  }
  else if(dim(gene_vars)[1]>300){
    s=13
  }
  else if (dim(gene_vars)[1]<250){
    s=9
  }
  else {
    s=11
  }
  
  p <- ggplot(vars_per_method, aes(x=Method, y=Numbers, fill=Prediction)) + 
    geom_bar(position="stack", stat="identity", colour = 'black', width=.7) + 
    theme_classic() +
    labs(
      y = 'Number of predicted missense variants',
      x= '',
      title=gene,
      subtitle = paste0(dim(gene_vars)[1], ' total missense variants')
    ) +
    scale_fill_manual(values = c("grey80", "skyblue2", "tomato"), labels = c("Missing", "Neutral", "Deleterious")) + 
    coord_cartesian(ylim = c(0, dim(gene_vars)[1]), # This focuses the x-axis on the range of interest
                    clip = 'off') +
    geom_text(data=subset(vars_per_method, Prediction=='.'), aes(label=Numbers, y=dim(gene_vars)[1]+(3*s), 
                                                                 fill=NULL), hjust = 0.5, size = 2.6) +
    geom_text(data=subset(vars_per_method, Prediction=='N'), aes(label=Numbers, y=dim(gene_vars)[1]+(2*s), 
                                                                 fill=NULL), hjust = 0.5, size = 2.6) +
    geom_text(data=subset(vars_per_method, Prediction=='D'), aes(label=Numbers, y=dim(gene_vars)[1]+s, 
                                                                 fill=NULL), hjust = 0.5, size = 2.6) +
    
    theme(plot.title = element_text(size = (12), face='bold', vjust = 8, hjust=0), 
          plot.subtitle = element_text(size = (10), vjust = 9.8, hjust=0, color="gray50", face='bold'), 
          legend.direction = "vertical",
          legend.position = c(1.07, 1.06),
          legend.key.size = unit(0.3, units = 'cm'),
          plot.margin = unit(c(3, 3, 0.5, 0.5), "cm"),
          axis.title = element_text(size = (11), face='bold'),
          axis.text = element_text(size = (10)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.text = element_text(size = 9),
          legend.title = element_text(size =10, face='bold'))
  return(p)
}

i=1
plots <- list()
for (gene in UGT_genes){
  plots[[i]] <- DNM_vars_per_method_per_gene(gene)
  i=i+1
}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], ncol=3)
ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method_UGT1_genes.pdf', width = 24, height = 19)

plot_grid(plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], ncol=5)
ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method_UGT2_genes.pdf', width = 39, height = 13)

plot_grid(plots[[20]], plots[[21]], ncol=2)
ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method_UGT3_genes.pdf', width = 16, height = 6)

plot_grid(plots[[22]])
ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method_UGT8_genes.pdf', width = 8, height = 5.8)

# --------------------------------------------------------------------------------------------------------------

## Plot allele frequencies of predicted D variants per method

colors = c('ADME'='mediumpurple2', 
           'AlphaMissense'='red3',
           'CADD'='darkorange3', 
           'ClinPred'='turquoise3',
           'DANN'='yellow3', 
           'Eigen-PC'='olivedrab3',
           'FATHMM'='lightsalmon2', 
           'fathmm-MKL'='lightcoral', 
           'LRT'='lightsteelblue3',
           'M-CAP'='yellow4', 
           'MetaLR'='steelblue2', 
           'MetaSVM'='dodgerblue3', 
           'MutationAssessor'='goldenrod', 
           'MutPred'='magenta2',
           'MVP'='blue2', 
           'Polyphen2 HDIV'='darkseagreen3', 
           'Polyphen2 HVAR'='mediumseagreen', 
           'PrimateAI'='darkred',
           'PROVEAN'='darkorchid3', 
           'REVEL'='cadetblue3', 
           'SIFT'='peachpuff3', 
           'VEST4'='aquamarine4')

for (variant in new_variants_predictions$Variant_ID){
  allele_freq <- vector()
  ## Allele freq of variant in each gene dataset
  for (gene in UGT_genes){
    UGT_missense_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))
    if(variant %in% UGT_missense_vars$Variant_ID){
      allele_freq <- append(allele_freq, UGT_missense_vars[which(UGT_missense_vars$Variant_ID==variant), 'Allele_Frequency'])
    }
  }
  if (length(unique(allele_freq))==1) {allele_freq <- unique(allele_freq)}
  new_variants_predictions[which(new_variants_predictions$Variant_ID==variant), 'Allele_Frequency'] <- allele_freq
}

## Allele frequencies of D variants per method
data <- vector()
for(method in colnames(new_variants_predictions)[2:23]){
  allele_freq_method <-  new_variants_predictions[which(new_variants_predictions[,method]=='D'), c('Variant_ID', 'Allele_Frequency')]
  method <- rep(method, length(allele_freq_method$Allele_Frequency))
  method_data <- cbind(allele_freq_method$Variant_ID, allele_freq_method$Allele_Frequency, method)
  data <- rbind(data, method_data)
}
data <- as.data.frame(data)
colnames(data) <- c('Variant_ID', 'Allele_Frequency', 'Method')
data$Allele_Frequency<- as.numeric(data$Allele_Frequency)
data$Method <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred', '', data$Method))))
data$Method <- factor(data$Method, levels=names(numD))

## Identify variants with allele freq >0.5
data[which(data$Allele_Frequency>0.5),]
#       Variant_ID  Allele_Frequency         Method
#   4-69795626-C-T         0.7568806 Polyphen2 HDIV
#   4-69795626-C-T         0.7568806 Polyphen2 HVAR
#   4-69795626-C-T         0.7568806           CADD
#   4-69795626-C-T         0.7568806           DANN
#  4-115589302-A-G         0.9953603            LRT

shapes <- c('4-69795626-C-T'=17,
            '4-115589302-A-G'=15)
data$label <- apply(data, 1, function(x){if (as.numeric(x['Allele_Frequency'])>0.5){x['Variant_ID']} else {NA}})
## Number of D vars per method
num_per_method <- as.data.frame(table(data$Method))
colnames(num_per_method) <- c('Method', 'n')

ggplot(data = data, mapping = aes(x = Method, y = Allele_Frequency, color = Method, shape=label)) +
  geom_jitter(data=subset(data, is.na(label)), shape=16, width = 0.1, height = 0, alpha = 0.7, size = 1.5) +
  geom_point(data=subset(data, !is.na(label)), aes(shape=label), alpha = 0.7, size = 1.3, stroke = 1) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  scale_color_manual(values = colors) +
  guides(color = 'none') + 
  geom_text(data = num_per_method, aes(x=Method, label=n,  y=-0.05, shape=NULL, color=NULL), size=2) +
  labs(x='', y='GMAF of missense variants predicted as deleterious', shape='Variant ID (GMAF>0.5)', 
       subtitle = paste0(dim(new_variants_predictions)[1], ' total missense variants across all UGT genes')) +
  theme(title = element_text(size = (9), face='bold'),
        plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
        axis.title = element_text(size = (8.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7.5))

ggsave(filename='plots/03_Anno_functional_impact/GMAF_allDvars_perMethod.pdf', width = 8, height = 5)

# ----------------------------------------------------------------------------------------------------

## Plot GMAF of variants predicted as D by each methods in each gene
i=1
plots <- list()
for (gene in UGT_genes){
  gene_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))$Variant_ID
  gene_D_vars <- data[which(data$Variant_ID %in% gene_vars),]
  gene_D_vars$Method <- factor(gene_D_vars$Method, levels=names(numD))
  ## Number of D vars per method
  num_per_method <- as.data.frame(table(gene_D_vars$Method))
  colnames(num_per_method) <- c('Method', 'n')
  num_per_method$Method <- factor(x = num_per_method$Method, levels=names(numD))
  
  ## No predictions by AlphaMissense for UGT8 tx
  if (gene=='UGT8'){
    gene_D_vars <- subset(gene_D_vars, Method!='AlphaMissense')
    num_per_method <- subset(num_per_method, Method!='AlphaMissense')
  }
  
  plots[[i]] <- ggplot(data = gene_D_vars, mapping = aes(x = Method, y = Allele_Frequency, color = Method, shape=label)) +
    geom_text(data = num_per_method, aes(label=n,  y=-0.03, shape=NULL, color=NULL), size=2) +
    geom_jitter(data=subset(gene_D_vars, is.na(label)), shape=16, width = 0.2, height = 0.0, alpha = 0.7, size = 1) +
    geom_point(data=subset(gene_D_vars, !is.na(label)), mapping=aes(shape=label), alpha = 0.7, size = 1, stroke = 1) +
    scale_shape_manual(values = shapes) +
    scale_y_continuous(limits = c(-0.03, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_color_manual(values = colors) +
    theme_bw() +
    guides(color = 'none') + 
    labs(x='', y='GMAF of missense variants predicted as deleterious', shape='Variant ID (GMAF>0.5)', 
         title=gene, subtitle = paste0(length(unique(gene_D_vars$Variant_ID)), ' missense variants predicted as deleterious by at least one method')) +
    theme(title = element_text(size = (9), face='bold'),
          plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
          axis.title = element_text(size = (8.5), face='bold'),
          axis.text = element_text(size = (8)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
          legend.title = element_text(size=6.5), 
          legend.text = element_text(size=6), 
          legend.key = element_blank(),,
          legend.background=element_blank(),
          legend.key.size = unit(0, 'lines'),
          legend.justification = c(0.9,0.9), legend.position = c(0.9,0.95))
  i=i+1
}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], ncol=3)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_UGT1_Dvars_perMethod.pdf', width = 16, height = 14)

plot_grid(plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], ncol=5)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_UGT2_Dvars_perMethod.pdf', width = 24, height = 8)

plot_grid(plots[[20]], plots[[21]], ncol=2)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_UGT3_Dvars_perMethod.pdf', width = 11, height = 5)

plot_grid(plots[[22]])
ggsave(filename='plots/03_Anno_functional_impact/GMAF_UGT8_Dvars_perMethod.pdf', width = 5.7, height = 4.2)

# ---------------------------------------------------------------------------------------------------------

## Plot cumulative MAF per method
ggplot() +
  sapply(names(numD), function(method){eval(parse_expr(paste0('geom_line(data=data[which(data$Method==\'', method, '\'),], aes(x=1:dim(data[which(data$Method==\'', method, '\'),])[1], y = cumsum(Allele_Frequency), color=Method), size=1, alpha=0.75)')))}) +
  theme_bw() +
  scale_color_manual(values = colors[names(numD)], breaks = names(numD)) +
  labs(x='Number of missense variants predicted as deleterious', y='Accumulated GMAF of variants') +
  theme(axis.title = element_text(size = (9), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename='plots/03_Anno_functional_impact/cumGMAF_Dvars_perMethod.pdf', width = 11, height = 7)  

# ----------------------------------------------------------------------------------------------------

## Plot carrier frequencies of predicted D variants per method

data$Carrier_Frequency <- (2*data$Allele_Frequency*(1-data$Allele_Frequency)) + (data$Allele_Frequency)**2
num_per_method <- as.data.frame(table(data$Method))
colnames(num_per_method) <- c('Method', 'n')

## Identify variants with allele freq >0.5
data[which(data$Carrier_Frequency>0.5),]
# Variant_ID      Allele_Frequency         Method           label Carrier_Frequency
# 2-234602202-A-C        0.3443850           SIFT            <NA>         0.5701690
#  4-69795626-C-T        0.7568806 Polyphen2 HDIV  4-69795626-C-T         0.9408930
#  4-69795626-C-T        0.7568806 Polyphen2 HVAR  4-69795626-C-T         0.9408930
#  4-69795626-C-T        0.7568806           CADD  4-69795626-C-T         0.9408930
#  4-69795626-C-T        0.7568806           DANN  4-69795626-C-T         0.9408930
#  4-70156313-T-A        0.4690027            LRT            <NA>         0.7180419
#  4-70160309-C-G        0.4894400            LRT            <NA>         0.7393285
# 4-115589302-A-G        0.9953603            LRT 4-115589302-A-G         0.9999785

shapes <- c('4-69795626-C-T'=17,
            '4-115589302-A-G'=15, 
            '2-234602202-A-C'=4,
            '4-70156313-T-A'=0,
            '4-70160309-C-G'=18)

data$label <- apply(data, 1, function(x){if (as.numeric(x['Carrier_Frequency'])>0.5){x['Variant_ID']} else {NA}})

## Number of D vars per method
num_per_method <- as.data.frame(table(data$Method))
colnames(num_per_method) <- c('Method', 'n')

ggplot(data = data, mapping = aes(x = Method, y = Carrier_Frequency, color = Method)) +
  geom_jitter(data=subset(data, is.na(label)), shape=16, width = 0.1, height = 0, alpha = 0.7, size = 1.5) +
  geom_point(data=subset(data, !is.na(label)), aes(shape=label), alpha = 0.7, size = 1.1, stroke = 1) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  scale_color_manual(values = colors) +
  guides(color = 'none') + 
  geom_text(data = num_per_method, aes(x=Method, label=n,  y=-0.05, shape=NULL, color=NULL), size=2) +
  labs(x='', y='Carrier frequency of missense variants predicted as deleterious', shape=paste0('Variant ID', '\n', '(Carrier frequency >0.5)'),
       subtitle = paste0(dim(new_variants_predictions)[1], ' total missense variants across all UGT genes')) +
  theme(title = element_text(size = (9), face='bold'),
        plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
        axis.title = element_text(size = (8.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7.5))

ggsave(filename='plots/03_Anno_functional_impact/CarrierFreq_Dvars_perMethod.pdf', width = 8, height = 5)



## Evaluate prediction accuracy of methods

## ClinVar variants for each gene (no variants in UGT2A3)
benchmark_data <- vector()
benchmark_whole_data <- vector()
for (gene in UGT_genes[which(UGT_genes!='UGT2A3')]){
  data <- read_delim(paste0('~/Desktop/UGT_genetic_profiling_KI/raw-data/ClinVar_data/clinvar_variants_', gene, '.txt'), delim='\t', show_col_types = FALSE)
  ## Add variant ID
  data$Ref <- sapply(data$Name, function(x){substr(strsplit(x, '>')[[1]][1], nchar(strsplit(x, '>')[[1]][1]), nchar(strsplit(x, '>')[[1]][1]))})
  data$Alt <- sapply(data$Name, function(x){substr(strsplit(x, '>')[[1]][2], 1,1)})
  data$Variant_ID <- paste(data$GRCh37Chromosome, data$GRCh37Location, data$Ref, data$Alt, sep='-')
  ## Add D/N effect
  data$effect <- sapply(data$`Clinical significance (Last reviewed)`, function(x){
    if (length(grep('Pathogenic|pathogenic', x))!=0) {'D'}
    else if(length(grep('Benign|benign', x))!=0) {'N'}
  })
  assign(paste0('clinvar_variants_', gene), data)
  print(paste0(dim(data)[1], ' benchmark variants in ', gene))
  ## Generate input file to run predictions in ANNOVAR for these variants
  benchmark_whole_data <-rbind(benchmark_whole_data, data[,c('Variant_ID', 'effect')])
  benchmark_data <- rbind(benchmark_data, data[, c('GRCh37Chromosome', 'GRCh37Location', 'GRCh37Location', 'Ref', 'Alt')])
}

## "18 benchmark variants in UGT1A1"
## "20 benchmark variants in UGT1A3"
## "24 benchmark variants in UGT1A4"
## "27 benchmark variants in UGT1A5"
## "30 benchmark variants in UGT1A6"
## "37 benchmark variants in UGT1A7"
## "38 benchmark variants in UGT1A8"
## "37 benchmark variants in UGT1A9"
## "3 benchmark variants in UGT2A1"
## "2 benchmark variants in UGT2A2"
## "1 benchmark variants in UGT2B7"
## "1 benchmark variants in UGT2B10"
## "1 benchmark variants in UGT2B11"
## "3 benchmark variants in UGT2B15"
## "5 benchmark variants in UGT2B17"
## "8 benchmark variants in UGT2B28"
## "1 benchmark variants in UGT3A1"
## "2 benchmark variants in UGT3A2"
## "3 benchmark variants in UGT8"

## ANNOVAR input 
benchmark_data <- unique(as.data.frame(benchmark_data))
colnames(benchmark_data) <- c('Chromosome', 'Start', 'End', 'Ref', 'Obs')
save(benchmark_data, file ='processed-data/03_Anno_functional_impact/benchmark_data.Rdata')
write.table(benchmark_data, file ='processed-data/03_Anno_functional_impact/benchmark_data.txt', row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(benchmark_data, file ='processed-data/03_Anno_functional_impact/benchmark_data.csv', row.names = FALSE, col.names = FALSE, sep = '\t')

## ANNOVAR output
benchmark_scores <- as.data.frame(read_csv('processed-data/03_Anno_functional_impact/ANNOVAR_output/myanno_benchmark.hg19_multianno.csv'))
colnames(benchmark_scores) <- lapply(strsplit(colnames(benchmark_scores), '\\.\\.\\.'), function(x){x[[1]]})

## Scores of interest
scores_algorithms <- c('SIFT_score', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'LRT_score','MutationAssessor_score', 'FATHMM_score', 
                       'fathmm-MKL_coding_score', 'PROVEAN_score', 'VEST3_score', 'VEST4_score', 'CADD_phred', 'DANN_score', 'MetaSVM_score', 'MetaLR_score', 
                       'REVEL_score', 'PrimateAI_score', 'M-CAP_score', 'ClinPred_score', 'Eigen-PC-raw_coding', 'MutPred_score', 'MVP_score')
benchmark_scores <- benchmark_scores[,scores_algorithms]
## Add variant ID and effect 
benchmark_whole_data <- unique(as.data.frame(benchmark_whole_data))
benchmark_scores <- cbind(benchmark_whole_data, benchmark_scores)
colnames(benchmark_scores)[c(9, 13, 19, 21)] <- c('fathmm.MKL_score', 'CADD_phred_score', 'M.CAP_score', 'Eigen.PC_score')

## Add ADME scores
ADME_pred<- data.frame(matrix(nrow=dim(benchmark_scores)[1], ncol=length(ADME_thresholds)+1))
colnames(ADME_pred) <- c('Variant_ID', paste0(names(ADME_thresholds), '_pred'))
ADME_pred$Variant_ID <- benchmark_scores$Variant_ID

for(algorithm in names(ADME_thresholds)){
  ## Evaluate if the algorithm score of each variant passes cutoff (1) or not (0)
  ADME_pred[paste0(algorithm, '_pred')] <- apply(benchmark_scores, 1, 
                                                                    function(x){if (x[paste0(algorithm, '_score')]=='.'){'.'}
                                                                      else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', ADME_thresholds[[algorithm]]))) ){1}
                                                                      else{0} })
}

## Add global ADME scores for each variant 
benchmark_scores$ADME_score <- signif(apply(ADME_pred[,-1], 1, function(x){mean(as.numeric(x[which(x!='.')]))}), digits=3)
benchmark_scores$ADME_score[which(is.nan(benchmark_scores$ADME_score))] <- '.'

## AlphaMissense scores (searching in canonical txs of UGT genes)
benchmark_scores$AlphaMissense_score <- unlist(sapply(benchmark_scores$Variant_ID, function(x){AM_score_pred(x)[1]}))

## Categorical predictions 
benchmark_pred <- data.frame(matrix(nrow=dim(benchmark_scores)[1], ncol=length(names(algorithms_thresholds))+1))
colnames(benchmark_pred) <- c('Variant_ID', paste0(names(algorithms_thresholds), '_pred'))
benchmark_pred$Variant_ID <- benchmark_scores$Variant_ID

for(algorithm in names(algorithms_thresholds)){
  benchmark_pred[paste0(algorithm, '_pred')] <- apply(benchmark_scores, 1, 
                                                               function(x){if (x[paste0(algorithm, '_score')]=='.' | x[paste0(algorithm, '_score')]=='-'){'.'}
                                                                 else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', algorithms_thresholds[[algorithm]]))) ){'D'}
                                                                 else{'N'} })
}

benchmark_scores_preds <- cbind(benchmark_scores, benchmark_pred[,-1])


## Determine if predictions are TP, TN, FP or FN
positive_negative_predictions <- data.frame(matrix(ncol=22, nrow=dim(benchmark_scores_preds)[1]))
colnames(positive_negative_predictions) <- colnames(benchmark_pred)[-1]
for (algorithm in colnames(benchmark_pred)[-1]){
  positive_negative_predictions[,algorithm] <- unlist(apply(benchmark_scores_preds, 1, function(x){if(x[algorithm]=='D' & x['effect']=='D'){'TP'}
    else if(x[algorithm]=='N' & x['effect']=='D'){'FN'}
    else if(x[algorithm]=='D' & x['effect']=='N'){'FP'}
    else if(x[algorithm]=='N' & x['effect']=='N'){'TN'}
    else if(x[algorithm]=='.'){'NA'}}))
}

## Specificity and sensitivity (at the given thresholds) of each method
sensitivities <- apply(positive_negative_predictions, 2, function(x){
                                      length(which(x=='TP'))/(length(which(x=='TP')) + length(which(x=='FN')))})
specificities <- apply(positive_negative_predictions, 2, function(x){
  length(which(x=='TN'))/(length(which(x=='TN')) + length(which(x=='FP')))})

## TP, TN, FP, FN with ADME optimized thresholds
ADME_pred$effect <- benchmark_scores_preds$effect
ADME_predictions <- data.frame(matrix(ncol=5, nrow=dim(ADME_pred)[1]))
colnames(ADME_predictions) <- colnames(ADME_pred)[2:6]
for (algorithm in colnames(ADME_pred)[2:6]){
  ADME_predictions[,algorithm] <- unlist(apply(ADME_pred, 1, function(x){if(x[algorithm]==1 & x['effect']=='D'){'TP'}
    else if(x[algorithm]==0 & x['effect']=='D'){'FN'}
    else if(x[algorithm]==1 & x['effect']=='N'){'FP'}
    else if(x[algorithm]==0 & x['effect']=='N'){'TN'}
    else if(x[algorithm]=='.'){'NA'}}))
}

ADME_sensitivities <- apply(ADME_predictions, 2, function(x){
  length(which(x=='TP'))/(length(which(x=='TP')) + length(which(x=='FN')))})
ADME_specificities <- apply(ADME_predictions, 2, function(x){
  length(which(x=='TN'))/(length(which(x=='TN')) + length(which(x=='FP')))})
names(ADME_sensitivities) <- names(ADME_specificities) <- sapply(names(ADME_sensitivities), function(x){strsplit(x, '_')[[1]][1]})
names(ADME_sensitivities)[4] <- names(ADME_specificities)[4] <- 'VEST4'


## ROC curves
r <- list()
for (algorithm in paste0(names(algorithms_thresholds), '_score')){
  r[[algorithm]] <- roc(response=as.factor(benchmark_scores_preds$effect), 
                        predictor=as.numeric(benchmark_scores_preds[,algorithm]), 
                        levels=c('N', 'D'), na.rm=TRUE)
}
## Add AUC per method
data <- as.data.frame(cbind('AUC'=paste0('AUC = ', signif(as.numeric(lapply(r, function(x){x$auc})), digits=3))))
## Add number of D and N variants used to evaluate each method
data$num_vars <- sapply(colnames(benchmark_pred)[2:23], function(x){
                       DN_num <- table(benchmark_scores_preds[which(benchmark_scores_preds[,x]!='.'), 'effect'])
                       paste0('n = ', DN_num['D'], ' D; ', DN_num['N'], ' N')
                       })
## Locate coordinate corresponding to the used threshold for each method
data$sensitivity <- sensitivities
data$specificity <- specificities

names(r) <- gsub(' phred', '', gsub('_', ' ', gsub('\\.', '-', gsub('_score', '', names(r)))))
data$name <- names(r)

## Add coordinates for ADME thresholds
data$ADME_sensitivity <- as.numeric(sapply(data$name, function(x){if(x %in% names(ADME_sensitivities)){ADME_sensitivities[x]} else {'NA'}}))
data$ADME_specificity <- as.numeric(sapply(data$name, function(x){if(x %in% names(ADME_specificities)){ADME_specificities[x]} else {'NA'}}))

ggroc(r) + 
  facet_wrap(~name) +
  theme_bw() + theme(legend.position = "none") + 
  geom_text(data = data, aes(0, 0.19, label= AUC, hjust = 1), size=3.2, fontface='bold') +
  theme(strip.background = element_rect(fill="gray95", size=1, color="gray60"),
        strip.text = element_text(face="bold"),
        axis.text = element_text( size = 6)) +
  geom_text(data = data, aes(0, 0.05, label= num_vars, hjust = 1), size=2.5, color='black') +
  scale_color_manual(values=colors) +
  ## Point corresponding to used threshold
  geom_point(data=data, aes(x=specificity, y=sensitivity)) +
  ## Point for ADME optimized thresholds
  geom_point(data=data, aes(x=ADME_specificity, y=ADME_sensitivity), shape=5, color='red', size=1.3, stroke = 1)
ggsave(filename='plots/03_Anno_functional_impact/AUC_ROC_methods.pdf', width = 8, height = 8)



# ____________________________________________________________________________________________
#  3.1.6  Development of a UGT-optimized prediction framework
# ____________________________________________________________________________________________

## Function to compute Youden index (J) for all potential thresholds of an algorithm

Youden_indices <- function(method){
  
  method_score <- paste0(method, '_score')
  method_pred <- paste0(method, '_pred')
  method_name <- gsub(' phred', '', gsub('\\.', '-', gsub('_', ' ', gsub('_pred', '', method))))
  
  ## Method scores for benchmark variants
  method_data <- benchmark_scores_preds[which(benchmark_scores_preds[,method_score]!='.'), c(method_score, 'effect')]
  
  ## Conventional threshold
  conv_threshold <- as.numeric(gsub('[<, =, >]', '', algorithms_thresholds[[method]]))
  ## J at the conventional threshold
  J_conventional <- sensitivities[method_pred] + specificities[method_pred] -1 
  
  ## Define potential thresholds (scores)
  scores <- seq(from=min(as.numeric(method_data[,method_score])), to=max(as.numeric(method_data[,method_score])), length.out=100)
  ## Add conventional threshold, max and min scores if not included 
  if (!conv_threshold %in% scores){scores <- append(scores, conv_threshold)}
  if (!max(as.numeric(method_data[,method_score])) %in% scores){scores <- append(scores, max(as.numeric(method_data[,method_score])))}
  if (!min(as.numeric(method_data[,method_score])) %in% scores){scores <- append(scores, min(as.numeric(method_data[,method_score])))}
  
  ## Sensitivity, specificity and J at each new threshold
  scores_Js <- vector()
  Js <- vector()
  for (score in scores){
    ## New threshold given by score
    ## Same direction
    direction <- gsub('[0-9, ., =, -]', '', algorithms_thresholds[[method]])
    threshold <- paste0(direction, '(', score, ')')
    
    ## Categorize based on new threshold
    pred_at_threshold <- sapply(method_data[,method_score], function(x){
                      if (eval(parse_expr(paste0('as.numeric(x)', threshold))) ){'D'}
                          else{'N'} })
    pred_vs_real <- as.data.frame(cbind(pred_at_threshold, 'effect'=method_data$effect))
    
    ## Define TP, TN, FP and FN
    pred_vs_real$comparison <- apply(pred_vs_real, 1, function(x){if(x['pred_at_threshold']=='D' & x['effect']=='D'){'TP'}
                                                       else if(x['pred_at_threshold']=='N' & x['effect']=='D'){'FN'}
                                                       else if(x['pred_at_threshold']=='D' & x['effect']=='N'){'FP'}
                                                       else if(x['pred_at_threshold']=='N' & x['effect']=='N'){'TN'} })
    
    ## Sensitivity and specificity at the given threshold
    sensitivity <- length(which(pred_vs_real$comparison=='TP'))/(length(which(pred_vs_real$comparison=='TP')) + length(which(pred_vs_real$comparison=='FN')))
    specificity <- length(which(pred_vs_real$comparison=='TN'))/(length(which(pred_vs_real$comparison=='TN')) + length(which(pred_vs_real$comparison=='FP')))
    
    ## J at the threshold
    J = sensitivity + specificity -1
    
    ## Add if J hasn't been reached before with another threshold; add conventional score
    if (! J %in% Js | score==conv_threshold | score==max(as.numeric(method_data[,method_score]))){
      scores_Js <- rbind(scores_Js, c('score'=score, 'sensitivity'=sensitivity, 'specificity'=specificity, 'J'=J))
    }
    Js <- append(Js, J)
  }
  
  scores_Js <- as.data.frame(scores_Js)
  max_J_data <- scores_Js[which.max(scores_Js$J),]
  max_J <- scores_Js[which.max(scores_Js$J),'J']
  max_J_score <- scores_Js[which.max(scores_Js$J),'score'] ## optimal threshold
  max_J_sensitivity <- scores_Js[which.max(scores_Js$J),'sensitivity']
  max_J_specificity <- scores_Js[which.max(scores_Js$J),'specificity']

  ## J at optimal threshold - J at conventional threshold = delta(J)
  delta_J <- max_J -J_conventional
  
  if(method=='SIFT'){
    nudge_x_c = 0.13
    nudge_y_c = 0.01
    nudge_x_o = 0.1
    nudge_y_o = 0.1
  }
  else if(method=='MutationAssessor'){
    nudge_x_c = -0.2
    nudge_y_c = 0.01
    nudge_x_o = 0.45
    nudge_y_o = 0.15
  }
  else if(method=='FATHMM'){
    nudge_x_c = -0.45
    nudge_y_c = 0.1
    nudge_x_o = 0.45
    nudge_y_o = 0.1
  }
  else if(method=='PROVEAN'){
    nudge_x_c = 0.9
    nudge_y_c = 0.055
    nudge_x_o = -0.5
    nudge_y_o = 0.055
    scores_to_show <- c(1,3,5)
  }
  else if(method=='PrimateAI'){
    nudge_x_c = -0.1
    nudge_y_c = 0
    nudge_x_o = 0.1
    nudge_y_o = 0.055
    scores_to_show <- c(1,3,5)
  }
  else if(method=='MutPred'){
    nudge_x_c = -0.1
    nudge_y_c = 0.055
    nudge_x_o = 0.1
    nudge_y_o = 0.055
  }
  else if(method=='MetaSVM'){
    nudge_x_c = 0.07
    nudge_y_c = 0.07
    nudge_x_o = 0.1
    nudge_y_o = 0.2
  }
  else if(method=='ClinPred' | method=='CADD_phred'){
    nudge_x_c = 0.06
    nudge_y_c = -0.05
    nudge_x_o = 0.09
    nudge_y_o = 0.16
  }
  else if(method=='Eigen.PC'){
    nudge_x_c = -0.2
    nudge_y_c = 0.09
    nudge_x_o = 0.13
    nudge_y_o = 0.15
  }
  else if(method=='MVP'){
    nudge_x_c = 0.01
    nudge_y_c = 0.09
    nudge_x_o = 0.03
    nudge_y_o = 0.19
  }
  else if(method=='DANN'){
    nudge_x_c = -0.1
    nudge_y_c = 0.03
    nudge_x_o = -0.05
    nudge_y_o = 0.1
  }
  else if (method=='AlphaMissense'){
    nudge_x_c = 0.06
    nudge_y_c = 0.05
    nudge_x_o = -0.1
    nudge_y_o = 0.06
  }
  else{
    nudge_x_c = -0.1
    nudge_y_c = -0.1
    nudge_x_o = 0.15
    nudge_y_o = 0.1
    scores_to_show <- c(1,5)
  }
  

  ## Plot scores vs J 
  p <- ggplot(scores_Js, aes(x=score, y=J))+
    geom_line(color=colors[method_name], size=1) + 
    ## Max J with optimized threshold
    geom_point(x=max_J_score, y=max_J, color='red') + 
    geom_text_repel(data = subset(scores_Js, J==max_J & score!=conv_threshold), label=signif(max_J, digits=3), 
                    size=2.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=1, 
                    box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = nudge_x_o, nudge_y = nudge_y_o, fontface='bold') + 
    ## Line for score that yields the max J
    geom_segment(aes(x = max_J_score, y = 0, xend = max_J_score, yend = max_J), linetype=1, linewidth=0.3, color='grey60') +
    
    ## J with conventional threshold
    geom_point(x=conv_threshold, y=J_conventional, color=colors[method_name]) +
    geom_text_repel(data = subset(scores_Js, score==conv_threshold), label=signif(J, digits=3), 
                    size=2.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=0, vjust=1, 
                    box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = nudge_x_c, nudge_y = nudge_y_c, fontface='bold') + 
    ## Line for conventional score 
    geom_segment(aes(x = conv_threshold, y = 0, xend = conv_threshold, yend = J_conventional), linetype=3, linewidth=0.6, color='grey60') +
    theme_classic() +
    coord_cartesian(ylim=c(0, NA), expand = FALSE) +
    scale_x_continuous(breaks = sort(c(min(as.numeric(scores_Js$score)), max(as.numeric(scores_Js$score)),
                                       signif(c(seq(from=min(as.numeric(scores_Js$score)), 
                                                  to=max(as.numeric(scores_Js$score)), 
                                                  length.out=5)[-scores_to_show], conv_threshold, max_J_score), digits=2)))) +
    labs(x='Score', y='Youden index (J)', title=method_name) +
    ## Label for delta J
    geom_label(x=(max(as.numeric(method_data[,method_score]))+min(as.numeric(method_data[,method_score])))/2, y=0.11, 
               label=paste0('ΔJ = ', signif(delta_J, digits=2)), 
               size=2.8, color='grey30', label.size = NA, fontface='bold', label.padding = unit(0.1, "lines"))+
    theme(title = element_text(size = (9), face='bold'),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = (8.5), face='bold'),
          axis.text = element_text(size = (8)), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6))
  
  if(method=='DANN'){
    p <- p + geom_text_repel(data = subset(scores_Js, score==max_J_score), label=signif(max_J_score, digits=3), 
                             aes(x=max_J_score, y=0.1),
                             size=2.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=0, vjust=1, 
                             box.padding = 0.5, lineheight=unit(1, 'lines'), segment.curvature = 0.2)
  }
  
  if(method=='LRT'){
    p <- ggplot(scores_Js, aes(x=score, y=J))+
      geom_line(color=colors[method_name], size=1) + 
      geom_point(x=max_J_score, y=max_J, color='red') + 
      geom_text_repel(data = subset(scores_Js, J==max_J), label=signif(max_J, digits=3), 
                      size=2.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=1, 
                      box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = 0.1, nudge_y = 0, fontface='bold') + 
      geom_segment(aes(x = max_J_score, y = 0, xend = max_J_score, yend = max_J), linetype=1, linewidth=0.3, color='grey60') +
      theme_classic() +
      coord_cartesian(ylim=c(-0.01, max(as.numeric(scores_Js$J))+0.05), 
                      xlim=c(-0.01, max(as.numeric(scores_Js$score))),
                      expand = FALSE) +
      labs(x='Score', y='Youden index (J)', title=method_name) +
      ## Label for delta J
      geom_label(x=(max(as.numeric(method_data[,method_score]))+min(as.numeric(method_data[,method_score])))/2, y=0.11, 
                 label=paste0('ΔJ = ', signif(delta_J, digits=2)), 
                 size=2.8, color='grey30', label.size = NA, fontface='bold', label.padding = unit(0.1, "lines"))+
      geom_text_repel(data = subset(scores_Js, score==max_J_score), label=signif(max_J_score, digits=3), 
                      aes(x=max_J_score, y=0.1),
                      size=2.5, color='grey40', min.segment.length = unit(0, 'lines'), hjust=0.5, vjust=1, 
                      box.padding = 0.5, lineheight=unit(1.2, 'lines'), segment.curvature = -0.2) +
      theme(title = element_text(size = (9), face='bold'),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = (8.5), face='bold'),
            axis.text = element_text(size = (8)), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6))
    
  }
  
  return(list(p, max_J_data))
}

i=1
plots <- list()
for (method in names(algorithms_thresholds)){
  plots[[i]] <- Youden_indices(method)[[1]]
  i=i+1
}

## Don't include PolyPhen2 HVAR
plot_grid(plots[[1]], plots[[2]], plots[[4]], plots[[5]], plots[[6]],
          plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]],
          plots[[13]], plots[[14]], plots[[15]], plots[[16]], plots[[17]], plots[[18]], 
          plots[[19]], plots[[20]], plots[[21]], plots[[22]], ncol=7)
ggsave(filename='plots/03_Anno_functional_impact/Youden_Index_plots.pdf', width = 26, height = 10)


## MCAP new threshold in line, digits in intermediate x scores, delta symbol, red rhombus 



## Add coordinate for new thresholds in ROC curves



## Predict variant effect by all algorithms using these new thresholds


## Number of missing scores per variant
table(apply(new_variants_predictions[,24:45], 1, function(x){length(which(x=='.'))}))
#    0    1    2    3    4    5    7    8    9   10   21 
# 2827 1383 1543  346  180   46    2   17    6    1    1 

## Consensus on prediction 



## UGT-optimized score as the mean of all algorithms' predictions for each variant


## D, N and M variants with this framework



################################################################################
##          3.2  Annotate functional consequence of all UGT variants 
################################################################################

## Define categories of predicted effect of exonic variant types
deleterious <- c('splice_donor_variant', 'splice_acceptor',
                 'stop_gained', 'frameshift_variant', 'start_lost')

neutral <- c('synonymous_variant', 'inframe_deletion', 'inframe_insertion', 
             'stop_retained_variant', 'splice_region_variant', 'stop_lost')

## Take missense predictions by MVP (61 missing variants)
valid_MVP_preds <- new_variants_predictions[which(new_variants_predictions$MVP_pred!='.'),c('Variant_ID', 'MVP_pred')]
rownames(valid_MVP_preds) <- valid_MVP_preds$Variant_ID

## Annotate functional impact of all exonic (+ 3 promoter variants in UGT1A1) variants in each gene
for (gene in UGT_genes){
  exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  exonic_data$Functional_impact <- sapply(exonic_data$Variant_ID, function(x){if (exonic_data[x, 'VEP_Annotation'] %in% deleterious){'D'}
    else if (exonic_data[x, 'VEP_Annotation'] %in% neutral){'N'}
    ## Annotate MVP prediction for missense variants
    else if (x %in% valid_MVP_preds$Variant_ID){valid_MVP_preds[x, 'MVP_pred']}
    else{'.'} })
  assign(paste0(gene, '_exonic_data'), exonic_data)
}

## Annotate effect of 5' upstream UGT1A1 variants according to ClinVar
UGT1A1_exonic_data[which(UGT1A1_exonic_data$VEP_Annotation=='5\' upstream'), 'Functional_impact'] <- c('D', 'N', 'D')


## Plot MAF of all D and N variants per gene in each population

## Define MAF in each population
populations <- c('African_or_African_American',
                 'Latino_or_Admixed_American',
                 'Ashkenazi_Jewish',
                 'East_Asian',
                 'European_Finnish',
                 'European_non_Finnish',
                 'South_Asian',
                 'Other')

plots <- list()
i=1
for (gene in UGT_genes){

  ## Exonic variants per gene
  exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  ## Number of missense variants missed by MVP prediction
  num_missed <- length(which(exonic_data$Functional_impact=='.'))
  ## Subset to variants with valid functional consequence 
  exonic_data <- exonic_data[which(exonic_data$Functional_impact!='.'),]
  
  for (p in populations){
    exonic_data[,paste0('MAF_', p)] <- exonic_data[,paste0('Allele_Count_', p)]/exonic_data[,paste0('Allele_Number_', p)]
    ## Define MAF 
    exonic_data[which(exonic_data[,paste0('MAF_', p)]>=0.5), paste0('MAF_', p)] <- 1-exonic_data[which(exonic_data[,paste0('MAF_', p)]>=0.5), paste0('MAF_', p)]
    ## If no MAF available in the population
    nan <- which(is.nan(exonic_data[,paste0('MAF_', p)]))
    #print(c(gene, p, length(nan)))
    exonic_data[nan, paste0('MAF_', p)] <- NA
  }
  ## Global MAF
  exonic_data[which(exonic_data$Allele_Frequency>=0.5), 'Allele_Frequency'] <- 1-exonic_data[which(exonic_data$Allele_Frequency>=0.5), 'Allele_Frequency']
  
  groups <- c(populations, 'global')
  data <- exonic_data[,c('Variant_ID', 'Functional_impact', 'Allele_Frequency', paste0('MAF_', populations))]
  data <- melt(data, id.vars = c('Variant_ID', 'Functional_impact'))
  colnames(data) <- c('Variant_ID', 'Functional_impact', 'Group', 'MAF')
  
  ## Different shapes for D variants with MAF>=0.01
  shapes <- c('2-234668879-C-CAT'=8,
              '2-234668879-C-CATAT'=11,
              '2-234675779-A-G'=14,
              '2-234676872-C-T'=9,
              '2-234637917-C-T'=7,
              '2-234638282-G-GT'=25,
              '2-234622331-GC-G'=3,
              '2-234590935-G-T'=15,
              '2-234545998-G-A'=2,
              '4-70512787-A-T'=12,
              '4-69811110-A-C'=6,
              '4-69693141-GT-G'=5,
              '4-69693242-T-C'=4,
              '4-70070366-A-T'=10,
              '4-70078393-C-T'=13,
              '4-69512937-T-A'=17,
              '4-69528742-G-A'=18,
              '4-69536234-G-T'=25)
  
  ## Label those variants
  data$Label <- apply(data, 1, function(x){if (is.na(x['MAF'])){NA}
                                           else if(x['Functional_impact']=='D' & as.numeric(x['MAF'])>=0.01){x['Variant_ID']} else{NA}})
  data$Functional_impact <- factor(data$Functional_impact, levels=c('N', 'D'))
  data$Group <- factor(data$Group, levels=c(paste0('MAF_',populations), 'Allele_Frequency'))
  
  plots[[i]] <- ggplot(data = data, mapping = aes(x = Group, y = MAF, color = Functional_impact)) +
        geom_point(data=subset(data, is.na(Label)), alpha = 0.9, size = 1.3, position = position_jitterdodge(seed=2)) +
        geom_point(data=subset(data, !is.na(Label)), aes(shape=Label), size=1.5, color='tomato3', stroke = 0.7, position = position_jitterdodge(seed=2, jitter.width=0.1, dodge.width = 0.5)) +
        theme_bw() +
        scale_color_manual(values = c("skyblue2", "tomato"), labels = c("Neutral", "Deleterious")) +
        scale_shape_manual(values=shapes[subset(data, !is.na(Label))$Variant_ID]) + 
        scale_x_discrete(breaks=c(paste0('MAF_',populations), 'Allele_Frequency'),
                         labels=c("African/African American",
                                  "Latino/Admixed American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European Finnish",
                                  "European non Finnish",
                                  "South Asian",
                                  "Other",
                                  "Global")) +
        labs(title=gene, 
             subtitle=paste0(table(exonic_data$Functional_impact)['D'], ' D variants; ', 
                             table(exonic_data$Functional_impact)['N'], ' N variants; ',
                             num_missed, ' missense variants missed by MVP'), 
             x='', y='MAF of variants', color='Functional impact', shape=paste0('Deleterious variant ID', '\n', '(MAF>0.01)')) +
        theme(title = element_text(size = (9), face='bold'),
              plot.subtitle = element_text(size = 8.5, color = "gray50"),
              axis.title = element_text(size = (8.5), face='bold'),
              axis.text = element_text(size = (8)),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
              legend.title = element_text(size=8.5), 
              legend.text = element_text(size=8))
  i=i+1

}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
          plots[[6]], plots[[7]], plots[[8]], plots[[9]], ncol=3)
ggsave(filename='plots/03_Anno_functional_impact/MAF_pop_vars_per_UGT1.pdf', width = 20, height = 13)

plot_grid(plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]],
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], ncol=2)
ggsave(filename='plots/03_Anno_functional_impact/MAF_pop_vars_per_UGT2.pdf', width = 15, height = 23)

plot_grid(plots[[20]], plots[[21]], ncol=2)
ggsave(filename='plots/03_Anno_functional_impact/MAF_pop_vars_per_UGT3.pdf', width = 15, height = 5)

plot_grid(plots[[22]])
ggsave(filename='plots/03_Anno_functional_impact/MAF_pop_vars_per_UGT8.pdf', width = 8, height = 5)



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
# date     2023-10-09
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
# glue          1.6.2   2022-02-24 [1] CRAN (R 4.3.0)
# gtable        0.3.3   2023-03-21 [1] CRAN (R 4.3.0)
# here        * 1.0.1   2020-12-13 [1] CRAN (R 4.3.0)
# hms           1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
# labeling      0.4.2   2020-10-20 [1] CRAN (R 4.3.0)
# lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.3.0)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
# munsell       0.5.0   2018-06-12 [1] CRAN (R 4.3.0)
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
