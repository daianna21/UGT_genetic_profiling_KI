
library(here)
library(rlang)
library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(corrplot)
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

## Load exonic data for each gene 
for (gene in UGT_genes){
  load(here(paste0('~/Desktop/UGT_genetic_profiling_KI/processed-data/01_Data_Processing/', gene, '_exonic_data.Rdata')),
       verbose=TRUE)
  assign( paste0(gene, '_exonic_data'), exonic_vars)
}

## Define categories of predicted effect of exonic variant types
deleterious <- c('splice_donor', 'splice_acceptor', 'stop_lost', 'stop_gained', 'frameshift')
neutral <- c('synonymous')
## missense functions with ANOVA
# deletereous?


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

############################ 1. Confirm all variants are annotated as exonic  ############################# 

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


################################ 2.  Check all variants are non-synonymous  ################################

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
                       'fathmm.MKL_coding_score', 'PROVEAN_score', 'VEST3_score', 'VEST4_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'MetaSVM_score', 'MetaLR_score', 
                       'REVEL_score', 'PrimateAI_score', 'M.CAP_score', 'ClinPred_score', 'Eigen.PC.raw_coding', 'MutPred_score', 'MVP_score')
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
      ## For shared variants, assume the predictions are the same in all genes and take the ones reported in the first one
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
colnames(variants_predictions)[c(9, 13, 14, 22, 30)] <- c('fathmm.MKL_score', 'CADD_score', 'CADD_phred_score', 'Eigen.PC_score', 'fathmm.MKL_pred')


## Add ADME-optimized model scores

## ADME-optimized algorithm thresholds to categorize variants 
ADME_thresholds <- list('SIFT'='<0.0376', 
                        'Polyphen2_HDIV'='>0.3841',
                        'LRT'='<0.0025',
                        'MutationAssessor'='>2.0566',
                        'FATHMM'='<0.486',
                        'fathmm.MKL'='>0.3982',
                        'PROVEAN'='< -3.286',
                        'VEST3'='>0.4534',
                        'CADD_phred'='>19.19',
                        'DANN'='>0.9688',
                        'MetaSVM'='> -0.3371',
                        'MetaLR'='>0.4039')

## Categorize variants in D/N by each algorithm using these thresholds
ADME_categorical_predictions <- data.frame(matrix(nrow=dim(variants_predictions)[1], ncol=13))
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

# ___________________________________________________________________________________________
#  3.1.4  Comparison and evaluation of predictive algorithms 
# ___________________________________________________________________________________________

####################  3.1.4.1 Compare predictions of different algorithms  ####################

#############################################################################
##   Compare binary categorical predictions of variants in all UGT genes
#############################################################################

## Function to create density plot of raw scores for variants in the different functional predicted categories

score_density_plot <- function(algorithm, predicted_cat_type){
  
  algorithm_score <- paste0(algorithm, '_score')
  algorithm_pred <- paste0(algorithm, '_pred')
  
  ## Define colors for predicted effect categories 
  colors <- list('D'='tomato2', 'T'='skyblue1', 'N'='skyblue1', 'B'='skyblue1', 'P'='lightsalmon', 
                 'H'= 'red4', 'M'='red3', 'L'='dodgerblue3', 'U'='grey90')
  
  if (predicted_cat_type=='new'){
    
    threshold <- algorithms_thresholds[[algorithm]]
    numeric_threshold <- as.numeric(gsub('[<, =, >]', '', threshold))
    
    ## Subset to variants with algorithm scores/predictions
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
    
    p1 <- ggplot(data = df, aes(x = x, ymin = 0, ymax = y, fill = pred)) +
      geom_ribbon(alpha=0.7) +
      theme_bw() +
      scale_fill_manual(values=colors[names(table(df$pred))]) +
      labs(x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm)), ' raw scores'), y= 'Density', fill='Predicted effect',
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
                'fathmm.MKL', 'PROVEAN', 'MetaSVM', 'MetaLR', 'M.CAP', 'ClinPred', 'LRT', 'PrimateAI')

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
          plots[[23]], plots[[24]], plots[[25]], plots[[26]],
          ncol=6, rel_widths = c(rep(c(0.74,1), 13)))

ggsave(filename='plots/03_Anno_functional_impact/Returned_RawScores_density_plots.pdf', width = 30, height = 15)



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
                              'CADD'='>0.73',              
                              'DANN'='>0.96', 
                              'REVEL'='>0.5',               
                              'Eigen.PC'='>=0', 
                              'MVP'= '>0.75',
                              'LRT'='<0.001',
                              'MutPred'='>0.5',
                              'PrimateAI'='>=0.803',
                              'VEST4'='>0.5',
                              'ADME'='>0.5')

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
## Data frame with scores and new binary predictions per algorithm 
new_variants_predictions <- cbind(categorical_predictions, variants_predictions[,paste0(names(algorithms_thresholds), '_score')])
new_variants_predictions$MutPred_score <- replace(new_variants_predictions$MutPred_score, which(new_variants_predictions$MutPred_score=='-'), '.')


## Distribution of raw scores in D and N variants defined by cutoffs
plots <- list()
for (i in 1:length(names(algorithms_thresholds))){
  plots[[i]] <- score_density_plot(names(algorithms_thresholds)[i], 'new')
}

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]],
          ncol=7)

ggsave(filename='plots/03_Anno_functional_impact/New_RawScores_density_plots.pdf', width = 30, height = 9)


# ------------------------------------------------------------------------------
## Correlation between raw scores from each pair of methods
raw_scores <- as.data.frame(apply(new_variants_predictions[,paste0(names(algorithms_thresholds), '_score')], 2, as.numeric))
corr <- matrix(nrow=21, ncol = 21)
colnames(corr) <- rownames(corr) <- colnames(raw_scores)

for (i in 1:length(colnames(raw_scores))){
  for (j in 1:length(colnames(raw_scores))){
    ## Subset to variants with valid scores in both algorithms
    raw_scores_subset <- raw_scores[which(!is.na(raw_scores[,colnames(raw_scores)[i]]) & !is.na(raw_scores[,colnames(raw_scores)[j]])),]
    corr[i, j] <- signif(cor(raw_scores_subset[,colnames(raw_scores)[i]], raw_scores_subset[,colnames(raw_scores)[j]], method = 'pearson'), digits=2)
  }
}


whole_corr <- corr
colnames(corr) <- rownames(corr) <- gsub('\\.','-', gsub('_', ' ', gsub('_score', '', colnames(corr))))
corr[lower.tri(corr)] <- NA
half_corr_data <- melt(as.data.frame(corr), na.rm = TRUE)
half_corr_data$value <- signif(as.numeric(half_corr_data$value), digits = 3)

## Mean corr coeff
unique_half_corr_data <- half_corr_data[which(half_corr_data$value!=1), ]
mean(unique_half_corr_data$value)
# [1] 0.3918721

## Highest corr coeffs
unique_half_corr_data[order(unique_half_corr_data$value, decreasing = TRUE), ][1:4,]
#               Var1             Var2    value
#     Polyphen2 HDIV   Polyphen2 HVAR    0.890
#            MetaSVM           MetaLR    0.744
#          fathmm-MKL         Eigen-PC   0.691
#            MetaSVM            REVEL    0.689

## Percentage of high coeffs (>0.5)
length(which(unique_half_corr_data$value>0.5))/dim(unique_half_corr_data)[1]*100
# [1] 27.20588

## Percentage of medium coeffs (0.3<r<0.5)
length(which(unique_half_corr_data$value>0.3 & unique_half_corr_data$value<0.5))/dim(unique_half_corr_data)[1]*100
# [1] 43.38235

## Percentage of low coeffs (<0.3)
length(which(unique_half_corr_data$value<0.3))/dim(unique_half_corr_data)[1]*100
# [1] 29.41176




# ------------------------------------------------------------------------------
## Agreement proportion between predictions from each pair of methods
predictions <- as.data.frame(new_variants_predictions[,paste0(names(algorithms_thresholds), '_pred')])
agreement_prop <- matrix(ncol=21, nrow=21)
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
colnames(agreement_prop) <- rownames(agreement_prop) <- gsub('\\.','-', gsub('_', ' ', gsub('_pred', '', colnames(agreement_prop))))
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
         col.lim = c(-1,1),
         col=colorRampPalette(c('dodgerblue4', 'deepskyblue3', 'lightskyblue1', "white",'mistyrose2',"tomato2", 'firebrick4'))(100)
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
# [1] 0.6993382

## Highest corr coeffs
unique_half_agreement_prop[order(unique_half_agreement_prop$value, decreasing = TRUE), ][1:4,]
#               Var1             Var2    value
# 80           VEST3         FATHMM      0.963
# 20  Polyphen2 HVAR Polyphen2 HDIV      0.944
# 128         MetaLR        MetaSVM      0.929
# 134          REVEL        MetaSVM      0.924

## Percentage of high prop (>0.5)
length(which(unique_half_agreement_prop$value>0.5))/dim(unique_half_agreement_prop)[1]*100
# [1] 92.64706


# ------------------------------------------------------------------------------
## Scatterplot of raw scores and categorical predictions of 2 methods

scatterplot_compare_2methods <- function(algorithm1, algorithm2){
  
  algorithm1_score <- paste0(algorithm1, '_score')
  algorithm2_score <- paste0(algorithm2, '_score')
  algorithm1_pred <- paste0(algorithm1, '_pred')
  algorithm2_pred <- paste0(algorithm2, '_pred')
  
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
                            else if(x[1]=='D' & x[2]=='N'){paste0('D in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[1])), '; ', 'N in ', 
                                                                  gsub('_',' ',gsub('\\.', '-', colnames(data)[2])))}
                            else if(x[1]=='N' & x[2]=='D'){paste0('D in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[2])), '; ', 'N in ', 
                                                                  gsub('_',' ',gsub('\\.', '-', colnames(data)[1])))}
                            else {'N in both'}
                      })
  colors <- list()
  colors[['D in both']]='indianred'
  colors[[paste0('D in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[1])), '; ', 'N in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[2])))]]='thistle2'
  colors[[paste0('D in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[2])), '; ', 'N in ', gsub('_',' ',gsub('\\.', '-', colnames(data)[1])))]]='lightpink1'
  colors[['N in both']]='lightblue3'
  
  ggplot(data, aes(x=eval(parse_expr(algorithm1)), y=eval(parse_expr(algorithm2)), color = categories)) +
    geom_point(size=2) +
    stat_smooth(geom = "line", alpha = 0.9, size = 1, method = lm, color = "orangered3", fullrange = FALSE) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      subtitle = paste0("Corr: ", correlation, '; Agreement: ', 100*agreement, '%'),
      x = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm1)), ' raw score'),
      y = paste0(gsub('_', ' ', gsub('\\.', '-', algorithm2)), ' raw score'),
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
scatterplot_compare_2methods('SIFT', 'Polyphen2_HDIV')
scatterplot_compare_2methods('Eigen.PC', 'CADD') **
scatterplot_compare_2methods('Polyphen2_HDIV', 'Polyphen2_HVAR')
scatterplot_compare_2methods('FATHMM', 'fathmm.MKL')
scatterplot_compare_2methods('MetaSVM', 'MetaLR')
scatterplot_compare_2methods('FATHMM', 'MutPred')
scatterplot_compare_2methods('FATHMM', 'CADD')
scatterplot_compare_2methods('FATHMM', 'DANN')
scatterplot_compare_2methods('FATHMM', 'VEST4') ***
scatterplot_compare_2methods('ADME', 'REVEL')
scatterplot_compare_2methods('ADME', 'fathmm.MKL')
scatterplot_compare_2methods('ADME', 'ClinPred')
scatterplot_compare_2methods('M.CAP', 'LRT')
scatterplot_compare_2methods('Eigen.PC', 'M.CAP') ** 






####################  3.1.4.2 Evaluate predictions of different algorithms  ####################

known_exonic_vars <- list('UGT1A1'=list('2-234669144-G-A'= 'reduced UGT1A1 expression'),
                          'UGT1A6'=list('2-234602191-A-G'= 'increased risk for severe neutropenia when treated with irinotecan ', 
                                        '2-234602202-A-C'= 'higher metabolic activity', 
                                        '2-234601669-T-G'= 'increase the likelihood of cardiotoxicity when treated with anticancer anthracyclines'),
                          'UGT1A7'=list('2-234591205-T-C'= 'increases risk of vomiting in colorectal cancer patients treated with anticancer drugs'),
                          'UGT1A8'=list('2-234526871-C-G'= 'increases chances of having diarrhea in kidney transplant patients treated with the immune suppressants'),
                          'UGT2B4'=list('4-70346565-A-T'= ''),
                          'UGT2B7'=list('4-69964338-T-C'= 'reduced response to oxycodone and improved response to oxcarbazepine', 
                                        '4-69962449-G-T'= 'LoF mutation; increases levels of valproic acid in the plasma of epilepsy patients'),
                          'UGT2B15'=list('4-69536084-A-C'= ''))

variant_data <- strsplit('2-234669144-G-A', '-')[[1]]
myanno_UGT1A1[which(myanno_UGT1A1$Chr==variant_data[1] & myanno_UGT1A1$Start==variant_data[2] & myanno_UGT1A1$End==variant_data[2] 
      & myanno_UGT1A1$Ref==variant_data[3] & myanno_UGT1A1$Alt==variant_data[4]), c('SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred',
                                                                                    'MutationAssessor_pred', 'FATHMM_pred', 'fathmm.MKL_coding_pred', 'PROVEAN_pred',
                                                                                    'VEST3_score', 'CADD_raw', 'DANN_score', 'MetaSVM_pred', 'MetaLR_pred')]
      
sapply(UGT_genes, function(gene){names(table(eval(parse_expr(paste0('myanno_', gene, '$Gene.refGene')))))})











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
# date     2023-09-23
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# ape                  * 5.7-1     2023-03-13 [1] CRAN (R 4.3.0)
# beachmat               2.16.0    2023-05-08 [1] Bioconductor
# Biobase              * 2.61.0    2023-06-02 [1] Bioconductor
# BiocGenerics         * 0.47.0    2023-06-02 [1] Bioconductor
# BiocParallel           1.35.3    2023-07-07 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# corrplot             * 0.92      2021-11-18 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6    2023-07-02 [1] Bioconductor
# DelayedMatrixStats     1.23.0    2023-04-25 [1] Bioconductor
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10    2023-05-28 [1] Bioconductor
# GenomicRanges        * 1.53.1    2023-06-02 [1] Bioconductor
# ggplot2              * 3.4.2     2023-04-03 [1] CRAN (R 4.3.0)
# glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# gtable                 0.3.3     2023-03-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# IRanges              * 2.35.2    2023-06-23 [1] Bioconductor
# labeling               0.4.2     2020-10-20 [1] CRAN (R 4.3.0)
# lattice                0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-0     2023-07-08 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                   3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
# pheatmap             * 1.0.12    2019-01-04 [1] CRAN (R 4.3.0)
# phylotools           * 0.2.4     2023-08-31 [1] Github (helixcn/phylotools@758d338)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.8     2022-11-11 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors            * 0.39.1    2023-06-02 [1] Bioconductor
# scales                 1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# scuttle              * 1.9.4     2023-01-23 [1] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.23.0    2023-04-25 [1] Bioconductor
# sparseMatrixStats      1.13.0    2023-05-20 [1] Bioconductor
# stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
# systemfonts            1.0.4     2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6     2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
# vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.3.0)
# withr                  2.5.0     2022-03-03 [1] CRAN (R 4.3.0)
# XVector                0.41.1    2023-06-02 [1] Bioconductor
# zlibbioc               1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────








