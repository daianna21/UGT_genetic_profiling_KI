
library(here)
library(rlang)
library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(corrplot)
library(pROC)
library(ggrepel)
library(ggExtra)
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
canonical_UGT1_txs <- list('UGT1A1'= 'ENST00000305208.5', 'UGT1A3'='ENST00000482026.1', 'UGT1A4'='ENST00000373409.3', 
                           'UGT1A5'='ENST00000373414.3', 'UGT1A6'='ENST00000305139.6', 'UGT1A7'='ENST00000373426.3', 
                           'UGT1A8'= 'ENST00000373450.4','UGT1A9'= 'ENST00000354728.4', 'UGT1A10'='ENST00000344644.5')

canonical_UGT2_txs <- list('UGT2A1'= 'ENST00000503640.1', 'UGT2A2'='ENST00000457664.2', 'UGT2A3'='ENST00000251566.4', 
                           'UGT2B4'='ENST00000305107.6', 'UGT2B7'='ENST00000305231.7', 'UGT2B10'='ENST00000265403.7', 
                           'UGT2B11'= 'ENST00000446444.1', 'UGT2B15'= 'ENST00000338206.5', 'UGT2B17'='ENST00000317746.2', 
                           'UGT2B28'='ENST00000335568.5')

canonical_UGT3_txs <- list('UGT3A1'= 'ENST00000274278.3', 'UGT3A2'='ENST00000282507.3')

canonical_UGT8_txs <- list('UGT8'= 'ENST00000310836.6')



## Load exonic data for each gene 
non_minor_alleles_allGenes <- vector()
for (gene in UGT_genes){
  exonic_vars <- eval(parse_expr(load(here(paste0('~/Documents/KI_projects/UGT_genetic_profiling_KI/processed-data/01_Data_Processing/', gene, '_exonic_data.Rdata')),
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
  write.table(missense_vars, file = paste0('processed-data/03_Anno_functional_impact/', 
                                                          gene, '_missense_variants.csv'), row.names = FALSE, col.names = TRUE, sep = '\t')
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
  save(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/ANNOVAR/input_data/', 
                                                                                      gene, '_missense_vars_ANNOVAR_format.Rdata'))
  write.table(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/ANNOVAR/input_data/', 
                                                   gene, '_missense_vars_ANNOVAR_format.txt'), row.names = FALSE, col.names = FALSE, sep = '\t')
  write.table(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/ANNOVAR/input_data/', 
                                                          gene, '_missense_vars_ANNOVAR_format.csv'), row.names = FALSE, col.names = FALSE, sep = '\t')
}
 

## --> ANNOVAR was run in 3.1_run_ANNOVAR.sh

# _______________________________________________________________________________
#  3.1.2  Examination of ANNOVAR gene-based annotation output 
# _______________________________________________________________________________

## Download ANNOVAR output for each gene
for (gene in UGT_genes){
  myanno <- read.csv(paste0('processed-data/03_Anno_functional_impact/ANNOVAR/output_data/myanno', gene, '.hg19_multianno.csv'))
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


## --> AlphaMissense data for UGT canonical txs were obtained in 3.2_extract_AlphaMissense_scores.sh

# _______________________________________________________________________________
#  3.1.4  Extract AlphaMissense (AM) predictions for all UGT variants
# _______________________________________________________________________________

## Retrieve AM scores and predictions for all variants of a gene (in its canonical tx)
for (gene in UGT_genes){
  tx <- canonical_txs[[gene]]
  data <- read_tsv(paste0('~/Documents/KI_projects/UGT_genetic_profiling_KI/raw-data/AlphaMissense_data/', tx, '_AlphaMissense_data'), show_col_types = FALSE)
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

## Real tool names
tool_names <- c('SIFT'='SIFT',               
               'Polyphen2_HDIV'='PolyPhen-2 HDIV',     
               'Polyphen2_HVAR'='PolyPhen-2 HVAR',     
               'MutationAssessor'='MutationAssessor',     
               'FATHMM'= 'FATHMM',           
               'fathmm.MKL'='FATHMM-MKL',   
               'PROVEAN'= 'PROVEAN',         
               'MetaSVM'='MetaSVM',               
               'MetaLR'='MetaLR',             
               'M.CAP'='M-CAP',            
               'ClinPred'='ClinPred',               
               'CADD_phred'='CADD',              
               'DANN'='DANN', 
               'REVEL'='REVEL',               
               'Eigen.PC'='Eigen-PC', 
               'MVP'= 'MVP',
               'LRT'='LRT',
               'MutPred'='MutPred',
               'PrimateAI'='PrimateAI',
               'VEST4'='VEST',
               'ADME'='ADME',
               'AlphaMissense'='AlphaMissense')

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
    
    if (algorithm=='CADD_phred'){score_type <- ' phred score'}
    else {score_type <- ' raw score'}
    
    p1 <- ggplot(data = df, aes(x = x, ymin = 0, ymax = y, fill = pred)) +
      geom_ribbon(alpha=0.7) +
      theme_bw() +
      scale_fill_manual(values=colors[names(table(df$pred))]) +
      labs(x = paste0(tool_names[algorithm], score_type), y= 'Density', fill='Predicted effect',
           subtitle=paste0('Missingness: ', signif(as.numeric(missingness), digits=3), '%; ',
                           num_vars, ' variants')) +
      geom_line(aes(y = y)) +
      geom_vline(xintercept = numeric_threshold, color = 'indianred3', linetype='dashed', linewidth=0.6) +
      geom_label(aes(x = numeric_threshold, y = max(df$y), color = 'indianred3', label = numeric_threshold), 
               hjust = hjust, vjust = 3, fontface = 2, fill = "white", show.legend = FALSE) 
    
    if(algorithm!='AlphaMissense'){
      p1 <- p1 + theme(legend.position='none',
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       plot.subtitle = element_text(size = 10, color = "gray30"), 
                       plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
                       axis.text = element_text(size = (10)),
                       legend.text = element_text(size = 10),
                       legend.title = element_text(size =11, face='bold'),
                       axis.title.x = element_text(size = (11.5), face='bold'),
                       axis.title.y = element_text(size = (11.5)))
    }
    else{
      p1 <- p1 + theme(legend.key = element_rect(fill = "white", colour = "black"),
                       plot.subtitle = element_text(size = 10, color = "gray30"),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
                       axis.text = element_text(size = (10)),
                       legend.text = element_text(size = 10),
                       legend.title = element_text(size =11, face='bold'),
                       axis.title.x = element_text(size = (11.5), face='bold'),
                       axis.title.y = element_text(size = (11.5)))
    }
    
    return(p1)
  }
  
  
  
  else if(predicted_cat_type=='without_Polyphen2_HVAR'){
    ## HEREE
    
    threshold <- algorithms_thresholds[[algorithm]]
    numeric_threshold <- as.numeric(gsub('[<, =, >]', '', threshold))
    
    data <- new_variants_predictions[which(eval(parse_expr(paste0('new_variants_predictions$',algorithm_score)))!='.'),]
    
    missingness <- apply(new_variants_predictions, 2, function(x){100*length(which(x=='.'))/dim(new_variants_predictions)[1]})[[algorithm_score]]
    num_vars <- apply(new_variants_predictions, 2, function(x){length(which(x!='.'))})[[algorithm_score]]
    
    ## Density function 
    density <- density(as.numeric(data[,algorithm_score]), 
                       from = min(as.numeric(data[,algorithm_score])), 
                       to = max(as.numeric(data[,algorithm_score]))) 
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
      labs(x = tool_names[algorithm], y= 'Density', fill='Predicted effect',
           subtitle=paste0(signif(as.numeric(missingness), digits=3), '% (n=',
                           num_vars, ')')) +
      geom_line(aes(y = y)) +
      geom_vline(xintercept = numeric_threshold, color = 'indianred3', linetype='dashed', linewidth=0.6) +
      geom_label(aes(x = numeric_threshold, y = max(df$y), color = 'indianred3', label = numeric_threshold), 
                 hjust = hjust, vjust = 3, fontface = 2, fill = "white", show.legend = FALSE) 
    
    if(algorithm!='AlphaMissense'){
      p1 <- p1 + theme(legend.position='none',
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       plot.subtitle = element_text(size = 10, color = "gray30"), 
                       plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
                       axis.text = element_text(size = (8)),
                       legend.text = element_text(size = 10),
                       legend.title = element_text(size =11, face='bold'),
                       axis.title.x = element_text(size = (11.5), face='bold'),
                       axis.title.y = element_text(size = (11.5)))
    }
    else{
      p1 <- p1 + theme(legend.key = element_rect(fill = "white", colour = "black"),
                       plot.subtitle = element_text(size = 10, color = "gray30"),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
                       axis.text = element_text(size = (8)),
                       legend.text = element_text(size = 10),
                       legend.title = element_text(size =11, face='bold'),
                       axis.title.x = element_text(size = (11.5), face='bold'),
                       axis.title.y = element_text(size = (11.5)))
    }
    
    return(p1)
  }
  
  
  
  else{
    
    data <- variants_predictions[which(eval(parse_expr(paste0('variants_predictions$',algorithm_score)))!='.'),]
    
    p2 <- ggplot(data = data, aes(x = as.numeric(eval(parse_expr(algorithm_score))))) +
      geom_density(alpha=0.6, fill='grey')+
      theme_bw() +
      labs(x = paste0(tool_names[algorithm], ' raw score'), y= 'Density') 
    
    p3 <- ggplot(data = data, aes(x = as.numeric(eval(parse_expr(algorithm_score))), 
                                  fill=eval(parse_expr(algorithm_pred))))+
      geom_density(alpha=0.6)+
      scale_fill_manual(values=colors[names(table(eval(parse_expr(paste0('data$', algorithm_pred)))))]) +
      theme_bw() +
      labs(x = paste0(tool_names[algorithm], ' raw score'), y= 'Density', fill='Predicted effect') 
    
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

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],<
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

## Shared legend
legend <- get_legend(
  # Create some space to the left of the legend
  plots[[22]] + theme(legend.box.margin = margin(0, 0, 0, 1))
)

plots[[22]] <-  plots[[22]] + theme(legend.position = 'none', 
                                    plot.subtitle = element_text(size = 10, color = "gray30"),
                                    plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'))

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]], plots[[22]], ncol=5, legend)

ggsave(filename='plots/03_Anno_functional_impact/New_RawScores_density_plots.pdf', width = 14.5, height = 13)


## Exclude Polyphen2-HVAR
no_Polyphen2_HVAR <- names(algorithms_thresholds)[names(algorithms_thresholds)!='Polyphen2_HVAR']
plots <- list()

for (i in 1:length(no_Polyphen2_HVAR)){
  plots[[i]] <- score_density_plot(no_Polyphen2_HVAR[i], 'without_Polyphen2_HVAR')
}

legend <- get_legend(
  plots[[21]] + theme(legend.position = 'right',
                      legend.box.margin = margin(0, 0, 0, 1))
)

plots[[21]] <-  plots[[21]] + theme(legend.position = 'none', 
                                    plot.subtitle = element_text(size = 10, color = "gray30"),
                                    plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'))

plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], 
          plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]], 
          ncol=7, legend, align = 'vh')
ggsave(filename='plots/03_Anno_functional_impact/New_RawScores_density_plots_without_Polyphen2_HVAR.pdf', width = 11.8, height = 7.7)


# ------------------------------------------------------------------------------
## Correlation between raw scores from each pair of methods

## Exclude Polyphen2_HVAR
new_variants_predictions_wP <- new_variants_predictions[,-4]
algorithms_thresholds_wP <- algorithms_thresholds[-3]
  
raw_scores <- as.data.frame(apply(new_variants_predictions_wP[,paste0(names(algorithms_thresholds_wP), '_score')], 2, as.numeric))
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
colnames(corr) <- rownames(corr) <- tool_names[-3]
## Half matrix
corr[lower.tri(corr)] <- NA
## Take absolute corr
corr <- abs(corr)
half_corr_data <- melt(corr, na.rm = TRUE)
half_corr_data$value <- signif(as.numeric(half_corr_data$value), digits = 3)

## Mean corr coeff
unique_half_corr_data <- half_corr_data[which(half_corr_data$value!=1), ]
mean(unique_half_corr_data$value)
# [1] 0.5612857

## Highest corr coeffs
unique_half_corr_data[order(unique_half_corr_data$value, decreasing = TRUE), ][1:4,]
#               Var1             Var2    value
#            MetaSVM           MetaLR     0.93
#         CADD phred         Eigen-PC     0.93
#            MetaSVM            REVEL     0.88
#           ClinPred             ADME     0.85

## Percentage of high coeffs (|r|>0.5)
length(which(unique_half_corr_data$value>0.5))/dim(unique_half_corr_data)[1]*100
# [1] 63.80952

## Percentage of medium coeffs (0.3=<|r|=<0.5)
length(which(unique_half_corr_data$value>=0.3 & unique_half_corr_data$value<=0.5))/dim(unique_half_corr_data)[1]*100
# [1] 27.14286

## Percentage of low coeffs (|r|<0.3)
length(which(unique_half_corr_data$value<0.3))/dim(unique_half_corr_data)[1]*100
# [1] 9.047619


# ------------------------------------------------------------------------------
## Agreement proportion between predictions from each pair of methods
predictions <- as.data.frame(new_variants_predictions_wP[,paste0(names(algorithms_thresholds_wP), '_pred')])
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
colnames(agreement_prop) <- rownames(agreement_prop) <- tool_names[-3]
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
         number.cex = 0.4,
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
         number.cex = 0.4,
         cl.cex = 0.4,
         col.lim = c(0,1),
         col=colorRampPalette(c("white","white", 'darkseagreen2', "darkgreen"))(100)
         
)

dev.off()

## Mean prop
unique_half_agreement_prop <- half_agreement_prop[which(half_agreement_prop$value!=1), ]
mean(unique_half_agreement_prop$value)
# [1] 0.6356

## Highest corr coeffs
unique_half_agreement_prop[order(unique_half_agreement_prop$value, decreasing = TRUE), ][1:4,]
#             Var1             Var2   value
#           MetaSVM          MetaLR   0.930
#           MetaSVM           REVEL   0.926
#            MetaLR           REVEL   0.905
#            FATHMM           REVEL   0.899

## Percentage of high prop (>0.5)
length(which(unique_half_agreement_prop$value>0.5))/dim(unique_half_agreement_prop)[1]*100
# [1] 81.42857


# ------------------------------------------------------------------------------
## Scatterplot of raw scores and categorical predictions of 2 methods

scatterplot_compare_2methods <- function(algorithm1, algorithm2){
  
  algorithm1_score <- paste0(algorithm1, '_score')
  algorithm2_score <- paste0(algorithm2, '_score')
  algorithm1_pred <- paste0(algorithm1, '_pred')
  algorithm2_pred <- paste0(algorithm2, '_pred')
  
  score_type1 <- ' raw score'
  if (algorithm2=='CADD_phred'){score_type2 <- ' phred score'}
  else {score_type2 <- ' raw score'}
  
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
                            else if(x[1]=='D' & x[2]=='N'){paste0('D in ', tool_names[colnames(data)[1]], '; ', '\n',
                                                                  'N in ', tool_names[colnames(data)[2]])}
                            else if(x[1]=='N' & x[2]=='D'){paste0('D in ', tool_names[colnames(data)[2]], '; ', '\n',
                                                                  'N in ', tool_names[colnames(data)[1]])}
                            else {'N in both'}
                      })
  data$categories <- factor(data$categories, levels=c('D in both', 
                                                      paste0('D in ', tool_names[colnames(data)[1]], '; ', '\n',
                                                                          'N in ', tool_names[colnames(data)[2]]),
                                                      paste0('D in ', tool_names[colnames(data)[2]], '; ', '\n',
                                                             'N in ', tool_names[colnames(data)[1]]),
                                                      'N in both'))
  colors <- list()
  colors[['D in both']]='indianred'
  colors[[paste0('D in ', tool_names[colnames(data)[1]], '; ', '\n',
                 'N in ', tool_names[colnames(data)[2]])]]='thistle2'
  colors[[paste0('D in ', tool_names[colnames(data)[2]], '; ', '\n',
                 'N in ', tool_names[colnames(data)[1]])]]='lightpink1'
  colors[['N in both']]='lightblue3'
  
  p <- ggplot(data, aes(x=eval(parse_expr(algorithm1)), y=eval(parse_expr(algorithm2)), color = categories)) +
    geom_point(size=2) +
    stat_smooth(geom = "line", alpha = 0.9, size = 1, method = lm, color = "orangered3", fullrange = FALSE) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      subtitle = paste0("Corr: ", correlation, '; Agreement: ', 100*agreement, '%'),
      x = paste0(tool_names[algorithm1], score_type1),
      y = paste0(tool_names[algorithm2], score_type2),
      color='Predicted effect'
    ) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.title = element_text(size = (11), face='bold'),
      axis.text = element_text(size = (10)),
      plot.subtitle = element_text(size = 9.6, color = "gray30"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size =10, face='bold'),
      legend.key.height = unit(0.95,"cm"), 
      legend.key.width = unit(0.25,"cm"),
      legend.key.size = unit(0.25,"cm"))
  
  return(p)

}

## Compare methods of interest
## Top; similar algorithms
p1 <- scatterplot_compare_2methods('Polyphen2_HDIV', 'Polyphen2_HVAR')
p2 <- scatterplot_compare_2methods('MetaSVM', 'MetaLR')
#p3 <- scatterplot_compare_2methods('FATHMM', 'fathmm.MKL')
## High corr, high agreement
p4 <- scatterplot_compare_2methods('Eigen.PC', 'CADD_phred')
## Top on agreement
p5 <- scatterplot_compare_2methods('MetaSVM', 'REVEL')
## Low corr, high agreement
# p6 <- scatterplot_compare_2methods('FATHMM', 'PrimateAI')
# p7 <- scatterplot_compare_2methods('FATHMM', 'VEST4') 
# p8 <- scatterplot_compare_2methods('Eigen.PC', 'M.CAP') 
## ADME
p9 <- scatterplot_compare_2methods('ADME', 'CADD_phred')
p10 <- scatterplot_compare_2methods('ADME', 'VEST4')
p11 <- scatterplot_compare_2methods('ADME', 'LRT')
p12 <- scatterplot_compare_2methods('ADME', 'MutationAssessor')
p13 <- scatterplot_compare_2methods('ADME', 'PROVEAN')
p14 <-scatterplot_compare_2methods('ADME', 'ClinPred')
## AlphaMissense
p15 <- scatterplot_compare_2methods('AlphaMissense', 'MetaSVM')
p16 <- scatterplot_compare_2methods('AlphaMissense', 'MetaLR')
p17 <- scatterplot_compare_2methods('AlphaMissense', 'REVEL')
p18 <- scatterplot_compare_2methods('AlphaMissense', 'VEST4')

plot_grid(p1, p2, p4, p5, ncol=2, align='vh', 
          labels = LETTERS[1:4])
ggsave('plots/03_Anno_functional_impact/Corr_agreement_tools.pdf', height = 6, width = 9.6)




####################  3.1.5.2 Evaluate predictions of different algorithms  ####################

## Plot the number of predicted D, N and missing variants per algorithm

vars_per_method <- melt(sapply(colnames(new_variants_predictions_wP)[2:22], function(x){table(new_variants_predictions_wP[, x])}))
colnames(vars_per_method) <- c('Prediction', 'Method', 'Numbers')
vars_per_method$Method <- tool_names[-3][gsub('_pred', '', vars_per_method$Method)]
# Order methods by number of D variants
numD<- sapply(colnames(new_variants_predictions_wP)[2:22], function(x){table(new_variants_predictions_wP[, x])['D']})
numD <- numD[order(numD, decreasing = TRUE)]
names(numD) <- tool_names[-3][gsub('_pred.D', '', names(numD))]
vars_per_method$Method <- factor(vars_per_method$Method, levels=names(numD))
## Order to have D variants first in each bar
cat_order <- c('.', 'N', 'D')
vars_per_method$Prediction <- factor(vars_per_method$Prediction, levels=cat_order)

ggplot(vars_per_method, aes(x=Method, y=Numbers, fill=Prediction)) + 
  geom_bar(position="stack", stat="identity", colour = 'black', width=.65) + 
  theme_classic() +
  labs(
    y = 'Number of predicted missense variants',
    x= '',
    subtitle = paste0(dim(new_variants_predictions)[1], ' total missense variants across all UGT genes')
  ) +
  scale_fill_manual(values = c("grey80", "skyblue2", "tomato"), labels = c("Missing", "Neutral", "Deleterious")) + 
  #scale_y_discrete(limits= c(0, 1500, 3000, 4500, dim(new_variants_predictions)[1])) +
  coord_cartesian(ylim = c(0, dim(new_variants_predictions)[1]),clip = 'off') +
  geom_text(data=subset(vars_per_method, Prediction=='.'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+810, 
                                                               fill=NULL), hjust = 0.5, size = 2) +
  geom_text(data=subset(vars_per_method, Prediction=='N'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+540, 
                                                               fill=NULL), hjust = 0.5, size = 2) +
  geom_text(data=subset(vars_per_method, Prediction=='D'), aes(label=Numbers, y=dim(new_variants_predictions)[1]+270, 
                                                               fill=NULL), hjust = 0.5, size = 2) +
  theme(plot.subtitle = element_text(size = (10), vjust = 9.8, hjust=0, color="gray50", face='bold'), 
        legend.direction = "vertical",
        legend.position = c(1.10, 1.07),
        legend.key.size = unit(0.17, units = 'cm'),
        plot.margin = unit(c(3, 3, 0.5, 0.5), "cm"),
        axis.title = element_text(size = (11), face='bold'),
        axis.text = element_text(size = (10)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size =9, face='bold'))

ggsave(filename='plots/03_Anno_functional_impact/D_N_M_vars_per_method.pdf', width = 6.5, height = 5)

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
  vars_per_method$Method <- tool_names[gsub('_pred', '', vars_per_method$Method)]
  
  # Order methods by number of D variants first, and N variants second
  DNM_numbers <- as.data.frame(t(DNM_numbers))
  DNM_numbers <- DNM_numbers[order(DNM_numbers$D, DNM_numbers$N,  decreasing = TRUE),]
  numD <- rownames(DNM_numbers)
  numD <- tool_names[gsub('_pred', '', numD)]
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
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
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
           'FATHMM-MKL'='lightcoral', 
           'LRT'='lightsteelblue3',
           'M-CAP'='yellow4', 
           'MetaLR'='steelblue2', 
           'MetaSVM'='dodgerblue3', 
           'MutationAssessor'='goldenrod', 
           'MutPred'='magenta2',
           'MVP'='blue2', 
           'PolyPhen-2 HDIV'='darkseagreen3', 
           'PolyPhen-2 HVAR'='mediumseagreen', 
           'PrimateAI'='darkred',
           'PROVEAN'='darkorchid3', 
           'REVEL'='cadetblue3', 
           'SIFT'='peachpuff3', 
           'VEST'='aquamarine4',
           'UGT-optimized'='palevioletred2')

for (variant in new_variants_predictions_wP$Variant_ID){
  allele_freq <- vector()
  ## Allele freq of variant in each gene dataset
  for (gene in UGT_genes){
    UGT_missense_vars <- eval(parse_expr(paste0(gene, '_missense_vars')))
    if(variant %in% UGT_missense_vars$Variant_ID){
      allele_freq <- append(allele_freq, UGT_missense_vars[which(UGT_missense_vars$Variant_ID==variant), 'Allele_Frequency'])
    }
  }
  if (length(unique(allele_freq))==1) {allele_freq <- unique(allele_freq)}
  new_variants_predictions_wP[which(new_variants_predictions_wP$Variant_ID==variant), 'Allele_Frequency'] <- allele_freq
}

## Allele frequencies of D variants per method
data <- vector()
for(method in colnames(new_variants_predictions_wP)[2:22]){
  allele_freq_method <-  new_variants_predictions_wP[which(new_variants_predictions_wP[,method]=='D'), c('Variant_ID', 'Allele_Frequency')]
  method <- rep(method, length(allele_freq_method$Allele_Frequency))
  method_data <- cbind(allele_freq_method$Variant_ID, allele_freq_method$Allele_Frequency, method)
  data <- rbind(data, method_data)
}
data <- as.data.frame(data)
colnames(data) <- c('Variant_ID', 'Allele_Frequency', 'Method')
data$Allele_Frequency<- as.numeric(data$Allele_Frequency)
data$Method <- tool_names[-3][gsub('_pred', '', data$Method)]
data$Method <- factor(data$Method, levels=names(numD))

## Identify variants with allele freq >0.5
data[which(data$Allele_Frequency>0.5),]
#       Variant_ID  Allele_Frequency          Method
#   4-69795626-C-T         0.7568806 PolyPhen-2 HDIV
#   4-69795626-C-T         0.7568806            CADD
#   4-69795626-C-T         0.7568806            DANN
#  4-115589302-A-G         0.9953603             LRT

## % of D variants with MAF>0.01 per method
sapply(levels(data$Method), function(x){dim(subset(data, Method==x & Allele_Frequency>0.01))[1] / dim(subset(data, Method==x))[1]*100})
#    MutPred             CADD          PROVEAN             DANN MutationAssessor             SIFT 
#  0.0000000        0.4688362        0.4182934        0.5136986        0.4184100        0.5382775 
#  PolyPhen-2 HDIV         ClinPred       FATHMM-MKL             ADME              LRT 
#        0.3984064        0.0000000        0.4063539        0.4219409        0.8574491 
#   Eigen-PC             VEST            M-CAP    AlphaMissense           MetaLR          MetaSVM 
#  0.3351955        0.2347418        0.0000000        0.4629630        0.0000000        0.1041667 
#       MVP            REVEL           FATHMM        PrimateAI 
# 0.0000000        0.1455604        2.1052632        0.0000000 


shapes <- c('4-69795626-C-T'=17,
            '4-115589302-A-G'=15)
data$label <- apply(data, 1, function(x){if (as.numeric(x['Allele_Frequency'])>0.5){x['Variant_ID']} else {NA}})
## Number of D vars per method
num_per_method <- as.data.frame(table(data$Method))
colnames(num_per_method) <- c('Method', 'n')

## Number of unique D variants that are common (across all methods)
length(unique(subset(data, Allele_Frequency>0.01)$Variant_ID))
# [1] 32

## Percentage of all D variants (across all methods) the common ones represent
32/length(unique(data$Variant_ID)) *100
# [1] 0.5873715

## Save data
write.table(data, file ='processed-data/03_Anno_functional_impact/GMAF_allDvars_perMethod.csv', row.names = FALSE, col.names = TRUE, sep = '\t')


ggplot(data = data, mapping = aes(x = Method, y = Allele_Frequency, color = Method, shape=label)) +
  geom_jitter(data=subset(data, is.na(label)), shape=16, width = 0.1, height = 0, alpha = 0.7, size = 1.5) +
  geom_point(data=subset(data, !is.na(label)), aes(shape=label), alpha = 0.7, size = 1.3, stroke = 1) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  scale_color_manual(values = colors) +
  guides(color = 'none') + 
  #geom_text(data = num_per_method, aes(x=Method, label=n,  y=-0.05, shape=NULL, color=NULL), size=1.6) +
  labs(x='', y='MAF of missense variants predicted as deleterious', shape='Variant ID (MAF>0.5)', 
       subtitle = paste0(dim(new_variants_predictions)[1], ' total missense variants across all UGT genes')) +
  theme(title = element_text(size = (6), face='bold'),
        plot.subtitle = element_text(size = (7), color="gray50", face='bold'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = (6.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=5.5),
        legend.position = c(0.8,0.8),
        legend.key = element_blank(),
        legend.background = element_rect(fill=NA, color='black'),
        legend.key.size = unit(0.6, 'lines'))

ggsave(filename='plots/03_Anno_functional_impact/GMAF_allDvars_perMethod.pdf', width = 4, height = 3.3)

# ----------------------------------------------------------------------------------------------------

## Plot GMAF of variants predicted as D by each method in each gene 
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
    labs(x='', y='MAF of missense variants predicted as deleterious', shape='Variant ID (MAF>0.5)', 
         title=gene, subtitle = paste0(length(unique(gene_D_vars$Variant_ID)), ' missense variants predicted as deleterious by at least one method')) +
    theme(title = element_text(size = (9), face='bold'),
          plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = (8.5), face='bold'),
          axis.text = element_text(size = (8)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
          legend.title = element_text(size=8), 
          legend.text = element_text(size=7.5),
          legend.position = c(0.9,0.95),
          legend.key = element_blank(),
          legend.background = element_rect(fill=NA, color='black'),
          legend.key.size = unit(0, 'lines'),
          legend.justification = c(0.9,0.9))
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

## Carrier frequency is q**2 + 2pq with q the frequency of the alternate allele
data$Carrier_Frequency <- (2*data$Allele_Frequency*(1-data$Allele_Frequency)) + (data$Allele_Frequency)**2
num_per_method <- as.data.frame(table(data$Method))
colnames(num_per_method) <- c('Method', 'n')

## Identify variants with allele freq >0.5
data[which(data$Carrier_Frequency>0.5),]
# Variant_ID      Allele_Frequency         Method           label Carrier_Frequency
# 2-234602202-A-C        0.3443850           SIFT            <NA>         0.5701690
#  4-69795626-C-T        0.7568806 Polyphen2 HDIV  4-69795626-C-T         0.9408930
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
  data <- read_delim(paste0('~/Documents/KI_projects/UGT_genetic_profiling_KI/raw-data/ClinVar_data/clinvar_variants_', gene, '.txt'), delim='\t', show_col_types = FALSE)
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
  print(paste0(dim(data)[1], ' benchmark variants in ', gene, ':', table(data$effect)))
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
save(benchmark_data, file ='processed-data/03_Anno_functional_impact/ANNOVAR/input_data/benchmark_data.Rdata')
write.table(benchmark_data, file ='processed-data/03_Anno_functional_impact/ANNOVAR/input_data/benchmark_data.txt', row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(benchmark_data, file ='processed-data/03_Anno_functional_impact/ANNOVAR/input_data/benchmark_data.csv', row.names = FALSE, col.names = FALSE, sep = '\t')

## --> ANNOVAR was run in 3.1_run_ANNOVAR.sh

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

colors = c('ADME'='mediumpurple2', 
           'AlphaMissense'='red3',
           'CADD'='darkorange3', 
           'ClinPred'='turquoise3',
           'DANN'='yellow3', 
           'Eigen-PC'='olivedrab3',
           'FATHMM'='lightsalmon2', 
           'FATHMM-MKL'='lightcoral', 
           'LRT'='lightsteelblue3',
           'M-CAP'='yellow4', 
           'MetaLR'='steelblue2', 
           'MetaSVM'='dodgerblue3', 
           'MutationAssessor'='goldenrod', 
           'MutPred'='magenta2',
           'MVP'='blue2', 
           'PolyPhen-2 HDIV'='darkseagreen3', 
           'PolyPhen-2 HVAR'='mediumseagreen', 
           'PrimateAI'='darkred',
           'PROVEAN'='darkorchid3', 
           'REVEL'='cadetblue3', 
           'SIFT'='peachpuff3', 
           'VEST'='aquamarine4',
           'UGT-optimized'='palevioletred2')

r <- list()
for (algorithm in paste0(names(algorithms_thresholds), '_score')){
  r[[algorithm]] <- roc(response=benchmark_scores_preds$effect, 
                        predictor=as.numeric(benchmark_scores_preds[,algorithm]), 
                        levels=c('N', 'D'), na.rm=TRUE)
}
## Add AUC per method
data <- as.data.frame(cbind('AUC'=paste0('AUC = ', signif(as.numeric(lapply(r, function(x){x$auc})), digits=3))))
data$auc_number <-  signif(as.numeric(lapply(r, function(x){x$auc})), digits=3)
## Add number of D and N variants used to evaluate each method
data$num_vars <- sapply(colnames(benchmark_pred)[2:23], function(x){
                       DN_num <- table(benchmark_scores_preds[which(benchmark_scores_preds[,x]!='.'), 'effect'])
                       paste0('n = ', DN_num['D'], ' D; ', DN_num['N'], ' N')
                       })
## Locate coordinate corresponding to the used threshold for each method
data$sensitivity <- sensitivities
data$specificity <- specificities

data$name <- gsub(' phred', '', gsub('_', ' ', gsub('\\.', '-', gsub('_score', '', names(r)))))
data$method_name <-  gsub('_score', '', names(r))
names(r) <- tool_names[gsub('_score', '', names(r))]

## Add coordinates for ADME thresholds
data$ADME_sensitivity <- as.numeric(sapply(data$name, function(x){if(x %in% names(ADME_sensitivities)){ADME_sensitivities[x]} else {'NA'}}))
data$ADME_specificity <- as.numeric(sapply(data$name, function(x){if(x %in% names(ADME_specificities)){ADME_specificities[x]} else {'NA'}}))

data$name <- names(r)

## Order by AUC 
data <- data[order(data$auc_number, decreasing = TRUE),]
r <- r[data$name]

ggroc(r) + 
  facet_wrap(~factor(name, levels=data$name)) +
  theme_bw() + 
  geom_text(data = data, aes(0, 0.19, label= AUC, hjust = 1), size=3.2, fontface='bold') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="gray95", size=1, color="gray60"),
        strip.text = element_text(face="bold"),
        axis.text = element_text( size = 8)) +
  geom_text(data = data, aes(0, 0.05, label= num_vars, hjust = 1), size=2.5, color='black') +
  scale_color_manual(values=colors) +
  guides(color='none') +
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
  method_name <- tool_names[method]
  
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
    
    ## Add if J hasn't been reached before with another threshold; add conventional and max score; 
    ## Don't include scores that yield J<0
    if ((! J %in% Js | score==conv_threshold | score==max(as.numeric(method_data[,method_score]))) & J>=0 ){
      scores_Js <- rbind(scores_Js, c('score'=score, 'sensitivity'=sensitivity, 'specificity'=specificity, 'J'=J))
    }
    Js <- append(Js, J)
  }
  
  scores_Js <- as.data.frame(scores_Js)
  ## New threshold
  max_J <- scores_Js[which.max(scores_Js$J),'J']
  max_J_score <- scores_Js[which.max(scores_Js$J),'score'] ## optimal threshold
  max_J_sensitivity <- scores_Js[which.max(scores_Js$J),'sensitivity']
  max_J_specificity <- scores_Js[which.max(scores_Js$J),'specificity']
  new_threshold <- paste0(direction, '(', max_J_score, ')')
  max_J_data <- c(scores_Js[which.max(scores_Js$J),], 'new_threshold'=new_threshold)

  ## J at optimal threshold - J at conventional threshold = delta(J)
  delta_J <- max_J -J_conventional
  
  ## df for plot arguments
  args <- list()
  args[['SIFT']] <- c('nudge_x_c' = 0.13, 'nudge_y_c' = 0.01, 'nudge_x_o' =-0.03, 'nudge_y_o' = 0.15, 'curv'=0.35, 'hjust'=1, 'vjust'=0)
  args[['Polyphen2_HDIV']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.2, 'nudge_y_o' = 0.15, 'curv'=0.25, 'hjust'=0, 'vjust'=-1)
  args[['Polyphen2_HVAR']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['MutationAssessor']] <- c('nudge_x_c' = -0.2, 'nudge_y_c' = 0.01, 'nudge_x_o' = 0.45, 'nudge_y_o' = 0.15, 'curv'=0.35, 'hjust'=1, 'vjust'=1)
  args[['FATHMM']] <- c('nudge_x_c' = -0.85, 'nudge_y_c' = 0.1, 'nudge_x_o' = -0.03, 'nudge_y_o' = 0.12, 'curv'=0.35, 'hjust'=1, 'vjust'=0)
  args[['fathmm.MKL']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.2, 'nudge_y_o' = 0.15, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['PROVEAN']] <- c('nudge_x_c' = 0.9, 'nudge_y_c' = 0.025, 'nudge_x_o' = -0.9, 'nudge_y_o' = 0.0, 'curv'=0.15, 'hjust'=0.1, 'vjust'=0.1)
  args[['PrimateAI']] <- c('nudge_x_c' = 0.0, 'nudge_y_c' = 0.2, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.02, 'curv'=0.35, 'hjust'=1, 'vjust'=0)
  args[['MutPred']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = 0.055, 'nudge_x_o' = 0.1, 'nudge_y_o' = 0.055, 'curv'=0.35, 'hjust'=1, 'vjust'=0)
  args[['MetaSVM']] <- c('nudge_x_c' = 0.07, 'nudge_y_c' = 0.07, 'nudge_x_o' = 0.1, 'nudge_y_o' = 0.2, 'curv'=0.35, 'hjust'=-0.63, 'vjust'=0.7)
  args[['MetaLR']] <- c('nudge_x_c' = -0.15, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['M.CAP']] <- c('nudge_x_c' = 0.05, 'nudge_y_c' = -0.17, 'nudge_x_o' = 0.14, 'nudge_y_o' = 0, 'curv'=-0.05, 'hjust'=1, 'vjust'=0.3)
  args[['ClinPred']] <- c('nudge_x_c' = 0.06, 'nudge_y_c' = -0.05, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.2, 'curv'=-0.3, 'hjust'=0, 'vjust'=0)
  args[['CADD_phred']] <- c('nudge_x_c' = -1, 'nudge_y_c' = 0.5, 'nudge_x_o' = 2, 'nudge_y_o' = 0, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['Eigen.PC']] <- c('nudge_x_c' = -0.4, 'nudge_y_c' = 0.1, 'nudge_x_o' = 0.8, 'nudge_y_o' = 0.15, 'curv'=-0.35, 'hjust'=0.8, 'vjust'=0)
  args[['MVP']] <- c('nudge_x_c' = 0.1, 'nudge_y_c' = 0, 'nudge_x_o' = -0.18, 'nudge_y_o' = 0, 'curv'=0.35, 'hjust'=0, 'vjust'=0.1)
  args[['DANN']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.02, 'nudge_x_o' = -0.05, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['REVEL']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['LRT']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)
  args[['VEST4']] <- c('nudge_x_c' = -0.15, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=0.35, 'hjust'=1, 'vjust'=0)
  args[['ADME']] <- c('nudge_x_c' = -0.1, 'nudge_y_c' = -0.1, 'nudge_x_o' = 0.15, 'nudge_y_o' = 0.1, 'curv'=-0.35, 'hjust'=1, 'vjust'=0)
  args[['AlphaMissense']] <- c('nudge_x_c' = 0.06, 'nudge_y_c' = 0.05, 'nudge_x_o' = -0.1, 'nudge_y_o' = 0.06, 'curv'=-0.35, 'hjust'=0, 'vjust'=0)

  if(method=='PROVEAN' | method=='PrimateAI' | method=='MetaSVM' | method=='MetaLR' | method=='REVEL' | method=='AlphaMissense'){
    scores_to_show <- c(1,3,5)
  }
  else if(method=='Eigen.PC'){
    scores_to_show <- c(1,4,5)
  }
  else{
    scores_to_show <- c(1,5)
  }
  
  ## Plot scores vs J 
  p <- ggplot(scores_Js, aes(x=score, y=J))+
    geom_line(color=colors[method_name], size=1.2) + 
    coord_cartesian(ylim=c(0, NA), expand = FALSE) +
    
    ## Line for score that yields the max J
    geom_segment(aes(x = max_J_score, y = 0, xend = max_J_score, yend = max_J), linetype=1, linewidth=0.5, color='grey60') +
    ## Max J with optimized threshold
    geom_point(x=max_J_score, y=max_J, color='red', size=2) + 
    geom_text_repel(data = subset(scores_Js, J==max_J & score!=conv_threshold), label=signif(max_J, digits=3), 
                    size=3.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=1, 
                    box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = args[[method]]['nudge_x_o'], 
                    nudge_y = args[[method]]['nudge_y_o'], fontface='bold') + 
    
    ## J with conventional threshold
    geom_point(x=conv_threshold, y=J_conventional, color=colors[method_name], size=2) +
    geom_text_repel(data = subset(scores_Js, score==conv_threshold), label=signif(J, digits=3), 
                    size=3.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=0, vjust=1, 
                    box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = args[[method]]['nudge_x_c'], 
                    nudge_y = args[[method]]['nudge_y_c'], fontface='bold') + 
    ## Line for conventional score 
    geom_segment(aes(x = conv_threshold, y = 0, xend = conv_threshold, yend = J_conventional), linetype=3, linewidth=0.8, color='grey60') +
    
    theme_classic() +
    labs(x='Score', y='Youden index (J)', title=method_name) +
    scale_x_continuous(breaks = sort(c(min(as.numeric(scores_Js$score)), max(as.numeric(scores_Js$score)),
                                       signif(c(seq(from=min(as.numeric(scores_Js$score)),
                                                  to=max(as.numeric(scores_Js$score)),
                                                  length.out=5)[-scores_to_show], conv_threshold), digits=2)))) +
  
    ## Label for delta J
    geom_label(x=(max(as.numeric(method_data[,method_score]))+min(as.numeric(method_data[,method_score])))/2, y=0.14, 
               label=paste0('ΔJ = ', signif(delta_J, digits=2)), 
               size=3.8, color='grey30', label.size = NA, fontface='bold', label.padding = unit(0.1, "lines")) +
    
    theme(title = element_text(size = (11), face='bold'),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = (11), face='bold'),
          axis.text = element_text(size = (11)), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=10))
  
  
  if(method!='MutationAssessor'){
    ## Text for new score
    p <- p + geom_text_repel(data = subset(scores_Js, score==max_J_score & score!=conv_threshold), label=signif(max_J_score, digits=3), 
                        aes(x=max_J_score, y=0.05), size=3.1, color='grey40', min.segment.length = unit(0, 'lines'), 
                        hjust=args[[method]]['hjust'], vjust=args[[method]]['vjust'], 
                        box.padding = unit(0.5, "lines"), lineheight=unit(2, 'lines'), linewidth=unit(2, 'cm'), segment.curvature = args[[method]]['curv']) 
  }
  
  
  if(method=='LRT'){
    p <- ggplot(scores_Js, aes(x=score, y=J))+
      geom_line(color=colors[method_name], size=1.2) + 
      geom_point(x=max_J_score, y=max_J, color='red', size=2) + 
      geom_text_repel(data = subset(scores_Js, J==max_J), label=signif(max_J, digits=3), 
                      size=3.7, color='grey40', min.segment.length = unit(0, 'lines'), hjust=1, 
                      box.padding = 0.5, lineheight=unit(1, 'lines'), nudge_x = 0.1, nudge_y = 0, fontface='bold') + 
      geom_segment(aes(x = max_J_score, y = 0, xend = max_J_score, yend = max_J), linetype=1, linewidth=0.5, color='grey60') +
      geom_text_repel(data = subset(scores_Js, score==max_J_score & score!=conv_threshold), label=signif(max_J_score, digits=3), 
                      aes(x=max_J_score, y=0.05), size=3.1, color='grey40', min.segment.length = unit(0, 'lines'), 
                      hjust=args[[method]]['hjust'], vjust=args[[method]]['vjust'], 
                      box.padding = unit(0.5, "lines"), lineheight=unit(1, 'lines'), segment.curvature = args[[method]]['curv']) +
      theme_classic() +
      coord_cartesian(ylim=c(0, 1), 
                      xlim=c(-0.01, max(as.numeric(scores_Js$score))),
                      expand = FALSE) +
      labs(x='Score', y='Youden index (J)', title=method_name) +
      ## Label for delta J
      geom_label(x=(max(as.numeric(method_data[,method_score]))+min(as.numeric(method_data[,method_score])))/2, y=0.14, 
                 label=paste0('ΔJ = ', signif(delta_J, digits=2)), 
                 size=3.8, color='grey30', label.size = NA, fontface='bold', label.padding = unit(0.1, "lines")) +
      theme(title = element_text(size = (11), face='bold'),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = (11), face='bold'),
            axis.text = element_text(size = (11)), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=10))
  }
  
  return(list(p, max_J_data))
}

i=1
plots <- list()
new_thresholds <- vector()
for (method in data$method_name){
  output <- Youden_indices(method)
  new_thresholds <- rbind(new_thresholds, output[[2]])
  plots[[i]] <- output[[1]]
  i=i+1
}
new_thresholds <- as.data.frame(new_thresholds)
rownames(new_thresholds) <- data$method_name

## Don't include PolyPhen2 HVAR
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],
          plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]],
          plots[[13]], plots[[14]], plots[[15]], plots[[16]], plots[[18]], 
          plots[[19]], plots[[20]], plots[[21]], plots[[22]], ncol=7)
ggsave(filename='plots/03_Anno_functional_impact/Youden_Index_plots.png', width = 24, height = 9)
ggsave(filename='plots/03_Anno_functional_impact/Youden_Index_plots.pdf', width = 24, height = 9)
## Smaller version
ggsave(filename='plots/03_Anno_functional_impact/Youden_Index_plots_smaller.pdf', width = 15.3, height = 7.5)


## Add coordinate for new thresholds in ROC curves
data$new_threshold_sensitivity <- unlist(new_thresholds$sensitivity)
data$new_threshold_specificity <- unlist(new_thresholds$specificity)

ggroc(r[-17]) + 
  facet_wrap(~factor(name, levels=names(r[-17])), ncol=7) +
  theme_bw() + theme(legend.position='none') +
  geom_text(data = subset(data, name!='PolyPhen-2 HVAR'), aes(0, 0.19, label= AUC, hjust = 1), size=3.2, fontface='bold') +
  theme(strip.background = element_rect(fill="gray95", size=1, color="gray60"),
        strip.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text( size = 6)) +
  geom_text(data = subset(data, name!='PolyPhen-2 HVAR'), aes(0, 0.05, label= num_vars, hjust = 1), size=2.5, color='black') +
  scale_color_manual(values=colors) +
  ## Point corresponding to conventional threshold
  geom_point(data=subset(data, name!='PolyPhen-2 HVAR'), aes(x=specificity, y=sensitivity)) +
  ## Point for UGT optimized thresholds
  geom_point(data=subset(data, name!='PolyPhen-2 HVAR'), aes(x=new_threshold_specificity, y=new_threshold_sensitivity), shape=5, color='red', size=1.3, stroke = 1) 

ggsave(filename='plots/03_Anno_functional_impact/AUC_ROC_new_thresholds_methods.pdf', width = 9, height = 4.5)


## Create table with results
results <- data.frame(cbind('algorithm'= as.vector(data$name), 'conv_threshold' = as.vector(algorithms_thresholds[data$method_name]), 
                               'conv_threshold_sensitivity'= as.vector(data$sensitivity), 'conv_threshold_specificity'= as.vector(data$specificity),
                               'conv_threshold_J'= as.vector(data$sensitivity + data$specificity -1), 
                               'new_threshold'= as.vector(new_thresholds$new_threshold), 'new_threshold_sensitivity'= as.vector(new_thresholds$sensitivity),
                               'new_threshold_specificity'= as.vector(new_thresholds$specificity), 'new_threshold_J'= as.vector(new_thresholds$J)) ) 
results <- apply(results, 2, as.character)
write.table(results, file = 'processed-data/03_Anno_functional_impact/UGT_optimization_results.csv', row.names = FALSE, col.names = TRUE, sep = '\t')


# ---------------------------------------------------------------------------------------------------------

## Predictions of benchmark variants with these cutoffs

## Don't include PolyPhen2 HVAR
new_algorithms_thresholds <- new_thresholds[which(rownames(new_thresholds)!='Polyphen2_HVAR'),'new_threshold']
names(new_algorithms_thresholds) <- data$method_name[-17]

new_benchmark_pred <- data.frame(matrix(nrow=dim(benchmark_scores)[1], ncol=length(names(new_algorithms_thresholds))+1))
colnames(new_benchmark_pred) <- c('Variant_ID', paste0(names(new_algorithms_thresholds), '_pred'))
new_benchmark_pred$Variant_ID <- benchmark_scores$Variant_ID

for(algorithm in names(new_algorithms_thresholds)){
  new_benchmark_pred[paste0(algorithm, '_pred')] <- apply(benchmark_scores, 1, 
                                                      function(x){if (x[paste0(algorithm, '_score')]=='.'){NA}
                                                        else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', new_algorithms_thresholds[[algorithm]]))) ){1}
                                                        else{0} })
}

## Global UGT-optimized score as the mean of all predictions for each variant
UGT_optimized_score <- apply(new_benchmark_pred[,-1], 1, function(x){mean(x[which(x!='.')])})
UGT_optimized_score[which(is.nan(UGT_optimized_score))] <- NA

r <- list()
r[['UGT-optimized']] <- roc(response=as.factor(benchmark_scores$effect), 
                        predictor=as.numeric(UGT_optimized_score), 
                        levels=c('N', 'D'), na.rm=TRUE)
## Add AUC 
AUC=paste0('AUC = ', r$`UGT-optimized`$auc)
## Add number of D and N variants used
DN_num <- table(benchmark_scores[which(!is.na(UGT_optimized_score)), 'effect'])
num_vars <- paste0('n = ', DN_num['D'], ' D; ', DN_num['N'], ' N')
## Compute J and find best threshold
Js <- as.data.frame(cbind('threshold'=signif(r$`UGT-optimized`$thresholds, digits=3), 
                          'specificity'=r$`UGT-optimized`$specificities, 
                          'sensitivity'=r$`UGT-optimized`$sensitivities, 
                          'J'=(r$`UGT-optimized`$specificities + r$`UGT-optimized`$sensitivities -1)))
ggroc(r) + 
  facet_wrap(~ name) +
  theme_bw() + theme(legend.position = "none") + 
  geom_text(aes(0, 0.19), label= AUC, hjust = 1, size=3.2, fontface='bold') +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_rect(fill="gray95", size=1, color="gray60"),
         strip.text = element_text(face="bold"),
         axis.text = element_text( size = 6)) +
  geom_text(x=0, y=0.05, label= num_vars, hjust = 1, size=2.5, color='black') +
  scale_color_manual(values=colors) +
  ## Point corresponding to best threshold
  annotate('point', x=Js[which.max(Js$J),'specificity'], y=Js[which.max(Js$J),'sensitivity'], color=colors[['UGT-optimized']]) +
  geom_text_repel(data=Js[which.max(Js$J),], aes(x=specificity, y=sensitivity), 
                  label=paste0('Threshold of ', Js[which.max(Js$J),'threshold']),
                  size=2.4, color='grey20', min.segment.length = unit(0, 'lines'),
                  box.padding = unit(0.8, "lines"), lineheight=unit(6, 'lines'), segment.curvature = 0)

ggsave(filename='plots/03_Anno_functional_impact/AUC_ROC_UGT-optimizedMethod.pdf', width = 2.1, height = 2)

# ---------------------------------------------------------------------------------------------------------

## Use framework for all missense variants

## New categories with optimized thresholds
new_categorical_predictions <- data.frame(matrix(nrow=dim(new_variants_predictions)[1], ncol=length(names(new_algorithms_thresholds))+1))
colnames(new_categorical_predictions) <- c('Variant_ID', paste0(names(new_algorithms_thresholds), '_pred'))
new_categorical_predictions$Variant_ID <- new_variants_predictions$Variant_ID

for(algorithm in names(new_algorithms_thresholds)){
  new_categorical_predictions[paste0(algorithm, '_pred')] <- apply(new_variants_predictions, 1, 
                                                               function(x){if (x[paste0(algorithm, '_score')]=='.'){NA}
                                                                 else if (eval(parse_expr(paste0('as.numeric(x[paste0(algorithm, \'_score\')])', new_algorithms_thresholds[[algorithm]]))) ){1}
                                                                 else{0} })
}

new_variants_predictions$UGT_optimized_score <- apply(new_categorical_predictions[,-1], 1, function(x){mean(x[which(!is.na(x))])}) 
## Predictions
new_variants_predictions$UGT_optimized_pred <- sapply(new_variants_predictions$UGT_optimized_score, function(x){if(x>0.821){'D'}
                                                                                                                else{'N'}})
save(new_variants_predictions, file='processed-data/03_Anno_functional_impact/new_variants_predictions.Rdata')

## Number of D and N variants predicted
table(new_variants_predictions$UGT_optimized_pred)
#   D     N 
# 448  5904 

# ---------------------------------------------------------------------------------------------------------

## Plot GMAF of D variants per gene

GMAFs_genes <- vector()
for (gene in UGT_genes){
  ## D variants in the gene
  gene_data <- eval(parse_expr(paste0(gene, '_missense_vars')))
  gene_vars <- gene_data$Variant_ID
  gene_D_vars <- new_variants_predictions[which(new_variants_predictions$Variant_ID %in% gene_vars & 
                                                new_variants_predictions$UGT_optimized_pred=='D'), 'Variant_ID']
  GMAFs <- gene_data[which(gene_data$Variant_ID %in% gene_D_vars), c('Variant_ID', 'Allele_Frequency')]
  
  ## Number of D vars
  num_D <- length(gene_D_vars)
  
  GMAFs_genes <- rbind(GMAFs_genes, cbind('Variant_ID'=GMAFs$Variant_ID, 'GMAFs'=GMAFs$Allele_Frequency, 'gene'=rep(gene, num_D)))
}
  
GMAFs_genes <- as.data.frame(GMAFs_genes)
GMAFs_genes$GMAFs <- as.numeric(GMAFs_genes$GMAFs)
GMAFs_genes$gene <- factor(GMAFs_genes$gene, levels=UGT_genes)

## Number of D variants per gene
num_D_per_gene <- as.data.frame(table(GMAFs_genes$gene))
colnames(num_D_per_gene) <- c('gene', 'number')

## Colors for genes
genes_colors <- c('UGT1A1'='thistle2',
                  'UGT1A3'='plum2',
                  'UGT1A4'='plum3',
                  'UGT1A5'='pink',
                  'UGT1A6'='lightpink2',
                  'UGT1A7'='lightpink3',
                  'UGT1A8'='hotpink3',
                  'UGT1A9'='orchid2',
                  'UGT1A10'='violetred',
                  'UGT1A[1-10]'='hotpink1',
                  'UGT2A1'='peachpuff2',
                  'UGT2A2'='sandybrown',
                  'UGT2A[1-2]'='chocolate2',
                  'UGT2A3'='tan3',
                  'UGT2B4'='lightsalmon1',
                  'UGT2B7'='bisque3',
                  'UGT2B10'='lightsalmon3',
                  'UGT2B11'='coral',
                  'UGT2B15'='sienna3',
                  'UGT2B17'='lightcoral',
                  'UGT2B28'='chocolate4',
                  'UGT3A1'='thistle3',
                  'UGT3A2'='plum4',
                  'UGT8'='lightskyblue3')

## Shapes for variants with GMAF>0.01
GMAFs_genes$label <- apply(GMAFs_genes, 1, function(x){if(as.numeric(x['GMAFs'])>0.01){x['Variant_ID']} else{NA}})
GMAFs_genes[which(GMAFs_genes$GMAFs>0.01),]
#     Variant_ID       GMAFs    gene           label
# 4-70462042-C-T  0.07988111  UGT2A1  4-70462042-C-T
# 4-70462042-C-T  0.07988111  UGT2A2  4-70462042-C-T
shapes <- c('4-70462042-C-T'=25)

## Compact all variants with GMAF<=1e-05 to the single category '<=1e-05'
GMAFs_genes$GMAFs <- sapply(GMAFs_genes$GMAFs, function(x){if(x<=1e-05){1e-05} else{x}})

ggplot(data = GMAFs_genes, mapping = aes(x = gene, y = GMAFs, color = gene)) +
    geom_jitter(data=subset(GMAFs_genes, is.na(label)), shape=16, width = 0.3, height = 0, alpha = 0.8, size = 1.2) +
    geom_point(data=subset(GMAFs_genes, !is.na(label)), aes(shape=label), alpha = 0.8, size = 1.2, stroke = 1) +
    scale_y_continuous(trans='log10', breaks=c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05),
                       labels=c('-1', '-2', '-3', '-4', '≤-5')) +
    scale_color_manual(values = genes_colors) +
    scale_shape_manual(values = shapes) +
    theme_bw() +
    guides(color = 'none') + 
    geom_text(data = num_D_per_gene, aes(x=gene, label=number,  y=(1e-05)/4*3, color=NULL), size=2) +
    labs(x='', y='log10(MAF) of missense variants predicted as deleterious', 
        subtitle = paste0(table(new_variants_predictions$UGT_optimized_pred)['D'], 
                          ' missense variants predicted as deleterious by this UGT-optimized framework'),
        shape='Variant ID (MAF>0.01)') +
    theme(plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
          axis.title = element_text(size = (8.5), face='bold'),
          axis.text = element_text(size = (8)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='italic'),
          legend.title = element_text(size=8, face='bold'), 
          legend.text = element_text(size=7.5, face='bold'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())

ggsave(filename='plots/03_Anno_functional_impact/GMAF_Dvars_byUGTopMethod.png', width = 8, height = 4.5)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_Dvars_byUGTopMethod.pdf', width = 8, height = 4.5)





################################################################################
##          3.2  Annotate functional consequence of all UGT variants 
################################################################################

## Define categories of predicted effect of exonic variant types
deleterious <- c('splice_donor_variant', 'splice_acceptor',
                 'stop_gained', 'frameshift_variant', 'start_lost', '5\' upstream')

neutral <- c('synonymous_variant', 'inframe_deletion', 'inframe_insertion', 
             'stop_retained_variant', 'splice_region_variant', 'stop_lost')

## Annotate functional impact of all exonic (+ 2 promoter variants in UGT1A1) variants in each gene
for (gene in UGT_genes){
  exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  exonic_data$Functional_impact <- sapply(exonic_data$Variant_ID, function(x){if (exonic_data[x, 'VEP_Annotation'] %in% deleterious){'D'}
    else if (exonic_data[x, 'VEP_Annotation'] %in% neutral){'N'}
    ## Take missense predictions by the UGT-optimized framework
    else if (exonic_data[x, 'VEP_Annotation']=='missense_variant'){new_variants_predictions[which(new_variants_predictions$Variant_ID==x), 'UGT_optimized_pred']}
    else{'NA'} })
  assign(paste0(gene, '_exonic_data'), exonic_data)
  save(exonic_data, file = paste0('processed-data/03_Anno_functional_impact/', gene, '_exonic_data.Rdata'))
}


## Plot GMAF of all D variants per gene 

GMAFs_Dvars_genes <- vector()
for (gene in UGT_genes){
  ## All variants in the gene
  gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  ## D variants only
  GMAFs_Dvars <- gene_data[which(gene_data$Functional_impact=='D'), c('Variant_ID', 'Allele_Frequency', 'Location_in_txs')]
  
  ## Number of D vars per gene
  num_D <- dim(GMAFs_Dvars)[1]
  
  GMAFs_Dvars_genes <- rbind(GMAFs_Dvars_genes, cbind(GMAFs_Dvars, 'gene'=rep(gene, num_D)))
}

GMAFs_Dvars_genes <- as.data.frame(GMAFs_Dvars_genes)
GMAFs_Dvars_genes$gene <- factor(GMAFs_Dvars_genes$gene, levels=UGT_genes)

## New category of UGT1 and UGT2 variants: shared (in Exon 2-5/6) or unique (in Exon 1 or promoter region)
GMAFs_Dvars_genes$shared_or_unique <- apply(GMAFs_Dvars_genes, 1, function(x){if(x['gene'] %in% UGT1_genes & !x['Location_in_txs'] %in% c('Exon 1', '5\' upstream')){'UGT1A[1-10]'} else{x['gene']}})

GMAFs_Dvars_genes$shared_or_unique <- apply(GMAFs_Dvars_genes, 1, function(x){if(x['gene'] %in% c('UGT2A1', 'UGT2A2') & !x['Location_in_txs'] %in% c('Exon 1')){'UGT2A[1-2]'} else{x['shared_or_unique']}})

## Collapse shared variants 
GMAFs_Dvars_shared_or_unique <- rbind(unique(subset(GMAFs_Dvars_genes, shared_or_unique %in% c('UGT1A[1-10]', 'UGT2A[1-2]'))[c(1,2,5)]),
                                      subset(GMAFs_Dvars_genes, !shared_or_unique %in% c('UGT1A[1-10]', 'UGT2A[1-2]'))[c(1,2,5)])

## Number of D variants per gene/gene locus
num_D_per_gene_or_locus <- as.data.frame(table(GMAFs_Dvars_shared_or_unique$shared_or_unique))
colnames(num_D_per_gene_or_locus) <- c('gene_or_locus', 'number')
# UGT1A[1-10]      UGT1A1     UGT1A10      UGT1A3      UGT1A4      UGT1A5      UGT1A6      UGT1A7      UGT1A8      UGT1A9  UGT2A[1-2]      UGT2A1 
#         78          50          61          51          31          34          35          49          56          50          95          37 
# UGT2A2      UGT2A3     UGT2B10     UGT2B11     UGT2B15     UGT2B17     UGT2B28      UGT2B4      UGT2B7      UGT3A1      UGT3A2        UGT8 
#     28          56          46          61          69          40          59          48          58          34          51          29 

## Confirm number of D unique and shared variants in UGT1 and UGT2 genes
for (gene in c(UGT1_genes, 'UGT2A1', 'UGT2A2')){
  gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  ## Unique variants are in Exon 1 (and 5' upstream for UGT1A1)
  print(paste0('Unique D vars in ', gene, ': ',  sum(table(subset(gene_data, Functional_impact=='D' &
                                                              Location_in_txs %in% c('Exon 1', '5\' upstream'))$Location_in_txs))))
  print(paste0('Shared D vars in ', gene, ': ',  sum(table(subset(gene_data, Functional_impact=='D'  &
                                                                  ! Location_in_txs %in% c('Exon 1', '5\' upstream'))$Location_in_txs))))
}

# [1] "Unique D vars in UGT1A1: 50"
# [1] "Shared D vars in UGT1A1: 78"
# [1] "Unique D vars in UGT1A3: 51"
# [1] "Shared D vars in UGT1A3: 78"
# [1] "Unique D vars in UGT1A4: 31"
# [1] "Shared D vars in UGT1A4: 78"
# [1] "Unique D vars in UGT1A5: 34"
# [1] "Shared D vars in UGT1A5: 78"
# [1] "Unique D vars in UGT1A6: 35"
# [1] "Shared D vars in UGT1A6: 78"
# [1] "Unique D vars in UGT1A7: 49"
# [1] "Shared D vars in UGT1A7: 78"
# [1] "Unique D vars in UGT1A8: 56"
# [1] "Shared D vars in UGT1A8: 78"
# [1] "Unique D vars in UGT1A9: 50"
# [1] "Shared D vars in UGT1A9: 78"
# [1] "Unique D vars in UGT1A10: 61"
# [1] "Shared D vars in UGT1A10: 78"
# [1] "Unique D vars in UGT2A1: 37"
# [1] "Shared D vars in UGT2A1: 95"
# [1] "Unique D vars in UGT2A2: 28"
# [1] "Shared D vars in UGT2A2: 95"

## Order 
GMAFs_Dvars_shared_or_unique$shared_or_unique <- factor(GMAFs_Dvars_shared_or_unique$shared_or_unique, levels=c(UGT1_genes, "UGT1A[1-10]", 'UGT2A1', 'UGT2A2', "UGT2A[1-2]",  UGT2_genes[3:10], UGT3_genes, UGT8_genes))

## Reduce all GMAF<=1e-05 to '<=1e-05'
GMAFs_Dvars_shared_or_unique$Allele_Frequency <- sapply(GMAFs_Dvars_shared_or_unique$Allele_Frequency, function(x){if(x<=1e-05){1e-05} else{x}})

## Define shapes for D variants with GMAF>0.01
GMAFs_Dvars_shared_or_unique$label <- apply(GMAFs_Dvars_shared_or_unique, 1, function(x){if(as.numeric(x['Allele_Frequency'])>0.01){x['Variant_ID']} else{NA}})
GMAFs_Dvars_shared_or_unique[which(GMAFs_Dvars_shared_or_unique$Allele_Frequency>0.01),]
#            Variant_ID   Allele_Frequency  shared_or_unique
#        4-70462042-C-T         0.07988111        UGT2A[1-2]
#   2-234668879-C-CATAT         0.01559741            UGT1A1
#     2-234668879-C-CAT         0.34657642            UGT1A1
#        4-70512787-A-T         0.03228072            UGT2A1

shapes <- c('2-234668879-C-CAT'=8,
            '2-234668879-C-CATAT'=11,
            '4-70462042-C-T'=25,
            '4-70512787-A-T'=3)

## Face for x axis labels
italic.labels <- ifelse(levels(GMAFs_Dvars_shared_or_unique$shared_or_unique) %in% c('UGT1A[1-10]' ,'UGT2A[1-2]'), yes = "bold", no = "bold.italic")

ggplot(data = GMAFs_Dvars_shared_or_unique, mapping = aes(x = shared_or_unique, y = Allele_Frequency, color = shared_or_unique)) +
  geom_jitter(data=subset(GMAFs_Dvars_shared_or_unique, is.na(label)), shape=16, width = 0.3, height = 0, alpha = 0.8, size = 1.2) +
  geom_point(data=subset(GMAFs_Dvars_shared_or_unique, !is.na(label)), aes(shape=label), alpha = 0.8, size = 1.2, stroke = 1) +
  scale_y_continuous(trans='log10', breaks=c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05),
                     labels=c('-1', '-2', '-3', '-4', '≤-5')) +
  scale_color_manual(values = genes_colors) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  guides(color = 'none', fill='none', shape='none') + 
  geom_text(data = num_D_per_gene_or_locus, aes(x=gene_or_locus, label=number,  y=(1e-05)/4*3, color=NULL), size=2.3) +
  geom_text_repel(aes(label=label), color='black', fontface='bold', size=2, min.segment.length=0.9)+
  labs(x='', y='log10(MAF) of unique and shared deleterious variants in UGT genes', 
       subtitle = paste0(dim(GMAFs_Dvars_shared_or_unique)[1], ' deleterious variants across all UGT genes'),
       shape='Variant ID (MAF>0.01)') +
  theme(plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
        axis.title = element_text(size = (8.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face=italic.labels),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(filename='plots/03_Anno_functional_impact/GMAF_totalDvars_per_gene.png', width = 6.2, height = 5)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_totalDvars_per_gene.pdf', width = 6.2, height = 5)
save(GMAFs_Dvars_genes, file='processed-data/03_Anno_functional_impact/GMAFs_Dvars_genes.Rdata')
save(GMAFs_Dvars_shared_or_unique, file='processed-data/03_Anno_functional_impact/GMAFs_Dvars_shared_or_unique.Rdata')


## Plot GMAF of all D variants per gene family

GMAFs_Dvars_gene_fam <- vector()

for (fam in gene_families){
  ## Genes 
  UGT_genes <- eval(parse_expr(paste0(fam, '_genes')))
  ## Subset to unique gene family variants
  fam_vars <- unique(subset(GMAFs_Dvars_genes, gene %in% UGT_genes)[,-4])
  ## Add gene fam info
  fam_vars <- cbind(fam_vars, 'gene_fam'=rep(fam, dim(fam_vars)[1]))
  
  GMAFs_Dvars_gene_fam <- rbind(GMAFs_Dvars_gene_fam, fam_vars)
  
}
GMAFs_Dvars_gene_fam$label <- apply(GMAFs_Dvars_gene_fam, 1, function(x){if(as.numeric(x['Allele_Frequency'])>0.01){x['Variant_ID']} else{NA}})

GMAFs_Dvars_gene_fam$Allele_Frequency <- sapply(GMAFs_Dvars_gene_fam$Allele_Frequency, function(x){if(x<=1e-05){1e-05} else{x}})

## Number of D vars per family
num_D_per_gene_fam <- as.data.frame(table(GMAFs_Dvars_gene_fam$gene_fam))
colnames(num_D_per_gene_fam) <- c('gene_fam', 'number')

gene_fams_colors <- c('UGT1'='hotpink', 'UGT2'='tan1', 'UGT3'='purple1', 'UGT8'='steelblue2')

ggplot(data = GMAFs_Dvars_gene_fam, mapping = aes(x = gene_fam, y = Allele_Frequency, color = gene_fam)) +
  geom_jitter(data=subset(GMAFs_Dvars_gene_fam, is.na(label)), shape=16, width = 0.3, height = 0, alpha = 0.8, size = 1) +
  geom_point(data=subset(GMAFs_Dvars_gene_fam, !is.na(label)), aes(shape=label), alpha = 0.8, size = 1, stroke = 1) +
  scale_y_continuous(trans='log10', breaks=c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05),
                     labels=c('-1', '-2', '-3', '-4', '≤-5')) +
  scale_color_manual(values = gene_fams_colors) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  guides(color = 'none', shape='none') + 
  geom_text(data = num_D_per_gene_fam, aes(x=gene_fam, label=number,  y=(1e-05)/4*3, color=NULL), size=2.4) +
  geom_text_repel(aes(label=label), color='black', fontface='bold', size=2, min.segment.length=0.9, force_pull=2) +
  labs(x='', y='log10(GMAF) of deleterious variants per gene family', 
       subtitle = paste0(length(unique(GMAFs_Dvars_gene_fam$Variant_ID)), 
                         ' deleterious variants in UGT genes of all families'),
       shape='Variant ID (GMAF>0.01)') +
  theme(plot.subtitle = element_text(size = (9), color="gray50", face='bold'), 
        axis.title = element_text(size = (8.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(face='bold'),
        legend.title = element_text(size=8, face='bold'), 
        legend.text = element_text(size=7.5, face='bold'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(filename='plots/03_Anno_functional_impact/GMAF_totalDvars_per_gene_fam.png', width = 3.8, height = 4)
ggsave(filename='plots/03_Anno_functional_impact/GMAF_totalDvars_per_gene_fam.pdf', width = 3.8, height = 4)
save(GMAFs_Dvars_gene_fam, file='processed-data/03_Anno_functional_impact/GMAFs_Dvars_gene_fam.Rdata')







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
# date     2023-10-24
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# ! package        * version date (UTC) lib source
# bit              4.0.5   2022-11-15 [1] CRAN (R 4.3.0)
# bit64            4.0.5   2020-08-30 [1] CRAN (R 4.3.0)
# cachem           1.0.8   2023-05-01 [1] CRAN (R 4.3.0)
# callr            3.7.3   2022-11-02 [1] CRAN (R 4.3.0)
# cli              3.6.1   2023-03-23 [1] CRAN (R 4.3.0)
# colorspace       2.1-0   2023-01-23 [1] CRAN (R 4.3.0)
# corrplot       * 0.92    2021-11-18 [1] CRAN (R 4.3.0)
# cowplot        * 1.1.1   2020-12-30 [1] CRAN (R 4.3.0)
# crayon           1.5.2   2022-09-29 [1] CRAN (R 4.3.0)
# curl             5.0.1   2023-06-07 [1] CRAN (R 4.3.0)
# desc             1.4.2   2022-09-08 [1] CRAN (R 4.3.0)
# devtools       * 2.4.5   2022-10-11 [1] CRAN (R 4.3.0)
# digest           0.6.33  2023-07-07 [1] CRAN (R 4.3.0)
# distributional   0.3.2   2023-03-22 [1] CRAN (R 4.3.0)
# dplyr            1.1.2   2023-04-20 [1] CRAN (R 4.3.0)
# ellipsis         0.3.2   2021-04-29 [1] CRAN (R 4.3.0)
# V fansi            1.0.4   2023-10-08 [1] CRAN (R 4.3.1) (on disk 1.0.5)
# farver           2.1.1   2022-07-06 [1] CRAN (R 4.3.0)
# fastmap          1.1.1   2023-02-24 [1] CRAN (R 4.3.0)
# fs               1.6.3   2023-07-20 [1] CRAN (R 4.3.0)
# generics         0.1.3   2022-07-05 [1] CRAN (R 4.3.0)
# ggdist         * 3.3.0   2023-05-13 [1] CRAN (R 4.3.0)
# ggExtra        * 0.10.1  2023-08-21 [1] CRAN (R 4.3.0)
# V ggplot2        * 3.4.2   2023-10-12 [1] CRAN (R 4.3.1) (on disk 3.4.4)
# ggrepel        * 0.9.3   2023-02-03 [1] CRAN (R 4.3.0)
# ggside         * 0.2.2   2023-10-24 [1] Github (jtlandis/ggside@83002b3)
# glue             1.6.2   2022-02-24 [1] CRAN (R 4.3.0)
# V gtable           0.3.3   2023-08-21 [1] CRAN (R 4.3.0) (on disk 0.3.4)
# here           * 1.0.1   2020-12-13 [1] CRAN (R 4.3.0)
# hms              1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
# htmltools        0.5.5   2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets      1.6.2   2023-03-17 [1] CRAN (R 4.3.0)
# httpuv           1.6.11  2023-05-11 [1] CRAN (R 4.3.0)
# insight          0.19.6  2023-10-12 [1] CRAN (R 4.3.1)
# V labeling         0.4.2   2023-08-29 [1] CRAN (R 4.3.0) (on disk 0.4.3)
# later            1.3.1   2023-05-02 [1] CRAN (R 4.3.0)
# lifecycle        1.0.3   2022-10-07 [1] CRAN (R 4.3.0)
# magrittr         2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
# memoise          2.0.1   2021-11-26 [1] CRAN (R 4.3.0)
# mime             0.12    2021-09-28 [1] CRAN (R 4.3.0)
# miniUI           0.1.1.1 2018-05-18 [1] CRAN (R 4.3.0)
# munsell          0.5.0   2018-06-12 [1] CRAN (R 4.3.0)
# pillar           1.9.0   2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild         1.4.2   2023-06-26 [1] CRAN (R 4.3.0)
# pkgconfig        2.0.3   2019-09-22 [1] CRAN (R 4.3.0)
# pkgload          1.3.2.1 2023-07-08 [1] CRAN (R 4.3.0)
# plyr             1.8.8   2022-11-11 [1] CRAN (R 4.3.0)
# prettyunits      1.1.1   2020-01-24 [1] CRAN (R 4.3.0)
# pROC           * 1.18.4  2023-07-06 [1] CRAN (R 4.3.0)
# processx         3.8.2   2023-06-30 [1] CRAN (R 4.3.0)
# profvis          0.3.8   2023-05-02 [1] CRAN (R 4.3.0)
# promises         1.2.0.1 2021-02-11 [1] CRAN (R 4.3.0)
# ps               1.7.5   2023-04-18 [1] CRAN (R 4.3.0)
# purrr            1.0.1   2023-01-10 [1] CRAN (R 4.3.0)
# R6               2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
# ragg             1.2.5   2023-01-12 [1] CRAN (R 4.3.0)
# Rcpp             1.0.11  2023-07-06 [1] CRAN (R 4.3.0)
# readr          * 2.1.4   2023-02-10 [1] CRAN (R 4.3.0)
# remotes          2.4.2.1 2023-07-18 [1] CRAN (R 4.3.0)
# reshape2       * 1.4.4   2020-04-09 [1] CRAN (R 4.3.0)
# rlang          * 1.1.1   2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot        2.0.3   2022-04-02 [1] CRAN (R 4.3.0)
# rstudioapi       0.15.0  2023-07-07 [1] CRAN (R 4.3.0)
# scales           1.2.1   2022-08-20 [1] CRAN (R 4.3.0)
# see            * 0.8.0   2023-06-05 [1] CRAN (R 4.3.0)
# sessioninfo    * 1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
# shiny            1.7.4.1 2023-07-06 [1] CRAN (R 4.3.0)
# stringi          1.7.12  2023-01-11 [1] CRAN (R 4.3.0)
# stringr          1.5.0   2022-12-02 [1] CRAN (R 4.3.0)
# systemfonts      1.0.4   2022-02-11 [1] CRAN (R 4.3.0)
# textshaping      0.3.6   2021-10-13 [1] CRAN (R 4.3.0)
# tibble           3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect       1.2.0   2022-10-10 [1] CRAN (R 4.3.0)
# tzdb             0.4.0   2023-05-12 [1] CRAN (R 4.3.0)
# urlchecker       1.0.1   2021-11-30 [1] CRAN (R 4.3.0)
# usethis        * 2.2.2   2023-07-06 [1] CRAN (R 4.3.0)
# V utf8             1.2.3   2023-10-22 [1] CRAN (R 4.3.1) (on disk 1.2.4)
# V vctrs            0.6.3   2023-10-12 [1] CRAN (R 4.3.1) (on disk 0.6.4)
# vroom            1.6.3   2023-04-28 [1] CRAN (R 4.3.0)
# V withr            2.5.0   2023-09-26 [1] CRAN (R 4.3.1) (on disk 2.5.1)
# xtable           1.8-4   2019-04-21 [1] CRAN (R 4.3.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# V ── Loaded and on-disk version mismatch.