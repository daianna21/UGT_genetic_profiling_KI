
library(here)
library(rlang)
library(ggplot2)
library(cowplot)

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

## Define effect of variant types




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
    if (length(strsplit(UGT_missense_vars$Reference[i], '')[[1]])==1 & UGT_missense_vars$Reference[i] %in% c('A', 'T', 'C', 'G')){TRUE}
    else {FALSE}
  }))), ' ',
  
  names(table(
    sapply(1:dim(UGT_missense_vars)[1], function(i){
      if (length(strsplit(UGT_missense_vars$Alternate[i], '')[[1]])==1 & UGT_missense_vars$Alternate[i] %in% c('A', 'T', 'C', 'G')){TRUE}
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
  save(missense_vars_ANNOVAR_format, file = paste0('processed-data/03_Anno_functional_impact/', 
                                                   gene, '_missense_vars_ANNOVAR_format.avinput'))
}


