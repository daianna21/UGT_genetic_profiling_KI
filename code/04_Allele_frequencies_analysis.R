

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