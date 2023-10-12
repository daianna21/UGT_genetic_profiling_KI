
library(here)
library(rlang)
library(ggplot2)
library(cowplot)
library(phylotools)
library(sessioninfo)


####################################################################################################
##                                2. Gene constraint analysis
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
  exonic_vars <- eval(parse_expr(load(here(paste0('~/Desktop/UGT_genetic_profiling_KI/processed-data/01_Data_Processing/', gene, '_exonic_data.Rdata')),
       verbose=TRUE)))
  assign( paste0(gene, '_exonic_data'), exonic_vars)
}

## Load info of canonical txs 
for (gene_family in gene_families){
  load(here(paste0('~/Desktop/UGT_genetic_profiling_KI/processed-data/01_Data_Processing/canonical_', gene_family, '_txs.Rdata')),
       verbose=TRUE)
}


#####################################################################################
##  2.1 Compare expected vs observed proportions of missense and synonymous variants
#####################################################################################

################# Observed proportions in each gene #################

## Proportions given by # missense variants / # synonymous variants
for (gene_family in gene_families){
  
  genes<- eval(parse_expr(paste0(gene_family, '_genes')))
  
  for (gene in genes){
    exonic_gene_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
    print(paste0(gene, ': ', signif(table(exonic_gene_data$VEP_Annotation)['missense_variant'] 
                                    / table(exonic_gene_data$VEP_Annotation)['synonymous_variant'], digits = 4) ))
    
  }
}

######################
####  UGT1 variants
######################

# [1] "UGT1A1: 2.384"
# [1] "UGT1A3: 2.335"
# [1] "UGT1A4: 2.275"
# [1] "UGT1A5: 2.085"
# [1] "UGT1A6: 2.406"
# [1] "UGT1A7: 2.188"
# [1] "UGT1A8: 2.107"
# [1] "UGT1A9: 2.241"
# [1] "UGT1A10: 2.23"

######################
####  UGT2 variants
######################

# [1] "UGT2A1: 3.118"
# [1] "UGT2A2: 3.056"
# [1] "UGT2A3: 2.725"
# [1] "UGT2B4: 2.784"
# [1] "UGT2B7: 2.927"
# [1] "UGT2B10: 3.211"
# [1] "UGT2B11: 3.099"
# [1] "UGT2B15: 2.798"
# [1] "UGT2B17: 2.924"
# [1] "UGT2B28: 2.831"

######################
####  UGT3 variants
######################

# [1] "UGT3A1: 2.382"
# [1] "UGT3A2: 2.424"

######################
####  UGT8 variants
######################

# [1] "UGT8: 1.606"


################# Expected proportions in each gene #################

for (gene_family in gene_families){
  
  genes <- eval(parse_expr(paste0(gene_family, '_genes')))
  txs <- eval(parse_expr(paste0('canonical_', gene_family, '_txs')))
  
  for (gene in genes){
    
    tx <- txs[[gene]]
    
    ## Read coding sequence (CDS) and protein sequence of each gene tx
    fastaFile = read.fasta(paste0("raw-data/CDS_seq_data/Homo_sapiens_", tx, "_sequence.fasta"))
    ## Extract CDS and peptide sequence
    cds_sequence = strsplit(fastaFile$seq.text[1], '')[[1]] 
    ## Remove initial 3 nts (to avoid start codon variants)
    cds_sequence = cds_sequence[-c(1:3)]
    ## Remove last 3 nts (avoid stop codon variants)
    cds_sequence = cds_sequence[-(c(length(cds_sequence)-2):length(cds_sequence))]
    pep_sequence = strsplit(fastaFile$seq.text[2], '')[[1]]
    ## Remove initial aa (Met)
    pep_sequence = pep_sequence[-1]
    
    ## Verify # nts in CDS = # amino acids x 3
    if (length(cds_sequence) == length(pep_sequence) *3) {
      
      nts <- c('A', 'G', 'T', 'C')
      effects <- vector()
      i=1
      j=1
      while (i< length(cds_sequence)-1) {
        
        ## Each codon
        codon = cds_sequence[i:(i+2)]
        ## Codified amino acid 
        aa = pep_sequence[j]
        
        ## 9 possible point mutations for codon
        mutations <- vector()
        for (nt in 1:3){
          mutations <- append(mutations, unlist(sapply(nts, function(nt2){
            if (codon[nt]!=nt2){
              ## Replace nt in codon 
              codon[nt] <- nt2
              paste(codon, collapse = '')} 
          })))
        }
        
        ## Evaluate if each mutation changes aa or not (i.e. is missense or synonymous)
        ## Exclude stop-gained mutations (that result in a stop codon identified as '*')
        effects <- append(effects, sapply(mutations, function(mutation){if (GENETIC_CODE[mutation] != "*"){
          if (GENETIC_CODE[mutation] != aa){'missense'} else{'synonymous'}
        } else {'stop gained'}
        }))
        
        ## Next codon and aa
        i=i+3
        j=j+1
      }
      
      ## Proportion 
      assign(paste0(gene, '_expected_variant_num'), effects)
      print(paste0(gene, ': ', signif(table(effects)['missense']/table(effects)['synonymous'], digits=4)))
      
    }
  }
}

######################
####  UGT1 variants
######################

# [1] "UGT1A1: 3.176"
# [1] "UGT1A3: 3.226"
# [1] "UGT1A4: 3.149"
# [1] "UGT1A5: 3.163"
# [1] "UGT1A6: 3.292"
# [1] "UGT1A7: 3.256"
# [1] "UGT1A8: 3.263"
# [1] "UGT1A9: 3.286"
# [1] "UGT1A10: 3.309"

######################
####  UGT2 variants
######################

# [1] "UGT2A1: 3.381"
# [1] "UGT2A2: 3.371"
# [1] "UGT2A3: 3.384"
# [1] "UGT2B4: 3.41"
# [1] "UGT2B7: 3.502"
# [1] "UGT2B10: 3.496"
# [1] "UGT2B11: 3.493"
# [1] "UGT2B15: 3.484"
# [1] "UGT2B17: 3.549"
# [1] "UGT2B28: 3.5"

######################
####  UGT3 variants
######################

# [1] "UGT3A1: 3.254"
# [1] "UGT3A2: 3.34"

######################
####  UGT8 variants
######################
# [1] "UGT8: 3.254"


## Create table with constraint metrics for each gene

constraint_metrics <- data.frame(matrix(ncol = 15, nrow = 22))
colnames(constraint_metrics) <- c('Gene', 'Gene_family', 'Num_miss_ob', 'Num_syn_ob', 'Prop_ob', 'Num_miss_exp', 'Num_syn_exp', 'Prop_exp', 'diff_Prop_ob_exp',
                                  'Num_miss_ob_gnomAD', 'Num_syn_ob_gnomAD', 'Num_miss_exp_gnomAD', 'Num_syn_exp_gnomAD', 'oe_miss', 'oe_syn')
constraint_metrics$Gene <- UGT_genes
constraint_metrics$Gene_family <- rep(c('UGT1', 'UGT2', 'UGT3', 'UGT8'), c(9,10,2,1))
constraint_metrics$Num_miss_exp <- sapply(paste0('table(', UGT_genes, '_expected_variant_num)[\'missense\']'), function(x){eval(parse_expr(x))})
constraint_metrics$Num_syn_exp <-  sapply(paste0('table(', UGT_genes, '_expected_variant_num)[\'synonymous\']'), function(x){eval(parse_expr(x))})
constraint_metrics$Prop_exp <- sapply(paste0('signif(table(', UGT_genes, '_expected_variant_num)[\'missense\']/table(', UGT_genes,
                                             '_expected_variant_num)[\'synonymous\'], digits=4)'), function(x){eval(parse_expr(x))})
constraint_metrics$diff_Prop_ob_exp <- constraint_metrics$Prop_ob - constraint_metrics$Prop_exp
constraint_metrics$Num_miss_ob <- sapply(paste0('table(', UGT_genes, '_exonic_data$VEP_Annotation)[\'missense_variant\']'), function(x){eval(parse_expr(x))})
constraint_metrics$Num_syn_ob <- sapply(paste0('table(', UGT_genes, '_exonic_data$VEP_Annotation)[\'synonymous_variant\']'), function(x){eval(parse_expr(x))})
constraint_metrics$Prop_ob <- sapply(paste0('signif(table(', UGT_genes, '_exonic_data$VEP_Annotation)[\'missense_variant\']/table(', UGT_genes,
                                            '_exonic_data$VEP_Annotation)[\'synonymous_variant\'], digits=4)'), function(x){eval(parse_expr(x))})

## Metrics from gnomAD:
constraint_metrics$Num_miss_ob_gnomAD <- c(309, 360, 378, 314, 294, 314, 360, 346, 366, 340, 358, 354, 331, 338, 397, 459, 313, 238, 373, 290, 277, 210)
constraint_metrics$Num_syn_ob_gnomAD <- c(126, 157, 169, 151, 120, 150, 157, 149, 166, 107, 118, 131, 123, 104, 122, 144, 117, 86, 125, 120, 110, 124)
constraint_metrics$Num_miss_exp_gnomAD <- c(297, 294.9, 298, 299.6, 295.4, 294.1, 294.9, 300.3, 297.9, 274.7, 278.1, 272.5, 277, 275.9, 265.8, 268.5, 272.4, 
                                            269.9, 271.9, 283.1, 284.7, 289.5)
constraint_metrics$Num_syn_exp_gnomAD <- c(122.9, 118.1, 120.3, 120.5, 114.6, 115.1, 118.1, 116.3, 117, 97, 96.8, 96.1, 98.6, 96.6, 93.6, 93.3, 97.5, 95.6, 
                                           93.6, 106.6, 108.4, 111.5)
constraint_metrics$oe_miss <- c(1.04, 1.22, 1.27, 1.05, 1, 1.07, 1.22, 1.15, 1.23, 1.24, 1.29, 1.3, 1.19, 1.22, 1.49, 1.71, 1.15, 0.88, 1.37, 1.02, 0.97, 0.73)
constraint_metrics$oe_syn <- c(1.03, 1.33, 1.4, 1.25, 1.05, 1.3, 1.33, 1.28, 1.42, 1.1, 1.22, 1.36, 1.25, 1.08, 1.3, 1.54, 1.2, 0.9, 1.34, 1.13, 1.01, 1.11)

save(constraint_metrics, file = paste0('processed-data/02_Constraint_analysis/constraint_metrics.Rdata'))


## Correlation plots

colors_gene_fam <- list('UGT1'='hotpink', 'UGT2'='tan1', 'UGT3'='purple1', 'UGT8'='steelblue2')

###############  Correlation between manually computed observed and expected miss/syn proportions  ###############

p1 <- ggplot(constraint_metrics, aes(x=Prop_exp, y=Prop_ob, color = Gene_family)) +
  ## Add scatterplot
  geom_point(size=2) +
  ## Add regression line
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  ## Colors
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  ## Add Pearson correlation coefficient as subtitle
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$Prop_exp, constraint_metrics$Prop_ob, method = "pearson"), digits = 3)),
    ## Add axis labels
    x = 'Expected proportion of missense/synonymous variants',
    y = 'Observed proportion of missense/synonymous variants',
    color='Gene family'
  ) +
  ## Plot margins and text size
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between manually computed observed miss/syn proportion and missense oe ratio  ###############  

p2 <- ggplot(constraint_metrics, aes(x=Prop_ob, y=oe_miss, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$Prop_ob, constraint_metrics$oe_miss, method = "pearson"), digits = 3)),
    x = 'Observed proportion of missense/synonymous variants',
    y = 'oe ratio for number of missense variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between manually computed expected miss/syn proportion and missense oe ratio  ###############  

p3 <- ggplot(constraint_metrics, aes(x=Prop_exp, y=oe_miss, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$Prop_exp, constraint_metrics$oe_miss, method = "pearson"), digits = 3)),
    x = 'Expected proportion of missense/synonymous variants',
    y = 'oe ratio for number of missense variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between manually computed observed miss/syn proportion and synonymous oe ratio  ###############  

p4 <- ggplot(constraint_metrics, aes(x=Prop_ob, y=oe_syn, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$Prop_ob, constraint_metrics$oe_syn, method = "pearson"), digits = 3)),
    x = 'Observed proportion of missense/synonymous variants',
    y = 'oe ratio for number of synonymous variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between manually computed expected miss/syn proportion and synonymous oe ratio  ###############  

p5 <- ggplot(constraint_metrics, aes(x=Prop_exp, y=oe_syn, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$Prop_exp, constraint_metrics$oe_syn, method = "pearson"), digits = 3)),
    x = 'Expected proportion of missense/synonymous variants',
    y = 'oe ratio for number of synonymous variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between missense oe ratio and observed - expected miss/syn proportion ###############  

p6 <- ggplot(constraint_metrics, aes(x=diff_Prop_ob_exp, y=oe_miss, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$diff_Prop_ob_exp, constraint_metrics$oe_miss, method = "pearson"), digits = 3)),
    x = 'Observed - Expected proportion of missense/synonymous variants',
    y = 'oe ratio for number of missense variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


###############  Correlation between synonymous oe ratio and observed - expected miss/syn proportion ###############  

p7 <- ggplot(constraint_metrics, aes(x=diff_Prop_ob_exp, y=oe_syn, color = Gene_family)) +
  geom_point(size=2) +
  stat_smooth(geom = "line", alpha = 0.6, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
  scale_color_manual(values = colors_gene_fam) +
  theme_bw() +
  labs(
    subtitle = paste0("Corr: ", signif(cor(constraint_metrics$diff_Prop_ob_exp, constraint_metrics$oe_syn, method = "pearson"), digits = 3)),
    x = 'Observed - Expected proportion of missense/synonymous variants',
    y = 'oe ratio for number of synonymous variants (gnomAD)',
    color='Gene family'
  ) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title = element_text(size = (11)),
    axis.text = element_text(size = (10)),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size =10))


plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow=2, rel_widths = c(1,1,1,1,1,1.5,1.5))
ggsave('plots/02_Constraint_analysis/Corr_plots.pdf', width = 25, height = 8.5)







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
# date     2023-09-22
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
# phylotools           * 0.2.4     2023-08-31 [1] Github (helixcn/phylotools@758d338)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
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

