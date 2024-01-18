
library(here)
library(ggplot2)
library(scales)
library(rlang)
library(reshape2)
library(sessioninfo)


####################################################################################################
##                    4. Population-scale analysis of deleterious UGT variants
####################################################################################################

gene_families <- c('UGT1', 'UGT2', 'UGT3', 'UGT8')

UGT1_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10')
UGT2_genes <- c('UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28')
UGT3_genes <- c('UGT3A1', 'UGT3A2')
UGT8_genes <- c('UGT8')

UGT_genes <- c('UGT1A1', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT1A10',
               'UGT2A1', 'UGT2A2', 'UGT2A3', 'UGT2B4', 'UGT2B7', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28',
               'UGT3A1', 'UGT3A2', 
               'UGT8')

## Load exonic vars per gene
for (gene in UGT_genes){
  exonic_vars <- eval(parse_expr(load(here(paste0('~/Desktop/UGT_genetic_profiling_KI/processed-data/03_Anno_functional_impact/', 
                                                  gene, '_exonic_data.Rdata')) )))
  assign(paste0(gene, '_exonic_data'), exonic_vars)
}
## Load Variant - GMAF info - Gene/Locus
load(here('processed-data/03_Anno_functional_impact/GMAFs_Dvars_shared_or_unique.Rdata'))
## Add gene family info
GMAFs_Dvars_shared_or_unique$gene_fam <- substring(GMAFs_Dvars_shared_or_unique$shared_or_unique, 1, 4)



# _______________________________________________________________________________
#  4.1   Examine MAF of deleterious UGT variants within each population
# _______________________________________________________________________________

## Plot 
populations <- c('African_or_African_American',
                 'Latino_or_Admixed_American',
                 'East_Asian',
                 'South_Asian',
                 'European_Finnish',
                 'European_non_Finnish',
                 'Ashkenazi_Jewish')

allD_vars_MAF_in_pops <- vector()
for (gene in UGT_genes){
  
  ## Exonic variants per gene
  exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  
  ## Obtain MAF of variants within each population 
  for (p in populations){
    exonic_data[,paste0('MAF_', p)] <- exonic_data[,paste0('Allele_Count_', p)]/exonic_data[,paste0('Allele_Number_', p)]
    ## If no MAF available in the population
    replace(exonic_data[,paste0('MAF_', p)], which(is.nan(exonic_data[,paste0('MAF_', p)])), NA)
  }
  ## Subset to D vars 
  data <- exonic_data[which(exonic_data$Functional_impact=='D'),c('Variant_ID', 'Allele_Frequency', paste0('MAF_', populations))]
  data <- melt(data, id.vars = c('Variant_ID'))
  colnames(data) <- c('Variant_ID', 'Group', 'MAF')
  allD_vars_MAF_in_pops <- rbind(allD_vars_MAF_in_pops, data)
  
}

allD_vars_MAF_in_pops <- unique(allD_vars_MAF_in_pops)
#allD_vars_MAF_in_pops$Group <- factor(allD_vars_MAF_in_pops$Group, levels=c(paste0('MAF_',populations), 'Allele_Frequency'))

## Number of variants with missing MAF per population
table(allD_vars_MAF_in_pops[which(is.na(allD_vars_MAF_in_pops$MAF)), 'Group'])
# MAF_African_or_African_American  MAF_Latino_or_Admixed_American                  MAF_East_Asian                 MAF_South_Asian 
#                               0                               0                               0                              75 
# MAF_European_Finnish        MAF_European_non_Finnish            MAF_Ashkenazi_Jewish                Allele_Frequency 
#                    0                               0                               0                               0 

## Manually add MAF of 2−234668879−C−CAT (UGT1A1*28) and 2-234668879-C-CATAT (UGT1A1*37) variants in South Asians 
allD_vars_MAF_in_pops[which(allD_vars_MAF_in_pops$Group=='MAF_South_Asian' & allD_vars_MAF_in_pops$Variant_ID=='2-234668879-C-CAT'), 'MAF'] <- 0.4557
allD_vars_MAF_in_pops[which(allD_vars_MAF_in_pops$Group=='MAF_South_Asian' & allD_vars_MAF_in_pops$Variant_ID=='2-234668879-C-CATAT'), 'MAF'] <- 0

## Label variants with MAF>0.01 or if they are the UGT1A1*28 allele; ignore missing MAFs
allD_vars_MAF_in_pops$Label <- apply(allD_vars_MAF_in_pops, 1, function(x){if (is.na(x['MAF'])){NA}
                                                                           else if(as.numeric(x['MAF'])>=0.01){x['Variant_ID']}
                                                                           else if(x['Variant_ID']=='2-234668879-C-CAT'){x['Variant_ID']}
                                                                           else{NA}})
## Order
allD_vars_MAF_in_pops$Label <- factor(allD_vars_MAF_in_pops$Label, levels=c("2-234668879-C-CAT", "2-234668879-C-CATAT",
                                                                            "2-234638282-G-GT", "2-234622331-GC-G",
                                                                            "2-234676872-C-T", "4-70512787-A-T",
                                                                            "4-70462042-C-T", "4-69811110-A-C", 
                                                                            "4-69693141-GT-G", "4-70078393-C-T",
                                                                            "4-69536234-G-T", "4-69512937-T-A", NA))

## Add in which genes those variants are
commonDvars_genes <- sapply(levels(allD_vars_MAF_in_pops$Label), function(x){GMAFs_Dvars_shared_or_unique[which(GMAFs_Dvars_shared_or_unique$Variant_ID==x), 'shared_or_unique']})
variant_labels_withGene <- paste0(levels(allD_vars_MAF_in_pops$Label), ' (', commonDvars_genes, ')')

## Different shapes for D variants with MAF>=0.01 in at least one population
shapes <- c('2-234668879-C-CAT'=8,
            '2-234668879-C-CATAT'=11,
            '4-70462042-C-T'=25,
            '4-70512787-A-T'=3,
            '2-234622331-GC-G'=14,
            '2-234638282-G-GT'=7,
            '2-234676872-C-T'=9,
            '4-69512937-T-A'=1,
            '4-69536234-G-T'=12,
            '4-69693141-GT-G'=5,
            '4-69811110-A-C'=4,
            '4-70078393-C-T'=0)


# ----------------------------------------------------------------------------------------------------------------
#  4.1.1  Expected number of deleterious UGT variants per individual in each population                 


## Number of D variants per individual in each population
num_pop <- list()

for (group in names(table(allD_vars_MAF_in_pops$Group))){
  
  ## Aggregated frequency of all D UGT variants (Qk) = number of D UGT variants per haploid genome 
  allDvars_agg <- c('all_vars'=sum(subset(allD_vars_MAF_in_pops, !is.na(MAF) & Group==group)$MAF))
  ## Total D UGT variants per diploid genome (2Qk)
  allDvars_num <- 2*allDvars_agg
  
  ## Number of unique and shared D variants in UGT genes in a diploid genome (2Qg)
  group_agg <- vector()
  group_num <- vector()
  for (gene_locus_group in names(table(GMAFs_Dvars_shared_or_unique$shared_or_unique))){
    ## Group variants
    group_vars <- subset(GMAFs_Dvars_shared_or_unique, shared_or_unique==gene_locus_group)$Variant_ID
    ## Aggregated MAF of variants in such genetic group (Qg)
    group_agg <- sum(subset(allD_vars_MAF_in_pops, !is.na(MAF) & Group==group & Variant_ID %in% group_vars)$MAF)
    ## Variants per individual 
    group_num[gene_locus_group] <- 2*group_agg
  }

  ## Number of D UGT variants in each gene family per diploid genome (2Qf) (UGT8 is the same as the gene) 
  gene_fam_num <- vector()
  for(fam in gene_families[-4]){
    ## Family variants
    fam_vars <- subset(GMAFs_Dvars_shared_or_unique, gene_fam==fam)$Variant_ID
    ## Aggregated frequency of family variants (Qf)
    fam_agg <- sum(subset(allD_vars_MAF_in_pops, !is.na(MAF) & Group==group & Variant_ID %in% fam_vars)$MAF)
    ## Number
    gene_fam_num[fam] <- (2*fam_agg)
  }
  
  ## Concatenate
  num_pop[[group]] <- c(allDvars_num, gene_fam_num, group_num)

}

## As data frame
num_pop <- as.data.frame(num_pop)
num_pop$genetic_group <- rownames(num_pop)

## Confirm 2Q = sum(2Qf's)
apply(num_pop[,-9], 2, function(x){c(signif(x['all_vars'])==signif(sum(x[gene_families])))})
# MAF_African_or_African_American  MAF_Latino_or_Admixed_American                  MAF_East_Asian                 MAF_South_Asian
#                            TRUE                            TRUE                            TRUE                            TRUE
#            MAF_European_Finnish        MAF_European_non_Finnish            MAF_Ashkenazi_Jewish                Allele_Frequency
#                            TRUE                            TRUE                            TRUE                            TRUE

## Confirm 2Qf = sum(2Qg) taking unique variants per gene and shared variants only once
apply(num_pop[,-9], 2, function(x){c(signif(x['UGT1'])==signif(sum(x[5:14])),
                                         signif(x['UGT2'])==signif(sum(x[15:25])),
                                         signif(x['UGT3'])==signif(sum(x[26:27])),
                                         signif(x['UGT8'])==signif(sum(x[28])))})
#      MAF_African_or_African_American MAF_Latino_or_Admixed_American MAF_East_Asian MAF_South_Asian MAF_European_Finnish
# UGT1                            TRUE                           TRUE           TRUE            TRUE                 TRUE
# UGT2                            TRUE                           TRUE           TRUE            TRUE                 TRUE
# UGT3                            TRUE                           TRUE           TRUE            TRUE                 TRUE
# UGT8                            TRUE                           TRUE           TRUE            TRUE                 TRUE
#      MAF_European_non_Finnish MAF_Ashkenazi_Jewish Allele_Frequency
# UGT1                     TRUE                 TRUE             TRUE
# UGT2                     TRUE                 TRUE             TRUE
# UGT3                     TRUE                 TRUE             TRUE
# UGT8                     TRUE                 TRUE             TRUE


## Subset to genetic groups
num_pop_in_groups <- subset(num_pop, ! genetic_group %in% num_pop$genetic_group[1:4])
## Number of variants in each genetic and human group
num_pop_in_groups <- melt(num_pop_in_groups)
num_pop_in_groups$genetic_group <- factor(num_pop_in_groups$genetic_group, 
                                                levels=num_pop$genetic_group[-c(1:4)])
num_pop_in_groups$gene_fam <- substring(num_pop_in_groups$genetic_group, 1,4)
colnames(num_pop_in_groups) <- c('genetic_group',  'group', 'num', 'gene_fam')

# Qf in each population
Qf_pop <- num_pop[gene_families, ]
Qf_pop$gene_fam <- Qf_pop$genetic_group
Qf_pop <- melt(Qf_pop)
Qf_pop$label <- signif(Qf_pop$value, digits=2)

## Stacked barplot 
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

## Order populations by Qf of UGT1s, global at the end
pop_Qf_UGT1_order <- c(as.vector(subset(Qf_pop, gene_fam=="UGT1" & variable!='Allele_Frequency')[order(subset(Qf_pop, gene_fam=="UGT1" & variable!='Allele_Frequency')$label,
                                                                                           decreasing = TRUE), 'variable']), 'Allele_Frequency')
num_pop_in_groups$group <- factor(num_pop_in_groups$group, levels=pop_Qf_UGT1_order)
italic.labels <- ifelse(levels(num_pop_in_groups$genetic_group) %in% c('UGT1A[1-10]' ,'UGT2A[1-2]'), yes = "bold", no = "italic")

ggplot(num_pop_in_groups) +
  geom_bar(aes(x = group, y = num, fill = genetic_group), 
           position = "stack", stat = "identity") +
  geom_text(data=Qf_pop, aes(label=label, y=value, x=variable, fill=NULL),
            vjust=-0.25, size=2) +
  facet_wrap(~ gene_fam, scales="free_y", ncol=4) + 
  ## Colors by 2nd variable
  scale_fill_manual(values = genes_colors) +
  theme_bw() +
  labs(x='', y = 'Averange number of deleterious UGT variants per individual', fill='Unique gene region/ Overlapping region') +
  scale_x_discrete(breaks=c(paste0('MAF_',populations), 'Allele_Frequency'),
                   labels=c("African/African American",
                            "Latino/Admixed American",
                            "East Asian",
                            "South Asian",
                            "Finnish",
                            "European non Finnish",
                            "Ashkenazi Jewish",
                            "Global")) +
  theme(axis.title = element_text(size = (11), face='bold'),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold', size=10),
        legend.title = element_text(size=11, face='bold'), 
        legend.text = element_text(size=10, face=italic.labels), 
        strip.background = element_rect(fill="gray95", size=1, color="gray60"),
        strip.text = element_text(face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(filename='plots/04_Population_scale_analysis/Num_Dvars_per_individual.pdf', width = 14, height = 7)   

# ----------------------------------------------------------------------------------------------------------------

## Plot
allD_vars_MAF_in_pops$Group <- factor(allD_vars_MAF_in_pops$Group, levels=pop_Qf_UGT1_order)

ggplot(data = allD_vars_MAF_in_pops, mapping = aes(x = Group, y = MAF, color = Group)) +
  geom_point(data=subset(allD_vars_MAF_in_pops, is.na(Label)), alpha = 0.65, size = 1.3, 
             position = position_jitter(width = 0.1, height = 0), color="tomato") +
  geom_point(data=subset(allD_vars_MAF_in_pops, !is.na(Label)), aes(shape=Label), size=1.5, 
             position = position_jitter(width = 0.3, height = 0), color='tomato3', stroke = 0.9) +
  theme_bw() +
  scale_shape_manual(values=shapes[subset(allD_vars_MAF_in_pops, !is.na(Label))$Variant_ID], 
                     labels = variant_labels_withGene) + 
  scale_x_discrete(breaks=c(paste0('MAF_',populations), 'Allele_Frequency'),
                   labels=c("African/African American",
                            "Latino/Admixed American",
                            "East Asian",
                            "South Asian",
                            "Finnish",
                            "European non Finnish",
                            "Ashkenazi Jewish",
                            "Global")) +
  labs(title='Deleterious variants in all UGT genes per population', 
      # subtitle=paste0(length(unique(allD_vars_MAF_in_pops$Variant_ID)), ' deleterious variants in total'), 
       x='', y='MAF of deleterious UGT variants in each human group', shape=paste0('Variant ID (MAF>0.01) & Gene(s)')) +
  ## Add aggregated frequency of all D variants per population
  coord_cartesian(ylim = c(0, max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)), clip = 'off') +
  geom_text(data=melt(num_pop['all_vars',]), aes(x=variable, y=max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)+0.035, shape=NULL, label=signif(value, 3)), size=3, color='grey20', fontface='bold') +
  
  theme(plot.title = element_text(size = (11), face='bold', vjust = 7.1, hjust=0),
        #plot.subtitle = element_text(size = 9, color = "gray50", vjust = 7, hjust=0, face='bold'),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "cm"),
        axis.title = element_text(size = (11), face='bold'),
        axis.text = element_text(size = (10)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.title = element_text(size=11, face='bold'), 
        legend.text = element_text(size=10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.2))
        # panel.grid.minor = element_blank(), 
        # panel.grid.major = element_blank())

ggsave(filename='plots/04_Population_scale_analysis/MAF_totalDvars_per_population.pdf', width = 9, height = 6.5)   



# _______________________________________________________________________________
#  4.2   Compare MAF of the regulatory UGT1A1 variants in each population
# _______________________________________________________________________________

## Variants are annotated in UGT1A8 dataset                                                                                           
var_data <- UGT1A8_data[UGT1A8_data$rsIDs=='rs34983651',]   

MAF_var_data <- data.frame(matrix(ncol=8, nrow=5))
colnames(MAF_var_data) <- c(paste0('MAF_',populations), 'Allele_Frequency')
rownames(MAF_var_data) <- rownames(var_data)

## MAF in each population
for (p in populations){
  MAF_var_data[,paste0('MAF_', p)] <- var_data[,paste0('Allele_Count_', p)]/var_data[,paste0('Allele_Number_', p)]
}

MAF_var_data$Allele_Frequency <- var_data$Allele_Frequency

## Add MAF of variants in South Asians 
MAF_var_data['2-234668879-C-CAT', 'MAF_South_Asian'] <- 0.4557
MAF_var_data['2-234668879-C-CATAT', 'MAF_South_Asian'] <- 0
MAF_var_data['2-234668879-CAT-C', 'MAF_South_Asian'] <- 0
MAF_var_data['2-234668879-C-CATATAT', 'MAF_South_Asian'] <- 0
MAF_var_data['2-234668879-C-CATATATAT', 'MAF_South_Asian'] <- 0

## Add MAF of reference allele (1 - MAF of all alternate alleles)
MAF_var_data['2-234668879-C-C', ] <- 1-apply(MAF_var_data, 2, sum)
MAF_var_data$variant <- factor(rownames(MAF_var_data), levels=c('2-234668879-CAT-C', '2-234668879-C-C', '2-234668879-C-CAT',
                               '2-234668879-C-CATAT', '2-234668879-C-CATATAT', '2-234668879-C-CATATATAT')) 
MAF_var_data <- as.data.frame(melt(MAF_var_data))
MAF_var_data$label <- sapply(MAF_var_data$value, function(x){if(x==0){NA}else{100*signif(x, digits=3)}})

## Pie charts
ggplot(MAF_var_data, aes(x="", y=value, fill=variant)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  facet_wrap(~ variable, labeller = labeller(variable = c("MAF_African_or_African_American" = "African/African American",
                                                          "MAF_Latino_or_Admixed_American" = "Latino/Admixed American",
                                                          "MAF_East_Asian" = "East Asian",
                                                          "MAF_South_Asian"="South Asian",
                                                          "MAF_European_Finnish"="Finnish",
                                                          "MAF_European_non_Finnish"="European non Finnish",
                                                          "MAF_Ashkenazi_Jewish"="Ashkenazi Jewish",
                                                          "Allele_Frequency"="Global")), ncol=3) + 
  labs(fill='UGT1A1 promoter variant')+
  theme_void() + 
  geom_label_repel(aes(label = label, fill=variant), position=position_stack(vjust=0.5),
                   force_pull=3, size=3, label.padding=0.1, show.legend=FALSE) +
  scale_fill_manual(values = c('2-234668879-C-CATAT'='mediumaquamarine',
                               '2-234668879-C-CATATATAT'='gold2',
                               '2-234668879-C-CATATAT'='lightskyblue',
                               '2-234668879-CAT-C'='pink2',
                               '2-234668879-C-CAT'='salmon',
                               '2-234668879-C-C'='peachpuff2'), 
                    labels=c('2-234668879-C-CATAT'='(TA)8 (UGT1A1*37)',
                             '2-234668879-C-CATATATAT'='(TA)10',
                             '2-234668879-C-CATATAT'='(TA)9',
                             '2-234668879-CAT-C'='(TA)5 (UGT1A1*36)',
                             '2-234668879-C-CAT'='(TA)7 (UGT1A1*28)',
                             '2-234668879-C-C'='(TA)6 (UGT1A1*1)')) +
  theme(legend.title = element_text(size=8.5, face='bold'), 
        legend.text = element_text(size=8),
        strip.text = element_text(face="bold"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))

ggsave(filename='plots/04_Population_scale_analysis/UGT1A1_promoter_variants.pdf', width = 8, height = 8)  







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
# date     2023-11-12
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
# RColorBrewer     1.1-3   2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp             1.0.11  2023-07-06 [1] CRAN (R 4.3.0)
# readr          * 2.1.4   2023-02-10 [1] CRAN (R 4.3.0)
# remotes          2.4.2.1 2023-07-18 [1] CRAN (R 4.3.0)
# reshape2       * 1.4.4   2020-04-09 [1] CRAN (R 4.3.0)
# rlang          * 1.1.1   2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot        2.0.3   2022-04-02 [1] CRAN (R 4.3.0)
# rstudioapi       0.15.0  2023-07-07 [1] CRAN (R 4.3.0)
# scales         * 1.2.1   2022-08-20 [1] CRAN (R 4.3.0)
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
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────