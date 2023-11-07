
library(here)
library(ggplot2)
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

## Plot MAF of all D variants within each population

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
allD_vars_MAF_in_pops$Group <- factor(allD_vars_MAF_in_pops$Group, levels=c(paste0('MAF_',populations), 'Allele_Frequency'))

## Number of variants with missing MAF per population
table(allD_vars_MAF_in_pops[which(is.na(allD_vars_MAF_in_pops$MAF)), 'Group'])
# MAF_African_or_African_American  MAF_Latino_or_Admixed_American                  MAF_East_Asian                 MAF_South_Asian 
#                               0                               0                               0                              75 
# MAF_European_Finnish        MAF_European_non_Finnish            MAF_Ashkenazi_Jewish                Allele_Frequency 
#                    0                               0                               0                               0 

## Manually add MAF of 2−234668879−C−CAT (UGT1A1*28) variant in South Asians (just taking Indians from )
allD_vars_MAF_in_pops[which(allD_vars_MAF_in_pops$Group=='MAF_South_Asian' & allD_vars_MAF_in_pops$Variant_ID=='2-234668879-C-CAT'), 'MAF'] <- 0.4557

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
#                  Expected number of deleterious UGT variants per individual in each population                 
# ----------------------------------------------------------------------------------------------------------------

## Number of D variants per individual in each population
num_pop <- list()

for (group in names(table(allD_vars_MAF_in_pops$Group))){
  
  ## Aggregated frequency of all D UGT variants (Qk) = number of D UGT variants per haploid genome 
  allDvars_agg <- c('all_vars'=sum(subset(allD_vars_MAF_in_pops, !is.na(MAF) & Group==group)$MAF))
  ## Total D UGT variants per diploid genome (2Qk)
  allDvars_num <- 2*allDvars_agg
  
  ## Number of unique and shared D variants in UGT genes in a diploid genome (2Qg)
  group_agg <- vector()
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
apply(num_pop[,-9], 2, function(x){c(signif(x['UGT1'])==signif(sum(x[6:15])),
                                         signif(x['UGT2'])==signif(sum(x[16:26])),
                                         signif(x['UGT3'])==signif(sum(x[27:28])),
                                         signif(x['UGT8'])==signif(sum(x[29])))})
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
num_pop_in_groups <- subset(num_pop, ! genetic_group %in% num_pop$genetic_group[1:5])
## Number od variants in each genetic and human group
num_pop_in_groups <- melt(num_pop_in_groups)
num_pop_in_groups$genetic_group <- factor(num_pop_in_groups$genetic_group, 
                                                levels=num_pop$genetic_group[-c(1:5)])
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
              

ggplot(num_pop_in_groups) +
  geom_bar(aes(x = group, y = num, fill = genetic_group), 
           position = "stack", stat = "identity") +
  geom_text(data=Qf_pop, aes(label=label, y=value, x=variable, fill=NULL),
            vjust=-0.25, size=2) +
  facet_wrap(~ gene_fam, scales="free_y", ncol=4) + 
  ## Colors by 2nd variable
  scale_fill_manual(values = genes_colors) +
  theme_bw() +
  labs(x='', y = 'Number of deleterious UGT variants', fill='Unique gene region/ Overlapping region') +
  scale_x_discrete(breaks=c(paste0('MAF_',populations), 'Allele_Frequency'),
                   labels=c("African/African American",
                            "Latino/Admixed American",
                            "East Asian",
                            "South Asian",
                            "Finnish",
                            "European non Finnish",
                            "Ashkenazi Jewish",
                            "Global")) +
  theme(axis.title = element_text(size = (7), face='bold'),
        axis.text = element_text(size = (6)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        plot.subtitle = element_text(size = 7, color = "gray40"),
        legend.title = element_text(size=9, face='bold'), 
        legend.text = element_text(size=7.5, face='bold'), 
        strip.background = element_rect(fill="gray95", size=1, color="gray60"),
        strip.text = element_text(face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# ----------------------------------------------------------------------------------------------------------------

## Plot
ggplot(data = allD_vars_MAF_in_pops, mapping = aes(x = Group, y = MAF, color = Group)) +
  geom_point(data=subset(allD_vars_MAF_in_pops, is.na(Label)), alpha = 0.65, size = 1.3, 
             position = position_jitter(width = 0.1, height = 0), color="tomato") +
  geom_point(data=subset(allD_vars_MAF_in_pops, !is.na(Label)), aes(shape=Label), size=1.5, 
             position = position_jitter(width = 0.25, height = 0), color='tomato3', stroke = 0.9) +
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
       subtitle=paste0(length(unique(allD_vars_MAF_in_pops$Variant_ID)), ' deleterious variants in total'), 
       x='', y='MAF of deleterious UGT variants in each human group', shape=paste0('Variant ID (MAF>0.01) & Gene(s)')) +
  coord_cartesian(ylim = c(0, max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)), clip = 'off') +
  ## Add aggregated frequency of all D variants per population
  geom_text(data=melt(carr_freq_pop['all_vars',]), aes(x=variable, y=max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)+0.03, shape=NULL, label=signif(value, 3)), size=2.4, color='grey20', fontface='bold') +
  
  theme(plot.title = element_text(size = (9), face='bold', vjust = 7.1, hjust=0),
        plot.subtitle = element_text(size = 8.5, color = "gray50", vjust = 7, hjust=0, face='bold'),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "cm"),
        axis.title = element_text(size = (8.5), face='bold'),
        axis.text = element_text(size = (8)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'),
        legend.title = element_text(size=8.5, face='bold'), 
        legend.text = element_text(size=8),
        panel.border = element_rect(colour = "black", fill = NA,
                                    size = 0.2))

ggsave(filename='plots/03_Anno_functional_impact/MAF_totalDvars_per_population.pdf', width = 8.5, height = 6)   


