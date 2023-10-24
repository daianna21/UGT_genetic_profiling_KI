

## Plot allele frequencies of all D variants within each population

populations <- c('African_or_African_American',
                 'Latino_or_Admixed_American',
                 'Ashkenazi_Jewish',
                 'East_Asian',
                 'European_Finnish',
                 'European_non_Finnish',
                 'South_Asian',
                 'Other')

allD_vars_MAF_in_pops <- vector()
for (gene in UGT_genes){
  
  ## Exonic variants per gene
  exonic_data <- eval(parse_expr(paste0(gene, '_exonic_data')))
  
  ## Compute MAF within each population 
  for (p in populations){
    exonic_data[,paste0('MAF_', p)] <- exonic_data[,paste0('Allele_Count_', p)]/exonic_data[,paste0('Allele_Number_', p)]
    ## If no MAF available in the population
    nan <- which(is.nan(exonic_data[,paste0('MAF_', p)]))
    exonic_data[nan, paste0('MAF_', p)] <- NA
  }
  ## MAF of D vars in each population
  data <- exonic_data[which(exonic_data$Functional_impact=='D'),c('Variant_ID', 'Allele_Frequency', paste0('MAF_', populations))]
  data <- melt(data, id.vars = c('Variant_ID'))
  colnames(data) <- c('Variant_ID', 'Group', 'MAF')
  allD_vars_MAF_in_pops <- rbind(allD_vars_MAF_in_pops, data)
  
}

allD_vars_MAF_in_pops <- unique(allD_vars_MAF_in_pops)
allD_vars_MAF_in_pops$Group <- factor(allD_vars_MAF_in_pops$Group, levels=c(paste0('MAF_',populations), 'Allele_Frequency'))

## Number of variants with missing MAF per population
table(allD_vars_MAF_in_pops[which(!is.na(allD_vars_MAF_in_pops$MAF)), 'Group'])
# Allele_Frequency MAF_African_or_African_American  MAF_Latino_or_Admixed_American            MAF_Ashkenazi_Jewish 
#             1206                            1206                            1206                            1206 
# MAF_East_Asian            MAF_European_Finnish        MAF_European_non_Finnish                 MAF_South_Asian 
#           1206                            1206                            1206                            1131 
# MAF_Other 
#      1206 

## Label variants with MAF>0.01; ignore missing MAFs
allD_vars_MAF_in_pops$Label <- apply(allD_vars_MAF_in_pops, 1, function(x){if (is.na(x['MAF'])){NA}
  else if(as.numeric(x['MAF'])>=0.01){x['Variant_ID']} else{NA}})
## Order
allD_vars_MAF_in_pops$Label <- factor(allD_vars_MAF_in_pops$Label, levels=c("2-234668879-C-CAT", "2-234668879-C-CATAT",
                                                                            "2-234638282-G-GT", "2-234622331-GC-G",
                                                                            "2-234676872-C-T", "4-70512787-A-T",
                                                                            "4-70462042-C-T", "4-69811110-A-C", 
                                                                            "4-69693141-GT-G", "4-70078393-C-T",
                                                                            "4-69536234-G-T", "4-69512937-T-A", NA))

## Add in which genes those variants are
commonDvars_genes <- sapply(levels(allD_vars_MAF_in_pops$Label), function(x){as.vector(GMAFs_Dvars_genes[which(GMAFs_Dvars_genes$Variant_ID==x), 'gene'])})
commonDvars_genes[5] <- 'UGT1A[1-10]'
commonDvars_genes[7] <- 'UGT2A1, UGT2A2'
variant_labels_withGene <- paste0(levels(allD_vars_MAF_in_pops$Label), ' (', commonDvars_genes, ')')

## Aggregated frequencies of D variants in each population
aggregated_MAFs_per_pop <- melt(sapply(names(table(allD_vars_MAF_in_pops$Group)), function(x){sum(allD_vars_MAF_in_pops[which(!is.na(allD_vars_MAF_in_pops$MAF) & allD_vars_MAF_in_pops$Group==x), 'MAF'])}))
aggregated_MAFs_per_pop$Group <- rownames(aggregated_MAFs_per_pop)
colnames(aggregated_MAFs_per_pop)[1] <- 'aggregated_MAF'
aggregated_MAFs_per_pop$aggregated_MAF <- signif(aggregated_MAFs_per_pop$aggregated_MAF, digits=3)
aggregated_MAFs_per_pop$max_MAF <- sapply(names(table(allD_vars_MAF_in_pops$Group)), function(x){max(allD_vars_MAF_in_pops[which(!is.na(allD_vars_MAF_in_pops$MAF) &  allD_vars_MAF_in_pops$Group==x),'MAF'])})+0.01

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
                            "Ashkenazi Jewish",
                            "East Asian",
                            "European Finnish",
                            "European non Finnish",
                            "South Asian",
                            "Other",
                            "Global")) +
  labs(title='Deleterious variants in all UGT genes', 
       subtitle=paste0(length(unique(allD_vars_MAF_in_pops$Variant_ID)), ' deleterious variants in total'), 
       x='', y='MAF of deleterious variants in each group', shape=paste0('Variant ID (MAF>0.01) & Gene(s)')) +
  coord_cartesian(ylim = c(0, max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)), clip = 'off') +
  geom_text(data=aggregated_MAFs_per_pop, aes(x=Group, y=max(subset(allD_vars_MAF_in_pops, !is.na(MAF))$MAF)+0.03, shape=NULL, label=aggregated_MAF), size=2.4, color='grey20', fontface='bold') +
  
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


## Aggregated MAF of D variants per gene and gene family in each population
for (UGTgene in UGT_genes){
  gene_vars <- subset(GMAFs_Dvars_genes, gene==UGTgene)$Variant_ID
  vars_freqs <- subset(allD_vars_MAF_in_pops, Variant_ID %in% gene_vars)
  ## Carrier frequencies: p**2 + 2pq 
  vars_freqs$carr_freqs <- sapply(vars_freqs$MAF, function(p){(p**2)+(2*p*(1-p))})
}