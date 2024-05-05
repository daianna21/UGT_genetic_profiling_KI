# 1. Data processing and filtering
  * **1\.1  Obtain variants within each gene feature**: total and exonic variants reported in [gnomAD](https://gnomad.broadinstitute.org/) per gene and gene familiy were obtained. 
  
    * *1\.1\.1 Detect overlapping variants between different genes of the same family*: the shared variants between all *UGT1*s and in *UGT2A*[1-2] were identified.
      *  1\.1\.1\.1 Subset to variants present in the canonical transcripts of UGT genes: only variants present in the canonical transcripts of the genes were considered. The % of the initial variants captured in canonical gene transcripts was obtained for each gene.  
      *  1\.1\.1\.2 Find variants shared in canonical txs of the genes in a family: we obtained which variants are common in multiple UGT transcripts and verified no variants were shared between genes from different families. The % and number of variants kept in canonical gene transcripts were computed for each gene family and for the whole UGT superfamily. 

    * *1\.1\.2 Manual annotation of variants based on boundaries of canonical txs*: we annotated the location of all given variants based on the boundaries of the canonical transcripts of UGT genes. 
      *  1\.1\.2\.1  Check manual annotation of all and shared variants in each gene family: we first examined the manually obtained transcript boundaries-based locations of all and shared UGT variants in each gene family. We confirmed the position of shared variants in overlapping regions between transcript features of different genes and corroborated the presence of the shared exonic variants in all overlapping genes encompassing such variants. Total variants were quantified per gene transcriptomic feature (and per kbs comprised by such variants).
      * 1\.1\.2\.2 Evaluate [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) annotation of exonic variants in each gene: exonic variants were extracted per gene and their exonic transcript location was validated comparing against VEP annotation; numbers of exonic variants per gene canonical transcript (and per CDS kb) were calculated. 
      
    * *1\.1\.3 Variant quantification in each gene family*: total and exonic variants were quantified per gene family. 
    * *1\.1\.4 Examination of the minor allele frequency of all UGT exonic variants*: the global minor allele frequency (GMAF) of all exonic variants across all 22 UGTs was examined to study the prevalence of rare variation in this gene superfamily.
    
  * **1\.2  Integration of punctual variants of pharmacogenomic relevance**: append in gene variant datasets or confirm the presence of missing regulatory and exonic variants of interest that have been suggested/demonstrated to be deleterious in the literature. 


# 2. Gene constraint analysis
  * **2\.1 Compare expected vs observed ratios of missense to synonymous variants**: the observed ratios of missense to synonymous variants per gene were computed based on the gnomAD reported variants. The expected ratios were calculated interrogating the consequence of all potential single nucleotide variants in the CDS of each gene canonical tx. These and additional gene constraint estimates from [gnomAD](https://gnomad.broadinstitute.org/) were compared. 



# 3. Annotation of the functional impact of variants
The functional effect (deleterious/neutral) of the exonic variants across all 22 human UGT genes was annotated. 

  * **3\.1  Functional prediction of missense variants**: the functional consequence of UGT missense variants was predicted by multiple functionality prediction algorithms through ANNOVAR, as well as implementing an ADME-optimized framework and extracting the predictions from AlphaMissense.
  
    * *3\.1\.1  ANNOVAR input file preparation*: variant input files were prepared to annotate them and predict their effects running ANNOVAR. 
    * *3\.1\.2  Examination of ANNOVAR gene-based annotation output*: we confirmed all input missense variants are annotated as exonic and non-synonymous, evaluated point cases of those variants additionally annotated as something other than exonic/non-synonymous, and corroborated in which gene(s) these unique and shared UGT variants are found. 
    
    * *3\.1\.3  Apply ADME-optimized functionality prediction framework*: a prediction framework optimized for variant effect assessment in pharmacogenes developed in [Zhou, Y. et al, 2018](https://www.nature.com/articles/s41397-018-0044-2) was implemented.
    
    * *3\.1\.4  Extract AlphaMissense (AM) predictions for all UGT variants*: the AM scores and predictions for UGT missense variants were retrieved from the AM catalog. 
    
    * *3\.1\.5  Comparison and evaluation of predictive algorithms*: the predictive outputs of all methods were contrasted and their predictive performance evaluated using [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) variants as benchmark data.
    
    * *3\.1\.6  Development of a UGT-optimized prediction framework*: a new score was developed integrating the predictions of all the previous algorithms by optimizing their score thresholds to define deleteriousness for UGT missense variation. 


  * **3\.2  Annotate functional consequence of all UGT variants**: all exonic (and promoter TA *UGT1A1*) variants were functionally annotated per gene and gene family and deleterious variants were quantified for each; the global minor allele frequency (GMAF) of deleterious variants was examined. 
  
  
# 4. Population-scale analysis of deleterious UGT variants



