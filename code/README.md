# 1. Data processing and filtering
  * **1.1  Obtain variants within each gene feature**: 
  * 




# 2. Gene constraint analysis




# 3. Annotation of the functional impact of variants
The functional effect (deleterious/neutral) of the exonic variants across all 22 human UGT genes was annotated. 

  * **3\.1  Functional prediction of missense variants**: the functional consequence of UGT missense variants was predicted by multiple functionality prediction algorithms through ANNOVAR, as well as implementing an ADME-optimized framework and extracting the predictions from AlphaMissense.
  
    * *3\.1\.1  ANNOVAR input file preparation*: variant input files were prepared to annotate them and predict their effects running ANNOVAR. 
    * *3\.1\.2  Examination of ANNOVAR gene-based annotation output*: we confirmed all input missense variants are annotated as exonic and non-synonymous, evaluated point cases of those variants additionally annotated as something other than exonic/non-synonymous, and corroborated in which gene(s) these unique and shared UGT variants are found. 
    
    * *3\.1\.3  Apply ADME-optimized functionality prediction framework*: a prediction framework optimized for variant effect assessment in pharmacogenes developed in [Zhou, Y. et al, 2018](https://www.nature.com/articles/s41397-018-0044-2) was implemented.
    
    * *3\.1\.4  Extract AlphaMissense (AM) predictions for all UGT variants*: the AM scores and predictions for UGT missense variants were retrieved from the AM catalog. 
    
    * *3\.1\.5  Comparison and evaluation of predictive algorithms*: the predictive outputs of all methods were contrasted and their predictive performance evaluated using ClinVar variants as benchmark data.
    
    * *3\.1\.6  Development of a UGT-optimized prediction framework*: a new method was developed integrating the predictions of all the previous algorithms by optimizing their score thresholds to define deleteriousness for UGT missense variation. 


  * **3\.2  Annotate functional consequence of all UGT variants**: all exonic (and promoter TA *UGT1A1*) variants were functionally annotated per gene and gene family and deleterious variants were quantified for each; the global minor allele frequencies (GMAF) of deleterious variants were examined. 
  
  
# 4. Population-scale analysis of deleterious UGT variants



