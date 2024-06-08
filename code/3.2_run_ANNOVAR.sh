
#!/bin/bash

##############################################################################################
##          Run ANNOVAR to annotate and functionally predict UGT missense variants
##############################################################################################

#-------------------------------------------------------------------------------
#               gnomAD v2.1.1 missense variants (Ghr37)
#-------------------------------------------------------------------------------

perl table_annovar.pl ANNOVAR_input/UGT3A2_input.txt humandb/ -buildver hg19 -out myannoUGT3A2 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout
