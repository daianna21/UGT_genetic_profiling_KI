
#!/bin/bash

##############################################################################################
##          Run ANNOVAR to annotate and functionally predict UGT missense variants
##############################################################################################

#-------------------------------------------------------------------------------
#               gnomAD v2.1.1 missense variants (GRCh37)
#-------------------------------------------------------------------------------

## Delete " " from Ref and Obs columns
sed 's/'\"'/''/g' input_data/UGT1A1_input.txt > input_data/UGT1A1_input.txt
## Run ANNOVAR
perl table_annovar.pl input_data/UGT1A1_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A1 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A3_input.txt > input_data/UGT1A3_input.txt
perl table_annovar.pl input_data/UGT1A3_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A3 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A4_input.txt > input_data/UGT1A4_input.txt
perl table_annovar.pl input_data/UGT1A4_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A4 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A5_input.txt > input_data/UGT1A5_input.txt
perl table_annovar.pl input_data/UGT1A5_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A5 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A6_input.txt > input_data/UGT1A6_input.txt
perl table_annovar.pl input_data/UGT1A6_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A6 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A7_input.txt > input_data/UGT1A7_input.txt
perl table_annovar.pl input_data/UGT1A7_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A7 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A8_input.txt > input_data/UGT1A8_input.txt
perl table_annovar.pl input_data/UGT1A8_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A8 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A9_input.txt > input_data/UGT1A9_input.txt
perl table_annovar.pl input_data/UGT1A9_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A9 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT1A10_input.txt > input_data/UGT1A10_input.txt
perl table_annovar.pl input_data/UGT1A10_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT1A10 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2A1_input.txt > input_data/UGT2A1_input.txt
perl table_annovar.pl input_data/UGT2A1_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2A1 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2A2_input.txt > input_data/UGT2A2_input.txt
perl table_annovar.pl input_data/UGT2A2_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2A2 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2A3_input.txt > input_data/UGT2A3_input.txt
perl table_annovar.pl input_data/UGT2A3_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2A3 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B4_input.txt > input_data/UGT2B4_input.txt
perl table_annovar.pl input_data/UGT2B4_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B4 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B7_input.txt > input_data/UGT2B7_input.txt
perl table_annovar.pl input_data/UGT2B7_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B7 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B10_input.txt > input_data/UGT2B10_input.txt
perl table_annovar.pl input_data/UGT2B10_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B10 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B11_input.txt > input_data/UGT2B11_input.txt
perl table_annovar.pl input_data/UGT2B11_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B11 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B15_input.txt > input_data/UGT2B15_input.txt
perl table_annovar.pl input_data/UGT2B15_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B15 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B17_input.txt > input_data/UGT2B17_input.txt
perl table_annovar.pl input_data/UGT2B17_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B17 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT2B28_input.txt > input_data/UGT2B28_input.txt
perl table_annovar.pl input_data/UGT2B28_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT2B28 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT3A1_input.txt > input_data/UGT3A1_input.txt
perl table_annovar.pl input_data/UGT3A1_input.txt humandb/ -buildver hg19 -out output_data/myannoUGT3A1 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT3A2_input.txt > input_data/UGT3A2_input.txt
perl table_annovar.pl input_data/UGT3A2_input.txt humandb/ -buildver hg19 -out output_data/output_data/myannoUGT3A2 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout

sed 's/'\"'/''/g' input_data/UGT8_input.txt > input_data/UGT8_input.txt
perl table_annovar.pl input_data/UGT8_input.txt humandb/ -buildver hg19 -out output_data/output_data/myannoUGT8 -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout


#-------------------------------------------------------------------------------
#                      ClinVar missense variants (GRCh37)
#-------------------------------------------------------------------------------

sed 's/'\"'/''/g' input_data/benchmark.txt > input_data/benchmark.txt
perl table_annovar.pl input_data/benchmark.txt  humandb/ -buildver hg19 -out output_data/output_data/myanno_benchmark -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout


#-------------------------------------------------------------------------------
#                  Non-minor variants for all genes (GRCh37)
#-------------------------------------------------------------------------------

sed 's/'\"'/''/g' input_data/non_minor_alleles_allGenes_ANNOVAR_format.txt > input_data/non_minor_alleles_allGenes_ANNOVAR_format.txt
perl table_annovar.pl input_data/non_minor_alleles_allGenes_ANNOVAR_format.txt  humandb/ -buildver hg19 -out output_data/output_data/myanno_non_minor -remove -protocol refGene,dbnsfp42a,dbnsfp30a -operation g,f,f -nastring . -csvout


