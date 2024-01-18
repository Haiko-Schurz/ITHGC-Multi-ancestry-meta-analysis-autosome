# ITHGC-Multi-ancestry-meta-analysis-autosome

These scripts correspond to the analysis done for the following publication: https://doi.org/10.7554/eLife.84394
The International Tuberculosis Host Genetics Consortium is a large scale collaboration to investigate genetic susceptibility to TB. It is the largest collection of TB GWAS data to date. 
The meta-analysis was done using GWAMA and MR-MEGA, as described in the publication. 

**HLA analysis**
The HLA analysis was done by first imputing the HLA alleles for each dataset and then doing association testing analysis for each HLA class I and class II gene. Finally the results from the association for each dataset were meta-analysed. 
1. HiBag.R: Script to impute HLA alleles
2. HiBag_assoc.R: Association testing for each HLA class I and class II gene
3. HLA_meta.R: Meta-analysis script for HLA association results across diffrent datasets.

**ITHGC_concordance_effect_direction.R**: Script used to do the concordance in direction of effects analysis between the different datasets. 

**ITHGC_forest_plot.R**: Script to produce forest plots for the SNPs of interest.

**Manhattan and QQ plot.R**: Script to produce the Manhattan and QQ-plot for the meta-analysis results. 
