## Predicting GTEx - comparing to observered - scripts
#####Place all the IDs, genes, and genotypes in correct order/format
    genotypescript.R

#####Run prediction on genotypes    
    SNP2GReX.pl


#####Format correctly 

    formatfiles.R
will run:

    ensemblids.pl

##### and compare predicted and observed gene expression (print out correlation & pvalue tables)

    pred_obs.R

#####Read in pvalues and correlation values to make pdf plots

    plots.R 




#Can format and predict with one script: 

    pred_obs_formatfiles.R
    
###For Polygenic model:
######Run to add reference and effect allele to eQTL data to make dosage files.
######an mix and match cis and trans files, >> onto one file.

    addalleles.pl
    
