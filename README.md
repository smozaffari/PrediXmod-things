## Predicting GTEx - comparing to observered - scripts
1. Place all the IDs, genes, and genotypes in correct order/format
    genotypescript.R

2. Run prediction on genotypes    
    SNP2GReX.pl

3. Format & compare predicted and observed
    1. One script:
        pred_obs_formatfiles.R

    2. Two scripts:
        * Format correctly 

        formatfiles.R
        
            * which will run

                ensemblids.pl

        * and compare predicted and observed gene expression (print out correlation & pvalue tables)

        pred_obs.R

3. Read in pvalues and correlation values to make pdf plots

    plots.R 


    
###For Polygenic model:
######Run to add reference and effect allele to eQTL data to make dosage files.
######an mix and match cis and trans files, >> onto one file.

    addalleles.pl
    
