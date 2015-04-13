## Predicting GTEx - comparing to observered - scripts
####1. Place all the IDs, genes, and genotypes in correct order/format

````
genotypescript.R
````

####2. Run prediction on genotypes    
````
SNP2GReX.pl
````

####3. Format & compare predicted and observed gene expression:
##### * One script (if formatted observed gene expression file doesn't exist): formats and compares in one script


````
pred_obs_formatfiles.R
````


##### * Two scripts (if formatted observed gene expression file doesn't exist)

  1st script Format observed gene expression file correctly: 
  which will run `ensemblids.pl` to convert ensembl ids to gene names
  

````
formatfiles.R
````

  2nd script: Compare predicted and observed gene expression (if formatted observed gene expression file already exists)


````
pred_obs.R
````

####4. Read in pvalues and correlation values to make pdf plots

````plots.R
````


###For Polygenic model:
######Run to add reference and effect allele to eQTL data to make dosage files.
######an mix and match cis and trans files, >> onto one file.

    `addalleles.pl`
    
