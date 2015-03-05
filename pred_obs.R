####### Rscript scripts/pred_obs.R  /group/im-lab/nas40t2/smozaffari/Lasso/GTEx_pilot_predicted_alpha1 /group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot WB
####### Rscript scripts/pred_obs.R  /group/im-lab/nas40t2/smozaffari/Polygenic_score/GTEx_pilot_predicted_cis0.05_trans10e-05 /group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot WB

args <- commandArgs(TRUE)

predicted <- args[1]  #predicted file output from SNP2GReX.pl
#Lasso: /group/im-lab/nas40t2/smozaffari/Lasso/GTEx_pilot_predicted_alpha1
#Elastic Net: /group/im-lab/nas40t2/smozaffari/Elastic_Net/GTEx_pilot_predicted
#polygenic:  /group/im-lab/nas40t2/smozaffari/Polygenic_score/GTEx_pilot_predicted_cis0.05_trans10e-05

observed <- args[2]
#/group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot_ensembl

tissue <- args[3]

obs <- read.table(observed, header = T)
pred <- read.table(predicted, header = T)

rownames(pred)<- pred$gene	       #put genes in rownames
pred$gene <- NULL		       #remove column with gene names - now in rows

pred_1 <- pred[,sort(colnames(pred))]          #sort individuals in columns
pred_2 <- pred_1[sort(rownames(pred_1)),]      #sort genes in rows

obs_1 <- obs[, sort(colnames(obs))]	#sort individuals in columns
obs_2 <- obs_1[order(obs_1$Gene),]	#sort genes in rows

obs_3 <- obs_2[which(obs_2$Gene%in%rownames(pred_2)),]		#keep overlap of predicted genes in new observed table
pred_3 <- pred_2[which(rownames(pred_2)%in%obs_2$Gene),]	#keep overlap of observed genes in new predicted table

duplicates <- c()  #could be duplicate genes in observed file
duplicates <- anyDuplicated(obs_3$Gene)  #count duplicate genes in observed file

vars <- c()

#if duplicate genes exist, make two observed tables with each- will be labeled with genename
if (duplicates > 0) {
   files <- c()
   duplicates <- c( duplicates, anyDuplicated(obs_3$Gene, fromLast=T))					#count the other duplicate
   for (i in 1:length(duplicates)) {
       x <- obs_3[-duplicates[i],]								        #remove duplicate gene
       filename <- paste(observed, obs_3$Gene[duplicates[i]], i, duplicates[i], sep="_");			#new filename to output table
       print (paste("There are 2 Ensembl IDs for 1 gene: ", obs_3$Gene[duplicates[i]], sep = "")) ; 	#notify which gene is duplicated
       var_name <- paste(obs_3$Gene[duplicates[i]], i, duplicates[i],  sep="_");	       	                #new name 
       vars <- c(vars, var_name)
       files <- c(files, filename);		                                                        #put new names of table to call on later
       rownames(x) <- x$Gene						                                #assign genenames to row
       x$Gene <- NULL 								                        #remove column with gene names
       x <- x[,-c(which(!colnames(x)%in%colnames(pred_3)))]						#only include individuals in predicted file
       write.table(x, filename, quote = F, row.names = T)			                        #write table to have for later
    }
} else {
  rownames(obs_3) <- obs_3$Gene
  obs_3$Gene <- NULL
  obs_4 <- obs_3[,-c(which(!colnames(obs_3)%in%colnames(pred_3)))]
    new_name <- paste(observed, "unique", sep = "_")
    write.table(obs_4, new_name , row.names = T, quote = F)
}

pred_4 <- pred_3[,-c(which(!colnames(pred_3)%in%colnames(obs_3)))]
new_pname <- paste(predicted, tissue, "observed_overlap", sep="_")
write.table(pred_4, new_pname, row.names = T, quote = F)

dim(pred_4)

for (j in 1:length(files)) {
    mypvals <- c()
    mycorvec <- c()
    obs <- read.table(files[j], header = T)
    for (i in 1:dim(pred_4)[1]) {
        pred <- pred_4
	genesum <- summary(lm(as.numeric(obs[i,]) ~ as.numeric(pred[i,])))
        ctest <- cor.test(as.numeric(obs[i,]), as.numeric(pred[i,]))
        cor <- c(ctest$estimate[1])
        mycorvec <- c(mycorvec, cor)
        pval <- genesum$coefficients[8];
    	mypvals <- c(mypvals, pval)
    }
    names(mypvals) <- rownames(pred)
    names(mycorvec) <- rownames(pred)

    write.table(mypvals, paste(predicted, tissue, vars[j],  "Pvals", sep = "_"), quote = F)
    write.table(mycorvec, paste(predicted, tissue, vars[j],  "Corvec", sep = "_"), quote = F)
}
