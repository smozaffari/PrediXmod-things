#Compare pvalues and correlations with qqplots and r^2 plots
#must load qqunif already

ENpvals <- read.table("Elastic_Net/GTEx_pilot_predicted_GCOM1_3980_Pvals")
ENcorvec <- read.table("Elastic_Net/GTEx_pilot_predicted_GCOM1_3980_Corvec")
LApvals <- read.table("Lasso/GTEx_pilot_predicted_alpha1_GCOM1_3980_Pvals")
LAcorvec <- read.table("Lasso/GTEx_pilot_predicted_alpha1_GCOM1_3980_Corvec")
#read in all the pvalues and correlation tables for the two different methods.
#here looking at only one version of GCOM1 gene though there are two in the observed file

pdf("Lasso_GTEX_WB.pdf")
qqunif(LApvals$x, main = "Lasso qqplot")  #qqunif.REDO script  
m <- dim(LApvals)[1]
n <- 167
nullcorvec = tanh(rnorm(m)/sqrt(n-3))
qqplot(nullcorvec^2,LAcorvec^2); abline(0,1); grid()
dev.off()

pdf("ElasticNet_GTEX_WB.pdf")
m <- dim(ENpvals)[1]
n <- 167
nullcorvec = tanh(rnorm(m)/sqrt(n-3))
qqunif(ENpvals$x, main = "Elastic Net qqplot")    #qqunif.REDO script 
qqplot(nullcorvec^2,LAcorvec^2); abline(0,1); grid()
dev.off()

#change tables into vectors with names
LApval <- LApvals$x
names(LApval) <- rownames(LApvals)
LAcor <- LAcorvec$x
names(LAcor) <- rownames(LAcorvec)
ENpval <- ENpvals$x
names(ENpval) <- rownames(ENpvals)
ENcor <- ENcorvec$x
names(ENcor) <- rownames(ENcorvec)

#a and b are values that are in both ENcorvec and LAcorvec for Encorvex and LAcorvec respectively
a <- which(rownames(ENcorvec)%in%rownames(LAcorvec))
ENcorboth <- ENcor[a]
b <- which(rownames(LAcorvec)%in%rownames(ENcorvec))
LAcorboth <- LAcor[b]
#keep those that are in both

#keep those that are in each and don't overlap
LApvalonly <- LApval[-b]
ENpvalonly <- ENpval[-a]
lb <- length(LApvalonly)
la <- length(ENpvalonly)

pdf("Lasso_EN_GTEX_WB.pdf")
plot(LAcorboth^2, ENcorboth^2, xlab = "Lasso R^2", ylab = "Elastic Net R^2", main = paste("Comparing overlapping genes R^2 ", length(a),sep=""))

qqunif(ENpvalonly, main = paste("Elastic Net only genes -", la, sep=""))
qqunif(LApvalonly, main = paste("Lasso only genes -", lb, sep=""))
dev.off()
