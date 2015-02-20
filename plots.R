ENpvals <- read.table("Elastic_Net/GTEx_pilot_predicted_GCOM1_3980_Pvals")
ENcorvec <- read.table("Elastic_Net/GTEx_pilot_predicted_GCOM1_3980_Corvec")
LApvals <- read.table("Lasso/GTEx_pilot_predicted_alpha1_GCOM1_3980_Pvals")
LAcorvec <- read.table("Lasso/GTEx_pilot_predicted_alpha1_GCOM1_3980_Corvec")


pdf("Lasso_GTEX_WB.pdf")
qqunif(LApvals$x, main = "Lasso qqplot") #from 
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

LApval <- LApvals$x
names(LApval) <- rownames(LApvals)
LAcor <- LAcorvec$x
names(LAcor) <- rownames(LAcorvec)
ENpval <- ENpvals$x
names(ENpval) <- rownames(ENpvals)
ENcor <- ENcorvec$x
names(ENcor) <- rownames(ENcorvec)

a <- which(rownames(ENcorvec)%in%rownames(LAcorvec))
ENcorboth <- ENcor[a]
b <- which(rownames(LAcorvec)%in%rownames(ENcorvec))
LAcorboth <- LAcor[b]

LApvalonly <- LApval[-b]
ENpvalonly <- ENpval[-a]
