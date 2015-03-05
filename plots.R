#### 
##args <- args <- commandArgs(TRUE)

tissue <- args[1]

ENpvals <- read.table("Elastic_Net/GTEx_pilot_predicted_WB_GCOM1_3980_WB_Pvals")
ENcorvec <- read.table("Elastic_Net/GTEx_pilot_predicted_WB_GCOM1_3980_WB_Corvec")
LApvals <- read.table("Lasso/GTEx_pilot_predicted_alpha1_WB_GCOM1_3981_WB_Pvals")
LAcorvec <- read.table("Lasso/GTEx_pilot_predicted_alpha1_WB_GCOM1_3981_WB_Corvec")

source("/group/im-lab/nas40t2/smozaffari/scripts/qqplot.R")

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
qqunif(ENpvals$x, main = "Elastic Net qqplot")
qqplot(nullcorvec^2,ENcorvec^2); abline(0,1); grid()
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

pdf("Lasso_EN_GTEX_WB.pdf")
plot(ENcorboth^2, LAcorboth^2, xlab = "ENcor^2", ylab = "LASSOcor^2", main =paste( "Comparing overlapping genes R^2 ", length(b), sep = ""))
qqunif(ENpvalonly, main= paste("Elastic Net only genes -", length(ENpvalonly), sep=""))
qqunif(LApvalonly, main = paste("LASSO only genes -", length(LApvalonly), sep = ""))
