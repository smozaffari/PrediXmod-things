#Sahar Mozaffari
#formats original observed expression file to have ensembl ids & to be compared to predited
####### Rscript scripts/pred_obs.R /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt /group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot

args <- commandArgs(TRUE)

original <- args[1]  #original file of observed gene expression, corrected for PCs
#/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt

# ID ENSG00000237613.2 ENSG00000186092.4 ENSG00000235249.1
#1 GTEX.PWOO       -0.16662285       0.032003635       -0.14115172
#2 GTEX.PX3G       -0.20667492       0.009758687       -0.08671254
#3 GTEX.PLZ5        0.05458398      -0.052856978        0.08435900
#4 GTEX.OXRK       -0.14251255      -0.116945977        0.15908600
#5 GTEX.QVJO        0.02352867      -0.176810495        0.04375056

observed <- args[2]
#observed <- read.table("/group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot", header = T)

origin <- read.table(original, header = T)
#origin <- read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt", header = T)

rownames(origin) <- origin$ID
origin$ID <- NULL
torigin <- t(origin)
write.table(torigin, observed, quote = F)
ensembl <- paste(observed, "ensembl", sep="_")
system(paste("perl /group/im-lab/nas40t2/smozaffari/scripts/ensemblids.pl", observed,  ensembl, sep=" "))
