print("Reading files")
bim <- read.table("/group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim")
SNPxID <- read.table("/group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNPxID")
SNPlist <- read.table("/group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNP.list")
ID <- read.table("/group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list")

#reads in all the different files

print("Merging files")
colnames(SNPxID) <- ID$V1
#assign individuals to columns

tab <- cbind(SNPlist, SNPlist, bim$V4, bim$V5, bim$V6, "NA", SNPxID)
#combine lists for SNP2GReX.pl
#rsid rsid position refallele dosallele MAF(NA) geneexpressions

colnames(tab)[1] <- "SNP"
colnames(tab)[2] <- "SNP"
colnames(tab)[3] <- "pos"
colnames(tab)[4] <- "Ref_allele"
colnames(tab)[5] <- "Dos_allele"
colnames(tab)[6] <- "MAF"
#labels for first 6 columns added onto gene expression

spt <- split(tab, bim$V1)
#split genotype file into tables by chromosome
print("Writing Files")
lapply(names(spt), function(x){write.table(spt[[x]], file = paste("GTEx_chr", x,".dos", sep = ""), quote = F, row.names = F, col.names = F)})

#write each table into a new file by chromosome name (for SNP2GReX.pl)
write.table(ID, "GTExsamples.txt", col.names = F, row.names = F, quote = F)
#write list of individual IDs for SNP2GReX.pl

print("Gunzipping files")
system("gzip GTEx_chr*")
#gzip all the files for SNP2GReX.pl
