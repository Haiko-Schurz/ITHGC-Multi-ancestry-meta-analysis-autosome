if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install("HIBAG")
library(HIBAG)

#Load HLA allele data from the prevoius script HiBag.R
hla_a<-read.table("HLA_A.txt",header=T)
#Load the .fam file for the samples to be analysed 
fam<-read.table("SAC_MEGA_chr6_forHLA.fam",header=F)
#Load any covariated to be used for the logistic regression
pca<-read.table("hla_pca.eigenvec", header = F)
pca<-pca[,c(2,3,4,5)]
#Get phenotype and genetic sex and sample ID from .fam file
pheno<-fam[,c(2,5,6)]

#Make the dataframe with alleles and covariates for the analysis
pheno_cov<-merge(pheno,pca, by.x="V2", by.y="V2")
colnames(pheno_cov)<-c("ID","sex","pheno","PC1","PC2","PC3")
hla<-merge(hla_a,pheno_cov, by.x="sample.id", by.y="ID")


#Analysis without PCA correction
#Set the hla.id to the gene to be analysed
hla.id <- "A"
hla_in <- hlaAllele(hla$sample.id,H1 = hla[, paste("allele", "1", sep="")],H2 = hla[, paste("allele", "2", sep="")],locus=hla.id, assembly="hg19") 

#Do associaiton testing, multiple models avaialble 
hlaAssocTest(hla_in, case ~ 1, data=dat, model="additive")
hlaAssocTest(hla_in, case ~ h , data=dat, showOR=TRUE)
hlaAssocTest(hla_in, case ~ h, data=dat, model="additive", showOR=TRUE)
hlaAssocTest(hla_in, case ~ h, data=dat, model="recessive", showOR=TRUE)
hlaAssocTest(hla_in, case ~ h, data=dat, model="genotype")

#Make a sumamry of the hla_in dataframe and write to file
hla_in$a1<-1
hla_in$a2<-2

for (i in 1:nrow(hla_in)){
  hla_in$a1[i]<- unlist(strsplit(hla_in$allele1,":"))[i]
  hla_in$a2[i]<-unlist(strsplit(hla_in$allele2,":"))[i]
}

unlist(strsplit(hla_in$allele1,":"))[1]

dat <- data.frame(case = c(rep(0, n/2), rep(1, n/2)), y = rnorm(n),pc1 = rnorm(n))
n<-nrow(hla)
write.table(dat, file = "dat.txt", col.names = T, row.names = F, quote = F)


############################################################################################
# Analysis with PCA correction

#make dataframe
dat<-hla[,c(1,5,6,7,8,9)]
dat<-hla[,c(1,5,6,7,8,9,10,11,12,13,14)]
colnames(dat)<-c("sample.id","sex","pheno","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
cont<-which(dat$pheno==1)
dat$pheno[cont]<-0
case<-which(dat$pheno==2)
dat$pheno[case]<-1

#Run additive model
add<-hlaAssocTest(hla_in, pheno ~ h + PC1+PC2 + PC3+PC4+ sex, data=dat, model="additive", showOR=T)

#Write out results
write.table(add, file = "SAC_HLA_DRB1_add.txt", col.names = T, row.names = T, quote = F)

#Make a QQ plot of th eassociation testing results 
obspval <- sort(add$h.pval)
r<-which(obspval[]==0)
#obspval<-obspval[-r]
logobspval <- -(log10(obspval))
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- trunc(max(logobspval))+1
expmax <- trunc(max(logexppval))+1
plot(c(0,expmax), c(0,expmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l")
points(logexppval, logobspval, pch=23, cex=.4, bg="black")
abline(0,1,lwd=3)

