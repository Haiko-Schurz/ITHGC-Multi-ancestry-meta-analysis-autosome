install.packages("metafor")
library(metafor)

#Load results off additive model from the HiBag_assoc.R script
china1<-read.table("hla.assoc1",header=T,stringsAsFactors=FALSE)
china1$SE<-(((log10(as.numeric(china1$h.2.5._OR)))-(log10(as.numeric(china1$h.97.5._OR))))/(2*1.96))

china2<-read.table("hla.assoc2",header=T,stringsAsFactors=FALSE)
china2$SE<-(((log10(as.numeric(china2$h.2.5._OR)))-(log10(as.numeric(china2$h.97.5._OR))))/(2*1.96))

germany<-read.table("hla.assoc3",header=T,stringsAsFactors=FALSE) #convert to OR

russia<-read.table("hla.assoc4",header=T,stringsAsFactors=FALSE) #convert to OR

gambia<-read.table("hla.assoc5",header=T,stringsAsFactors=FALSE) 
gambia$SE<-(((log10(as.numeric(gambia$h.2.5._OR)))-(log10(as.numeric(gambia$h.97.5._OR))))/(2*1.96))

ghana<-read.table("hla.assoc6",header=T,stringsAsFactors=FALSE)
ghana$SE<-(((log10(as.numeric(ghana$h.2.5._OR)))-(log10(as.numeric(ghana$h.97.5._OR))))/(2*1.96))

sac_a<-read.table("hla.assoc7",header=T,stringsAsFactors=FALSE)
sac_a$SE<-(((log10(as.numeric(sac_a$h.2.5._OR)))-(log10(as.numeric(sac_a$h.97.5._OR))))/(2*1.96))

sac_m<-read.table("hla.assoc8",header=T,stringsAsFactors=FALSE)
sac_m$SE<-(((log10(as.numeric(sac_m$h.2.5._OR)))-(log10(as.numeric(sac_m$h.97.5._OR))))/(2*1.96))


#Make a set of unique alleles 
allele<-c(germany$allele, russia$allele)
allele.u<-unique(allele)
out<-rbind(germany,russia)
out<-out[c(1:length(allele.u)),c(1:7)]
colnames(out)<-c("allele","Beta","CIL","CIH","SE","p.value","n.stud")

#Get the odds ratios and standard errors for each SNP from each dataset
for (i in 1:length(allele.u)){
  c1<-which(china1$allele==allele.u[i])
  c2<-which(china2$allele==allele.u[i])
  ger<-which(germany$allele==allele.u[i])
  rus<-which(russia$allele==allele.u[i])
  gam<-which(gambia$allele==allele.u[i])
  gha<-which(ghana$allele==allele.u[i])
  sa<-which(sac_a$allele==allele.u[i])
  sm<-which(sac_m$allele==allele.u[i])
  yi<-c(china1$h.est_OR[c1],china2$h.est_OR[c2],germany$h.est[ger],russia$h.est[rus],gambia$h.est_OR[gam]
        ,ghana$h.est_OR[gha],sac_a$h.est_OR[sa],sac_m$h.est_OR[sm])
  sei<-c(china1$SE[c1],china2$SE[c2],germany$SE[ger],russia$SE[rus],gambia$SE[gam]
         ,ghana$SE[gha],sac_a$SE[sa],sac_m$SE[sm])
  print(i)
  
  #Run the meta-analysis for each HLA allele
  meta<-try(rma(yi= yi,sei = sei,method="REML"),TRUE)
  if(isTRUE(class(meta)=="try-error")) { 
    out$allele[i]<-allele.u[i]
    out$Beta[i]<-NA
    out$CIL[i]<- NA
    out$CIH[i]<- NA
    out$SE[i]<-NA
    out$p.value[i]<-NA
    out$n.stud[i]<-length(yi)
    next } 
  else { 
    
    out$allele[i]<-allele.u[i]
    out$Beta[i]<-meta$beta
    out$CIL[i]<-meta$ci.lb
    out$CIH[i]<-meta$ci.ub
    out$SE[i]<-meta$se
    out$p.value[i]<-meta$pval
    out$n.stud[i]<-length(yi) } 
  
}

#Write meta-analysis results to file
out$p.adj<-out$p.value*197
write.table(out, file = "HLA_A_meta.txt", col.names = T, row.names = F, quote = F)

#Frequency analysis
data<-read.table("add_hla-A",header=T,stringsAsFactors=FALSE)
allele.u<-unique(data$allele)
out<-data[c(1:length(allele.u)), c(1:2)]
colnames(out)<-c("allele","freq")

for(i in 1:nrow(out)){
  al<-which(data$allele==allele.u[i])
  out$allele[i]<-allele.u[i]
  out$freq[i]<-(data$X.h.[al]/data$X...[al])*100
}

write.table(out, file = "HLA_A_freq.txt", col.names = T, row.names = F, quote = F)


