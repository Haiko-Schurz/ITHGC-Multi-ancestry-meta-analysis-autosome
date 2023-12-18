install.packages("forestplot")
library(forestplot)


snp<-read.table("Data_input",stringsAsFactors=FALSE,header=TRUE,sep = "\t")
#Data structure example
#Cohort   MARKERNAME               EA  NEA   OR     OR_95L  OR_95U  EAF          N         CHROMOSOME POSITION   
#China1   rs28383206:32575167:A:G  G   A     0.53   0.34    0.82    0.09527269   1069.999  6          32575167

row_names<-snp$Cohort
#forestplot(row_names,snp$OR,snp$OR_95L,snp$OR_95U,zero = 1,cex  = 2,lineheight = "auto",xlab = "EXP(beta)=OR", title  = "rs10485155")

#Set up table for plot input
tabletext<-cbind(
  c("Cohort", c(snp$Cohort[c(1:12)]), NA, "Meta-analysis"),
  c("OR", round(c(snp$OR[c(1:12)]),digits=2), NA, round(snp$OR[13],digits=2)),
  c("95%CI-L",round(c(snp$OR_95L[c(1:12)]),digits=2), NA, round(snp$OR_95L[13], digits=2)),
  c("95%CI-H",round(c(snp$OR_95U[c(1:12)]),digits=2), NA, round(snp$OR_95U[13],digits=2)))

in.data <- 
  structure(list(
    OR  = c(NA,c(snp$OR[c(1:12)]) , NA, snp$OR[13]), 
    OR_95L = c(NA, c(snp$OR_95L[c(1:12)]), NA, snp$OR_95L[13]),
    OR_95U = c(NA, c(snp$OR_95U[c(1:12)]), NA, snp$OR_95U[13])),
    .Names = c("OR", "OR_95L", "OR_95U"), 
    row.names = c(NA, -15L), 
    class = "data.frame")

#Plot forest plot
tiff(file = "plot_file_name.tif",height=850,width=750, compression = "lzw")
forestplot(tabletext, 
           in.data,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,12),TRUE),
           graphwidth = unit(120,"mm"), 
           xlog=TRUE, lwd.xaxis = 1.5, cex = 2,
           title = "Chr 6: rs28383206 G allele (MAF=0.17)",
           txt_gp = fpTxtGp(ticks=gpar(cex=1)),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
dev.off()

