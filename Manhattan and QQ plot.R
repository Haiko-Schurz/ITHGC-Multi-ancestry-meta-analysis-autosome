#This script is an updated version of the script available with the MRMEGA software, available here: https://genomics.ut.ee/en/tools

data<-read.table("MRMEGA_output_file",stringsAsFactors=FALSE,header=TRUE,sep = "\t")

#Make QQ plot
tiff(file = "QQplot.tif",height=600,width=600, compression = "lzw")
obschi <- (data$chisq_association)
obsndf <- (data$ndf_association)
obspval <- pchisq(obschi, obsndf, lower.tail=FALSE)
obspval <- ifelse (obspval < 1e-325, 1e-325, obspval)
obspval <- sort(obspval[complete.cases(obspval)])
logobspval <- -(log10(obspval))
exppval <- c(1:length(obspval))
logexppval <- -(log10( (exppval-0.5)/length(exppval)))
obsmax <- trunc(max(logobspval))+1
expmax <- trunc(max(logexppval))+1
plot(c(0,expmax), c(0,expmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l")
points(logexppval, logobspval, pch=23, cex=.4, bg="black")
dev.off()

#Make manhattan plot
tiff(file = "manhat.plot.tif",height=600,width=800,compression = "lzw")

obschi <- (data$chisq_association)
obsndf <- (data$ndf_association)
obspval <- pchisq(obschi, obsndf, lower.tail=FALSE)
obspval <- ifelse (obspval < 1e-325, 1e-325, obspval)
chr <- (data$Chromosome)
pos <- (data$Position)
obsmax <- trunc(max(-log10(obspval[complete.cases(obspval)])))+1

sort.ind <- order(chr, pos) 
chr <- chr[sort.ind]
pos <- pos[sort.ind]
obspval <- obspval[sort.ind]

x <- 1:22
x2<- 1:22

for (i in 1:22)
{
  curchr=which(chr==i)
  x[i] <- trunc((max(pos[curchr]))/100) +100000
  x2[i] <- trunc((min(pos[curchr]))/100) -100000
}

x[1]=x[1]-x2[1]
x2[1]=0-x2[1]

for (i in 2:24)
{
  x[i] <- x[i-1]-x2[i]+x[i]
  x2[i] <- x[i-1]-x2[i]
  
}
locX = trunc(pos/100) + x2[chr]
locY = -log10(obspval)
col1=rgb(0,0,108,maxColorValue=255)
col2=rgb(100,149,237,maxColorValue=255)
col3=rgb(0,205,102,maxColorValue=255)
col4 <- ifelse (chr%%2==0, col1, col2)
curcol <- ifelse (obspval<5e-8, col3, col4) 
plot(locX,locY,pch=20,col=curcol,axes=F,ylab="-log10 p-value",xlab="",bty="n",ylim=c(0,obsmax),cex=0.8)
axis(2,las=1)
for (i in 1:22)
{
  labpos = (x[i] + x2[i]) / 2
  mtext(i,1,at=labpos,cex=0.8,line=0)
}
mtext("Chromosome",1,at=x[22]/2,cex=1,line=1)
dev.off()

