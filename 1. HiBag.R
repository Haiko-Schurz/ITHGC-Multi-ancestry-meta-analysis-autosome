source("http://bioconductor.org/biocLite.R")
BiocManager::install(c("biocLite.R"))
biocLite("HIBAG")
library(HIBAG)

#Load genetic information files in .bed, .fam, and .nim files for chromosome 6
yourgeno <- hlaBED2Geno(bed.fn=".bed", fam.fn=".fam", bim.fn=".bim")
#Load the base pair and SNP list from the genotyping chip used for the data loaded above
model.list <- get(load("InfiniumMEGA-Broad-HLA4-hg19.RData"))

#Run imputation models for the different HLA class 1 and HLA class 2 genes
hla.id <- "A"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/HLA_A.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "B"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/HLA_B.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "C"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/HLA_C.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "DRB1"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/DRB1.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "DQA1"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/DQA1.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "DQB1"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/DQB1.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

hla.id <- "DPB1"
model <- hlaModelFromObj(model.list[[hla.id]])
sink ("output/DPB1.txt")
summary(model)
summary(yourgeno)
pred.guess <- predict(model, yourgeno, type="response", match.type="Position")
print(pred.guess)
summary(pred.guess)
sink ()

