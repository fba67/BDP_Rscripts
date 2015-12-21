args <- commandArgs(trailingOnly=T)
fullList <- read.table(args[1])
geneList <- read.table(args[2])#first column is plus, 2nd column is minus. I only intersect the fullList with plus genes, as the plus and minus are supposed to be linked, therefore the same row indices would be returned
geneList <- geneList[-1,]
idx <- sapply(seq(length(geneList[,1])),function(i)which(as.character(geneList[i,1])==as.character(fullList[,2])))
writeLines(as.character(idx),args[3])
