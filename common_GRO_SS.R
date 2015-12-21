args <- commandArgs(trailingOnly = T)
#args <- c('plus_TSSregions_genes_SS.txt.tssGeneBodyIntersected','minus_TSSregions_genes_SS.txt.tssGeneBodyIntersected','k562_SS_common.txt')
plus <- read.table(args[1],colClasses = 'character')
minus <- read.table(args[2],colClasses = 'character')
common <- as.numeric(readLines(args[3]))

plus.common <- plus[common,c(1,7,8)]
plus.common <- cbind(plus.common,cbind(rep('M2',length(common)),rep('0',length(common)),rep('+',length(common))),plus[common,5:6])
minus.common <- minus[common,c(1,7,8)]
minus.common <- cbind(minus.common,cbind(rep('M2',length(common)),rep('0',length(common)),rep('-',length(common))),minus[common,5:6])
idx <- as.numeric(readLines(args[3]))
common.plus <- sapply(seq(length(idx)),function(i)which(as.numeric(plus[,9])==idx[i])[1])
common.minus <- sapply(seq(length(idx)),function(i)which(as.numeric(minus[,9])==idx[i])[1])
plus.common <- plus[common.plus,c(1,7,8)]
plus.common <- cbind(plus.common,cbind(rep('M2',length(idx)),rep('0',length(idx)),rep('+',length(idx))),plus[common.plus,5:6])
minus.common <- minus[common.minus,c(1,7,8)]
minus.common <- cbind(minus.common,cbind(rep('M2',length(idx)),rep('0',length(idx)),rep('-',length(idx))),minus[common.minus,5:6])
pp <- which(plus.common[,8]=='protein_coding'&minus.common[,8]=='protein_coding')
np <- which(plus.common[,8]=='protein_coding'&minus.common[,8]!='protein_coding')idx <- as.numeric(readLines(args[3]))
write.table(row.names=F,col.names=F,quote=F,plus.common[pp,c(1:6)],paste(args[1],'_PP.common',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,minus.common[pp,c(1:6)],paste(args[2],'_PP.common',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,cbind(plus.common[pp,]),paste(args[1],'_PP.common_withID',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,cbind(minus.common[pp,]),paste(args[2],'_PP.commonwithID',sep=''),sep="\t")

write.table(row.names=F,col.names=F,quote=F,plus.common[np,c(1:6)],paste(args[1],'_NP.common',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,minus.common[np,c(1:6)],paste(args[2],'_NP.common',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,cbind(plus.common[np,]),paste(args[1],'_NP.common_withID',sep=''),sep="\t")
write.table(row.names=F,col.names=F,quote=F,cbind(minus.common[np,]),paste(args[2],'_NP.commonwithID',sep=''),sep="\t")
