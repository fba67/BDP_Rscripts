# to prepare the TSS regions for overlapping with CAGE and GRO expressions in order to find the gene product associated to each transcription site
args <- commandArgs(trailingOnly=T)
ens <- read.table('/MMCI/MS/DEEP-liver/work/Preprocessing/references/gencode.v19.annotation_trimmed.gtf')
print('done reading ens file')
window <- as.numeric(args[1])

plus <- ens[which(ens[,5]=='+'),]
minus <- ens[which(ens[,5]=='-'),]
minus[,3] <- minus[,4]

TSS <- as.numeric(plus[,3])
plus[,3] <- TSS - window
plus[,4] <- TSS + window
TSS <- as.numeric(minus[,3])
minus[,3] <- TSS - window
minus[,4] <- TSS + window
write.table(plus[,c(1,3,4,2,5,6,8)],paste(args[2],'plus_TSSregions.txt',sep='/'),quote=F,row.names=F,col.names=F)
write.table(minus[,c(1,3,4,2,5,6,8)],paste(args[2],'minus_TSSregions.txt',sep='/'),quote=F,row.names=F,col.names=F)
