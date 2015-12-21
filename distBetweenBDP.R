args <- commandArgs(trailingOnly=T)
plus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/41/LiHe/di_plus_mRNA_LiHe_41_Hf03.txt')
minus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/41/LiHe/di_minus_mRNA_LiHe_41_Hf03.txt')
idx <- as.numeric(readLines(args[1]))
plus.tss <- as.numeric(plus[idx,5])
minus.tss <- as.numeric(minus[idx,5])

dif <- plus.tss - minus.tss
pdf(args[2])
barplot(dif,main='dist. between + and - TSS',horiz=T)
dev.off()
