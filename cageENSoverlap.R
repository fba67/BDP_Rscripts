###ENS annotated divergent genes
plus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/di_plus_mRNA_HepG2.txt')
minus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/di_minus_mRNA_HepG2.txt')
i
overlapping.thresh <- 200
###CAGE TSSs
#cage.plus <- read.table('/MMCI/MS/ExpRegulation/work/data/ENC_CAGE/wgEncodeRikenCageK562CellPapPlusSignalRep1.wig.noComment.sorted')
#cage.minus <- read.table('/MMCI/MS/ExpRegulation/work/data/ENC_CAGE/wgEncodeRikenCageK562CellPapMinusSignalRep1.wig.noComment.sorted')
args <- commandArgs(trailingOnly=T)
cage.plus <- read.table(args[1])
cage.minus <- read.table(args[2])
chrs <- paste('chr',as.character(plus[,4]),sep='')
plus.cage.overlap <- sapply(seq(nrow(plus)),function(i){matched.chroms <- which(chrs[i]==as.character(cage.plus[,1]));
							intval <- findInterval(plus[i,5],cage.plus[matched.chroms,2])
							if(intval==0)return(0)
							diff.upper <- abs(cage.plus[matched.chroms[intval],2]-plus[i,5])
							diff.lower <- abs(cage.plus[matched.chroms[(intval+1)],2]-plus[i,5])
							if(diff.upper>diff.lower){
                if(diff.lower>overlapping.thresh){
                   return(0)}else
                     return(matched.chroms[1]-1+intval+1)}
              else{
                if(diff.upper>overlapping.thresh){
                  return(0)}else 
                    return(matched.chroms[1]-1+intval)
              }
							})
length(which(plus.cage.overlap==0))



minus.cage.overlap <- sapply(seq(nrow(minus)),function(i){matched.chroms <- which(chrs[i]==as.character(cage.minus[,1]));
							intval <- findInterval(minus[i,5],cage.minus[matched.chroms,3])
							if(intval==0)return(0)
							diff.upper <- abs(cage.minus[matched.chroms[intval],3]-minus[i,5])
							diff.lower <- abs(cage.minus[matched.chroms[(intval+1)],3]-minus[i,5])
							if(diff.upper>diff.lower){
							  if(diff.lower>overlapping.thresh){
							    return(0)}else
							      return(matched.chroms[1]-1+intval+1)}
							else{
							  if(diff.upper>overlapping.thresh){
							    return(0)}else 
							      return(matched.chroms[1]-1+intval)
							}
})
length(which(minus.cage.overlap==0))

common.idx <- seq(nrow(plus))[which((plus.cage.overlap & minus.cage.overlap)==T)]

diff.plus <- abs(plus$V5[common.idx] - cage.plus[plus.cage.overlap[common.idx],2])
diff.minus <- abs(minus$V5[common.idx] - cage.minus[minus.cage.overlap[common.idx],3])

cage.plus.normalized <- cage.plus
cage.minus.normalized <- cage.minus

cage.plus.normalized[plus.cage.overlap[common.idx],4] <- cage.plus[plus.cage.overlap[common.idx],4]/(cage.plus[plus.cage.overlap[common.idx],3]-cage.plus[plus.cage.overlap[common.idx],2]+1)
cage.minus.normalized[minus.cage.overlap[common.idx],4] <- cage.minus[minus.cage.overlap[common.idx],4]/(cage.minus[minus.cage.overlap[common.idx],3]-cage.minus[minus.cage.overlap[common.idx],2]+1)
write.table(cbind(plus[common.idx,c(1,2)],rep('+',length(common.idx)),cage.plus.normalized[plus.cage.overlap[common.idx],c(1,2,4)]),file='CAGE_K562_commonInEns_plus.bed',quote=F,row.names=F,col.names = F,sep='\t')
write.table(cbind(minus[common.idx,c(1,2)],rep('-',length(common.idx)),cage.minus.normalized[minus.cage.overlap[common.idx],c(1,3,4)]),file='CAGE_K562_commonInEns_minus.bed',quote=F,row.names=F,col.names = F,sep='\t')


write.table(plus[common.idx,c(1,2,3,4,5)],file='geneInfo_K562_commonInEns_plus.txt',quote=F,row.names=F,col.names = F)
write.table(minus[common.idx,c(1,2,3,4,5)],file='geneInfo_K562_commonInEns_minus.txt',quote=F,row.names=F,col.names = F)
