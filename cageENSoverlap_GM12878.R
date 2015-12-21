###ENS annotated divergent genes
#plus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/di_plus_mRNA_HepG2.txt')
#minus <- read.table('/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/di_minus_mRNA_HepG2.txt')
ens <- read.table('/MMCI/MS/DEEP-liver/work/Preprocessing/references/gencode.v19.annotation_trimmed.gtf')
print('done reading ens file')
ens[,3] <- as.numeric(ens[,3])
ens[,4] <- as.numeric(ens[,4])
plus <- ens[which(as.character(ens[,5])=="+"),]
print('done separating plus')
minus <- ens[which(as.character(ens[,5])=="-"),]
print('done separating minus')
temp <- minus[,3]
minus[,3] <- minus[,4]
minus[,3] <- temp

overlapping.thresh <- 50
###CAGE TSSs
#cage.plus <- read.table('/MMCI/MS/ExpRegulation/work/data/ENC_CAGE/wgEncodeRikenCageK562CellPapPlusSignalRep1.wig.noComment.sorted')
#cage.minus <- read.table('/MMCI/MS/ExpRegulation/work/data/ENC_CAGE/wgEncodeRikenCageK562CellPapMinusSignalRep1.wig.noComment.sorted')
args <- commandArgs(trailingOnly=T)
cage.plus <- read.table(args[1])
cage.minus <- read.table(args[2])
chrs <- plus[,1]
plus.cage.overlap <- sapply(seq(nrow(plus)),function(i){matched.chroms <- which(chrs[i]==as.character(cage.plus[,1]));
							intval <- findInterval(plus[i,3],cage.plus[matched.chroms,2])
							if(intval==0)return(0)
							diff.upper <- abs(cage.plus[matched.chroms[intval],2]-plus[i,3])
							diff.lower <- abs(cage.plus[matched.chroms[(intval+1)],2]-plus[i,3])
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
print('done with plus.cage.overlap')

chrs <- minus[,1]
minus.cage.overlap <- sapply(seq(nrow(minus)),function(i){matched.chroms <- which(chrs[i]==as.character(cage.minus[,1]));
							intval <- findInterval(minus[i,3],cage.minus[matched.chroms,3])
							if(intval==0)return(0)
							diff.upper <- abs(cage.minus[matched.chroms[intval],3]-minus[i,3])
							diff.lower <- abs(cage.minus[matched.chroms[(intval+1)],3]-minus[i,3])
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

diff.plus <- abs(plus[common.idx,3] - cage.plus[plus.cage.overlap[common.idx],2])
diff.minus <- abs(minus[common.idx,3] - cage.minus[minus.cage.overlap[common.idx],3])

cage.plus.normalized <- cage.plus
cage.minus.normalized <- cage.minus

cage.plus.normalized[plus.cage.overlap[common.idx],4] <- cage.plus[plus.cage.overlap[common.idx],4]#/(cage.plus[plus.cage.overlap[common.idx],3]-cage.plus[plus.cage.overlap[common.idx],2]+1)
cage.minus.normalized[minus.cage.overlap[common.idx],4] <- cage.minus[minus.cage.overlap[common.idx],4]#/(cage.minus[minus.cage.overlap[common.idx],3]-cage.minus[minus.cage.overlap[common.idx],2]+1)
#write.table(cbind(plus[common.idx,c(1,2)],rep('+',length(common.idx)),cage.plus.normalized[plus.cage.overlap[common.idx],c(1,2,4)]),file='CAGE_K562_commonInEns_plus.bed',quote=F,row.names=F,col.names = F,sep='\t')
#write.table(cbind(minus[common.idx,c(1,2)],rep('-',length(common.idx)),cage.minus.normalized[minus.cage.overlap[common.idx],c(1,3,4)]),file='CAGE_K562_commonInEns_minus.bed',quote=F,row.names=F,col.names = F,sep='\t')


write.table(plus[common.idx,c(1,2,3,4,5,6,8)],file='geneInfo_commonInEns_plus.txt',quote=F,row.names=F,col.names = F)
write.table(minus[common.idx,c(1,2,3,4,5,6,8)],file='geneInfo_commonInEns_minus.txt',quote=F,row.names=F,col.names = F)

match <- sapply(seq(nrow(plus)),function(i)which(as.character(plus[i,1])==as.character(ens[,6]))[1])
match.minus <- sapply(seq(nrow(minus)),function(i)which(as.character(minus[i,1])==as.character(ens[,6]))[1])
pp <- which(as.character(plus[common.idx,8])=='protein_coding'&as.character(minus[common.idx,8])=='protein_coding')
np <- which(as.character(plus[common.idx,8])=='protein_coding'&as.character(minus[common.idx,8])!='protein_coding')
pn <- which(as.character(plus[common.idx,8])!='protein_coding'&as.character(minus[common.idx,8])=='protein_coding')
nn <- which(as.character(plus[common.idx,8])!='protein_coding'&as.character(minus[common.idx,8])!='protein_coding')
write.table(plus[pp,],'PP_plus.txt',quote=F,row.names=F,col.names=F);write.table(minus[pp,],'PP_minus.txt',quote=F,row.names=F,col.names=F)
write.table(plus[pn,],'PN_plus.txt',quote=F,row.names=F,col.names=F);write.table(minus[pn,],'PN_minus.txt',quote=F,row.names=F,col.names=F)
write.table(plus[np,],'NP_plus.txt',quote=F,row.names=F,col.names=F);write.table(minus[np,],'NP_minus.txt',quote=F,row.names=F,col.names=F)
write.table(plus[nn,],'NN_plus.txt',quote=F,row.names=F,col.names=F);write.table(minus[nn,],'NN_minus.txt',quote=F,row.names=F,col.names=F)
