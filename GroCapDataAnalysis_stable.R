makeSubregionBins <- function(tss,region.size,bin.size){
		options(scipen=999)
		bins.start <- sapply(seq(length(tss$start)),function(i) return(seq((tss$start[i]-(floor(region.size/2))),tss$start[i]+(floor(region.size/2)-bin.size),bin.size)))
		bins.end <- sapply(seq(length(tss$start)),function(i) return(seq(tss$start[i]-(floor(region.size/2))+bin.size,tss$start[i]+(floor(region.size/2)),bin.size)))
		chrs <- gsub('chr',' ',tss$chr)
		chrs.rep <- rep(chrs,each=floor(region.size/bin.size))
		final.res <- cbind(chrs.rep,as.vector(bins.start),as.vector(bins.end))
		return(final.res)
}

args <- commandArgs(trailingOnly = T)
tss.plus <- read.table(args[1],col.names = c('chr','start','end','peaktype','val','strand'))
tss.minus <- read.table(args[2],col.names = c('chr','start','end','peaktype','val','strand'))
tss.plus$start <- as.numeric(tss.plus$start)
tss.plus$end <- as.numeric(tss.plus$end)
tss.minus$start <- as.numeric(tss.minus$start)
tss.minus$end <- as.numeric(tss.minus$end)
print(tss.plus$peaktype[1:10])
#filter all the M1 genes out, keeping only M2
if(length(which(tss.plus$peaktype=="M1"))!=0)
  tss.plus <- tss.plus[-which(tss.plus$peaktype=="M1"),]
if(length(which(tss.minus$peaktype=="M1"))!=0)
  tss.minus <- tss.minus[-which(tss.minus$peaktype=="M1"),]
print(tss.plus[1:10,])
print(tss.minus[1:10,])
region.size <- as.numeric(args[3])
#Correct for the minus strand TSS, in which the TTS is actually the TSS
if(tss.plus$strand[1]=='-'){
  tss.start <- tss.plus$start
  tss.plus$start <- tss.plus$end
  tss.plus$end <- tss.start
}
if(tss.minus$strand[1]=='-'){
  tss.start <- tss.minus$start
  tss.minus$start <- tss.minus$end
  tss.minus$end <- tss.start
}
#discard all the nearby genes that are not more than region.size bp apart (due to the reader tool limitatin)
tss.dist <- sapply(seq((length(tss.plus$start)-1)),function(i)tss.plus$start[(i+1)]-tss.plus$start[i])
tss.plus.unfiltered <- tss.plus
tss.minus.unfiltered <- tss.minus
print(dim(tss.plus))
tss.plus <- tss.plus.unfiltered[which(tss.dist>region.size),]
tss.minus <- tss.minus.unfiltered[which(tss.dist>region.size),]
print(dim(tss.plus))
print("Done filtering the nearby genes out!")
#make the binning coordinates to provide the bed file required for] the reader tool
plus.bins <- makeSubregionBins(tss.plus,region.size,bin.size)
minus.bins <- makeSubregionBins(tss.minus,region.size,bin.size)


###============================================================================###
write.table(file=args[4],plus.bins,row.names=F,col.names=F,quote=F)
write.table(file=args[5],minus.bins,row.names=F,col.names=F,quote=F)
#write.table(tss.pairs.promoters,file=args[3],quote=F,row.names=F,col.names = F)
#diff <- tss.pairs.promoters$V3-tss.pairs.promoters$V2
