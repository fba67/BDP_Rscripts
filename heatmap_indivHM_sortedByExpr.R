library('gplots')
library('ggplot2')
hm.heatmap <- function(x,y,bin.cnt,k){
  h.name <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac')
    par(mfrow=c(2,1))
#	y.neg <- which(y<0)
#	y <- abs(y)
	y.ord <- order(y)
	print('ordered the PNR')
	for(hm in seq(length(h.name))){
	  x.hm  <- x[,seq((hm-1)*bin.cnt+1,hm*bin.cnt,1)]
#	  x.hm[y.neg,] <- x.hm[y.neg,seq(ncol(x.hm),1)]
	  bestKM <- NULL
	  bestTotss <- Inf
	  for(trial in seq(100)){
	  	km <- kmeans(x.hm,k,iter.max=2000)
	    if(bestTotss>km$totss){
				bestTotss <- km$totss
				bestKM <- km
		}
	  }
	  km <- bestKM
	  medians <- sapply(seq(k),function(i)median(y[km$cluster==i]))
	  medians.ord <- order(medians)
	  new.cluster <- km$cluster
	  for(i in seq(k)){
			  new.cluster[km$cluster==medians.ord[i]] <- i
	  }
	  km$cluster <- new.cluster
	  km$centers <- km$centers[medians.ord,]
	  heatmap.annotation.cols <- vector(mode='character',length=nrow(x))
	  cols <- rainbow(k)
	  heatmap.annotation.cols <- c(rep(cols[1],times=length(which(km$cluster==1))),
								   rep(cols[2],times=length(which(km$cluster==2))),
								   rep(cols[3],times=length(which(km$cluster==3))))

	  plot(km$centers[1,],col=cols[1],pch=20,type='l',ylim=c(min(min(km$centers)),max(max(km$centers))),lwd=5,ylab='avg. hist. abundance',xlab='dist. from TSS',main=h.name[hm],xaxt='n')
	  axis(side=1, at=seq(40),labels=seq(-19,20,1),las=2)
	    for(i in 2:k){
				    points(km$centers[i,],col=cols[i],pch=20,type='l',lwd=4)
	    }
	    boxes <- lapply(seq(k),function(i)y[km$cluster==i])
		  boxplot(boxes,col=cols,ylab='log2 expr.',main='PNR',lwd=4)
	
	  colnames(x.hm) <- paste(round(seq(-1.9,2.0,.1),2),'k',sep='')
	  colnames(x.hm)[bin.cnt/2] <- 'TSS'

	  heatmap.2(x.hm[y.ord,],dendrogram='none',main=h.name[hm],trace="none",density.info="none",labRow="",Rowv=FALSE,Colv=FALSE,RowSideColors=as.character(heatmap.annotation.cols))
	}
#	  y.lim <- lapply(1:6,function(i)c(min(min(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)])),max(max(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))))
}
args <- commandArgs(trailingOnly = T)
x <- as.matrix(read.table(args[1]))
y <- as.numeric(readLines(args[2]))
y.minus <- as.numeric(readLines(args[3]))
bin.cnt <- 40
h.name <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac')
if(F){
idx <- which(y.minus > y)
temp <- y[idx]
y[idx] <- y.minus[idx]
y.minus[idx] <- temp
inverse.idx <- NULL
for(hm in seq(length(h.name)))
	inverse.idx  <- c(inverse.idx,seq((hm)*bin.cnt,(hm-1)*bin.cnt+1,-1))
x[idx,] <- x[idx,inverse.idx]
}
k <- as.numeric(args[4])
print(summary(log2((1+y)/(1+y.minus))))
pdf(args[5])
print('calling the hm.heatmap function')
boxplot(list(N=log2(1+y.minus),P=log2(1+y)))
hm.heatmap(log2(1+x),log2((1+y)/(1+y.minus)),bin.cnt,k)
dev.off()
