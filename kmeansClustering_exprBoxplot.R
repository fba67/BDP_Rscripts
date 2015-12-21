kmeans.boxplot <- function(x,y,bin.cnt,k){
  h.name <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac')
  km <- kmeans(x,k)
    par(mfrow=c(2,1))
    colors <- rainbow(k)
	  y.lim <- lapply(1:6,function(i)c(min(min(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)])),max(max(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))))
	  plot(km$centers[1,],col=colors[1],pch=20,type='l',ylim=c(min(min(km$centers)),max(max(km$centers))),ylab='avg. hist. abundance',xlab='dist. from TSS',xaxt='n')
	  axis(side=1,at=seq(40),labels=seq(-19,20,1))
	    for(i in 2:k){
				    points(km$centers[i,],col=colors[i],pch=20,type='l')
	    }
	    boxes <- lapply(seq(k),function(i)y[km$cluster==i])
		  boxplot(boxes,col=colors,ylab='log2 expr.')

	for(i in seq(length(h.name))){
		plot(km$centers[1,seq((i-1)*bin.cnt+1,i*bin.cnt,1)],col=colors[1],pch=20,type='l',ylim=c(min(min(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)])),max(max(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))),ylab='avg. hist. abundance',xlab='',main=h.name[i])
	for(kk in 2:k)
			points(km$centers[kk,seq((i-1)*bin.cnt+1,i*bin.cnt,1)],col=colors[kk],pch=20,type='l')
	boxes <- lapply(seq(k),function(ii)y[km$cluster==ii])
	boxplot(boxes,col=colors,ylab='log2 expr.',main='SAR')
	}
}
args <- commandArgs(trailingOnly = T)
x <- as.matrix(read.table(args[1]))
y <- as.numeric(readLines(args[2]))
y.minus <- as.numeric(readLines(args[3]))
k <- as.numeric(args[4])

bin.cnt <- 40
pdf(args[5])
kmeans.boxplot(log2(1+x),log2(1+y),bin.cnt,k)
dev.off()
