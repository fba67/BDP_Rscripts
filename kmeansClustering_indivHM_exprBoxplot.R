#library('ggplot2')
kmeans.boxplot <- function(x,y,bin.cnt,k){
  h.name <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac')
    par(mfrow=c(2,1))
    colors <- rainbow(k)
	for(hm in seq(length(h.name))){
	  km <- kmeans(x[,seq((hm-1)*40+1,hm*40,1)],k)
#	  y.lim <- lapply(1:6,function(i)c(min(min(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)])),max(max(km$centers[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))))
	  plot(km$centers[1,],col=colors[1],pch=20,type='l',ylim=c(min(min(km$centers)),max(max(km$centers))),ylab='avg. hist. abundance',xlab='dist. from TSS',main=h.name[hm],xaxt='n')
	  axis(side=1, at=seq(40),labels=seq(-19,20,1),las=2)
	    for(i in 2:k){
				    points(km$centers[i,],col=colors[i],pch=20,type='l')
	    }
	    boxes <- lapply(seq(k),function(i)y[km$cluster==i])
		  boxplot(boxes,col=colors,ylab='log2 expr.',main='SAR')
	}
}
args <- commandArgs(trailingOnly = T)
x <- as.matrix(read.table(args[1]))
y <- as.numeric(readLines(args[2]))
y.minus <- as.numeric(readLines(args[3]))
k <- as.numeric(args[4])
print(summary(log2((1+y)/(1+y.minus))))
bin.cnt <- 40
pdf(args[5])
kmeans.boxplot(log2(1+x),log2((1+y)/(1+y.minus)),bin.cnt,k)
dev.off()
