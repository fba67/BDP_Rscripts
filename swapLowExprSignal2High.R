args <- commandArgs(trailingOnly=T)
bin.cnt <- 40
x.plus <- as.matrix(read.table(args[1]))
x.minus <- as.matrix(read.table(args[2]))
plus <- as.numeric(readLines(args[3]))
minus <- as.numeric(readLines(args[4]))
idx <- which(plus < minus)
temp <- x.plus[idx,]
x.plus[idx,] <- x.minus[idx,]
x.minus[idx,] <- temp

temp <- plus[idx]
plus[idx] <- minus[idx]
minus[idx] <- temp

print(c(length(plus),length(minus),dim(x.plus)))
writeLines(text=as.character(plus),paste(args[3],'.swapped',sep=''))
writeLines(text=as.character(minus),paste(args[4],'.swapped',sep=''))
write.table(x.plus,paste(args[1],'.swapped',sep=''),row.names=F,col.names=F,quote=F)
write.table(x.minus,paste(args[2],'.swapped',sep=''),row.names=F,col.names=F,quote=F)
