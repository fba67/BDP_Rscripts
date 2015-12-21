#extracts the corresponding rows fro the indices provided in the idx (args[2]) file.
args <- commandArgs(trailingOnly=T)
x <- read.table(args[1])
idx <- as.numeric(readLines(args[2]))
write.table(x[idx,],args[3])
