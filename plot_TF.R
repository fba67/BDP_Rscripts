library(ggplot2)
args <- commandArgs(trailingOnly=T)
x <- as.matrix(read.table(args[1]))
tf.names <- c('Bcl3','Thap1','Nrf1','Bclaf1m33','Cmyc','Pol3','Atf3','Rad21','Ctcf','Ccnt2','Rpc155','E2f6h50V0416102','Cfos','Sirt6','Egr1V0416101','Chd21250','Smc3ab9263','Elf1sc631','Cjun','Tal1','Fosl1','Gtf2b','Tbp','Hdac2','Hmgn3','Tf3c110','Hey1','Irf1Ifna30','Usf2','Pol24h8','Jund','Xrcc4','Pu1','Mafkab50322','bGata1','Sin3ak20','Mxi1bhlh','bGata2','Sp1','Nelfe','bTr4','Sp2','Nfe2','bYy1','Srf','Nfya','bZnf274','Taf1','Nfyb') 

#pdf('TF_plots_errBars.pdf')
pdf(args[2])
par(mfrow=c(3,2))
for(i in seq(length(tf.names))){
    se <- apply(x[,seq((i-1)*40+1,i*40,1)],2,FUN=sd)
    print(qplot(x=seq(-19,20,1),y=colMeans(x[,seq((i-1)*40+1,i*40,1)]),main=tf.names[i]) + geom_errorbar(aes(ymin=colMeans(x[,seq((i-1)*40+1,i*40,1)])-se/2, ymax=colMeans(x[,seq((i-1)*40+1,i*40,1)])+se/2),size=.3,width=.2,position=position_dodge(.9)) +xlab("bins")+ylab("avg. TF")+theme(panel.background = element_rect(fill = 'white', colour = 'red')))
}
dev.off()

