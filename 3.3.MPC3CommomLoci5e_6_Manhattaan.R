library(CMplot)
library(data.table)
file_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/04.3cancersMPC/fastGWA_GLMM_final.fastGWA"
gwasResult <- fread(file_path,header = T,stringsAsFactors = F,data.table = F)
gwasResult.plot <- gwasResult[,c(2,1,3,13)]
colnames(gwasResult.plot) <- c("SNP","chr","pos","p")
setwd("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/04.3cancersMPC/")
CMplot(gwasResult.plot,plot.type = "m",LOG10 = T,threshold = c(5e-6), threshold.lty = c(2), 
       threshold.lwd = c(2), threshold.col=c("black"), signal.col=c("#B0282F"), highlight.col=c("#B0282F"), col=c("#d40f63","#e10b45","#f2c829","#9fc259","#47accc"), 
       signal.cex=c(1), signal.pch=c(19), cex = 0.7,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC3_5e_6_loci")

CMplot(gwasResult.plot,plot.type = "q",conf.int = F,box = F,threshold.lty = 2,threshold.col = "red",cex=0.6,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC3_5e_6_QQ")
