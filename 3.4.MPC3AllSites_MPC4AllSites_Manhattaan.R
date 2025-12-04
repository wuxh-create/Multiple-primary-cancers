# MPC3 ALLSites
rm(list=ls())
library(CMplot)
library(data.table)
file_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/04.3cancersMPCAllSite/fastGWA_GLMM_final.fastGWA"
gwasResult <- fread(file_path,header = T,stringsAsFactors = F,data.table = F)
gwasResult.plot <- gwasResult[,c(2,1,3,13)]
colnames(gwasResult.plot) <- c("SNP","chr","pos","p")
gwasResult.plot$p<-log10(gwasResult.plot$p)
head(gwasResult.plot)
setwd("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/04.3cancersMPCAllSite/")
CMplot(gwasResult.plot,plot.type = "m",LOG10 = F,threshold = c(5e-8), threshold.lty = c(2), 
       threshold.lwd = c(2), threshold.col=c("black"), signal.col=c("#B0282F"), highlight.col=c("#B0282F"), 
       col=c("#ec4e7c","#99bf67","#f1c23f","#45acc9","#f2a691","#f6ccdf"), 
       signal.cex=c(1), signal.pch=c(19), cex = 0.7,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC3_AllSites_loci")
CMplot(gwasResult.plot,plot.type = "q",conf.int = F,LOG10 = F,box = F,threshold.lty = 2,threshold.col = "red",cex=0.6,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC3_AllSites_loci")


library(CMplot)
library(data.table)
file_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/05.4cancersMPCAllSite/fastGWA_GLMM_final.fastGWA"
gwasResult <- fread(file_path,header = T,stringsAsFactors = F,data.table = F)
gwasResult.plot <- gwasResult[,c(2,1,3,13)]
colnames(gwasResult.plot) <- c("SNP","chr","pos","p")

gwasResult.plot_clean <- subset(gwasResult.plot, p > 0 & p < 1)


setwd("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/05.4cancersMPCAllSite/")
CMplot(gwasResult.plot_clean,plot.type = "m",LOG10 = T,threshold = c(5e-8), threshold.lty = c(2), 
       threshold.lwd = c(2), threshold.col=c("black"), signal.col=c("#B0282F"), highlight.col=c("#B0282F"), col=c("#d40f63","#e10b45","#f2c829","#9fc259","#47accc"), 
       signal.cex=c(1), signal.pch=c(19), cex = 0.7,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC4_AllSites_loci")
CMplot(gwasResult.plot_clean,plot.type = "q",conf.int = F,box = F,threshold.lty = 2,threshold.col = "red",cex=0.6,cex.axis=1.5,cex.lab=2,file="pdf",main = "MPC4_AllSites_loci")




p_value=gwasResult.plot$p
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
