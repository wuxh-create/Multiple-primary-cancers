# MPC3 ALLSites
rm(list=ls())
library(CMplot)
library(data.table)
# 定义文件夹列表
folders <- c(
  "01.breast", "02.colorectal", "03.prostate", "04.lung", "05.melanoma", "06.urinaryBladder",
  "07.lymphoidNeoplasms", "08.cervix", "09.kidney", "10.ovary", "11.otherDigestive", "12.uterus",
  "13.esophagus", "14.stomach", "15.myeloidNeoplasms", "16.headAndNeck", "17.pancreas", "18.liver",
  "19.otherFemaleGenital", "20.tCellAndNKCellNeoplasms", "21.otherRespiratory", "22.gallbladderAndBiliaryTract",
  "23.softTissueSarcoma", "24.otherUrinaryOrgans", "25.thyroid", "26.smallIntestine", "27.brain", "28.anus",
  "29.eyeAndOrbit", "30.bone", "31.testis", "32.mesothelioma", "33.otherEndocrine"
)

# 基础路径
base_path <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/03.SpecificMPC/"

# 遍历每个文件夹
for (folder in folders) {
  # 构建文件路径
  file_path <- paste0(base_path, folder, "/fastGWA_GLMM_final.fastGWA")
  
  # 读取文件
  gwasResult <- fread(file_path, header = T, stringsAsFactors = F, data.table = F)
  
  # 数据处理
  gwasResult.plot <- gwasResult[, c(2, 1, 3, 13)]
  colnames(gwasResult.plot) <- c("SNP", "chr", "pos", "p")
  
  # 设置工作目录
  setwd(paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/03.SpecificMPC/",folder))
  
  # 绘制曼哈顿图
  CMplot(gwasResult.plot, plot.type = "m", LOG10 = T, threshold = c(5e-8), threshold.lty = c(2),
         threshold.lwd = c(2), threshold.col = c("black"), signal.col = c("#B0282F"),
         highlight.col = c("#B0282F"), col = c("#ec4e7c", "#99bf67", "#f1c23f", "#45acc9", "#f2a691", "#f6ccdf"),
         signal.cex = c(1), signal.pch = c(19), cex = 0.7, cex.axis = 1.5, cex.lab = 2, file = "pdf",
         main = paste0(folder, "_5e_8_loci"))
  
  # 绘制QQ图
  CMplot(gwasResult.plot, plot.type = "q", conf.int = F, LOG10 = T, box = F, threshold.lty = 2,
         threshold.col = "red", cex = 0.6, cex.axis = 1.5, cex.lab = 2, file = "pdf",
         main = paste0(folder, "_5e_8_loci_QQ"))
}
