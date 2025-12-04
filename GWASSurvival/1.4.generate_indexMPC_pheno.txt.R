rm(list=ls())
gc()

pheno<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/pheno.txt",header = T)


# 设置工作目录
setwd("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData")

# 定义子文件夹的列表
folders <- c("01.breast", "02.colorectal", "03.prostate", "04.lung", "05.melanoma",
             "06.urinaryBladder", "07.lymphoidNeoplasms", "08.cervix", "09.kidney", "10.ovary",
             "11.otherDigestive", "12.uterus", "13.esophagus", "14.stomach", "15.myeloidNeoplasms",
             "16.headAndNeck", "17.pancreas", "18.liver", "19.otherFemaleGenital", "20.tCellAndNKCellNeoplasms",
             "21.otherRespiratory", "22.gallbladderAndBiliaryTract", "23.softTissueSarcoma", "24.otherUrinaryOrgans",
             "25.thyroid", "26.smallIntestine", "27.brain", "28.anus", "29.eyeAndOrbit", "30.bone", "31.testis",
             "32.mesothelioma", "33.otherEndocrine")

# 循环遍历每个子文件夹
for (folder in folders) {
  # 构建完整的文件路径
  path <- file.path(getwd(), folder, "*.mpc.sampleID.txt")
  
  # 获取匹配的文件列表
  files <- list.files(path = folder, pattern = "*.mpc.sampleID.txt", full.names = TRUE)
  
  # 循环遍历每个文件
  for (file in files) {
    # 读取文件
    sample_data <- read.table(file, header = FALSE)
    # 设置列名为"sampleID"
    colnames(sample_data) <- "sampleID"
    
    # 与pheno数据框取交集
    intersected_data <- merge(pheno, sample_data, by = "sampleID", all.x = FALSE)
    
    # 构建输出文件名
    output_file <- paste0(sub("\\.txt$", "", basename(file)), ".pheno")
    output_path <- file.path(folder, output_file)
    
    # 保存交集后的文件
    write.table(intersected_data, file = output_path, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # 打印完成信息
    cat("Processed and saved:", output_path, "\n")
  }
}


