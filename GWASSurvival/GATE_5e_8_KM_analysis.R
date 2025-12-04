

# GATE分析得到的SNP，进一步进行KM生存分析


rm(list=ls())
library(dplyr)
library(survival)
library(tidyr)
library(data.table)
#install.packages("survminer")
library(survminer)

file="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/newGenotype/GATE/step2Results.GWAS_5e-8_snp_genotype.raw.raw"

genotype_data <- fread(file,header = TRUE)

survival_data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/clinic.tsv",
                            header = TRUE)

names(genotype_data)[1]<-names(survival_data)[1]
# 整合两个数据集
combined_data <- inner_join(survival_data,genotype_data)
combined_data$osDays<-as.numeric(combined_data$osDays)
combined_data$OSMonths <- combined_data$osDays / 30
# 添加这行代码将月份转换为年
combined_data$OSYears <- combined_data$OSMonths / 12  
#筛选15年生存期内的
combined_data <- combined_data %>% filter(OSYears <= 15)
colnames(combined_data) <- gsub("_.*", "", colnames(combined_data))
#combined_data中的SNP列名以"rs"开头
snp_columns <- grep("^rs", names(combined_data), value = TRUE)
gg=c("AA","Aa","aa")
names(gg)=c("0","1","2")
# 初始化一个数据框来存储每个SNP的P值和SNP名称
p_values_df <- data.frame(SNP = character(0), P_Value = numeric(0))


# 对每个SNP进行KM生存分析和绘图
for (snp in snp_columns) {
  if (!(snp %in% names(combined_data))) {
    warning(paste("Column", snp, "not found in combined_data"))
    next
  }
  combined_data[[snp]] <- ifelse(is.na(combined_data[[snp]]), 1, combined_data[[snp]])
  # 计算基因型的计数
  genotype_counts <- table(combined_data[[snp]])
  genotype_groups <- unique(combined_data[[snp]])
  length(genotype_groups)
  if (length(genotype_groups) < 2) {
    # 如果只有1个组别，跳过该SNP的分析
    print(paste(snp, "skipped: Only one genotype group present"))
    next
  }
  
  # 进行生存分析
  model1 = survdiff(Surv(OSYears, status == "dead") ~ combined_data[[snp]], 
                    data = combined_data, 
                    na.action = na.exclude)
  
  # 计算P值
  p_value <- format(1 - pchisq(model1$chisq, length(levels(factor(combined_data[[snp]])))-1), 
                    digits = 3)
  
  # 只在P值小于0.05时进行绘图
  fit1 = survfit(Surv(OSYears, status == "dead") ~ combined_data[[snp]], 
                 data = combined_data, 
                 na.action = na.exclude)
  
  # 创建PDF
  pdf(paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/newGenotype/GWASSurvival_5e_8_SNP_Survival/",
             snp,".pdf"),5, 5)
  
  # 绘制生存曲线
  plot(fit1,
       xlab = "Times(year)", 
       col = c("red", "blue", "orange"),
       mark.time = TRUE,
       ylab = "Probability of survival",
       main = "")       
  
  # 添加标题
  mtext(side = 3,
        at = max(combined_data$OSYears)/2, 
        paste("MPC, ", snp, ", KM P = ", p_value, sep = ""), 
        line = 1, 
        cex = 1, 
        las = 1)
  
  # 添加图例
  legend("bottomleft", 
         col = c("red", "blue", "orange"),  
         cex = 0.8, 
         lty = 1, 
         bty = "n",
         legend = paste(as.vector(gg[attributes(as.factor(combined_data[[snp]]))$levels]), 
                        ", N=", 
                        model1$n, 
                        sep = ""))
  dev.off()
  # 将SNP名称和P值保存到数据框中
  p_values_df <- rbind(p_values_df, data.frame(SNP = snp, P_Value = p_value))
  # 打印进度
  print(paste(snp, " P-Value:", p_value, sep = ""))
}
write.table(p_values_df,
            paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/newGenotype/GWASSurvival_5e_8_SNP_Survival/z_snp_p_values"), 
            row.names = F,col.names = T,sep = "\t",quote = FALSE)
