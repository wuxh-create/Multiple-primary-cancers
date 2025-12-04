
{
  rm(list=ls())
  gc()
  library(dplyr)
  library(survival)
  library(tidyr)
  library(data.table)
  library(survminer)
  library(ggplot2)  
}
{
  file <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/05.GWASSurvival_1e_5_SNP_Survival/01.data/step2Results.GWAS_5e-8_snp_genotype.raw"
  genotype_data <- fread(file, header = TRUE)
  
  
  {
    result <- apply(genotype_data[, 7:581], 2, function(x) {
      c(`0` = sum(x == 0),    # 元素名 = 值，必须加空格
        `1` = sum(x == 1),
        `2` = sum(x == 2))
    })
    result_df <- as.data.frame(t(result))
    
  }
  
  
  
  survival_data <- read.table(
    "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/clinic.tsv",header = TRUE)
  
  names(genotype_data)[1] <- names(survival_data)[1]
  combined_data <- inner_join(survival_data, genotype_data)
  combined_data$osDays <- as.numeric(combined_data$osDays)
  combined_data$OSYears <- combined_data$osDays / 365  # 更精确的年转换
  combined_data <- combined_data %>% filter(OSYears <= 15)
  colnames(combined_data) <- gsub("_.*", "", colnames(combined_data))
  
  snp_columns <- names(combined_data)
  excluded_columns <- c(
    "OSYears", "genotype_group", 
    "sampleID", "sex", "AgeOfEventFinal", "birthYear", 
    "population", "osDays", "status", "ageAtDiagnosis", 
    "stage", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"
  )
  
  snp_columns <- names(combined_data)[!names(combined_data) %in% excluded_columns]
  
  
  gg_colors <- c("Carries" = "red", "Non Carries" = "blue")  # 更新颜色映射
  
  # p_values_df <- data.frame(SNP = character(0), P_Value = numeric(0))
  
  for (snp in snp_columns) {
    if (!(snp %in% names(combined_data))) {
      warning(paste("Column", snp, "not found in combined_data"))
      next
    }
    combined_data <- combined_data %>% 
      mutate(genotype_group = factor(
        ifelse(combined_data[[snp]] == 2, "Non Carries","Carries"),
        levels = c("Non Carries","Carries")
      ))
    
    combined_data[[snp]] <- ifelse(is.na(combined_data[[snp]]), 1, combined_data[[snp]])
    # 检查基因型分组数量
    if (nlevels(combined_data$genotype_group) < 2) {
      print(paste0(snp, "skipped: Insufficient genotype groups after reclassification"))
      next
    }
    # 新增样本量检查
    genotype_counts <- table(combined_data$genotype_group)
    if (min(genotype_counts) == 0) {
      print(paste0(snp, "skipped: One or more genotype groups have zero samples"))
      next
    }
    # 生存分析模型
    model1 <- survdiff(Surv(OSYears, status == "dead") ~ genotype_group,
                       data = combined_data, na.action = na.exclude)
    
    
    cox_model <- coxph(Surv(OSYears, status == "dead") ~ genotype_group,
                       data = combined_data, na.action = na.exclude)
    summary(cox_model)
    
    # 计算生存曲线
    fit1 <- survfit(Surv(OSYears, status == "dead") ~ genotype_group,
                    data = combined_data, na.action = na.exclude)
    pdf(paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/select_pdf/",snp, "_pval_survival.pdf"),
        width = 5, height = 5)

    # ggsurvplot(fit1,
    #            pval = TRUE, conf.int = TRUE,
    #            risk.table = TRUE, # 添加风险表
    #            risk.table.col = "strata", # 根据分层更改风险表颜色
    #            linetype = "strata", # 根据分层更改线型
    #            surv.median.line = "hv", # 同时显示垂直和水平参考线
    #            ggtheme = theme_bw(), # 更改ggplot2的主题
    #            palette = c("#ec4e7c", "#45acc9"))#定义颜色
    
    
    ggsurvplot(fit1,
               pval = TRUE, 
               conf.int = FALSE,
               risk.table = TRUE, # 添加风险表
               risk.table.col = "strata", # 根据分层更改风险表颜色
               ggtheme = theme_bw(), # 更改ggplot2的主题
               palette = c("#ec4e7c", "#45acc9"))#定义颜色
    
    dev.off()
    pdf(paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/select_pdf/",snp, "_event_survival.pdf"),
        width = 5, height = 5)
    # ggsurvplot(fit1,
    #            conf.int = TRUE,
    #            risk.table.col = "strata",
    #            risk.table = TRUE, # 添加风险表
    #            ggtheme = theme_bw(), 
    #            palette = c("#ec4e7c", "#45acc9"),
    #            fun = "event")
    
    
    ggsurvplot(fit1,
               conf.int = FALSE,
               risk.table.col = "strata",
               risk.table = TRUE, # 添加风险表
               ggtheme = theme_bw(), 
               palette = c("#ec4e7c", "#45acc9"),
               fun = "event")
    
    dev.off()
    
    # 更新结果数据框
    # p_values_df <- rbind(p_values_df, 
    #                      data.frame(SNP = snp, 
    #                                 P_Value = p_value))
    # # 打印进度
    # print(paste(snp, " P-Value:", p_value, sep = ""))
  }
  # 查看结果
  print(p_values_df)
}






