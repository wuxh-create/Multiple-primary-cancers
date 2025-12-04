
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
  file <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/breast_colorectal/04.SNP_genotype/split_3480_extracted.raw"
  genotype_data <- fread(file, header = TRUE)
  
  survival_data <- read.table(
    "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/clinic.tsv",header = TRUE)
  
  names(genotype_data)[1] <- names(survival_data)[1]
  combined_data <- inner_join(survival_data, genotype_data)
  combined_data$osDays <- as.numeric(combined_data$osDays)
  combined_data$OSYears <- combined_data$osDays / 365
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
  
  # 创建结果数据框
  survival_results <- data.frame(
    SNP = character(0),
    AA_5yr_survival = numeric(0),
    Aa_5yr_survival = numeric(0),
    aa_5yr_survival = numeric(0),
    aa_vs_AA_reduction = numeric(0),
    aa_vs_Aa_reduction = numeric(0),
    aa_vs_combined_reduction = numeric(0),
    P_Value = numeric(0)
  )
  
  for (snp in snp_columns) {
    if (!(snp %in% names(combined_data))) {
      warning(paste("Column", snp, "not found in combined_data"))
      next
    }
    
    # 创建三种基因型分组
    combined_data <- combined_data %>% 
      mutate(genotype_group = factor(
        case_when(
          combined_data[[snp]] == 0 ~ "AA",
          combined_data[[snp]] == 1 ~ "Aa", 
          combined_data[[snp]] == 2 ~ "aa"
        ),
        levels = c("AA", "Aa", "aa")
      ))
    
    # 处理缺失值
    combined_data[[snp]] <- ifelse(is.na(combined_data[[snp]]), 1, combined_data[[snp]])
    
    # 检查基因型分组数量
    genotype_counts <- table(combined_data$genotype_group)
    if (length(genotype_counts) < 2 || min(genotype_counts) < 5) {
      print(paste0(snp, " skipped: Insufficient samples in genotype groups"))
      next
    }
    
    # 生存分析
    model1 <- survdiff(Surv(OSYears, status == "dead") ~ genotype_group,
                       data = combined_data, na.action = na.exclude)
    
    p_value <- 1 - pchisq(model1$chisq, length(model1$n) - 1)
    
    # 计算生存曲线
    fit1 <- survfit(Surv(OSYears, status == "dead") ~ genotype_group,
                    data = combined_data, na.action = na.exclude)
    
    # 提取5年生存率
    survival_summary <- summary(fit1, times = 5)
    
    # 创建基因型到生存率的映射
    genotype_survival <- data.frame(
      genotype = survival_summary$strata,
      survival_5yr = survival_summary$surv
    )
    
    # 提取各基因型的5年生存率
    AA_survival <- ifelse("genotype_group=AA" %in% genotype_survival$genotype,
                          genotype_survival$survival_5yr[genotype_survival$genotype == "genotype_group=AA"], NA)
    Aa_survival <- ifelse("genotype_group=Aa" %in% genotype_survival$genotype,
                          genotype_survival$survival_5yr[genotype_survival$genotype == "genotype_group=Aa"], NA)
    aa_survival <- ifelse("genotype_group=aa" %in% genotype_survival$genotype,
                          genotype_survival$survival_5yr[genotype_survival$genotype == "genotype_group=aa"], NA)
    
    # 计算aa相对于其他基因型的生存率降低百分比
    aa_vs_AA_reduction <- ifelse(!is.na(AA_survival) & !is.na(aa_survival),
                                 (AA_survival - aa_survival) / AA_survival * 100, NA)
    
    aa_vs_Aa_reduction <- ifelse(!is.na(Aa_survival) & !is.na(aa_survival),
                                 (Aa_survival - aa_survival) / Aa_survival * 100, NA)
    
    # 计算aa相对于AA+Aa合并组的生存率降低百分比
    combined_AA_Aa_survival <- ifelse(!is.na(AA_survival) & !is.na(Aa_survival),
                                      (AA_survival + Aa_survival) / 2, 
                                      ifelse(!is.na(AA_survival), AA_survival, Aa_survival))
    
    aa_vs_combined_reduction <- ifelse(!is.na(combined_AA_Aa_survival) & !is.na(aa_survival),
                                       (combined_AA_Aa_survival - aa_survival) / combined_AA_Aa_survival * 100, NA)
    
    # 保存结果
    survival_results <- rbind(survival_results, 
                              data.frame(
                                SNP = snp,
                                AA_5yr_survival = ifelse(is.na(AA_survival), 0, AA_survival),
                                Aa_5yr_survival = ifelse(is.na(Aa_survival), 0, Aa_survival),
                                aa_5yr_survival = ifelse(is.na(aa_survival), 0, aa_survival),
                                aa_vs_AA_reduction = ifelse(is.na(aa_vs_AA_reduction), 0, aa_vs_AA_reduction),
                                aa_vs_Aa_reduction = ifelse(is.na(aa_vs_Aa_reduction), 0, aa_vs_Aa_reduction),
                                aa_vs_combined_reduction = ifelse(is.na(aa_vs_combined_reduction), 0, aa_vs_combined_reduction),
                                P_Value = p_value
                              ))
    
    # # 生成生存曲线图
    # pdf(paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/select_pdf/",snp, "_three_genotypes_survival.pdf"),
    #     width = 6, height = 5)
    
    (ggsurvplot(fit1,
                     pval = TRUE, 
                     conf.int = FALSE,
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     ggtheme = theme_bw(),
                     palette = c("#2E9FDF", "#E7B800", "#FC4E07"), # 蓝色、黄色、红色
                     title = paste("Survival Analysis for", snp),
                     legend.title = "Genotype",
                     legend.labs = c("AA", "Aa", "aa")))
    
    # dev.off()
    
    # 打印结果
    (paste("SNP:", snp))
    (paste("AA 5-year survival:", round(ifelse(is.na(AA_survival), 0, AA_survival), 3)))
    (paste("Aa 5-year survival:", round(ifelse(is.na(Aa_survival), 0, Aa_survival), 3)))
    (paste("aa 5-year survival:", round(ifelse(is.na(aa_survival), 0, aa_survival), 3)))
    (paste("aa vs AA reduction:", round(ifelse(is.na(aa_vs_AA_reduction), 0, aa_vs_AA_reduction), 1), "%"))
    (paste("aa vs Aa reduction:", round(ifelse(is.na(aa_vs_Aa_reduction), 0, aa_vs_Aa_reduction), 1), "%"))
    (paste("aa vs combined reduction:", round(ifelse(is.na(aa_vs_combined_reduction), 0, aa_vs_combined_reduction), 1), "%"))
    print("---")
  }
  
  # 查看最终结果
  print(survival_results)
  
  # 保存结果到文件
  write.csv(survival_results, 
            "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/genotype_survival_analysis_results.csv",
            row.names = FALSE)
}