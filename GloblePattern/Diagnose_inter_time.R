

#绘制患者第一次患癌和第二次患癌的时间间隔平均值和四分位图
#输入文件：
  # 所有具有诊断时间的MPC样本
  # 所有病人的出生时间
  # 所有样本的人种信息
# 输出：/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_first_second.pdf
  #     /home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_Interval_years.pdf
# 2024.12.8
# 2025.10.3.将癌症按照诊断时间间隔从高到低排序
# 武肖红


#绘制第一次患癌和第二次患癌的时间四分位图和间隔年时间图
{
  rm(list=ls())
  gc()
  getwd()
  .libPaths()
  library(dplyr)
  library(tidyr)
  library(tools)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(survminer)
  library(ggplot2)
}

# 绘制每种癌症中按照MPC比例排序的柱形图
{
  ### 总癌症数量
  library(tidyverse)
  library(randomcoloR)
  
  ### 读取数据
  cancerData <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/2.indexCancer.NumAndProportion.data", header = T, stringsAsFactors = F)
  cancerData$cancerName <- factor(cancerData$cancerName, levels = cancerData$cancerName[order(cancerData$totalNumber, decreasing = T)])
  cancerData2 <- reshape2::melt(cancerData, id = c('cancerName', "totalNumber", "proportion"))
  colnames(cancerData2) <- c("CancerClass", "totalNumber", "proportion", "Class", "SampleNumber")
  cancerData2$Class <- factor(cancerData2$Class, rev(unique(cancerData2$Class)))
  
  ### 计算比例
  cancerData2$Proportion <- (cancerData2$SampleNumber / cancerData2$totalNumber) * 100
  
  ### 根据 Class 为 mpcNumber 的 Proportion 排序横坐标
  mpc_order <- cancerData2 %>%
    filter(Class == "mpcNumber") %>%
    arrange(desc(Proportion)) %>%
    pull(CancerClass)
  
  cancerData2$CancerClass <- factor(cancerData2$CancerClass, levels = mpc_order)
  
  # 将kaposiSarcoma这一列去掉
  {
    cancerData2_select<-cancerData2[-which(cancerData2$CancerClass=="kaposiSarcoma"),]
  }
  
  ### 数量柱形图 7 x 13
  number <- ggplot(cancerData2_select, aes(x = CancerClass, y = SampleNumber, fill = Class)) +
    geom_bar(stat = 'identity', width = 0.6) +
    theme_classic() +
    ylab(label = "Number of Individuals") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.0, color = 'black', size = 14),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = 'black', size = 20, face = 'bold'),
          legend.position = "",
          panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) +
    # geom_text(aes(label = totalNumber), size = 3.5, vjust = -0.5) +
    scale_y_continuous(expand = c(0, 300), limits = c(0, 22000)) +
    scale_fill_manual(values = c("#FF8D00", "#BB1839"))
  
  number
  
  ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.1.mpc.num_new_mpc_sorted_mpc_proportion.pdf", 
         plot = number, device = "pdf", width = 12, height = 5.2)
  
  ### 比例柱形图 4.5 x 13
  Proportion <- ggplot(cancerData2_select, aes(x = CancerClass, y = Proportion, fill = Class)) +
    geom_bar(stat = 'identity', width = 0.6) +
    theme_classic() +
    ylab(label = "Proportion (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.0, color = 'black', size = 14),
          axis.text.y = element_text(color = 'black', size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = 'black', size = 20, face = 'bold'),
          legend.position = "",
          panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) +
    scale_y_continuous(expand = c(0, 2)) +
    scale_fill_manual(values = c("#FF8D00", "#BB1839"))
  
  Proportion
  
  ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.2.mpc.proportion_new_mpc_sorted_mpc_proportion.pdf", 
         plot = Proportion, device = "pdf", width = 12, height = 5.2)
  
}

{
  # 读入所有有诊断时间的MPC样本
  {
    allMPC<-read.table("/home/wuxh/03.MultiPCancer/03.predict.MPC.recurrence/all_single_multiple_cancer_recurrence",header = T,sep = "\t")
    allMPC$V1<-as.integer(allMPC$V1)
    allMPC<-allMPC[-which(allMPC$V2=="V2"),]
    allMPC<-allMPC[,1:17]
    names(allMPC)[1]<-c("sampleID")
    allMPC<-separate(allMPC,"V3",into="first_diagnose",sep = "-")
    allMPC<-separate(allMPC,"V5",into="second_diagnose",sep = "-")
    allMPC<-separate(allMPC,"V7",into="third_diagnose",sep = "/")
    allMPC<-separate(allMPC,"V9",into="fourth_diagnose",sep = "/")
    allMPC<-separate(allMPC,"V11",into="fifth_diagnose",sep = "/")
    allMPC<-separate(allMPC,"V13",into="sixth_diagnose",sep = "/")
    allMPC<-separate(allMPC,"V15",into="seventh_diagnose",sep = "/")
    allMPC<-separate(allMPC,"V17",into="eighth_diagnose",sep = "/")
    names(allMPC)<-c("sampleID","first_cancer","first_diagnose","second_cancer","second_diagnose","third_cancer","third_diagnose","fourth_cancer","fourth_diagnose",
                     "fifth_cancer","fifth_diagnose","sixth_cancer","sixth_diagnose","seventh_cancer","seventh_diagnose","eighth_cancer","eighth_diagnose")
  }
  # 读入所有病人的出生时间
  {
    clinic<-read.table("/home/wuxh/UKB50wData/01.Download/04.PhenotypeData/clinic_processed.txt",sep="\t",quote="",na.strings = c("", "NA", "N/A", "NULL"), fill=TRUE,header = F)
    colnames(clinic)<-clinic[1,]
    clinic<-clinic[-1,]
    #40005是癌症诊断日期
    #40006是ICD10编码的癌症类型
    #40008是癌症诊断年龄
    #40009是自我报告癌症发生
    #40010是自爆导致的死亡时间
    #40011是癌症肿瘤的组织学信息
    #40012癌症肿瘤的行为
    a<-names(clinic)
    #将诊断时间缺失的数据去除
    clinic_select<-clinic[,c(1,3)]
    names(clinic_select)<-c("sampleID","Year_of_birth")
  }
  #为所有有诊断时间顺序的样本添加出生日期
  {
    allMPC_clinic<-merge(clinic_select,allMPC,by="sampleID")
    # 将时间相关的列转化为number类型
    allMPC_clinic <- allMPC_clinic %>%
      mutate(across(c(2, 4, 6, 8, 10,12,14,16,18), as.numeric))
    allMPC_clinic$sampleID<-as.integer(allMPC_clinic$sampleID)
  }
  # 筛选白人样本
  {
    # 读入所有白人MPC样本
    MPC<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/multiPCancer.gcta.sample",sep = "\t")
    MPC<-MPC[,1] %>% as.data.frame()
    names(MPC)<-names(allMPC)[1]
    # 筛选所有有诊断时间的白人样本
    allMPC_clinic_white<-inner_join(allMPC_clinic,MPC)
  }
  # 计算首次诊断和第二次诊断的年龄
  {
    allMPC_clinic_white$first_diagnose_age<-allMPC_clinic_white$first_diagnose-allMPC_clinic_white$Year_of_birth
    allMPC_clinic_white$second_diagnose_age<-allMPC_clinic_white$second_diagnose-allMPC_clinic_white$Year_of_birth
  }
  #绘制第一种癌症和第二种癌症诊断的时间差柱形图
  {
    allMPC_clinic_white$Interval_years<-allMPC_clinic_white$second_diagnose_age-allMPC_clinic_white$first_diagnose_age
  }
}

# 绘制两张图的boxplot图
{
  # Interval_years的boxplot图
  {
    
    # # 将first_cancer转换为因子并按照新顺序排序
    # allMPC_clinic_white_select$first_cancer <- factor(allMPC_clinic_white_select$first_cancer, levels = mpc_order)
    # 
    # # 绘制boxplot图
    # plot_interval <- ggplot(allMPC_clinic_white_select, aes(x = first_cancer, y = Interval_years, fill = first_cancer)) +
    #   geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +  # 隐藏离群点
    #   theme_classic() +
    #   labs(x = "Cancer Type", y = "Interval Years", fill = "Cancer Type") +
    #   theme(
    #     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    #     axis.text.y = element_text(size = 12),
    #     axis.title.x = element_text(size = 14, face = "bold"),
    #     axis.title.y = element_text(size = 14, face = "bold"),
    #     legend.position = "none",  # 隐藏图例
    #     panel.border = element_rect(color = "black", fill = NA, size = 1)
    #   ) +
    #   scale_fill_manual(values = rep("#D1B2FF", length(mpc_order)))  # 单一颜色填充
    # 
    # plot_interval
    # 
    
    
    # 修改画图横坐标顺序
    {
      # 使用热图中提取的癌症类型顺序（heatmap_cancer_order）
      # 确保所有图形使用相同的顺序
      median_order <- allMPC_clinic_white_select %>%
        group_by(first_cancer) %>%
        summarise(median_interval = median(Interval_years, na.rm = TRUE)) %>%
        arrange(desc(median_interval)) %>%  # 按中位数降序排列
        pull(first_cancer)
      
      # 将median_order赋给mpc_order
      mpc_order <- median_order[c(2, 1, 3:length(median_order))]
      # 将first_cancer转换为因子并按照中位数降序排序
      allMPC_clinic_white_select$first_cancer <- factor(allMPC_clinic_white_select$first_cancer, 
                                                        levels = mpc_order)
      
      # 绘制boxplot图
      plot_interval <- ggplot(allMPC_clinic_white_select, aes(x = first_cancer, y = Interval_years, fill = first_cancer)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +  # 隐藏离群点
        theme_classic() +
        labs(x = "Cancer Type", y = "Interval Years", fill = "Cancer Type") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.position = "none",  # 隐藏图例
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
        scale_fill_manual(values = rep("#D1B2FF", length(median_order)))  # 使用新顺序的长度设置颜色
      
      plot_interval
    }
    ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_Interval_years_ABCDE_order.pdf", 
           plot = plot_interval, device = "pdf", width = 12, height = 5.2)
  }
  # 第一次诊断和第二次诊断的boxplot图
  {
    # allMPC_clinic_white_select<-allMPC_clinic_white[,c(1,3,19,20,21)]
    # 
    # allMPC_clinic_white_select$first_cancer <- factor(allMPC_clinic_white_select$first_cancer, levels = mpc_order)
    # 
    # long_data <- allMPC_clinic_white_select %>%
    #   pivot_longer(cols = c("first_diagnose_age", "second_diagnose_age"), 
    #                names_to = "diagnosis", 
    #                values_to = "age") %>%
    #   mutate(diagnosis = ifelse(diagnosis == "first_diagnose_age", "First Diagnosis", "Second Diagnosis"))
    # # 将first_cancer转换为因子并按照中位数排序
    # long_data$first_cancer <- factor(long_data$first_cancer, levels = mpc_order)
    # 
    # 
    # plot <- ggplot(long_data, aes(x = first_cancer, y = age, fill = diagnosis)) +
    #   geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8,outlier.shape = NA) +
    #   theme_classic() +
    #   labs(x = "Cancer Type", y = "Age at Diagnosis", fill = "Diagnosis") +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    #         axis.text.y = element_text(size = 12),
    #         axis.title.x = element_text(size = 14, face = "bold"),
    #         axis.title.y = element_text(size = 14, face = "bold"),
    #         legend.position = "top",
    #         panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    #   scale_fill_manual(values = c("First Diagnosis" = "#6CB4EE", "Second Diagnosis" = "#FA5D63"))
    # plot 
    # 
    
    # 修改一下画图顺序
    {
      # 加载必要的库
      library(tidyverse)
      library(lubridate)
      
      # 第一部分：热图数据准备（保持原始顺序）
      time <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.first2Cancer.date")
      time_sep <- separate(time, "V2", c("MPC1", "time1"), sep = ":")
      time_sep2 <- separate(time_sep, "V3", c("MPC2", "time2"), sep = ":")
      
      time_sep2 <- time_sep2 %>% 
        mutate(time1 = ymd(time1),
               time2 = ymd(time2)) %>%
        mutate(date_diff_days = as.numeric(time2 - time1),
               date_diff_years = date_diff_days / 365.25,
               new_column = if_else(time1 == time2, "Y", "N"))
      
      time_sep2_filter <- time_sep2 %>% filter(new_column == 'N')
      
      time_sep2_filter_mean <- time_sep2_filter %>%
        group_by(MPC1, MPC2) %>%
        summarise(avg_date_diff_years = mean(date_diff_years, na.rm = TRUE), .groups = 'drop')
      
      # 获取热图的癌症类型顺序（按首次出现顺序或字母顺序）
      heatmap_cancer_order <- unique(time_sep2_filter_mean$MPC1)  # 使用热图中MPC1的实际顺序
      # 或者按字母排序：heatmap_cancer_order <- sort(unique(time_sep2_filter_mean$MPC1))
      
      # 第二部分：箱线图（横坐标顺序与热图一致）
      allMPC_clinic_white_select <- allMPC_clinic_white[, c(1, 3, 19, 20, 21)]
      
      # 强制使用热图的癌症顺序
      allMPC_clinic_white_select$first_cancer <- factor(allMPC_clinic_white_select$first_cancer, 
                                                        levels = heatmap_cancer_order)
      
      long_data <- allMPC_clinic_white_select %>%
        pivot_longer(cols = c("first_diagnose_age", "second_diagnose_age"), 
                     names_to = "diagnosis", 
                     values_to = "age") %>%
        mutate(diagnosis = ifelse(diagnosis == "first_diagnose_age", "First Diagnosis", "Second Diagnosis"))
      
      long_data$first_cancer <- factor(long_data$first_cancer, levels = heatmap_cancer_order)
      
      # 绘制箱线图
      plot_box <- ggplot(long_data, aes(x = first_cancer, y = age, fill = diagnosis)) +
        geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8, outlier.shape = NA) +
        theme_classic() +
        labs(x = "Cancer Type", y = "Age at Diagnosis", fill = "Diagnosis") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 14, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "black", fill = NA, size = 1)) +
        scale_fill_manual(values = c("First Diagnosis" = "#6CB4EE", "Second Diagnosis" = "#FA5D63"))
    }
    ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_first_second_order_ABCDEFG.pdf", 
           plot = plot_box, device = "pdf", width = 12, height = 5.2)
    
  }
  # 这里是将展示顺序修改为按照癌症的首字母进行排序后展示
  {
    # 第一次诊断和第二次诊断的四分位图
    {
      allMPC_clinic_white_select <- allMPC_clinic_white[, c(1, 3, 19, 20, 21)]
      allMPC_clinic_white_select$first_cancer <- factor(allMPC_clinic_white_select$first_cancer, levels = mpc_order)
      # 转换为长格式数据
      long_data <- allMPC_clinic_white_select %>%
        pivot_longer(cols = c("first_diagnose_age", "second_diagnose_age"), 
                     names_to = "diagnosis", 
                     values_to = "age") %>%
        mutate(diagnosis = ifelse(diagnosis == "first_diagnose_age", "First Diagnosis", "Second Diagnosis"))
      # 确保第一癌种按字母顺序排列
      long_data$first_cancer <- factor(long_data$first_cancer, levels = mpc_order)
      # 绘制箱线图
      plot <- ggplot(long_data, aes(x = first_cancer, y = age, fill = diagnosis)) +
        geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8, outlier.shape = NA) +
        theme_classic() +
        labs(x = "Cancer Type", y = "Age at Diagnosis", fill = "Diagnosis") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 14, face = "bold"),
              legend.position = "top",
              panel.border = element_rect(color = "black", fill = NA, size = 1)) +
        scale_fill_manual(values = c("First Diagnosis" = "#6CB4EE", "Second Diagnosis" = "#FA5D63"))
      plot
      ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_first_second_order_by_cancer.pdf", 
             plot = plot, device = "pdf", width = 12, height = 5.2)
    }
    # Interval_years的boxplot图
    {
      # 绘制boxplot图
      plot_interval <- ggplot(allMPC_clinic_white_select, aes(x = first_cancer, y = Interval_years, fill = first_cancer)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +  # 隐藏离群点
        theme_classic() +
        labs(x = "Cancer Type", y = "Interval Years", fill = "Cancer Type") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.position = "none",  # 隐藏图例
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
        scale_fill_manual(values = rep("#D1B2FF", length(mpc_order)))  # 单一颜色填充
      
      plot_interval
      ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Diagnose_Interval_years_order_by_cancer.pdf", 
             plot = plot_interval, device = "pdf", width = 12, height = 5.2)
    }
    
  }
  
  
  
  
}



# 提取MPC大于等于3种的发病起始年龄，发病类型，时间间隔。
{
  allMPC_clinic_white_MPC_more_than_4<-allMPC_clinic_white[which(allMPC_clinic_white$fourth_cancer !=""),]
  nrow(allMPC_clinic_white_MPC_more_than_4)
  # allMPC_clinic_white_MPC_more_than_4的样本有66个
  allMPC_clinic_white_MPC_more_than_3<-allMPC_clinic_white[which(allMPC_clinic_white$third_cancer !=""),]
  # allMPC_clinic_white_MPC_more_than_3的样本有694个
  nrow(allMPC_clinic_white_MPC_more_than_3)
  
  write.table(allMPC_clinic_white_MPC_more_than_3, file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/MPC3_diadnose_time.txt", 
              quote = FALSE, row.names = FALSE,sep = "\t")
  
  write.table(allMPC_clinic_white_MPC_more_than_4, file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/MPC4_diadnose_time.txt", 
              quote = FALSE, row.names = FALSE,sep = "\t")
  
  
  allMPC_clinic_white_MPC_more_than_3$third_diagnose_age<-allMPC_clinic_white_MPC_more_than_3$third_diagnose-allMPC_clinic_white_MPC_more_than_3$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_3$third_second_diagnose_interval<-allMPC_clinic_white_MPC_more_than_3$third_diagnose_age-allMPC_clinic_white_MPC_more_than_3$second_diagnose_age
  
  allMPC_clinic_white_MPC_more_than_3$fourth_diagnose_age<-allMPC_clinic_white_MPC_more_than_3$fourth_diagnose-allMPC_clinic_white_MPC_more_than_3$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_3$fifth_diagnose_age<-allMPC_clinic_white_MPC_more_than_3$fifth_diagnose-allMPC_clinic_white_MPC_more_than_3$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_3$sixth_diagnose_age<-allMPC_clinic_white_MPC_more_than_3$sixth_diagnose-allMPC_clinic_white_MPC_more_than_3$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_3<-allMPC_clinic_white_MPC_more_than_3[,c(3,19,20,22,24:26)]
  
  
  
  {
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    
    long_data <- allMPC_clinic_white_MPC_more_than_3 %>% 
      pivot_longer(cols = c("first_diagnose_age", "second_diagnose_age", "third_diagnose_age", 
                            "fourth_diagnose_age", "fifth_diagnose_age", "sixth_diagnose_age"),
                   names_to = "diagnosis",
                   values_to = "age") %>%
      mutate(diagnosis = case_when(
        diagnosis == "first_diagnose_age" ~ "Diagnosis Age 1",
        diagnosis == "second_diagnose_age" ~ "Diagnosis Age 2",
        diagnosis == "third_diagnose_age" ~ "Diagnosis Age 3",
        diagnosis == "fourth_diagnose_age" ~ "Diagnosis Age 4",
        diagnosis == "fifth_diagnose_age" ~ "Diagnosis Age 5",
        diagnosis == "sixth_diagnose_age" ~ "Diagnosis Age 6"
      )) %>%
      mutate(diagnosis = factor(diagnosis, levels = c("Diagnosis Age 1", "Diagnosis Age 2", 
                                                      "Diagnosis Age 3", "Diagnosis Age 4", 
                                                      "Diagnosis Age 5", "Diagnosis Age 6")))
    
    # 指定癌症顺序
    cancer_order <- c("breast", "prostate", "colorectal", "melanoma","cervix","lung", "urinaryBladder","lymphoidNeoplasms",
                      "uterus",  "kidney", "myeloidNeoplasms", 
                      "ovary", "headAndNeck", "pancreas","esophagus","stomach","brain", "liver", "otherDigestive", "tCellAndNKCellNeoplasms", 
                      "otherFemaleGenital","thyroid","otherRespiratory","testis", "softTissueSarcoma", "eyeAndOrbit",
                      "gallbladderAndBiliaryTract", "mesothelioma","anus","smallIntestine", 
                      "otherUrinaryOrgans","bone","otherEndocrine","otherMaleGenital", 
                      "otherNervousSystem")
    long_data$first_cancer <- factor(long_data$first_cancer, levels = cancer_order)
    
    # 创建分组变量，格式为 "癌症类型 - Diagnose Age X"
    long_data <- long_data %>%
      mutate(cancer_diagnosis = paste(first_cancer, diagnosis))
    
    # 创建 boxplot，显示所有癌症类型的数据在同一个图形中
    plot <- ggplot(long_data, aes(x = first_cancer, y = age, fill = diagnosis)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8, outlier.shape = NA) +
      theme_classic() +
      labs(x = "Cancer Type", y = "Age at Diagnosis", fill = "Diagnosis") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Cancer types on x-axis
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 1)) +
      scale_fill_manual(values = c("Diagnosis Age 1" = "#F78978", 
                                   "Diagnosis Age 2" = "#7B68EE", 
                                   "Diagnosis Age 3" = "#8CD0FF", 
                                   "Diagnosis Age 4" = "#FF6EC7",
                                   "Diagnosis Age 5" = "#DAA520",
                                   "Diagnosis Age 6" = "brown")) +
      scale_x_discrete(labels = unique(long_data$first_cancer))  # Only cancer types on x-axis
    
    plot
    
    ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/allMPC_clinic_white_MPC_more_than_3_Diagnose_Age.pdf", 
           plot = plot, device = "pdf", width = 15, height = 5.2)
    
    
  }
  
  
}


# 提取MPC大于等于4种的发病起始年龄，发病类型，时间间隔。
{
  allMPC_clinic_white_MPC_more_than_4<-allMPC_clinic_white[which(allMPC_clinic_white$fourth_cancer !=""),]
  nrow(allMPC_clinic_white_MPC_more_than_4)
  # allMPC_clinic_white_MPC_more_than_4的样本有66个
  
  allMPC_clinic_white_MPC_more_than_4$third_diagnose_age<-allMPC_clinic_white_MPC_more_than_4$third_diagnose-allMPC_clinic_white_MPC_more_than_4$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_4$fourth_diagnose_age<-allMPC_clinic_white_MPC_more_than_4$fourth_diagnose-allMPC_clinic_white_MPC_more_than_4$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_4$fifth_diagnose_age<-allMPC_clinic_white_MPC_more_than_4$fifth_diagnose-allMPC_clinic_white_MPC_more_than_4$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_4$sixth_diagnose_age<-allMPC_clinic_white_MPC_more_than_4$sixth_diagnose-allMPC_clinic_white_MPC_more_than_4$Year_of_birth
  
  allMPC_clinic_white_MPC_more_than_4<-allMPC_clinic_white_MPC_more_than_4[,c(3,19,20,22:25)]
  
  
  
  {
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    
    long_data <- allMPC_clinic_white_MPC_more_than_4 %>% 
      pivot_longer(cols = c("first_diagnose_age", "second_diagnose_age", "third_diagnose_age", 
                            "fourth_diagnose_age", "fifth_diagnose_age", "sixth_diagnose_age"),
                   names_to = "diagnosis",
                   values_to = "age") %>%
      mutate(diagnosis = case_when(
        diagnosis == "first_diagnose_age" ~ "Diagnosis Age 1",
        diagnosis == "second_diagnose_age" ~ "Diagnosis Age 2",
        diagnosis == "third_diagnose_age" ~ "Diagnosis Age 3",
        diagnosis == "fourth_diagnose_age" ~ "Diagnosis Age 4",
        diagnosis == "fifth_diagnose_age" ~ "Diagnosis Age 5",
        diagnosis == "sixth_diagnose_age" ~ "Diagnosis Age 6"
      )) %>%
      mutate(diagnosis = factor(diagnosis, levels = c("Diagnosis Age 1", "Diagnosis Age 2", 
                                                      "Diagnosis Age 3", "Diagnosis Age 4", 
                                                      "Diagnosis Age 5", "Diagnosis Age 6")))
    
    # 指定癌症顺序
    cancer_order <- c("breast", "prostate", "colorectal", "melanoma","cervix","lung", "urinaryBladder","lymphoidNeoplasms",
                      "uterus",  "kidney", "myeloidNeoplasms", 
                      "ovary", "headAndNeck", "pancreas","esophagus","stomach","brain", "liver", "otherDigestive", "tCellAndNKCellNeoplasms", 
                      "otherFemaleGenital","thyroid","otherRespiratory","testis", "softTissueSarcoma", "eyeAndOrbit",
                      "gallbladderAndBiliaryTract", "mesothelioma","anus","smallIntestine", 
                      "otherUrinaryOrgans","bone","otherEndocrine","otherMaleGenital", 
                      "otherNervousSystem")
    long_data$first_cancer <- factor(long_data$first_cancer, levels = cancer_order)
    
    # 创建分组变量，格式为 "癌症类型 - Diagnose Age X"
    long_data <- long_data %>%
      mutate(cancer_diagnosis = paste(first_cancer, diagnosis))
    
    # 创建 boxplot，显示所有癌症类型的数据在同一个图形中
    plot_MPC4 <- ggplot(long_data, aes(x = first_cancer, y = age, fill = diagnosis)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8, outlier.shape = NA) +
      theme_classic() +
      labs(x = "Cancer Type", y = "Age at Diagnosis", fill = "Diagnosis") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Cancer types on x-axis
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 1)) +
      scale_fill_manual(values = c("Diagnosis Age 1" = "#F78978", 
                                   "Diagnosis Age 2" = "#7B68EE", 
                                   "Diagnosis Age 3" = "#8CD0FF", 
                                   "Diagnosis Age 4" = "#FF6EC7",
                                   "Diagnosis Age 5" = "#DAA520",
                                   "Diagnosis Age 6" = "brown")) +
      scale_x_discrete(labels = unique(long_data$first_cancer))  # Only cancer types on x-axis
    
    plot_MPC4
    
    ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/allMPC_clinic_white_MPC_more_than_4_Diagnose_Age.pdf", 
           plot = plot_MPC4, device = "pdf", width = 15, height = 5.2)
    
    
  }
  
  
}


# 分析单癌和MPC生存时间的差异boxplot图
{
  # 读取MPC样本分析
  {
    # 读入生存数据
    MPC_index_sample<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/02.indexMPCSample/merged_samples.txt")
    names(MPC_index_sample)<-c("sampleID","index_cancer")
    
    MPC_clinic<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/clinic.tsv",header = T)
    
    MPC_clinic_index_sample<-inner_join(MPC_clinic,MPC_index_sample)
    
    # index MPC
    {
      
      head(MPC_clinic_index_sample)
      table(MPC_clinic_index_sample$sex)
      # 匹配每种MPC的第一种癌症是什么
      MPC_clinic_add_index_cancer<-inner_join(MPC_clinic_index_sample,allMPC)
      head(MPC_clinic_add_index_cancer)
      MPC_clinic_add_index_cancer_select<-MPC_clinic_add_index_cancer[,c(10,2,6,7)]
      head(MPC_clinic_add_index_cancer_select)
      MPC_clinic_add_index_cancer_select$osYears<-MPC_clinic_add_index_cancer_select$osDays/365
      head(MPC_clinic_add_index_cancer_select)
      nrow(MPC_clinic_add_index_cancer_select)
      MPC_clinic_add_index_cancer_select$type<-"MPC"
      head(MPC_clinic_add_index_cancer_select)
      
      MPC_clinic_add_index_cancer_select<-MPC_clinic_add_index_cancer_select[,c(1,5,6,4)]
      head(MPC_clinic_add_index_cancer_select)
    }
  }

  
  
  # 读取单癌样本分析
  {
    singlecancer<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/03.sampleData/02.singleCancerSample/merged_file.txt",
                             header = T)
    # 获取生存时间，计算每种癌症的生存时间（天）
    {
      clinic<-read.table("/home/wuxh/UKB50wData/01.Download/04.PhenotypeData/clinic_processed.txt",sep="\t",quote="",na.strings = c("", "NA", "N/A", "NULL"), fill=TRUE,header = F)
      colnames(clinic)<-clinic[1,]
      clinic<-clinic[-1,]
      #40005是癌症诊断日期
      #40006是ICD10编码的癌症类型
      #40008是癌症诊断年龄
      # 40007是死亡年龄
      #40009是自我报告癌症发生
      #40010是自爆导致的死亡时间
      #40011是癌症肿瘤的组织学信息
      #40012癌症肿瘤的行为
      a<-names(clinic)
      #将诊断时间缺失的数据去除
      clinic_select<-clinic[,c(1,10,8,9)]
      
      non_na_count <- sum(!is.na(clinic_select[, 3]))
      
      clinic_select[, 3] <- ifelse(
        is.na(clinic_select[, 3]) & !is.na(clinic_select[, 4]),  # 条件：第三列为空且第四列不为空
        clinic_select[, 4],  # 替换为第四列的值
        clinic_select[, 3]   # 否则保持原值
      )
      
      non_na_count <- sum(!is.na(clinic_select[, 3]))
      
      clinic_select<-clinic_select[,c(1,2,3)]
      
      names(clinic_select)<-c("sampleID","Date_of_cancer_diagnosis","Date_of_death")
      
      # 如果4000这个死亡日期为空，则生存状态为alive，否则生存状态为dead。
      clinic_select <- clinic_select %>%
        mutate(os.status = ifelse(is.na(clinic_select[, 3]), "alive", "dead"))
      
      # 生存时间=（2021-11-12）-癌症诊断日期（取每个月的1号）
      
      # 修改代码，计算 osDays
      clinic_select <- clinic_select %>%
        mutate(
          # 将出生年份转为日期格式，默认加上 "01-01"
          Date_of_cancer_diagnosis_date =  as.Date(Date_of_cancer_diagnosis, format = "%Y-%m-%d"),
          # 将死亡日期转为 Date 对象，假设其格式为 YYYY-MM-DD
          Date_of_death = as.Date(Date_of_death, format = "%Y-%m-%d"),
          # 添加 osDays 列，计算生存天数
          osDays = ifelse(
            !is.na(Date_of_death),
            as.numeric(Date_of_death - Date_of_cancer_diagnosis_date),
            as.numeric(as.Date("2021-11-12") - Date_of_cancer_diagnosis_date)
          )
        )
      clinic_select<-clinic_select[-which(clinic_select$osDays < 0),]
    }

    names(singlecancer)[1]<-names(clinic_select)[1]
    clinic_select$sampleID<-as.integer(clinic_select$sampleID)
    singlecancer_add_clinic<-inner_join(clinic_select,singlecancer)
    singlecancer_add_clinic$osYears<-singlecancer_add_clinic$osDays/365
    singlecancer_add_clinic$type<-"SingleCancer"
    
    singlecancer_add_clinic<-singlecancer_add_clinic[,c(7,8,9,4)]
    
    head(singlecancer_add_clinic)   
     
    names(singlecancer_add_clinic)<-names(MPC_clinic_add_index_cancer_select)
    SingleCancer_MPC<-rbind(MPC_clinic_add_index_cancer_select,singlecancer_add_clinic)
    
    SingleCancer_MPC<-rbind(MPC_clinic_add_index_cancer_select,singlecancer_add_clinic)
    
    SingleCancer_MPC<-SingleCancer_MPC[which(SingleCancer_MPC$index_cancer != "kaposiSarcoma"),]
    
    nrow(SingleCancer_MPC)
    
    # status==0是指定0是结局即死亡，~1是看总体生存的情况。（status，0=死亡；1=生存）
    
    SingleCancer_MPC$status[SingleCancer_MPC$status=="dead"]<-0
    SingleCancer_MPC$status[SingleCancer_MPC$status=="alive"]<-1
    
    SingleCancer_MPC$status<-as.numeric(SingleCancer_MPC$status)
    
    
    SingleCancer_MPC$type[SingleCancer_MPC$type=="MPC"]<-"Multiple primary cancers"
  }
  
  # index MPC
  {
    # 将 first_cancer 列设置为 factor，并按照指定顺序排序
    SingleCancer_MPC$index_cancer <- factor(SingleCancer_MPC$index_cancer, levels = mpc_order)
    
    {
      # 计算 p 值并生成生存曲线图
      p_values <- SingleCancer_MPC %>%
        group_by(index_cancer) %>%
        summarize(
          # 创建 Surv 对象并进行生存分析
          p_value = {
            # 创建生存分析对象
            surv_obj <- Surv(osYears, status)  # osYears 为生存时间, status 为生存状态 (1 = 事件, 0 = 截尾)
            
            # 使用 log-rank 检验比较组间生存差异
            survdiff(surv_obj ~ type, data = cur_data())$p
          }
        ) %>%
        mutate(significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "NS"
        ))
      
      # 绘制生存曲线图，包含风险表
      plots <- list()  # 存储所有癌症类型的图
      
      for (cancer_type in unique(SingleCancer_MPC$index_cancer)) {
        # 筛选单一癌症数据
        data_subset <- SingleCancer_MPC %>% filter(index_cancer == cancer_type)
        # 创建生存对象
        surv_obj <- Surv(data_subset$osYears, data_subset$status)
        # KM 生存曲线拟合
        fit <- survfit(surv_obj ~ type, data = data_subset)
        # 获取横轴刻度时间点
        time_points <- pretty(data_subset$osYears)
        # 分组计算样本数
        sample_counts_single <- sapply(time_points, function(t) {
          sum(data_subset$osYears >= t & data_subset$type == "SingleCancer")
        })
        sample_counts_multiple <- sapply(time_points, function(t) {
          sum(data_subset$osYears >= t & data_subset$type == "Multiple primary cancers")
        })
        # 生成风险表的标签
        risk_table_labels <- data.frame(
          time = rep(time_points, 2),
          group = rep(c("SingleCancer", "Multiple primary cancers"), each = length(time_points)),
          N = c(sample_counts_single, sample_counts_multiple)
        )
        # 调整 ggplot2 的主题设置
        plot <- ggsurvplot(
          fit,
          data = data_subset,
          title = paste("Survival Curve for", cancer_type),  # 添加动态图形标题
          pval = TRUE,                      # 显示 P 值
          conf.int = TRUE,                  # 显示置信区间
          conf.int.style = "step",          # 使用折线样式标记置信区间（与生存曲线一致）
          conf.int.alpha = 0.8,             # 置信区间透明度
          risk.table = TRUE,                # 显示风险表
          risk.table.title = "Number at risk", # 风险表标题
          risk.table.height = 0.25,         # 调整风险表高度
          risk.table.col = "strata",        # 风险表按组着色
          xlab = "Time (Years)",            # 横轴标题
          ylab = "Survival Probability",    # 纵轴标题
          legend.title = "Group",           # 图例标题
          legend.labs = c("Multiple primary cancers", "SingleCancer"),  # 分组名称
          palette = c("#E41A1C", "#377EB8"),# 自定义两条曲线的颜色
          break.time.by = diff(time_points)[1],  # 设置刻度间隔
          ggtheme = theme_minimal()         # 图形主题
        )
        # 修改主图主题：增强横纵坐标、字体大小
        plot$plot <- plot$plot +
          theme(
            axis.title.x = element_text(size = 14),  # 横坐标轴标题字体加大加粗
            axis.title.y = element_text(size = 14),  # 纵坐标轴标题字体加大加粗
            axis.text.x = element_text(size = 14),                 # 横坐标轴刻度字体加大
            axis.text.y = element_text(size = 14),                 # 纵坐标轴刻度字体加大
            panel.grid = element_blank(),                          # 去掉背景网格线
            panel.background = element_rect(fill = "white"),        # 背景设置为白色
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
          )
        # 修改风险表：调整字体大小，去掉背景虚线
        plot$table <- plot$table +
          theme(
            axis.text.x = element_text(size = 12),                 # 风险表横坐标刻度字体加大
            # axis.text.y = element_text(size = 12),                 # 风险表纵坐标刻度字体加大
            panel.grid = element_blank(),                          # 去掉背景网格线
            panel.background = element_rect(fill = "white", color = "black") # 添加黑色边框
          )
        filename <- paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Single_index_MPC_Survival/SingleCancer_MPC_SurvivalCurveFor_",cancer_type, ".pdf")  # 设置文件名
        ggsave(filename, plot = plot$plot, device = "pdf", width = 5.6, height = 5.4)
        # 保存图形到列表中
        plots[[cancer_type]] <- plot
      }
      # 打印所有图
      for (plot in plots) {
        print(plot)
      }
    }
    
    # 绘制单癌和MPC的生存差异boxplot图
    
    plot <- ggplot(SingleCancer_MPC, aes(x = index_cancer, y = osYears, fill = type)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8, outlier.shape = NA) +
      theme_classic() +
      labs(x = "Cancer type", y = "Survival times (years)", fill = "type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 1)) +
      scale_fill_manual(values = c("SingleCancer" = "#6CB4EE", "Multiple primary cancers" = "#FA5D63"))
    
    plot2 <- plot + 
      geom_text(
        data = p_values,
        aes(
          x = index_cancer,
          y = 62,  # 显著性标注位置
          label = significance,
        ),
        inherit.aes = FALSE,
        color = "black"
      ) +
      geom_segment(
        aes(
          x = as.numeric(factor(index_cancer)) - 0.4,  # 左短线位置
          xend = as.numeric(factor(index_cancer)) - 0.4,
          y = 60,  # 短线起始位置
          yend = 61
        ),
        inherit.aes = FALSE,
        color = "black"
      ) +
      geom_segment(
        aes(
          x = as.numeric(factor(index_cancer)) + 0.4,  # 右短线位置
          xend = as.numeric(factor(index_cancer)) + 0.4,
          y = 60,  # 短线起始位置
          yend = 61
        ),
        inherit.aes = FALSE,
        color = "black"
      ) +
      geom_segment(
        aes(
          x = as.numeric(factor(index_cancer)) - 0.4,
          xend = as.numeric(factor(index_cancer)) + 0.4,
          y = 61,  # 横线位置
          yend = 61
        ),
        inherit.aes = FALSE,
        color = "black"
      )
    
    plot2
    
    ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/Single_index_MPC_Survival_time_boxplot.pdf", 
           plot = plot2, device = "pdf", width = 15, height = 7.2)
  }
  
}


# 分析MPC不同癌症男女生生存时间的差异boxplot图
{
  MPC_index_sample<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/02.indexMPCSample/merged_samples.txt")
  names(MPC_index_sample)<-c("sampleID","index_cancer")
  
  MPC_clinic<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/clinic.tsv",header = T)
  
  MPC_clinic_index_sample<-inner_join(MPC_clinic,MPC_index_sample)
  
  head(MPC_clinic_index_sample)
  nrow(MPC_clinic_index_sample)
  table(MPC_clinic_index_sample$sex)
  # 匹配每种MPC的第一种癌症是什么
  MPC_clinic_add_index_cancer<-inner_join(MPC_clinic_index_sample,allMPC)
  
  head(MPC_clinic_add_index_cancer)
  
  MPC_clinic_add_index_cancer_select<-MPC_clinic_add_index_cancer[,c(10,2,6,7)]
  head(MPC_clinic_add_index_cancer_select)
  MPC_clinic_add_index_cancer_select$osYears<-MPC_clinic_add_index_cancer_select$osDays/365
  head(MPC_clinic_add_index_cancer_select)
  nrow(MPC_clinic_add_index_cancer_select)
  # MPC_clinic_add_index_cancer_select$type<-"MPC"
  MPC_clinic_add_index_cancer_select <- MPC_clinic_add_index_cancer_select[!(MPC_clinic_add_index_cancer_select[, 1] %in% c("prostate","breast","ovary",
                                                                                                                            "testis",
                                                                                                          "cervix", "uterus","otherFemaleGenital",
                                                                                                          "otherMaleGenital","otherNervousSystem",
                                                                                                          "kaposiSarcoma")), ]

  # 要去除的癌症类型
  exclude_types <- c(
    "prostate", "breast", "ovary", "testis", "cervix", 
    "uterus", "otherFemaleGenital", "otherMaleGenital", 
    "otherNervousSystem", "kaposiSarcoma"
  )
  
  # 去除指定癌症类型
  remaining_cancer_types <- setdiff(mpc_order, exclude_types)
  
  # 转换为因子
  remaining_cancer_factors <- factor(remaining_cancer_types)
  
  
  
  MPC_clinic_add_index_cancer_select$index_cancer <- factor(MPC_clinic_add_index_cancer_select$index_cancer, levels = remaining_cancer_factors)
  head(MPC_clinic_add_index_cancer_select)
  MPC_clinic_add_index_cancer_select$status[MPC_clinic_add_index_cancer_select$status=="dead"]<-0
  MPC_clinic_add_index_cancer_select$status[MPC_clinic_add_index_cancer_select$status=="alive"]<-1
  MPC_clinic_add_index_cancer_select$status<-as.numeric(MPC_clinic_add_index_cancer_select$status)
{
  p_values <- MPC_clinic_add_index_cancer_select %>%
    group_by(index_cancer) %>%
    filter(n_distinct(sex) == 2) %>%  # 判断是否有两种性别
    filter(all(table(sex) >= 2)) %>%  # 判断每个性别的样本量是否 >= 2
    summarize(
      # 使用 Kaplan-Meier 生存分析（log-rank 检验）
      p_value = {
        # 创建 Surv 对象
        surv_obj <- Surv(osYears, status)  # osYears 为生存时间, status 为生存状态 (1 = 事件, 0 = 截尾)
        
        # log-rank 检验
        surv_diff <- survdiff(surv_obj ~ sex, data = cur_data()) # 比较性别间生存差异
        pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)  # 计算 p 值
      }
    ) %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ))
  
  plots <- list()  # 存储所有癌症类型的图
  
  for (cancer_type in unique(MPC_clinic_add_index_cancer_select$index_cancer)) {
    # 筛选当前癌症类型的数据
    data_subset <- MPC_clinic_add_index_cancer_select %>% filter(index_cancer == cancer_type)
    
    # 创建生存对象
    surv_obj <- Surv(data_subset$osYears, data_subset$status)
    
    # KM 生存曲线拟合
    fit <- survfit(surv_obj ~ sex, data = data_subset)
    
    # 获取横轴刻度时间点
    time_points <- pretty(data_subset$osYears)
    
    # 绘制生存曲线
    plot <- ggsurvplot(
      fit,
      data = data_subset,
      title = paste("Survival Curve for", cancer_type),  # 添加动态图形标题
      pval = TRUE,                      # 显示 P 值
      conf.int = TRUE,                  # 显示置信区间
      conf.int.style = "step",          # 使用折线样式标记置信区间（与生存曲线一致）
      conf.int.alpha = 0.8,             # 置信区间透明度
      risk.table = TRUE,                # 显示风险表
      risk.table.title = "Number at risk", # 风险表标题
      risk.table.height = 0.25,         # 调整风险表高度
      risk.table.col = "strata",        # 风险表按组着色
      xlab = "Time (Years)",            # 横轴标题
      ylab = "Survival Probability",    # 纵轴标题
      legend.title = "Sex",             # 图例标题
      legend.labs = c("Female", "Male"),  # 分组名称
      palette = c("#E41A1C", "#377EB8"),# 自定义两条曲线的颜色
      break.time.by = diff(time_points)[1],  # 设置刻度间隔
      ggtheme = theme_minimal()         # 图形主题
    )
    
    # 修改主图主题：增强横纵坐标、字体大小
    plot$plot <- plot$plot +
      theme(
        axis.title.x = element_text(size = 14),  # 横坐标轴标题字体加大加粗
        axis.title.y = element_text(size = 14),  # 纵坐标轴标题字体加大加粗
        axis.text.x = element_text(size = 14),   # 横坐标轴刻度字体加大
        axis.text.y = element_text(size = 14),   # 纵坐标轴刻度字体加大
        panel.grid = element_blank(),            # 去掉背景网格线
        panel.background = element_rect(fill = "white"), # 背景设置为白色
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    # 修改风险表：调整字体大小，去掉背景虚线
    plot$table <- plot$table +
      theme(
        axis.text.x = element_text(size = 12),   # 风险表横坐标刻度字体加大
        panel.grid = element_blank(),            # 去掉背景网格线
        panel.background = element_rect(fill = "white", color = "black") # 添加黑色边框
      )
    
    plot
    
    
    # 保存为 PDF 文件
    filename <- paste0("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/MPC_Sex/MPC_Sex_SurvivalCurve_", cancer_type, ".pdf")  # 文件名
    ggsave(filename, plot = plot$plot, device = "pdf", width = 5.6, height = 5.4)
    
    # 保存图形到列表中
    plots[[cancer_type]] <- plot
  }
  
}  
  plot <- ggplot(MPC_clinic_add_index_cancer_select, aes(x = index_cancer, y = osYears, fill = sex)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8,outlier.shape = NA) +
    theme_classic() +
    labs(x = "Cancer type", y = "Survival times (years)", fill = "sex") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.position = "top",
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    scale_fill_manual(values = c("male" = "#6CB4EE", "female" = "#FA5D63"))
  
  plot2<-plot+ geom_text(
    data = p_values,
    aes(
      x = index_cancer,
      y = 30,  # 显著性标注位置
      label = significance
    ),
    inherit.aes = FALSE,
    color = "black"
  ) + geom_segment(
    aes(
      x = as.numeric(factor(index_cancer)) - 0.4,  # 左短线位置
      xend = as.numeric(factor(index_cancer)) - 0.4,
      y = 28,  # 短线起始位置
      yend = 29
    ),
    inherit.aes = FALSE,
    color = "black"
  ) +
    # 垂直短线 - 右侧
    geom_segment(
      aes(
        x = as.numeric(factor(index_cancer)) + 0.4,  # 右短线位置
        xend = as.numeric(factor(index_cancer)) + 0.4,
        y = 28,  # 短线起始位置
        yend = 29
      ),
      inherit.aes = FALSE,
      color = "black"
    ) +
    # 横线连接两短线
    geom_segment(
      aes(
        x = as.numeric(factor(index_cancer)) - 0.4,
        xend = as.numeric(factor(index_cancer)) + 0.4,
        y = 29,  # 横线位置
        yend = 29
      ),
      inherit.aes = FALSE,
      color = "black"
    ) 
  
  plot2
  
  ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/MPC_survival_between_female_and_male.pdf", 
         plot = plot2, device = "pdf", width = 12, height = 5.2)
  
}

