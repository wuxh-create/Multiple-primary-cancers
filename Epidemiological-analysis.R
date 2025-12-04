




# 连续变量流行病学分析，计算wilcox P值和logFC值

###连续变量流行病学分析
#数据的提取和手工划分类型


#----------------------------------------------------------------------------part1----------------------------------------------------------------------------------------------------------------------------------------------
#载入包
{
  rm(list=ls())
  library(boot)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(table1)
  library(stringr)
  library(vcd)
}
#----------------------------------------------------------------------------part2----------------------------------------------------------------------------------------------------------------------------------------------
#读入单癌，所有癌症，MPC≥3，MPC≥4的样本
{
  single<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/08.singleCancer/singleCacner.filtered.white.sample")
  names(single)<-"sampleID"
  all_MPC<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/multiPCancer.sample")
  names(all_MPC)<-"sampleID"
  MPC2<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.2cancers.sample")
  names(MPC2)<-"sampleID"
  MPC3<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.3cancers.sample")
  names(MPC3)<-"sampleID"
  MPC4<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.4cancers.sample")
  names(MPC4)<-"sampleID"
  MPC5<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.5cancers.sample")
  names(MPC5)<-"sampleID"
  MPC6<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.6cancers.sample")
  names(MPC6)<-"sampleID"
  
  MPC_more3<-rbind(MPC3,MPC4,MPC5,MPC6)
  MPC_more4<-rbind(MPC4,MPC5,MPC6)
}
#---------------------------------------------------------part2 描述性特征分析--------------------------------------
#----------------------------------------------------------------------------part3----------------------------------------------------------------------------------------------------------------------------------------------
#读入不同field数据，进行倾向性评分表格计算
{
  # 指定文件夹路径
  path <- "/home/wuxh/03.MultiPCancer/01.Epidemiology/UKB_metadata/Z-Continuous"
  # 使用正则表达式获取所有以'D'开头的文件名
  files <- list.files(path, pattern = "*")
  # 遍历每个文件，对每个文件执行操作
  for (file_name in files) {
    # 构建完整的文件路径
    file_path <- file.path(path, file_name)
    # 读取文件
    # 应用到读取的数据上
    field <- read.table(file_path, na.strings = c("", "NA", "N/A", "NULL"), fill = TRUE, header = T) %>%
      as.data.frame() 
    names(field)[1] <- "sampleID"
    cat('---------------------------------------------------------------------------------------------\t')
    cat(file_name)
    cat('---------------------------------------------------------------------------------------------\t')
    
    # 选择特定列并重命名（根据需要调整选择的列）
    # 注意：这里假设您之前想要处理并选择的列是"21001-0.0"作为BMI值
    field_select <- field[,c(1,2)]
    names(field_select)[2] <- file_name
    # 与其他数据框合并示例，需要根据实际情况调整
    field_select$sampleID <- gsub("[^0-9]", "", field_select$sampleID)
    field_select$sampleID<-as.integer(field_select$sampleID)
    
    single_field <- left_join(single, field_select, by = "sampleID")
    all_MPC_field <- left_join(all_MPC, field_select, by = "sampleID")
    MPC3_field <- left_join(MPC_more3, field_select, by = "sampleID")
    MPC4_field <- left_join(MPC_more4, field_select, by = "sampleID")
    
    {
      single_field$group<-rep(1,nrow(single_field))
      #73650个样本  
      df=single_field
      single_field_remove_na<-df[!is.na(df[,2]),]
      #73650个样本，没有减少样本
    }
    #all_MPC_1001_diagnose_bim_smok_alcohol_remove_na
    {
      all_MPC_field$group<-rep(0,nrow(all_MPC_field))
      #15379个样本
      df=all_MPC_field
      all_MPC_field_remove_na<-df[!is.na(df[,2]),]
      #15379个样本，没有减少样本
    }
    #MPC3_1001_diagnose_bim_smok_alcohol_remove_na
    {
      MPC3_field$group<-rep(0,nrow(MPC3_field))
      #2235个样本
      df=MPC3_field
      MPC3_field_remove_na<-df[!is.na(df[,2]),]
      #2235个样本，没有减少
    }
    #MPC4_1001_diagnose_bim_smok_alcohol_remove_na
    {
      MPC4_field$group<-rep(0,nrow(MPC4_field))
      #394个样本
      df=MPC4_field
      MPC4_field_remove_na<-df[!is.na(df[,2]),]
      #394个样本，没有减少
    }
    
    pvalue <- function(x, y) {
      # 检查输入是否为数值型向量且不为空
      if (!is.numeric(x) || !is.numeric(y) || length(x) == 0 || length(y) == 0) {
        stop("Error: Input data must be two sets of continuous variables")
      }
      
      # 使用 Wilcox 检验计算 P 值和 95% 置信区间
      wilcox_test <- wilcox.test(x, y, conf.int = TRUE, conf.level = 0.95)
      test_result <- list()
      test_result$type <- "Wilcoxon rank-sum test"
      test_result$p.value <- wilcox_test$p.value
      test_result$logFC <- log2(mean(y) / mean(x))  # 计算 logFC
      test_result$conf.int <- wilcox_test$conf.int  # 95% 置信区间
      
      return(test_result)
    }
    
    #倾向性评分
    {
      all_results <- data.frame()
      #single和所有的
      # single 和所有的
      {
        single_all <- rbind(single_field_remove_na, all_MPC_field_remove_na)
        single_all$group <- factor(single_all$group, levels = c(1, 0), labels = c("single", "MPC"))
        # 确保只传递数值型列到 pvalue 函数
        x_values <- single_all[single_all$group == "single", 2]
        y_values <- single_all[single_all$group == "MPC", 2]
        
        # 检查并处理 NA 和非数值数据
        x_values <- x_values[!is.na(x_values) & is.numeric(x_values)]
        y_values <- y_values[!is.na(y_values) & is.numeric(y_values)]
        
        # 确保 x_values 和 y_values 不为空
        if (length(x_values) == 0 || length(y_values) == 0) {
          stop("Error: After cleaning, one of the input datasets is empty.")
        }
        
        # 计算 P 值和 logFC
        p_value_result <- pvalue(x_values, y_values)
        # 原始数据统计信息，以single和MPC作为列，其他统计信息作为行
        mpc_column_name <- paste("P=", format(as.numeric(p_value_result$p.value), scientific = TRUE, digits = 2),
                                 "_logFC=", format(as.numeric(p_value_result$logFC), scientific = TRUE, digits = 2),
                                 "_CI=(", format(as.numeric(p_value_result$conf.int[1]), scientific = TRUE, digits = 2),
                                 "_", format(as.numeric(p_value_result$conf.int[2]), scientific = TRUE, digits = 2), ")", sep = "")
        summary_stats <- data.frame(
          Statistic = c("Mean", "Median", "Std_Dev", "Min", "Max", "Count"),
          single = c(mean(x_values, na.rm = TRUE),
                     median(x_values, na.rm = TRUE),
                     sd(x_values, na.rm = TRUE),
                     min(x_values, na.rm = TRUE),
                     max(x_values, na.rm = TRUE),
                     length(x_values)),
          # 使用之前生成的列名作为 MPC4 列的名字
          MPC = c(
            mean(y_values, na.rm = TRUE),
            median(y_values, na.rm = TRUE),
            sd(y_values, na.rm = TRUE),
            min(y_values, na.rm = TRUE),
            max(y_values, na.rm = TRUE),
            length(y_values)
          ),
          check.names = FALSE  # 防止将列名转换为合法的R变量名，保留括号等符号
        )
        
        # 保存到文件
        file_name1 <- paste0(file_name,"_single_vs_MPC_P=", 
                             format(p_value_result$p.value, scientific = TRUE, digits = 2), "_logFC=", format(p_value_result$logFC, scientific = TRUE, digits = 2), "_CI=", 
                             "(", format(p_value_result$conf.int[1], scientific = TRUE, digits = 2), ",", format(p_value_result$conf.int[2], scientific = TRUE, digits = 2), ")")
        write.table(summary_stats, file = paste0("/home/wuxh/03.MultiPCancer/01.Epidemiology/UKB_metadata/Z-Continuous.res/single/", file_name1, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
        
        # 将结果按行合并到所有结果的数据框中
        summary_stats$Comparison <- mpc_column_name
        all_results <- rbind(all_results, summary_stats)
      }
      #single和MPC≥3种的
      {
        single_MPC3 <- rbind(single_field_remove_na, MPC3_field_remove_na)
        single_MPC3$group <- factor(single_MPC3$group, levels = c(1, 0), labels = c("single", "MPC≥3"))
        
        # 确保只传递数值型列到 pvalue 函数
        x_values <- single_MPC3[single_MPC3$group == "single", 2]
        y_values <- single_MPC3[single_MPC3$group == "MPC≥3", 2]
        # 检查并处理 NA 和非数值数据
        x_values <- x_values[!is.na(x_values) & is.numeric(x_values)]
        y_values <- y_values[!is.na(y_values) & is.numeric(y_values)]
        
        # 确保 x_values 和 y_values 不为空
        if (length(x_values) == 0 || length(y_values) == 0) {
          stop("Error: After cleaning, one of the input datasets is empty.")
        }
        
        # 计算 P 值和 logFC
        p_value_result <- pvalue(x_values, y_values)
        
        
        mpc3_column_name <- paste("P=", format(as.numeric(p_value_result$p.value), scientific = TRUE, digits = 2),
                                  "_logFC=", format(as.numeric(p_value_result$logFC), scientific = TRUE, digits = 2),
                                  "_CI=(", format(as.numeric(p_value_result$conf.int[1]), scientific = TRUE, digits = 2),
                                  "_", format(as.numeric(p_value_result$conf.int[2]), scientific = TRUE, digits = 2), ")", sep = "")
        summary_stats <- data.frame(
          Statistic = c("Mean", "Median", "Std_Dev", "Min", "Max", "Count"),
          single = c(mean(x_values, na.rm = TRUE),
                     median(x_values, na.rm = TRUE),
                     sd(x_values, na.rm = TRUE),
                     min(x_values, na.rm = TRUE),
                     max(x_values, na.rm = TRUE),
                     length(x_values)),
          # 使用之前生成的列名作为 MPC3 列的名字
          MPC3 = c(
            mean(y_values, na.rm = TRUE),
            median(y_values, na.rm = TRUE),
            sd(y_values, na.rm = TRUE),
            min(y_values, na.rm = TRUE),
            max(y_values, na.rm = TRUE),
            length(y_values)
          ),
          check.names = FALSE  # 防止将列名转换为合法的R变量名，保留括号等符号
        )
        
        # 使用定义好的字符串替换列名
        # names(summary_stats)[3] <- mpc3_column_name
        
        # 保存到文件
        file_name2 <- paste0(file_name,"_single_vs_MPC3_P=", 
                             format(p_value_result$p.value, scientific = TRUE, digits = 2), "_logFC=", format(p_value_result$logFC, scientific = TRUE, digits = 2), "_CI=", 
                             "(", format(p_value_result$conf.int[1], scientific = TRUE, digits = 2), ",", format(p_value_result$conf.int[2], scientific = TRUE, digits = 2), ")")
        write.table(summary_stats, file = paste0("/home/wuxh/03.MultiPCancer/01.Epidemiology/UKB_metadata/Z-Continuous.res/single/", file_name2, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
        
        # 将结果按行合并到所有结果的数据框中
        summary_stats$Comparison <- mpc3_column_name
        all_results <- cbind(all_results, summary_stats)
      }
      #single和MPC≥4种的
      {
        single_MPC4 <- rbind(single_field_remove_na, MPC4_field_remove_na)
        single_MPC4$group <- factor(single_MPC4$group, levels = c(1, 0), labels = c("single", "MPC≥4"))
        
        # 确保只传递数值型列到 pvalue 函数
        x_values <- single_MPC4[single_MPC4$group == "single", 2]
        y_values <- single_MPC4[single_MPC4$group == "MPC≥4", 2]
        
        
        # 检查并处理 NA 和非数值数据
        x_values <- x_values[!is.na(x_values) & is.numeric(x_values)]
        y_values <- y_values[!is.na(y_values) & is.numeric(y_values)]
        
        # 确保 x_values 和 y_values 不为空
        if (length(x_values) == 0 || length(y_values) == 0) {
          stop("Error: After cleaning, one of the input datasets is empty.")
        }
        
        # 计算 P 值和 logFC
        p_value_result <- pvalue(x_values, y_values)
        
        
        # 原始数据统计信息，以single和MPC≥4作为列，其他统计信息作为行
        # 合并统计信息和 P 值、logFC 和 95% CI 信息
        mpc4_column_name <- paste("P=", format(as.numeric(p_value_result$p.value), scientific = TRUE, digits = 2),
                                  "_logFC=", format(as.numeric(p_value_result$logFC), scientific = TRUE, digits = 2),
                                  "_CI=(", format(as.numeric(p_value_result$conf.int[1]), scientific = TRUE, digits = 2),
                                  "_", format(as.numeric(p_value_result$conf.int[2]), scientific = TRUE, digits = 2), ")", sep = "")
        summary_stats <- data.frame(
          Statistic = c("Mean", "Median", "Std_Dev", "Min", "Max", "Count"),
          single = c(mean(x_values, na.rm = TRUE),
                     median(x_values, na.rm = TRUE),
                     sd(x_values, na.rm = TRUE),
                     min(x_values, na.rm = TRUE),
                     max(x_values, na.rm = TRUE),
                     length(x_values)),
          # 使用之前生成的列名作为 MPC4 列的名字
          MPC4 = c(
            mean(y_values, na.rm = TRUE),
            median(y_values, na.rm = TRUE),
            sd(y_values, na.rm = TRUE),
            min(y_values, na.rm = TRUE),
            max(y_values, na.rm = TRUE),
            length(y_values)
          ),
          check.names = FALSE  # 防止将列名转换为合法的R变量名，保留括号等符号
        )
        
        # 使用定义好的字符串替换列名
        # names(summary_stats)[3] <- mpc4_column_name
        
        # 保存到文件
        file_name3 <- paste0(file_name,"_single_vs_MPC4_P=", 
                             format(p_value_result$p.value, scientific = TRUE, digits = 2), "_logFC=", format(p_value_result$logFC, scientific = TRUE, digits = 2), "_CI=", 
                             "(", format(p_value_result$conf.int[1], scientific = TRUE, digits = 2), ",", format(p_value_result$conf.int[2], scientific = TRUE, digits = 2), ")")
        write.table(summary_stats, file = paste0("/home/wuxh/03.MultiPCancer/01.Epidemiology/UKB_metadata/Z-Continuous.res/single/", file_name3, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
        
        # 将结果按行合并到所有结果的数据框中
        summary_stats$Comparison <- mpc4_column_name
        all_results <- cbind(all_results, summary_stats)
      }
      
      all_results<-all_results[,-c(5,6,9,10)]
      # 将所有分组的统计结果保存到一个整合文件中（按行合并）
      integrated_file_name <- paste0("/home/wuxh/03.MultiPCancer/01.Epidemiology/UKB_metadata/Z-Continuous.res/inter/", file_name, "_integrated.csv")
      write.table(all_results, file = integrated_file_name, sep = ",", row.names = FALSE, quote = FALSE)
    }
  }
}
#后续结果整理
{
  # res<-read.csv("/home/wuxh/02.UKBB_SNP_epistasis/04.UKBB/phenodata/all_for_MPC_single_res/zz_res.csv")
  # res_sep<-separate(res,"Column9",into = c("feild","OR"),sep = "_OR=",remove = T)
  # 
  # res_sep<-separate(res_sep,"feild",into = c("feild2","P"),sep = "_P=",remove = T)
  # res_sep<-separate(res_sep,"Column11",into = c("CI2","log2fc"),sep = "_log2fc=",remove = T)
  # 
  # res_sep<-separate(res_sep,"log2fc",into = c("log2fc_new","CI_shang"),sep = "_CI=",remove = T)
  # 
  # res_sep1<-res_sep[grepl("^_log2fc",res_sep[[ncol(res_sep)]]),]
  # res_sep1_sep<-separate(res_sep1,"Column12",c("NA","group"),sep="_log2fc=NA_CI=NA_",remove = T)
  # res_sep1_sep<-separate(res_sep1_sep,"Column10",into=c("no","CI_shang"),sep="_CI=",remove = T)
  # res_sep1_sep_unit<-unite(res_sep1_sep,"CI",c("CI_shang","CI2"),sep = "",remove = T)
  # res_sep1_sep_unit <- res_sep1_sep_unit[which(res_sep1_sep_unit$OR != ""), ]
  # res_sep1_sep_unit_select<-res_sep1_sep_unit[,c(1,2,3,5,8)]
  # names(res_sep1_sep_unit_select)<-c("fields","P-value","OR","95%CI","group")
  # res_sep1_sep_unit_select$type<-rep("OR",nrow(res_sep1_sep_unit_select))
  # 
  # 
  # res_sep2<-res_sep[!grepl("^_log2fc",res_sep[[ncol(res_sep)]]),]
  # res_sep2_sep<-separate(res_sep2,"Column12",c("CI_xia","group1","group2"),sep = "_",remove=T)
  # res_sep2_sep_unit<-unite(res_sep2_sep,"CI",c("CI_shang","CI_xia"),sep = "",remove = T)
  # res_sep2_sep_unit2<-unite(res_sep2_sep_unit,"group",c("group1","group2"),sep = "_",remove = T)
  # res_sep2_sep_unit_select<-res_sep2_sep_unit2[,c(1,2,6,7,8)]
  # names(res_sep2_sep_unit_select)<-c("fields","P-value","OR","95%CI","group")
  # res_sep2_sep_unit_select$type<-rep("log2FC",nrow(res_sep2_sep_unit_select))
  # 
  # 
  # res_total<-rbind(res_sep1_sep_unit_select,res_sep2_sep_unit_select)
  # 
  # write.csv(res_total, file = "/home/wuxh/02.UKBB_SNP_epistasis/04.UKBB/phenodata/all_for_MPC_single_res/zz_res_R_process.csv", row.names = FALSE)
}


