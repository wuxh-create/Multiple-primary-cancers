
# ^-^
# 2025.6.19
# 原发癌globalpattern画热图
# 时空分析。两两转换比例，和时间间隔，除了桑吉图，是不是还可以用heatmap来展示。
# wu xiaohong
# ^-^


#-----------------------------------------part1------------------------------------
#载入包
{
  rm(list=ls())
  gc()
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(table1)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(scales)
  library(lubridate)
}
#-----------------------------------------part2 诊断顺序------------------------------------

# 第一原发癌到第二原发癌的点图
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # 读取和处理数据
  data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_order_times",
                     header = T, sep = "\t", fill = T)
  data_sep <- separate(data, "MPC1.to.MPC2", into = c("MPC1", "MPC2"), sep = "-", remove = T)
  
  # 获取所有唯一的癌症类型
  all_cancers <- sort(unique(c(data_sep$MPC1, data_sep$MPC2)))
  
  # 创建所有可能的癌症对组合
  all_combinations <- expand.grid(MPC1 = all_cancers, MPC2 = all_cancers, stringsAsFactors = FALSE)
  
  # 将原始数据与完整组合合并，缺失值填充为0
  data_complete <- all_combinations %>%
    left_join(data_sep, by = c("MPC1", "MPC2")) %>%
    mutate(Times = ifelse(is.na(Times), 0, Times))
  
  # 查看数据统计
  cat("总癌症类型数:", length(all_cancers), "\n")
  cat("总组合数:", nrow(data_complete), "\n")
  cat("有数据的组合数:", sum(data_complete$Times > 0), "\n")
  cat("Times范围:", min(data_complete$Times), "-", max(data_complete$Times), "\n")
  
  # 绘制自定义颜色渐变热图
  # 创建自定义颜色方案 - 从浅蓝色到白色到红色
  col <- colorRampPalette(c("#f2eeef","#f3bac0","#f48591","#f55162","#f7021b"))(20)
  
  ggplot(data_complete, aes(MPC1, MPC2, color = Times, size = Times)) +
    geom_point(alpha = 0.7) +
    
    # 修改颜色映射：使用自定义颜色方案
    scale_color_gradientn(colors = col, name = "Frequency") +
    
    # 设置大小范围，避免0值点太小看不见
    scale_size_continuous(name = "Frequency", range = c(0.5, 6)) +
    
    labs(title = "Frequency of First-to-Second Cancer Sequences",
         x = "First Cancer",
         y = "Second Cancer") +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black", size = 1),
      plot.background = element_rect(fill = "white", color = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
}
# 第一原发癌到第二原发癌的热图
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
  
  # 读取和处理数据
  data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_order_times",
                     header = T, sep = "\t", fill = T)
  
  data_sep <- separate(data, "MPC1.to.MPC2", into = c("MPC1", "MPC2"), sep = "-", remove = T)
  
  # 获取所有癌症类型
  all_cancers <- sort(unique(c(data_sep$MPC1, data_sep$MPC2)))
  
  # 创建完整的数据框，填充缺失的组合
  data_complete <- data_sep %>%
    complete(MPC1 = all_cancers, MPC2 = all_cancers, fill = list(Times = 0))
  
  # 方法3: 将数据分成8个离散区间
  data_complete <- data_complete %>%
    mutate(
      Times_category = case_when(
        Times == 0 ~ "0",
        Times %in% 1:50 ~ "1-50",
        Times %in% 51:100 ~ "51-100",
        Times %in% 101:150 ~ "101-150",
        Times %in% 151:200 ~ "151-200",
        Times %in% 201:250 ~ "201-250",
        Times > 250 ~ ">250",
        TRUE ~ "0"
      )
    )
  
  # 设置分类顺序
  data_complete$Times_category <- factor(data_complete$Times_category, 
                                         levels = c("0", "1-50", "51-100", "101-150", "151-200", "201-250",">250"))
  
  # 创建8色调色板
  discrete_colors <- c("#F0F0F0", "#FFEEEE", 
                       "#FF9999", "#FF6666", "#FF3333", "#CC0000", "#800000")
  
  
  # discrete_colors <- c("#F0F0F0", "#CCDDFF", "#99BBFF", "#6699FF", 
  #                      "#FF9999", "#FF6666", "#FF3333", "#CC0000")
  names(discrete_colors) <- levels(data_complete$Times_category)
  
  plot3 <- ggplot(data_complete, aes(MPC1, MPC2, fill = Times_category)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(
      values = discrete_colors,
      name = "Frequency",
      na.value = "lightgrey"
    ) +
    labs(title = "Frequency of First-to-Second Cancer Sequences",
         x = "First Cancer",
         y = "Second Cancer") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black", size = 1),
      plot.background = element_rect(fill = "white", color = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
  plot3
}
#-----------------------------------------part2 诊断顺序------------------------------------
# 第一原发癌到第二原发癌的时间间隔诊断点图
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(lubridate)
  
  # 读取和处理时间数据
  time <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.first2Cancer.date")
  
  time_sep <- separate(time, "V2", c("MPC1", "time1"), sep = ":")
  time_sep2 <- separate(time_sep, "V3", c("MPC2", "time2"), sep = ":")
  
  time_sep2 <- time_sep2 %>%
    mutate(
      time1 = ymd(time1),
      time2 = ymd(time2)
    )
  
  # 计算两列的差值（以天为单位，然后转换为年）
  time_sep2 <- time_sep2 %>%
    mutate(
      date_diff_days = as.numeric(time2 - time1),
      date_diff_years = date_diff_days / 365.25  # 转换为年（考虑闰年）
    )
  
  time_sep2 <- time_sep2 %>%
    mutate(new_column = if_else(time1 == time2, "Y", "N"))
  
  # 过滤掉同一天诊断的病例
  time_sep2_filter <- time_sep2 %>% filter(new_column == 'N')
  
  # 按癌症类型分组，计算平均时间间隔（年）
  time_sep2_filter_mean <- time_sep2_filter %>%
    group_by(MPC1, MPC2) %>%
    summarise(avg_date_diff_years = mean(date_diff_years, na.rm = TRUE), .groups = 'drop')
  
  # 获取所有唯一的癌症类型
  all_cancers <- sort(unique(c(time_sep2_filter_mean$MPC1, time_sep2_filter_mean$MPC2)))
  
  # 创建所有可能的癌症对组合
  all_combinations <- expand.grid(MPC1 = all_cancers, MPC2 = all_cancers, stringsAsFactors = FALSE)
  
  # 将原始数据与完整组合合并，缺失值填充为0
  time_complete <- all_combinations %>%
    left_join(time_sep2_filter_mean, by = c("MPC1", "MPC2")) %>%
    mutate(avg_date_diff_years = ifelse(is.na(avg_date_diff_years), 0, avg_date_diff_years))
  
  # 查看数据统计
  cat("总癌症类型数:", length(all_cancers), "\n")
  cat("总组合数:", nrow(time_complete), "\n")
  cat("有数据的组合数:", sum(time_complete$avg_date_diff_years > 0), "\n")
  cat("时间间隔统计（年）:\n")
  cat("最小值:", round(min(time_complete$avg_date_diff_years, na.rm = TRUE), 2), "年\n")
  cat("最大值:", round(max(time_complete$avg_date_diff_years, na.rm = TRUE), 2), "年\n")
  cat("平均值:", round(mean(time_complete$avg_date_diff_years[time_complete$avg_date_diff_years > 0], na.rm = TRUE), 2), "年\n")
  cat("中位数:", round(median(time_complete$avg_date_diff_years[time_complete$avg_date_diff_years > 0], na.rm = TRUE), 2), "年\n")
  
  # 创建自定义颜色方案 - 参考频次热图的配色
  col <- colorRampPalette(c("#f2eeef","#f3bac0","#f48591","#f55162","#f7021b"))(20)
  
  # 绘制时间间隔热图
  ggplot(time_complete, aes(MPC1, MPC2, color = avg_date_diff_years, size = avg_date_diff_years)) +
    geom_point(alpha = 0.7) +
    
    # 使用自定义颜色方案
    scale_color_gradientn(colors = col, name = "Time Interval\n(Years)") +
    
    # 设置大小范围，参考频次热图的设置
    scale_size_continuous(name = "Time Interval\n(Years)", range = c(0.5, 5.5)) +
    
    labs(title = "Average Time Interval Between First and Second Cancer Diagnosis",
         x = "First Cancer",
         y = "Second Cancer") +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black", size = 1),
      plot.background = element_rect(fill = "white", color = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
}