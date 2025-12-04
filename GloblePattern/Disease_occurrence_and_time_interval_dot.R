
# 2025.7.4.
# 绘制发病次数点图和发病时间间隔点图


{
  rm(List=ls())
  gc()
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(lubridate)
}
# 绘制发病时间间隔点图
{
  time <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.first2Cancer.date")
  time_sep <- separate(time, "V2", c("MPC1", "time1"), sep = ":")
  time_sep2 <- separate(time_sep, "V3", c("MPC2", "time2"), sep = ":")
  time_sep2 <- time_sep2 %>% 
    mutate(time1 = ymd(time1),
           time2 = ymd(time2)
    )
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

  
  
  # 创建自定义颜色方案 - 参考频次热图的配色
  col <- colorRampPalette(c("grey","#fb0c0c","#da0909","#fd5959","#18ceff","#1f9bff","#1932ff"))(20)
  # 绘制时间间隔热图
  ggplot(time_complete, aes(MPC1, MPC2, color = avg_date_diff_years, size = avg_date_diff_years)) +
    geom_point(alpha = 0.7) +
    scale_color_gradientn(colors = col, name = "Time Interval\n(Years)") +
    scale_size_continuous(name = "Time Interval\n(Years)", range = c(0.1, 5.5)) +
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

{
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  
  # 数据读取和处理部分（保持不变，但稍作优化）
  time <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.first2Cancer.date")
  
  # 数据处理管道优化
  time_processed <- time %>%
    separate(V2, c("MPC1", "time1"), sep = ":") %>%
    separate(V3, c("MPC2", "time2"), sep = ":") %>%
    mutate(
      time1 = ymd(time1),
      time2 = ymd(time2),
      date_diff_days = as.numeric(time2 - time1),
      date_diff_years = date_diff_days / 365.25,
      same_day = (time1 == time2)
    ) %>%
    filter(!same_day) %>%  # 过滤掉同一天诊断的病例
    group_by(MPC1, MPC2) %>%
    summarise(avg_date_diff_years = mean(date_diff_years, na.rm = TRUE), .groups = 'drop')
  
  # 创建完整的癌症对组合矩阵
  all_cancers <- sort(unique(c(time_processed$MPC1, time_processed$MPC2)))
  time_complete <- expand.grid(MPC1 = all_cancers, MPC2 = all_cancers, stringsAsFactors = FALSE) %>%
    left_join(time_processed, by = c("MPC1", "MPC2")) %>%
    mutate(
      avg_date_diff_years = ifelse(is.na(avg_date_diff_years), 0, avg_date_diff_years),
      # 创建分类变量用于颜色映射
      interval_category = case_when(
        avg_date_diff_years == 0 ~ "No Data",
        TRUE ~ "Has Data"
      ),
      # 创建用于大小映射的变量
      point_size = ifelse(avg_date_diff_years == 0, 0.5, avg_date_diff_years)
    )
  
  time_complete <- time_complete %>%
    mutate(
      MPC1 = factor(MPC1, levels = mpc_order),
      MPC2 = factor(MPC2, levels = mpc_order)
    )
  
  # 为非0值创建颜色渐变（深红→浅红→浅蓝→深蓝）
  non_zero_data <- time_complete %>% filter(avg_date_diff_years > 0)
  if(nrow(non_zero_data) > 0) {
    color_breaks <- seq(min(non_zero_data$avg_date_diff_years), 
                        max(non_zero_data$avg_date_diff_years), 
                        length.out = 6)
    
    # 自定义颜色：深红→浅红→浅蓝→深蓝
    custom_colors <- c("#b92308","#fb0c0c","#899bed", "#4169E1", "#000080")
  }
  
  # 绘制优化后的热图
  p <- ggplot(time_complete, aes(x = MPC1, y = MPC2)) +
    # 先绘制所有0值点（灰色小点）
    geom_point(data = filter(time_complete, avg_date_diff_years == 0),
               color = "grey70", size = 0.5, alpha = 0.6) +
    # 再绘制非0值点（彩色渐变）
    geom_point(data = filter(time_complete, avg_date_diff_years > 0),
               aes(color = avg_date_diff_years, size = avg_date_diff_years),
               alpha = 0.8) +
    # 颜色映射
    scale_color_gradientn(
      colors = custom_colors,
      name = "Time Interval\n(Years)",
      breaks = color_breaks,
      labels = round(color_breaks, 1)
    ) +
    # 大小映射
    scale_size_continuous(
      name = "Time Interval\n(Years)",
      range = c(1, 6),
      breaks = color_breaks,
      labels = round(color_breaks, 1)
    ) +
    # 图形设置
    labs(
      title = "Average Time Interval Between First and Second Cancer Diagnosis",
      subtitle = "Grey dots indicate no data available",
      x = "First Cancer",
      y = "Second Cancer"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey50"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right"
    )
  
  # 显示图形
  p
  
  # 如果需要保存图形
  ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/cancer_time_interval_heatmap.pdf", 
         plot = p, width = 11.4, height = 7.6)
}



# 发病次数点图
{
  data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_order_times",
                     header = T, sep = "\t", fill = T)
  data_sep <- separate(data, "MPC1.to.MPC2", into = c("MPC1", "MPC2"), sep = "-", remove = T)
  all_cancers <- sort(unique(c(data_sep$MPC1, data_sep$MPC2)))
  all_combinations <- expand.grid(MPC1 = all_cancers, MPC2 = all_cancers, stringsAsFactors = FALSE)
  # data_complete <- all_combinations %>%
  #   left_join(data_sep, by = c("MPC1", "MPC2")) %>%
  #   mutate(Times = ifelse(is.na(Times), 0, Times))
  
  data_complete <- all_combinations %>%
    left_join(data_sep, by = c("MPC1", "MPC2"))
  
  col <- colorRampPalette(c("#3ea3c0","#6baed6","#f48591","#f55162","#f7021b"))(20)
  ggplot(data_complete, aes(MPC1, MPC2, color = Times, size = Times)) +
    geom_point(alpha = 0.7) +
    scale_color_gradientn(colors = col, name = "Frequency") +
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


{
  library(tidyverse)
  library(ggplot2)
  
  # 数据读取和处理
  data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_order_times",
                     header = T, sep = "\t", fill = T)
  
  # 数据处理优化
  data_processed <- data %>%
    separate(MPC1.to.MPC2, into = c("MPC1", "MPC2"), sep = "-", remove = T)
  
  # 创建完整的癌症对组合矩阵
  all_cancers <- sort(unique(c(data_processed$MPC1, data_processed$MPC2)))
  data_complete <- expand.grid(MPC1 = all_cancers, MPC2 = all_cancers, stringsAsFactors = FALSE) %>%
    left_join(data_processed, by = c("MPC1", "MPC2")) %>%
    mutate(
      Times = ifelse(is.na(Times), 0, Times),
      # 创建分类变量用于颜色映射
      frequency_category = case_when(
        Times == 0 ~ "No Data",
        TRUE ~ "Has Data"
      )
    )
  
  # 为非0值创建50的整数倍断点（强制使用50的间隔）
  non_zero_data <- data_complete %>% filter(Times > 0)
  if(nrow(non_zero_data) > 0) {
    # 获取数据范围
    min_val <- min(non_zero_data$Times)
    max_val <- max(non_zero_data$Times)
    
    # 将最小值向下取整到最近的50的倍数
    min_break <- floor(min_val / 50) * 50
    # 将最大值向上取整到最近的50的倍数
    max_break <- ceiling(max_val / 50) * 50
    
    # 强制创建50的整数倍断点，不管数量多少
    color_breaks <- seq(min_break, max_break, by = 50)
    
    # 如果断点过多（超过12个），可以选择性地减少断点数量，但仍保持50的倍数
    if(length(color_breaks) > 12) {
      # 创建更稀疏的50倍数断点，比如每150或每200
      step_size <- ceiling(length(color_breaks) / 8) * 50  # 确保步长是50的倍数
      color_breaks <- seq(min_break, max_break, by = step_size)
    }
    
    # 确保包含实际的最小值和最大值附近的断点
    if(min_val > 0 && !any(color_breaks >= min_val)) {
      color_breaks <- c(color_breaks, ceiling(min_val / 50) * 50)
    }
    if(!any(color_breaks >= max_val)) {
      color_breaks <- c(color_breaks, ceiling(max_val / 50) * 50)
    }
    
    # 去重、排序并确保都是正数
    color_breaks <- sort(unique(color_breaks[color_breaks >= 0]))
    
    # 自定义颜色：深蓝→浅红色→橙红→橙→深红
    # 在深蓝和橙红之间插入浅红色
    
    "#b92308","#fb0c0c","#899bed", "#4169E1", "#000080"
    
    
    custom_colors <- c("#a7b4ef","#FF4500","#fb0c0c","#b92308")
    
    # 如果你想要更浅的红色，可以选择以下任一选项：
    # custom_colors <- c("#4169E1", "#FFA0A0", "#FF6347", "#FF4500", "#8B0000")  # 浅红色选项1
    # custom_colors <- c("#4169E1", "#FF9999", "#FF6347", "#FF4500", "#8B0000")  # 浅红色选项2
    # custom_colors <- c("#4169E1", "#FFCCCB", "#FF6347", "#FF4500", "#8B0000")  # 浅红色选项3（更浅）
  }
  
  # 绘制优化后的频次热图
  p <- ggplot(data_complete, aes(x = MPC1, y = MPC2)) +
    # 先绘制所有0值点（灰色小点）
    geom_point(data = filter(data_complete, Times == 0),
               color = "grey70", size = 0.5, alpha = 0.6) +
    # 再绘制非0值点（彩色渐变）
    geom_point(data = filter(data_complete, Times > 0),
               aes(color = Times, size = Times),
               alpha = 0.8) +
    # 颜色映射 - 使用50的整数倍标签
    scale_color_gradientn(
      colors = custom_colors,
      name = "Frequency",
      breaks = color_breaks,
      labels = color_breaks  # 直接使用断点值，都是50的倍数
    ) +
    # 大小映射 - 使用50的整数倍标签
    scale_size_continuous(
      name = "Frequency",
      range = c(1, 6),
      breaks = color_breaks,
      labels = color_breaks  # 直接使用断点值，都是50的倍数
    ) +
    # 图形设置
    labs(
      title = "Frequency of First-to-Second Cancer Sequences",
      subtitle = "Grey dots indicate no observed sequences",
      x = "First Cancer",
      y = "Second Cancer"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey50"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right"
    )
  
  # 显示图形
  p
  
  ggsave("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/cancer_diagnose_times_heatmap.pdf", 
         plot = p, width = 11.4, height = 7.6)
}




