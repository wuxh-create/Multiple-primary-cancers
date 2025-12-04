# 读入所有有诊断时间的MPC样本
{
  rm(list=ls())
  gc()
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggrepel) # 用于避免文字遮挡
}

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

{
  # 生成需要的文件类型
  {
    # 设置目标文件路径
    input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC/"
    
    # 获取所有以 .mpc.sample 结尾的文件
    file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample$", full.names = TRUE)
    
    # 遍历文件并处理
    for (file_path in file_list) {
      # 提取文件名
      file_name <- basename(file_path)
      
      # 提取癌症类型（即 .mpc.sample 之前的部分）
      cancer_type <- sub("\\.mpc\\.sample$", "", file_name)
      
      # 读取文件内容
      sample_ids <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
      colnames(sample_ids) <- c("SampleID") # 设置列名为 SampleID
      
      # 添加新列：癌症类型
      sample_ids$CancerType <- cancer_type
      
      # 生成新文件路径，添加 .cancer 后缀
      new_file_path <- file.path(input_dir, paste0(file_name, ".cancer"))
      
      # 写入新文件
      write.table(sample_ids, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    }
    
    # 提示完成
    cat("所有文件已处理并生成新的文件。\n")
  }
  # 创建颜色
  {
    all_cancer_types <- c(
      "anus", "bone", "brain", "breast", "cervix", "colorectal", "esophagus",
      "eyeAndOrbit", "gallbladderAndBiliaryTract", "headAndNeck",
      "kidney", "liver", "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma",
      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", "otherFemaleGenital",
      "otherMaleGenital", "otherNervousSystem", "otherRespiratory", "otherUrinaryOrgans",
      "ovary", "pancreas", "prostate", "smallIntestine", "softTissueSarcoma",
      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", "urinaryBladder",
      "uterus")
    # 创建35种独特的颜色
    color_palette <- c(
      # 粉红色系列
      "#FF69B4", "#FF1493", "#FF66CC", "#FF99CC", "#FF3366", "#FF99FF", "#FF0066",
      # 黄色系列
      "#FFD700", "#FFA500", "#FFFF00", "#FFE666", "#FFF68F", "#FFB366", "#FFDB4D",
      # 蓝色系列
      "#00BFFF", "#1E90FF", "#4169E1", "#87CEEB", "#40E0D0", "#00CED1", "#48D1CC",
      # 绿色系列 (偏黄)
      "#98FB98", "#90EE90", "#32CD32", "#9ACD32", "#7FFF00", "#ADFF2F", "#7CFF00",
      # 绿色系列 (偏蓝)
      "#00FA9A", "#00FF7F", "#3CB371", "#66CDAA", "#20B2AA", "#48FFBD", "#43CD80",
      # 额外补充的糖果色
      "#FF1744", "#FF4081", "#F50057", "#D500F9", "#651FFF"
    )
    # 创建颜色映射
    cancer_colors <- setNames(color_palette, all_cancer_types)
  }
  # Bone
  {
    # 数据预处理
    {
      input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
      # 读取 bone.mpc.sample.cancer 文件
      bone_file <- file.path(input_dir, "bone.mpc.sample.cancer")
      bone_data <- read.table(bone_file, header = TRUE, stringsAsFactors = FALSE)
      # 获取所有以 .mpc.sample.cancer 结尾的文件
      file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample\\.cancer$", full.names = TRUE)
      # 初始化结果存储
      result <- data.frame(FileName = character(), RowCount = numeric(), stringsAsFactors = FALSE)
      # 遍历其他文件，与 bone_file 文件依次取交集
      for (file_path in file_list) {
        # 跳过 bone 文件自身
        if (file_path == bone_file) {
          next
        }
        # 提取当前文件的癌症类型
        file_name <- basename(file_path)
        cancer_type <- sub("\\.mpc\\.sample\\.cancer$", "", file_name)
        # 读取当前文件内容
        current_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        # 取交集（根据 SampleID 列）
        intersect_data <- merge(bone_data, current_data, by = "SampleID")
        # 提取新文件名（格式为 bone 和当前文件的癌症类型组合）
        new_file_name <- paste0("bone_", cancer_type)
        new_file_path <- file.path(input_dir, new_file_name)
        # 写入新文件
        # write.table(intersect_data, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        # 获取交集文件的行数
        row_count <- nrow(intersect_data)
        # 将结果添加到 result 数据框中
        result <- rbind(result, data.frame(FileName = new_file_name, RowCount = row_count))
      }
      
      # 删除行数为 0 的记录
      result <- result[result$RowCount > 0, ]
      # 按 RowCount 列降序排序
      result <- result[order(-result$RowCount), ]
      result<-separate(result,"FileName",c("yi","er"),sep = "_")
      result<-result[,c(2,3)] %>% as.data.frame()
      names(result)<-c("FileName","RowCount")
      result$RowCount <- as.numeric(result$RowCount)
    }
    # 画图数据预处理
    {
      lab1 <- as.vector(result$FileName)
      lab1 <- paste(lab1, "(", round(result$RowCount/sum(result$RowCount)*100, 2), "%)", sep="")
      lab2 <- round(result$RowCount/sum(result$RowCount)*100, 1)
      # 构建数据框并排序
      A <- data.frame(lab2, lab1)
      A_sorted <- A[order(-A$lab2),]
      # 更新排序后的标签
      lab2_sorted <- A_sorted$lab2
      lab1_sorted <- A_sorted$lab1
      # 获取当前数据中癌症类型的原始名称（不包含百分比）
      current_cancer_types <- gsub("\\(.*\\)$", "", lab1_sorted)  # 去除括号及其中的内容
      current_cancer_types <- trimws(current_cancer_types)         # 去除可能的空格
      # 获取对应的颜色
      colors_sorted <- cancer_colors[current_cancer_types]
    }
    # 画图
    {
      # 设置绘图参数
      par(mar = c(5, 4, 4, 14))  # 显著增加右边距
      par(xpd = TRUE)  # 允许绘制超出图形区域
      
      pdf("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.3.Bone_related_Subsequent_Cancer_Pie_Chart.pdf",
          width = 5.8, height = 8.6)
      
      # 绘制空心圆环图
      pie(lab2_sorted, 
          labels = lab1_sorted,
          radius = 1.0, 
          clockwise = TRUE,
          col = colors_sorted,
          main = "2.3.Bone_related_Subsequent_Cancer_Pie_Chart.pdf",
          cex = 0.8)  # 稍微减小主图中的文字大小
      
      # 在中心添加白色圆形以创建环形效果
      par(new = TRUE)
      pie(1, 
          radius = 0.5,
          col = "white", 
          border = "white",
          labels = "")
      
      # 添加图例，调整位置和大小以确保完整显示
      legend("topright",
             legend = lab1_sorted,
             fill = colors_sorted,
             cex = 0.65,        # 减小图例文字大小
             xpd = TRUE,        # 允许图例绘制在图形区域外
             inset = c(-0.5, 0), # 增加向右的偏移量
             bty = "n",         # 移除图例边框
             y.intersp = 0.8)   # 减小图例项之间的间距
      dev.off()
    }
  }
  # gallbladderAndBiliaryTract
  {
    # 数据预处理
    {
      input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
      # 读取 gallbladderAndBiliaryTract.mpc.sample.cancer 文件
      gallbladderAndBiliaryTract_file <- file.path(input_dir, "gallbladderAndBiliaryTract.mpc.sample.cancer")
      gallbladderAndBiliaryTract_data <- read.table(gallbladderAndBiliaryTract_file, header = TRUE, stringsAsFactors = FALSE)
      # 获取所有以 .mpc.sample.cancer 结尾的文件
      file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample\\.cancer$", full.names = TRUE)
      # 初始化结果存储
      result <- data.frame(FileName = character(), RowCount = numeric(), stringsAsFactors = FALSE)
      # 遍历其他文件，与 gallbladderAndBiliaryTract_file 文件依次取交集
      for (file_path in file_list) {
        # 跳过 bone 文件自身
        if (file_path == gallbladderAndBiliaryTract_file) {
          next
        }
        # 提取当前文件的癌症类型
        file_name <- basename(file_path)
        cancer_type <- sub("\\.mpc\\.sample\\.cancer$", "", file_name)
        # 读取当前文件内容
        current_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        # 取交集（根据 SampleID 列）
        intersect_data <- merge(gallbladderAndBiliaryTract_data, current_data, by = "SampleID")
        # 提取新文件名（格式为 bone 和当前文件的癌症类型组合）
        new_file_name <- paste0("gallbladderAndBiliaryTract_", cancer_type)
        new_file_path <- file.path(input_dir, new_file_name)
        # 写入新文件
        # write.table(intersect_data, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        # 获取交集文件的行数
        row_count <- nrow(intersect_data)
        # 将结果添加到 result 数据框中
        result <- rbind(result, data.frame(FileName = new_file_name, RowCount = row_count))
      }
      # 删除行数为 0 的记录
      result <- result[result$RowCount > 0, ]
      # 按 RowCount 列降序排序
      result <- result[order(-result$RowCount), ]
      result<-separate(result,"FileName",c("yi","er"),sep = "_")
      result<-result[,c(2,3)] %>% as.data.frame()
      names(result)<-c("FileName","RowCount")
      result$RowCount <- as.numeric(result$RowCount)
    }
    # 画图数据预处理
    {
      lab1 <- as.vector(result$FileName)
      lab1 <- paste(lab1, "(", round(result$RowCount/sum(result$RowCount)*100, 2), "%)", sep="")
      lab2 <- round(result$RowCount/sum(result$RowCount)*100, 1)
      # 构建数据框并排序
      A <- data.frame(lab2, lab1)
      A_sorted <- A[order(-A$lab2),]
      # 更新排序后的标签
      lab2_sorted <- A_sorted$lab2
      lab1_sorted <- A_sorted$lab1
      # 获取当前数据中癌症类型的原始名称（不包含百分比）
      current_cancer_types <- gsub("\\(.*\\)$", "", lab1_sorted)  # 去除括号及其中的内容
      current_cancer_types <- trimws(current_cancer_types)         # 去除可能的空格
      # 获取对应的颜色
      colors_sorted <- cancer_colors[current_cancer_types]
    }
    # 画图
    {
      # 设置绘图参数
      par(mar = c(5, 4, 4, 14))  # 显著增加右边距
      par(xpd = TRUE)  # 允许绘制超出图形区域
      
      pdf("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.3.gallbladderAndBiliaryTract_related_Subsequent_Cancer_Pie_Chart.pdf",
          width = 5.8, height = 8.6)
      
      # 绘制空心圆环图
      pie(lab2_sorted, 
          labels = lab1_sorted,
          radius = 1.0, 
          clockwise = TRUE,
          col = colors_sorted,
          main = "2.3.gallbladderAndBiliaryTract_related_Subsequent_Cancer_Pie_Chart.pdf",
          cex = 0.8)  # 稍微减小主图中的文字大小
      
      # 在中心添加白色圆形以创建环形效果
      par(new = TRUE)
      pie(1, 
          radius = 0.5,
          col = "white", 
          border = "white",
          labels = "")
      
      # 添加图例，调整位置和大小以确保完整显示
      legend("topright",
             legend = lab1_sorted,
             fill = colors_sorted,
             cex = 0.65,        # 减小图例文字大小
             xpd = TRUE,        # 允许图例绘制在图形区域外
             inset = c(-0.5, 0), # 增加向右的偏移量
             bty = "n",         # 移除图例边框
             y.intersp = 0.8)   # 减小图例项之间的间距
      dev.off()
    }
  }
  # smallIntestine
  {
    # 数据预处理
    {
      input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
      # 读取 smallIntestine.mpc.sample.cancer 文件
      smallIntestine_file <- file.path(input_dir, "smallIntestine.mpc.sample.cancer")
      smallIntestine_data <- read.table(smallIntestine_file, header = TRUE, stringsAsFactors = FALSE)
      # 获取所有以 .mpc.sample.cancer 结尾的文件
      file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample\\.cancer$", full.names = TRUE)
      # 初始化结果存储
      result <- data.frame(FileName = character(), RowCount = numeric(), stringsAsFactors = FALSE)
      # 遍历其他文件，与 bone_file 文件依次取交集
      for (file_path in file_list) {
        # 跳过 bone 文件自身
        if (file_path == smallIntestine_file) {
          next
        }
        # 提取当前文件的癌症类型
        file_name <- basename(file_path)
        cancer_type <- sub("\\.mpc\\.sample\\.cancer$", "", file_name)
        # 读取当前文件内容
        current_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        # 取交集（根据 SampleID 列）
        intersect_data <- merge(smallIntestine_data, current_data, by = "SampleID")
        # 提取新文件名（格式为 smallIntestine 和当前文件的癌症类型组合）
        new_file_name <- paste0("smallIntestine_", cancer_type)
        new_file_path <- file.path(input_dir, new_file_name)
        # 写入新文件
        # write.table(intersect_data, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        # 获取交集文件的行数
        row_count <- nrow(intersect_data)
        # 将结果添加到 result 数据框中
        result <- rbind(result, data.frame(FileName = new_file_name, RowCount = row_count))
      }
      
      # 删除行数为 0 的记录
      result <- result[result$RowCount > 0, ]
      # 按 RowCount 列降序排序
      result <- result[order(-result$RowCount), ]
      result<-separate(result,"FileName",c("yi","er"),sep = "_")
      result<-result[,c(2,3)] %>% as.data.frame()
      names(result)<-c("FileName","RowCount")
      result$RowCount <- as.numeric(result$RowCount)
    }
    # 画图数据预处理
    {
      lab1 <- as.vector(result$FileName)
      lab1 <- paste(lab1, "(", round(result$RowCount/sum(result$RowCount)*100, 2), "%)", sep="")
      lab2 <- round(result$RowCount/sum(result$RowCount)*100, 1)
      # 构建数据框并排序
      A <- data.frame(lab2, lab1)
      A_sorted <- A[order(-A$lab2),]
      # 更新排序后的标签
      lab2_sorted <- A_sorted$lab2
      lab1_sorted <- A_sorted$lab1
      # 获取当前数据中癌症类型的原始名称（不包含百分比）
      current_cancer_types <- gsub("\\(.*\\)$", "", lab1_sorted)  # 去除括号及其中的内容
      current_cancer_types <- trimws(current_cancer_types)         # 去除可能的空格
      # 获取对应的颜色
      colors_sorted <- cancer_colors[current_cancer_types]
    }
    # 画图
    {
      # 设置绘图参数
      par(mar = c(5, 4, 4, 14))  # 显著增加右边距
      par(xpd = TRUE)  # 允许绘制超出图形区域
      
      pdf("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.3.smallIntestine_related_Subsequent_Cancer_Pie_Chart.pdf",
          width = 5.8, height = 8.6)
      
      # 绘制空心圆环图
      pie(lab2_sorted, 
          labels = lab1_sorted,
          radius = 1.0, 
          clockwise = TRUE,
          col = colors_sorted,
          main = "2.3.smallIntestine_related_Subsequent_Cancer_Pie_Chart.pdf",
          cex = 0.8)  # 稍微减小主图中的文字大小
      
      # 在中心添加白色圆形以创建环形效果
      par(new = TRUE)
      pie(1, 
          radius = 0.5,
          col = "white", 
          border = "white",
          labels = "")
      
      # 添加图例，调整位置和大小以确保完整显示
      legend("topright",
             legend = lab1_sorted,
             fill = colors_sorted,
             cex = 0.65,        # 减小图例文字大小
             xpd = TRUE,        # 允许图例绘制在图形区域外
             inset = c(-0.5, 0), # 增加向右的偏移量
             bty = "n",         # 移除图例边框
             y.intersp = 0.8)   # 减小图例项之间的间距
      dev.off()
    }
  }
  # otherUrinaryOrgans
  {
    # 数据预处理
    {
      input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
      # 读取 otherUrinaryOrgans.mpc.sample.cancer 文件
      otherUrinaryOrgans_file <- file.path(input_dir, "otherUrinaryOrgans.mpc.sample.cancer")
      otherUrinaryOrgans_data <- read.table(otherUrinaryOrgans_file, header = TRUE, stringsAsFactors = FALSE)
      # 获取所有以 .mpc.sample.cancer 结尾的文件
      file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample\\.cancer$", full.names = TRUE)
      # 初始化结果存储
      result <- data.frame(FileName = character(), RowCount = numeric(), stringsAsFactors = FALSE)
      # 遍历其他文件，与 otherUrinaryOrgans_file 文件依次取交集
      for (file_path in file_list) {
        # 跳过 bone 文件自身
        if (file_path == otherUrinaryOrgans_file) {
          next
        }
        # 提取当前文件的癌症类型
        file_name <- basename(file_path)
        cancer_type <- sub("\\.mpc\\.sample\\.cancer$", "", file_name)
        # 读取当前文件内容
        current_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        # 取交集（根据 SampleID 列）
        intersect_data <- merge(otherUrinaryOrgans_data, current_data, by = "SampleID")
        # 提取新文件名（格式为 bone 和当前文件的癌症类型组合）
        new_file_name <- paste0("otherUrinaryOrgans_", cancer_type)
        new_file_path <- file.path(input_dir, new_file_name)
        # 写入新文件
        # write.table(intersect_data, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        # 获取交集文件的行数
        row_count <- nrow(intersect_data)
        # 将结果添加到 result 数据框中
        result <- rbind(result, data.frame(FileName = new_file_name, RowCount = row_count))
      }
      # 删除行数为 0 的记录
      result <- result[result$RowCount > 0, ]
      # 按 RowCount 列降序排序
      result <- result[order(-result$RowCount), ]
      result<-separate(result,"FileName",c("yi","er"),sep = "_")
      result<-result[,c(2,3)] %>% as.data.frame()
      names(result)<-c("FileName","RowCount")
      result$RowCount <- as.numeric(result$RowCount)
    }
    # 画图数据预处理
    {
      lab1 <- as.vector(result$FileName)
      lab1 <- paste(lab1, "(", round(result$RowCount/sum(result$RowCount)*100, 2), "%)", sep="")
      lab2 <- round(result$RowCount/sum(result$RowCount)*100, 1)
      # 构建数据框并排序
      A <- data.frame(lab2, lab1)
      A_sorted <- A[order(-A$lab2),]
      # 更新排序后的标签
      lab2_sorted <- A_sorted$lab2
      lab1_sorted <- A_sorted$lab1
      # 获取当前数据中癌症类型的原始名称（不包含百分比）
      current_cancer_types <- gsub("\\(.*\\)$", "", lab1_sorted)  # 去除括号及其中的内容
      current_cancer_types <- trimws(current_cancer_types)         # 去除可能的空格
      # 获取对应的颜色
      colors_sorted <- cancer_colors[current_cancer_types]
    }
    # 画图
    {
      # 设置绘图参数
      par(mar = c(5, 4, 4, 14))  # 显著增加右边距
      par(xpd = TRUE)  # 允许绘制超出图形区域
      
      pdf("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.3.otherUrinaryOrgans_related_Subsequent_Cancer_Pie_Chart.pdf",
          width = 5.8, height = 8.6)
      
      # 绘制空心圆环图
      pie(lab2_sorted, 
          labels = lab1_sorted,
          radius = 1.0, 
          clockwise = TRUE,
          col = colors_sorted,
          main = "2.3.otherUrinaryOrgans_related_Subsequent_Cancer_Pie_Chart.pdf",
          cex = 0.8)  # 稍微减小主图中的文字大小
      
      # 在中心添加白色圆形以创建环形效果
      par(new = TRUE)
      pie(1, 
          radius = 0.5,
          col = "white", 
          border = "white",
          labels = "")
      
      # 添加图例，调整位置和大小以确保完整显示
      legend("topright",
             legend = lab1_sorted,
             fill = colors_sorted,
             cex = 0.65,        # 减小图例文字大小
             xpd = TRUE,        # 允许图例绘制在图形区域外
             inset = c(-0.5, 0), # 增加向右的偏移量
             bty = "n",         # 移除图例边框
             y.intersp = 0.8)   # 减小图例项之间的间距
      dev.off()
    }
  }
  # otherDigestive
  {
    # 数据预处理
    {
      input_dir <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
      # 读取 otherDigestive.mpc.sample.cancer 文件
      otherDigestive_file <- file.path(input_dir, "otherDigestive.mpc.sample.cancer")
      otherDigestive_data <- read.table(otherDigestive_file, header = TRUE, stringsAsFactors = FALSE)
      # 获取所有以 .mpc.sample.cancer 结尾的文件
      file_list <- list.files(input_dir, pattern = "\\.mpc\\.sample\\.cancer$", full.names = TRUE)
      # 初始化结果存储
      result <- data.frame(FileName = character(), RowCount = numeric(), stringsAsFactors = FALSE)
      # 遍历其他文件，与 otherDigestive_file 文件依次取交集
      for (file_path in file_list) {
        # 跳过 bone 文件自身
        if (file_path == otherDigestive_file) {
          next
        }
        # 提取当前文件的癌症类型
        file_name <- basename(file_path)
        cancer_type <- sub("\\.mpc\\.sample\\.cancer$", "", file_name)
        # 读取当前文件内容
        current_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
        # 取交集（根据 SampleID 列）
        intersect_data <- merge(otherDigestive_data, current_data, by = "SampleID")
        # 提取新文件名（格式为 bone 和当前文件的癌症类型组合）
        new_file_name <- paste0("otherDigestive_", cancer_type)
        new_file_path <- file.path(input_dir, new_file_name)
        # 写入新文件
        # write.table(intersect_data, file = new_file_path, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        # 获取交集文件的行数
        row_count <- nrow(intersect_data)
        # 将结果添加到 result 数据框中
        result <- rbind(result, data.frame(FileName = new_file_name, RowCount = row_count))
      }
      # 删除行数为 0 的记录
      result <- result[result$RowCount > 0, ]
      # 按 RowCount 列降序排序
      result <- result[order(-result$RowCount), ]
      result<-separate(result,"FileName",c("yi","er"),sep = "_")
      result<-result[,c(2,3)] %>% as.data.frame()
      names(result)<-c("FileName","RowCount")
      result$RowCount <- as.numeric(result$RowCount)
    }
    # 画图数据预处理
    {
      lab1 <- as.vector(result$FileName)
      lab1 <- paste(lab1, "(", round(result$RowCount/sum(result$RowCount)*100, 2), "%)", sep="")
      lab2 <- round(result$RowCount/sum(result$RowCount)*100, 1)
      # 构建数据框并排序
      A <- data.frame(lab2, lab1)
      A_sorted <- A[order(-A$lab2),]
      # 更新排序后的标签
      lab2_sorted <- A_sorted$lab2
      lab1_sorted <- A_sorted$lab1
      # 获取当前数据中癌症类型的原始名称（不包含百分比）
      current_cancer_types <- gsub("\\(.*\\)$", "", lab1_sorted)  # 去除括号及其中的内容
      current_cancer_types <- trimws(current_cancer_types)         # 去除可能的空格
      # 获取对应的颜色
      colors_sorted <- cancer_colors[current_cancer_types]
    }
    # 画图
    {
      # 设置绘图参数
      par(mar = c(5, 4, 4, 10))  # 显著增加右边距
      par(xpd = TRUE)  # 允许绘制超出图形区域
      
      pdf("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/04.plot/2.3.otherDigestive_related_Subsequent_Cancer_Pie_Chart.pdf",
          width = 5.8, height = 8.6)
      
      # 绘制空心圆环图
      pie(lab2_sorted, 
          labels = lab1_sorted,
          radius = 1.0, 
          clockwise = TRUE,
          col = colors_sorted,
          main = "2.3.otherDigestive_related_Subsequent_Cancer_Pie_Chart.pdf",
          cex = 0.8)  # 稍微减小主图中的文字大小
      
      # 在中心添加白色圆形以创建环形效果
      par(new = TRUE)
      pie(1, 
          radius = 0.5,
          col = "white", 
          border = "white",
          labels = "")
      
      # 添加图例，调整位置和大小以确保完整显示
      legend("topright",
             legend = lab1_sorted,
             fill = colors_sorted,
             cex = 0.65,        # 减小图例文字大小
             xpd = TRUE,        # 允许图例绘制在图形区域外
             inset = c(-0.5, 0), # 增加向右的偏移量
             bty = "n",         # 移除图例边框
             y.intersp = 0.8)   # 减小图例项之间的间距
      dev.off()
    }
  }
}














