

#单癌GWAS结果
#path:/home/wuxh/luohh/7.five_leadSNP_inner_join_single_cancer_GWAS_res/04.GWAS.res
#-----------------------------------------part1------------------------------------
#载入包
{
  rm(list=ls())
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
}
#-----------------------------------------part2------------------------------------
# 5个leadSNP和单癌样本GWAS结果取交集
{
  five_leadSNP<-read.table("/home/luohh/UKB50wMultiPcancer/03.result/01.MPC/sigGWA.5e-8.header.results",header = T)
  
  five_leadSNP$OR<-exp(five_leadSNP$BETA)
  
  # 设置工作目录为指定路径
  setwd("/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/01.data/")
  
  path <- "/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/01.data"
  
  subdirs <- list.dirs(path, recursive = TRUE)
  
  combined_df <- data.frame()
  
  for (dir in subdirs) {
    # 获取当前子文件夹中所有文件的路径
    files <- list.files(dir, full.names = TRUE)
    
    # 选择以'lala'命名开头的文件
    lala_files <- files[str_detect(basename(files), "^fastGWA_GLMM_final.fastGWA_")]
    
    # 读取并合并这些文件到一个数据框中
    if (length(lala_files) > 0) {
      # 读取第一个文件获取列名
      first_file <- read.table(lala_files[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      col_names <- paste0("V", 1:ncol(first_file))
      
      # 读取所有文件并按列合并
      files_df <- bind_cols(lapply(lala_files, function(file) {
        read.table(file, header = FALSE, sep = "\t", col.names = col_names, stringsAsFactors = FALSE)
      }))  
      
      # 合并到主数据框中
      combined_df <- bind_rows(combined_df, files_df)
    }
  }
  
  
  
  combined_df<-combined_df[,-16]
  names(combined_df)<-c("cancertype","CHR","SNP","POS","A1","A2","N","AF1","T","SE_T","P_noSPA","BETA","SE","P","CONVERGE","Case","Control")
  
  combined_df$OR<-exp(combined_df$BETA)
}


#读取case和control每个基因型数目文件
{
  case_control_count<-read.table("/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/five_SNP_012_count",header = T,sep = "\t")
  
  names(case_control_count)<-c("FullPath_ColumnName","count_0","count_1","count_2","count_NA")
  case_control_count_select_sep<-separate(case_control_count,"FullPath_ColumnName",into = c("yi","er","san","si","wu","liu","qi","ba"),sep = "/")
  case_control_count_select_sep<-case_control_count_select_sep[,6:12]
  case_control_count_select_sep<-separate(case_control_count_select_sep,"ba",into = c("yi","er","SNP"),sep = "_")
  case_control_count_select_sep<-case_control_count_select_sep[,-c(3,4)]
  names(case_control_count_select_sep)<-c("cancertype","group","SNP","count_0","count_1","count_2","count_NA")
}
#将每个SNP对应的case和control中012数目与上述merged_data合并
{
  merge_data<-inner_join(merged_data,case_control_count_select_sep)
  
  merge_data<-unite(merge_data,"pos",c("CHR","POS"),sep = ":")
  merge_data<-unite(merge_data,"alleles",c("A1","A2"),sep = "/")
  
  write.csv(merge_data, file = "/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/MPC_single_Cancer_res_all.csv", 
            quote = FALSE, row.names = FALSE)
  
  merge_data_select<-merge_data[,c(1:13)] %>% unique()
  
  write.csv(merge_data_select, file = "/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/MPC_single_Cancer_res_select.csv", 
            quote = FALSE, row.names = FALSE)
  
}
#计算每种癌症对应case和control中的AF
{
  AF_data <- data.frame()
  # 遍历所有子文件夹
  for (sub_dir in sub_dirs) {
    # 定义需要处理的子目录类型
    types <- c("cancer", "noncancer")
    # 遍历每种类型的子目录
    for (type in types) {
      # 构建完整的子目录路径
      type_sub_dir <- file.path(sub_dir, type)
      # 检查该子目录是否存在
      if (dir.exists(type_sub_dir)) {
        # 获取当前子目录下所有以"AF_count.afreq"结尾的文件
        files <- list.files(path = type_sub_dir, pattern = "AF_count\\.acount$", full.names = TRUE)
        # 遍历当前文件夹下的所有文件
        for (file in files) {
          # 检查文件是否为空
          if (file.size(file) == 0) {
            # 如果文件为空，则打印提示信息并跳过当前文件
            cat("文件为空，跳过文件:", file, "\n")
            next
          }
          # 读取文件内容
          file_data <- read.table(file, header = F, stringsAsFactors = FALSE) # 假设文件有表头
          # 添加来源信息列，以区分数据来自哪个子目录和类型
          file_data$sub_dir <- basename(sub_dir) # 添加子目录名称为列
          file_data$type <- type # 添加类型（cancer/noncancer）为列
          # 将文件内容添加到合并的数据框中
          AF_data <- rbind(AF_data, file_data)
        }
      }
    }
  }
  names(AF_data)<-c("CHROM","SNP","REF","ALT","PROVISIONAL_REF","ALT_CTS","OBS_CT","cancertype","group")
  AF_data<-AF_data[,1:9]
  AF_data$AF_alt<-AF_data$ALT_CTS/AF_data$OBS_CT
  AF_data$AF_ref<-(AF_data$OBS_CT-AF_data$ALT_CTS)/AF_data$OBS_CT
}
#将merge_data和AF_data结合
{
  merge_data_AF<-inner_join(merge_data,AF_data)
  
  merge_data_AF<-merge_data_AF[,c(1:20,25:26)]
  
  merge_data_AF_show<-merge_data_AF[,1:15]
  
  merge_data_AF_plot<-merge_data_AF[,c(1:5,14:24)]
  
  merge_data_AF_plot_select<-merge_data_AF_plot[,c(1,3,8:12)]
  
  #单性别疾病样本基因型count替换
  
  
  
  merge_data_AF_plot_select<-merge_data_AF_plot_select[,c(1:3,8:11)]
  names(merge_data_AF_plot_select)<-c("cancertype","SNP","group","count_0","count_1","count_2","count_NA")
  
  merge_data_AF_plot_long<-melt(merge_data_AF_plot_select,id.vars=c("cancertype","SNP","group"))
  
  levels(merge_data_AF_plot_long$variable) <- c("AA", "Aa", "aa","NA")
  merge_data_AF_plot_long$value<-as.numeric(merge_data_AF_plot_long$value)
}
#画图展示
#加性模型
{
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  # 假定 merge_data_AF 和 merge_data_AF_plot_long 已经被定义和加载
  
  Cond1 <- "cancer"
  Cond4 <- "noncancer"
  
  for (Cancertype in unique(merge_data_AF_plot_long$cancertype)) {
    for (snp in unique(merge_data_AF_plot_long$SNP[merge_data_AF_plot_long$cancertype == Cancertype])) {
      
      Position <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "pos"] # 这里假设你的位置列名是"PositionColumnName"
      snp_id <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "SNP"] # 这里假设你的SNP列名是"SNPColumnName"
      Alleles <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "alleles"] # 这里假设你的等位基因列名是"AllelesColumnName"
      Case <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Case"] # 这里假设你的等位基因列名是"AllelesColumnName"
      Control <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Control"] # 这里假设你的等位基因列名是"AllelesColumnName"
      
      AF_case <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond1, ][1, "ALT_FREQS"] # 这里假设你的等位基因列名是"AllelesColumnName"
      AF_case <- round(AF_case, 3)
      AF_control <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond4, ][1, "ALT_FREQS"] # 这里假设你的等位基因列名是"AllelesColumnName"
      AF_control <- round(AF_control, 3)
      snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond1)
      
      # 使用paste函数来创建标题，变量之间用逗号和空格分隔
      plot_title <- paste(Cancertype, snp_id, Position, Alleles,Case,Control,AF_case,AF_control, sep = ", ")
      
      #绘制Case
      p1 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
        labs(x = "Genotype", y = paste0("Number of ", Cond1, " samples")) +
        scale_fill_manual(values = c("#9BDFDF", "#99CCFF", "#BEBCDF","#2ECC71")) +
        theme_minimal() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
              axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
        scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
      
      #绘制Control
      snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond4)
      snp_data$value <- as.numeric(as.character(snp_data$value))
      
      p2 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
        labs(x = "Genotype", y = paste0("Number of ", Cond4, " samples")) +
        scale_fill_manual(values = c("#C6D6A6", "#FED7D7", "#FFE2CE","#F8B195")) +
        theme_minimal() +  
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
              axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
        scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
      
      title_grob <- textGrob(plot_title, gp = gpar(fontsize = 13), hjust = 0.5)
      
      p <- grid.arrange(p1, p2, ncol = 2, top = title_grob)
      
      # 设置保存路径
      save_path <- file.path("/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/plot", paste0(Cancertype, "_", snp_id))
      
      # ggsave(paste0(save_path, ".png"), plot = p, width = 12, height = 5.4) # 调整保存图像的尺寸
      ggsave(paste0(save_path, ".pdf"), plot = p, width = 6.7, height = 3.7)
    }
  }
  
}
#显性模型
#0，1和2
{
  #将merge_data和AF_data结合
  {
    merge_data_AF<-inner_join(merge_data,AF_data)
    merge_data_AF<-merge_data_AF[,c(1:20,25:28)]
    
    merge_data_AF_show<-merge_data_AF[,1:15]
    
    merge_data_AF_plot<-merge_data_AF[,c(1:5,14:24)]
    
    merge_data_AF_plot_select<-merge_data_AF_plot[,c(1,3,8:12)]
    
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        count_0 = ifelse(cancertype == "02.prostate" & SNP == "rs191674933" & group == "noncancer", 86458, count_0),
        count_1 = ifelse(cancertype == "02.prostate" & SNP == "rs191674933" & group == "noncancer", 54584, count_1),
        count_2 = ifelse(cancertype == "02.prostate" & SNP == "rs191674933" & group == "noncancer", 8534, count_2),
        count_NA = ifelse(cancertype == "02.prostate" & SNP == "rs191674933" & group == "noncancer", 3831, count_NA)
      )
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        count_0 = ifelse(cancertype == "02.prostate" & SNP == "rs3218020" & group == "noncancer", 65000, count_0),
        count_1 = ifelse(cancertype == "02.prostate" & SNP == "rs3218020" & group == "noncancer", 68745, count_1),
        count_2 = ifelse(cancertype == "02.prostate" & SNP == "rs3218020" & group == "noncancer", 18013, count_2),
        count_NA = ifelse(cancertype == "02.prostate" & SNP == "rs3218020" & group == "noncancer", 1649, count_NA)
      )
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        col4 = ifelse(cancertype == "02.prostate" & SNP == "rs11597399" & group == "noncancer", 126026, count_0),
        col5 = ifelse(cancertype == "02.prostate" & SNP == "rs11597399" & group == "noncancer", 25950, count_1),
        col6 = ifelse(cancertype == "02.prostate" & SNP == "rs11597399" & group == "noncancer", 1319, count_2),
        col7 = ifelse(cancertype == "02.prostate" & SNP == "rs11597399" & group == "noncancer", 112, count_NA)
      )
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        col4 = ifelse(cancertype == "02.prostate" & SNP == "rs269812" & group == "noncancer", 118012, count_0),
        col5 = ifelse(cancertype == "02.prostate" & SNP == "rs269812" & group == "noncancer", 31996, count_1),
        col6 = ifelse(cancertype == "02.prostate" & SNP == "rs269812" & group == "noncancer", 2180, count_2),
        col7 = ifelse(cancertype == "02.prostate" & SNP == "rs269812" & group == "noncancer", 1219, count_NA)
      )
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        col4 = ifelse(cancertype == "02.prostate" & SNP == "rs62237617" & group == "noncancer", 152033, count_0),
        col5 = ifelse(cancertype == "02.prostate" & SNP == "rs62237617" & group == "noncancer", 540, count_1),
        col6 = ifelse(cancertype == "02.prostate" & SNP == "rs62237617" & group == "noncancer", 1, count_2),
        col7 = ifelse(cancertype == "02.prostate" & SNP == "rs62237617" & group == "noncancer", 833, count_NA)
      )
    
    
    
    
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        col4 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs191674933" & group == "noncancer", 86458, count_0),
        col5 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs191674933" & group == "noncancer", 54584, count_1),
        col6 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs191674933" & group == "noncancer", 8534, count_2),
        col7 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs191674933" & group == "noncancer", 3831, count_NA)
      )
    merge_data_AF_plot_select <- merge_data_AF_plot_select %>% 
      mutate(
        col4 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs3218020" & group == "noncancer", 65000, count_0),
        col5 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs3218020" & group == "noncancer", 68745, count_1),
        col6 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs3218020" & group == "noncancer", 18013, count_2),
        col7 = ifelse(cancertype == "28.otherMaleGenital" & SNP == "rs3218020" & group == "noncancer", 1649, count_NA)
      )
    
    
    merge_data_AF_plot_select$count_12<-merge_data_AF_plot_select$count_1+merge_data_AF_plot_select$count_2
    
    merge_data_AF_plot_select<-merge_data_AF_plot_select[,c(1:4,8)]
    
    merge_data_AF_plot_long<-melt(merge_data_AF_plot_select,id.vars=c("cancertype","SNP","group"))
    
    levels(merge_data_AF_plot_long$variable) <- c("AA", "Aa and aa")
    merge_data_AF_plot_long$value<-as.numeric(merge_data_AF_plot_long$value)
  }
  #画图展示
  {
    library(ggplot2)
    library(gridExtra)
    library(grid)
    
    # 假定 merge_data_AF 和 merge_data_AF_plot_long 已经被定义和加载
    
    Cond1 <- "cancer"
    Cond4 <- "noncancer"
    
    for (Cancertype in unique(merge_data_AF_plot_long$cancertype)) {
      for (snp in unique(merge_data_AF_plot_long$SNP[merge_data_AF_plot_long$cancertype == Cancertype])) {
        
        Position <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "pos"] # 这里假设你的位置列名是"PositionColumnName"
        snp_id <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "SNP"] # 这里假设你的SNP列名是"SNPColumnName"
        Alleles <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "alleles"] # 这里假设你的等位基因列名是"AllelesColumnName"
        Case <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Case"] # 这里假设你的等位基因列名是"AllelesColumnName"
        Control <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Control"] # 这里假设你的等位基因列名是"AllelesColumnName"
        
        AF_case <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond1, ][1, "AF_alt"] # 这里假设你的等位基因列名是"AllelesColumnName"
        AF_case <- round(AF_case, 3)
        AF_control <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond4, ][1, "AF_alt"] # 这里假设你的等位基因列名是"AllelesColumnName"
        AF_control <- round(AF_control, 3)
        snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond1)
        
        # 使用paste函数来创建标题，变量之间用逗号和空格分隔
        plot_title <- paste(Cancertype, snp_id, Position, Alleles,Case,Control,AF_case,AF_control, sep = ", ")
        
        #绘制Case
        p1 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
          labs(x = "Genotype", y = paste0("Number of ", Cond1, " samples")) +
          scale_fill_manual(values = c("#9BDFDF", "#99CCFF", "#BEBCDF","#2ECC71")) +
          theme_minimal() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
                axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
          scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
        
        #绘制Control
        snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond4)
        snp_data$value <- as.numeric(as.character(snp_data$value))
        
        p2 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
          labs(x = "Genotype", y = paste0("Number of ", Cond4, " samples")) +
          scale_fill_manual(values = c("#C6D6A6", "#FED7D7", "#FFE2CE","#F8B195")) +
          theme_minimal() +  
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
                axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
          scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
        
        title_grob <- textGrob(plot_title, gp = gpar(fontsize = 13), hjust = 0.5)
        
        p <- grid.arrange(p1, p2, ncol = 2, top = title_grob)
        
        # 设置保存路径
        save_path <- file.path("/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/plot", paste0(Cancertype, "_", snp_id,"_dominant"))
        
        # ggsave(paste0(save_path, ".png"), plot = p, width = 12, height = 5.4) # 调整保存图像的尺寸
        ggsave(paste0(save_path, ".pdf"), plot = p, width = 6.7, height = 3.7)
      }
    }
    
  }
  
}
#隐性模型
#0，1和2
{
  #将merge_data和AF_data结合
  {
    merge_data_AF<-inner_join(merge_data,AF_data)
    merge_data_AF<-merge_data_AF[,c(1:20,25:28)]
    
    merge_data_AF_show<-merge_data_AF[,1:15]
    
    merge_data_AF_plot<-merge_data_AF[,c(1:5,14:24)]
    
    merge_data_AF_plot_select<-merge_data_AF_plot[,c(1,3,8:12)]
    
    merge_data_AF_plot_select$count_01<-merge_data_AF_plot_select$count_0+merge_data_AF_plot_select$count_1
    
    merge_data_AF_plot_select<-merge_data_AF_plot_select[,c(1:3,6,8)]
    
    merge_data_AF_plot_long<-melt(merge_data_AF_plot_select,id.vars=c("cancertype","SNP","group"))
    
    levels(merge_data_AF_plot_long$variable) <- c("aa","AA and Aa")
    merge_data_AF_plot_long$value<-as.numeric(merge_data_AF_plot_long$value)
  }
  #画图展示
  {
    library(ggplot2)
    library(gridExtra)
    library(grid)
    
    # 假定 merge_data_AF 和 merge_data_AF_plot_long 已经被定义和加载
    
    Cond1 <- "cancer"
    Cond4 <- "noncancer"
    
    for (Cancertype in unique(merge_data_AF_plot_long$cancertype)) {
      for (snp in unique(merge_data_AF_plot_long$SNP[merge_data_AF_plot_long$cancertype == Cancertype])) {
        
        Position <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "pos"] # 这里假设你的位置列名是"PositionColumnName"
        snp_id <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "SNP"] # 这里假设你的SNP列名是"SNPColumnName"
        Alleles <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "alleles"] # 这里假设你的等位基因列名是"AllelesColumnName"
        Case <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Case"] # 这里假设你的等位基因列名是"AllelesColumnName"
        Control <- merge_data_AF[merge_data_AF$SNP == snp, ][1, "Control"] # 这里假设你的等位基因列名是"AllelesColumnName"
        
        AF_case <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond1, ][1, "AF_ref"] # 这里假设你的等位基因列名是"AllelesColumnName"
        AF_case <- round(AF_case, 3)
        AF_control <- merge_data_AF[merge_data_AF$SNP == snp & merge_data_AF$group == Cond4, ][1, "AF_ref"] # 这里假设你的等位基因列名是"AllelesColumnName"
        AF_control <- round(AF_control, 3)
        snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond1)
        
        # 使用paste函数来创建标题，变量之间用逗号和空格分隔
        plot_title <- paste(Cancertype, snp_id, Position, Alleles,Case,Control,AF_case,AF_control, sep = ", ")
        
        #绘制Case
        p1 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
          labs(x = "Genotype", y = paste0("Number of ", Cond1, " samples")) +
          scale_fill_manual(values = c("#9BDFDF", "#99CCFF", "#BEBCDF","#2ECC71")) +
          theme_minimal() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
                axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
          scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
        
        #绘制Control
        snp_data <- subset(merge_data_AF_plot_long, cancertype == Cancertype & SNP == snp_id & group == Cond4)
        snp_data$value <- as.numeric(as.character(snp_data$value))
        
        p2 <- ggplot(data = snp_data, aes(x = variable, y = value, fill = variable)) +
          geom_bar(stat = "identity", position = position_dodge()) +
          geom_text(aes(label = value), vjust = -0.3, position = position_dodge(width = 0.9)) +
          labs(x = "Genotype", y = paste0("Number of ", Cond4, " samples")) +
          scale_fill_manual(values = c("#C6D6A6", "#FED7D7", "#FFE2CE","#F8B195")) +
          theme_minimal() +  
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12), panel.border = element_rect(colour = "black", fill = NA, size = 1),
                axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13)) +
          scale_y_continuous(limits = c(0, max(snp_data$value) * 1.2)) # 自动调整Y轴范围
        
        title_grob <- textGrob(plot_title, gp = gpar(fontsize = 13), hjust = 0.5)
        
        p <- grid.arrange(p1, p2, ncol = 2, top = title_grob)
        
        # 设置保存路径
        save_path <- file.path("/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/plot", paste0(Cancertype, "_", snp_id,"_recessive"))
        
        # ggsave(paste0(save_path, ".png"), plot = p, width = 12, height = 5.4) # 调整保存图像的尺寸
        ggsave(paste0(save_path, ".pdf"), plot = p, width = 6.7, height = 3.7)
      }
    }
    
  }
  
}

merge_data_select<-merge_data[,c(1:7,11,13,15:20)]
merge_data_select <- merge_data_select %>%
   mutate(
      count_0_per = ifelse(group == "cancer", count_0 / Case, count_0 / Control),
      count_1_per = ifelse(group == "cancer", count_1 / Case, count_1 / Control),
      count_2_per = ifelse(group == "cancer", count_2 / Case, count_2 / Control)
    )
write.csv(merge_data_select, file = "/home/wuxh/luohh/6.five_leadSNP_inner_join_singlePC_cancer_GWAS_res/MPC_single_Cancer_012_count_per.csv", 
             quote = FALSE, row.names = FALSE)
