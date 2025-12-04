#^-^
#2024.10.25
#原发癌globalpattern部分的图和写作细节描述数据
#wu xiaohong
#^-^


#数目统计
{
  rm(list=ls())
  gc()
  # 读取癌症对数据
  pcancerData <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/1.MPCpair.table", header = T, stringsAsFactors = F)
  pcancerData <- as.matrix(pcancerData)
  
  # 将邻接矩阵转换为长格式数据框
  library(reshape2)
  cancer_pairs <- melt(pcancerData)
  colnames(cancer_pairs) <- c("Cancer1", "Cancer2", "CaseCount")
  
  # 移除对角线和NA值（如果有）
  cancer_pairs <- cancer_pairs[cancer_pairs$Cancer1 != cancer_pairs$Cancer2, ]  # 移除对角线（即相同癌症对自己）
  cancer_pairs <- na.omit(cancer_pairs)  # 移除NA值
  
  # 按照病例数进行降序排序
  cancer_pairs_sorted <- cancer_pairs[order(-cancer_pairs$CaseCount), ]
  # 查看排名最高的前10个癌症对
  top_pairs <- head(cancer_pairs_sorted, 20)
  (top_pairs)
  
  # 查看有多少去重后的癌症对
  {
    df<-cancer_pairs_sorted[which(cancer_pairs_sorted$CaseCount > 0),] %>% as.data.frame()
    names(df)<-c("Col1", "Col2","Value")
    df_sorted <- t(apply(df[, c("Col1", "Col2")], 1, sort)) # 排序第一列和第二列
    df_sorted <- as.data.frame(df_sorted)                  # 转换为数据框
    colnames(df_sorted) <- c("Col1", "Col2")               # 恢复列名
    # 添加原始列（比如 Value）
    df_sorted$Value <- df$Value
    # 去重操作
    df_unique <- df_sorted[!duplicated(df_sorted[, c("Col1", "Col2")]), ]
  }
}
#确定上述统计出数目的顺序
{
  mpc2order<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date",
                        header = F)
  mpc2order_sep1<-separate(mpc2order,"V2",c("cancer1","data1"),sep = ":")
  
  mpc2order_sep2<-separate(mpc2order_sep1,"V3",c("cancer2","data2"),sep = ":")
  
  MPC1_position_times<-table(mpc2order_sep2$cancer1) %>% as.data.frame()
  names(MPC1_position_times)<-c("MPC1_occurenece_position","Times")
  write.table(MPC1_position_times, 
              file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_MPC1_occurence_MPC1_position_times", 
              quote = FALSE, row.names = FALSE,sep = "\t")
  
  mpc2order_sep2_unite<-unite(mpc2order_sep2,"unite_cancer",c("cancer1","cancer2"),sep ="-")
  
  result <- subset(mpc2order_sep2, cancer1 == "urinaryBladder" & cancer2 == "prostate")
  
  result_mpc2_order_times<-table(mpc2order_sep2_unite$unite_cancer) %>% as.data.frame()
  
  names(result_mpc2_order_times)<-c("MPC1-to-MPC2","Times")
  
  result_mpc2_order_times_more_25<-result_mpc2_order_times[which(result_mpc2_order_times$Times > 25),]
  
  write.table(result_mpc2_order_times, 
            file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/01.data/03.diagnosisOrder/cancerDiagnose.allMPC.2mpc.date_order_times", 
            quote = FALSE, row.names = FALSE,sep = "\t")
  
  
}
#circos画圈图
{
  ---
    title: "MPC circo"
  author: "luohh"
  date: "2023-12-20"
  output: html_document
  ---
    
    
    
    ```{r}
  ### CircoPlot
  ### 示意图
  library(circlize)
  ####################################################################################################################################
  ### 10 x 10
  
  library(circlize)
  library(reshape2)
  
  pcancerData <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/1.MPCpair.table",header = T,stringsAsFactors = F)
  
  pcancerData <- as.matrix(pcancerData)
  
  pcancerData_frame<-as.data.frame(pcancerData)
  rownames(pcancerData_frame)<-rownames(pcancerData)
  colnames(pcancerData_frame)<-colnames(pcancerData)
  
  #circos.initialize()
  circos.par(gap.after = c(rep(1, nrow(pcancerData_frame)-1), 1))
  
  
  
  {
    pcancerData <- structure(
      list(
        colorectal = c(186, 32, 23, 21, 20, 18, 11, 9, 6, 6, 6, 5, 4, 3),
        breast = c(11, 441, 504, 388, 102, 418, 22, 195, 115, 91, 66, 148, 47, 95),
        cervix = c(30, 43, 37, 19, 1, 25, 48, 21, 12, 4, 2, 17, 0, 3),
        lung = c(462, 504, 119, 119, 101, 96, 0, 54, 33, 22, 19, 33, 18, 28),
        otherFemaleGenital = c(108, 108, 37, 257, 37, 231, 613, 152, 131, 79, 511, 218, 61, 78),
        melanoma = c(75, 462, 108, 257, 37, 231, 613, 152, 131, 79, 511, 218, 61, 78),
        prostate = c(87, 56, 22, 62, 6, 24, 83, 27, 23, 13, 16, 35, 11, 59),
        lymphoidNeoplasms = c(14, 53, 106, 28, 4, 51, 12, 12, 11, 8, 1, 7, 1, 2),
        kidney = c(31, 34, 14, 11, 1, 9, 26, 9, 7, 2, 13, 10, 240, 2),
        myeloidNeoplasms = c(78, 95, 28, 109, 5, 40, 88, 30, 15, 13, 7, 29, 14, 12)
      ),
      row.names = c("colorectal", "breast", "cervix", "lung", "otherFemaleGenital", "melanoma", "prostate", "lymphoidNeoplasms", "kidney", "myeloidNeoplasms"),
      class = "data.frame"
    )
    pcancerData$id <- rownames(pcancerData)
    long_pcancerData <- melt(pcancerData, id.vars = "id")
    colnames(long_pcancerData) <- c("from", "to", "value")
    arc_lengths <- aggregate(value ~ from, data = long_pcancerData, FUN = sum)
   
    sorted_order <- arc_lengths[order(arc_lengths$value), ]$from
    arc_lengths <- arc_lengths[order(arc_lengths$value, decreasing = TRUE), ]
    # new_order <- c(arc_lengths$from)
    
    new_order <- as.character(arc_lengths$from)
    num_sectors <- length(unique(long_pcancerData$from))
    gap_degree <- rep(1, num_sectors)
    circos.par(gap.after = gap_degree)
    chordDiagram(long_pcancerData, order = new_order)
    
  }
  
  # 添加id列
  pcancerData_frame$id <- rownames(pcancerData_frame)
  
  # 将数据转换为长格式
  long_pcancerData <- melt(pcancerData_frame, id.vars = "id")
  colnames(long_pcancerData) <- c("from", "to", "value")
  
  # 计算每个节点的总连接数
  arc_lengths <- aggregate(value ~ from, data = long_pcancerData, FUN = sum)
  
  # 按照总连接数降序排列并获取排序后的顺序
  sorted_order <- arc_lengths[order(arc_lengths$value, decreasing = TRUE), ]$from
  
  # 设置circos参数
  num_sectors <- length(unique(long_pcancerData$from))
  gap_degree <- rep(1, num_sectors)
  circos.par(gap.after = gap_degree,start.degree = 90)
  
  # 绘制chord图，按新顺序排列
  chordDiagram(long_pcancerData, order = sorted_order,annotationTrack = c("name","grid"))
  # annotationTrack = c("name","grid")去掉刻度线
  
  
  
  
  circos.clear()
  ```
  
  
  {
    pcancerData <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/1.MPCpair.table", header = TRUE, stringsAsFactors = FALSE)
    pcancerData <- as.matrix(pcancerData)
    # 检查并替换零值或 NA
    pcancerData[pcancerData == 0 | is.na(pcancerData)] <- 1  # 避免空白扇区
    # 确保 sector_names 包含所有行列名称
    sector_names <- unique(c(rownames(pcancerData), colnames(pcancerData)))
    # 动态生成与扇区数量一致的糖果色调
    candy_colors <- colorRampPalette(c(
      "#D2B4DE", "#FF69B4", "#77DD77", "#32CD32", "#00FF00",
      "#FFB347", "#F7DC6F", "#FFB6C1","#B0E0E6", "#FF69B4",
      "#00FA9A", "#1E90FF", "#19D3DA", "#00BFFF", "#8CD0FF",
      "#9932CC", "#8A2BE2", "#FF1493", "#FF00FF", "#FF1493",
      "#FF5722", "#EE82EE", "#FFC0CB", "#FA5D63", "#EE82EE"
    ))(length(sector_names))  # 动态生成足够的颜色
    
    candy_colors <-c("#00FF00","#FF5722","#00BFFF","#1E90FF","#77DD77","#FFB347","#F3A3EAFF","#32CD32",
                       "#98EAE8FF","#FA5D63","#1504C6FF","#A4FAC6FF","#4DD8DCFF","#0AD0CCFF","#08669BFF","#FDC1D9FF",
                       "#17E9B4FF","#500DEBFF","#7CE2FCFF","#4A7204FF","#DE5779FF","#A4F198FF","#FF7F11","#F02882FF",
                       "#17C65DFF","#060106FF","#FD085AFF","#71BAEBFF","#EE82EE","#AFB8BDFF","#3FEE0BFF","#D2FC89FF",
                       "#92F5D7FF","#AAF2FDFF","#3115D4FF","#146ADCFF")
    
    
    
    # 创建一个带有扇区名称的颜色向量
    grid_col <- setNames(candy_colors, sector_names)
    # circos 参数设置
    circos.par(gap.after = c(rep(1, length(sector_names) - 1), 1))
    # 绘制 chord diagram
    chordDiagram(
      pcancerData,
      order = sector_names,  # 按指定顺序
      grid.col = grid_col,   # 使用带命名的颜色向量
      transparency = 0.05    # 较低透明度以保持颜色鲜亮
    )
  }
  
  
}

#绘制圈图中转移数目top 20样本的转移柱形图
{
  top20<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/3.1.number.Top20.cancerpair")
}

#统计每种癌症中单癌样本数量和MPC样本数量以及占比
{
  rm(list=ls())
  NumandProportion<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/2.indexCancer.NumAndProportion.data",header = T)
}
#alluvial图，什么使用灰色表示？什么使用蓝色表示？
{
  ---
    title: "4.2.Rplot-alluvial"
  author: "luohh"
  date: "2023-12-21"
  output: html_document
  ---
    
    
    
    ```{r}
  ### alluvial 包 （炜炜师兄）
  remotes::install_github("mbojan/alluvial", build_vignettes=TRUE)
  library(alluvial)
  library(data.table)
  library(randomcoloR)
  
  ### 导入数据
  data <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/4.1.alluvial.2mpc.gt25.data",header = T,stringsAsFactors = F)
  
  cols <- hsv(h = sample(1:8/10), s = sample(3:8)/8, v = sample(3:8)/8)
  cols <- randomColor(count = 68) 
  cols <- rep("#C2C2C2",68)
  cols[data$value>100] <- "#D1B2FF"
  ord <- list(NULL, with(data, order(cancer1, cancer2)))
  
  ### 冲击图 
  alluvial(
    data[,1:2], alpha=1,
    freq = data$value,
    blocks = T,
    gap.width=0.2, # 块的间隔
    xw = 0.3, # 线的曲度
    cw = 0.01, # 字块大小
    # col=c('#BDBDBD','#4EAE4A','#984EA3')
    # axis_labels='',
    axes=T,
    ann=T,
    #col = cols[match(data$cancer1, unique(data$cancer2))],
    col = cols,
    ordering=ord
  )
  ```
  
  
}

#癌症器官的特异性分析中，多原发癌样本之间是否存在重叠
{
  breast<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/03.SpecificMPC/01.breast/breast.mpc.sample")
  
  
}


#读取有诊断时间的样本信息
{
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
  # 筛选白人样本
  {
    # 读入所有白人MPC样本
    MPC<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/multiPCancer.gcta.sample",sep = "\t")
    MPC<-MPC[,1] %>% as.data.frame()
    names(MPC)<-names(allMPC)[1]
    # 筛选所有有诊断时间的白人样本
    allMPC_clinic_white<-inner_join(allMPC,MPC)
  }
  
  
  # 将每种癌症中MPC的占比排序
  {
    MPC_order<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/01.globlePattern/03.result/2.indexCancer.NumAndProportion.data",
                          header = T)
  }
  #筛选排名前5的，统计发病规律
  {
    df_filtered_gallbladderAndBiliaryTract <- allMPC_clinic_white %>%
      filter(
        allMPC_clinic_white[[2]] == "gallbladderAndBiliaryTract" |
          allMPC_clinic_white[[4]] == "gallbladderAndBiliaryTract" |
          allMPC_clinic_white[[6]] == "gallbladderAndBiliaryTract" |
          allMPC_clinic_white[[8]] == "gallbladderAndBiliaryTract"
      )
    df_filtered_bone <- allMPC_clinic_white %>%
      filter(
        allMPC_clinic_white[[2]] == "bone" |
          allMPC_clinic_white[[4]] == "bone" |
          allMPC_clinic_white[[6]] == "bone" |
          allMPC_clinic_white[[8]] == "bone"
      )
    df_filtered_smallIntestine <- allMPC_clinic_white %>%
      filter(
        allMPC_clinic_white[[2]] == "smallIntestine" |
          allMPC_clinic_white[[4]] == "smallIntestine" |
          allMPC_clinic_white[[6]] == "smallIntestine" |
          allMPC_clinic_white[[8]] == "smallIntestine"
      )
    df_filtered_stomach <- allMPC_clinic_white %>%
      filter(
        allMPC_clinic_white[[2]] == "stomach" |
          allMPC_clinic_white[[4]] == "stomach" |
          allMPC_clinic_white[[6]] == "stomach" |
          allMPC_clinic_white[[8]] == "stomach"
      )
  }
}

#TWAS结果画图
{
  
  
  
  #这个脚本的目的是绘制MPC TWAS结果的柱形图，方便直观的展示一个基因在多少个组织中发挥作用。
  
  #part1-------------------------------------------load R packages-------------------------------------------------------------
  {
    rm(list=ls())
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(table1)
    library(stringr)
    library(ggplot2)
  }
  #part1-------------------------------------------generate a data frame with genes and their counts-------------------------------------------------------------
  {
    string_column <- c("AC078795.2", "AL133355.1", "AS3MT", "BORCS7","CALHM2","CALHM3","CCDC159","CDKN2B","CNNM2","ELOF1","GSTO1","GSTO2","ITPRIP",
                       "KLHL9","LINC01239","LRRC34","MYNN","NT5C2","PCGF6","PDCD11","RNASEH2A","SH3PXD2A","SH3PXD2A-AS1","SLC7A14","SLK","STN1",
                       "ZNF20","ZNF441","ZNF491","ZNF625","ZNF627","ZNF700","ZNF709","ZNF763")
    numeric_column <- c(1, 4, 3, 3,1,1,1,1,2,1,3,1,1,
                        1,1,2,1,1,1,3,1,1,1,1,3,11,
                        3,1,1,1,1,1,1,1)
    
    # 创建数据框
    df <- data.frame(Gene = string_column, Count = numeric_column, stringsAsFactors = FALSE)
    
    #控制画图的顺序根据count值从大到小排序
    df<-df[order(df$Count,decreasing = TRUE),]
    df$Gene<-factor(df$Gene,levels = df$Gene)
    
    ggplot(df, aes(x = Gene, y = Count,fill = Gene)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 75, hjust = 1,vjust = 1),
            axis.ticks.x = element_line(color = "black"),
            axis.ticks.length.x = unit(0, "cm"),
            legend.position = "none",
            panel.background = element_blank(),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(title = "Gene Counts", x = "Gene", y = "Count")+
      coord_cartesian(ylim = c(0, max(df$Count)))
    
    
    
    
    p<-ggplot(df, aes(x = Gene, y = Count, fill = Gene)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 12),
            axis.ticks.x = element_line(color = "black"),
            axis.ticks.length.x = unit(0.3, "cm"),
            legend.position = "none",
            panel.background = element_blank(),
            axis.line.x = element_line(color = "black"), # 添加黑色横坐标轴
            axis.line.y = element_line(color = "black"), # 添加黑色纵坐标轴
            axis.title.x = element_text(size = 15), # 调整横坐标轴标题字体大小
            axis.title.y = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      labs(x = "Gene", y = "The number of tissues")
    
    
    
    p + scale_y_continuous(expand = c(0,0),limits = c(0,12)) +
      ggtitle("TWAS results for significant SNP loci in MPC") +
      # 调整标题居中
      theme(plot.title = element_text(hjust = 0.5,size = 16))
  }
  
  
  
  
  
}

#绘制GWAS结果的曼哈顿图
{
  library(qqman)
  library(CMplot)
  library(data.table)
  file_path <- "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/03.result/02.colorectal/fastGWA_GLMM_final.fastGWA"
  
  gwasResult <- fread(file_path,header = T,stringsAsFactors = F,data.table = F)
  gwasResult.plot <- gwasResult[,c(2,1,3,13)]
  colnames(gwasResult.plot) <- c("SNP","chr","pos","p")
  
  gwasResult.plot<-gwasResult.plot[1:10000,]
  
  CMplot(gwasResult.plot,plot.type = "m",LOG10 = T,threshold = c(5e-6), threshold.lty = c(1), threshold.lwd = c(1), \
         threshold.col=c("black"), signal.col=c("#B0282F"), highlight.col=c("#B0282F"), col=c("#2A5888","#4396B1"), 
         signal.cex=c(1), signal.pch=c(19), cex = 0.7,cex.axis=1.5,cex.lab=2,file="pdf")
  
  CMplot(gwasResult.plot,plot.type="m",c(5e-6),threshold.col=c('grey','black'), 
         threshold.lty=c(1,2),threshold.lwd=c(1,1),amplify=T, signal.cex=c(1,1),signal.pch=c(20,20),signal.col=c("red","orange"))
  CMplot(gwasResult.plot,plot.type = "q",conf.int = F,box = FALSE,threshold.lty = 2,threshold.col = "red",cex=0.6,cex.axis=1.5,cex.lab=2,file="pdf")
}

#绘制GWAS结果的QQ plot图
{
  # 安装并加载所需的R包
  # install.packages("qqman")
  library(qqman)
  library(ggplot2)
  
  # 假设您的GWAS结果数据有一列包含p值，名为 'P'
  # 导入您的数据，例如data.csv
  # 替换'path/to/your/file.csv'为您的文件路径
  gwas_results <- read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/fastGWA_GLMM_final.fastGWA",header = T)
  
  # QQ-plot 使用 ggplot2
  qq_plot <- function(pvalues) {
    # 计算预期的 -log10(p) 值
    n <- length(pvalues)
    expected <- -log10(ppoints(n))
    
    # 计算观察的 -log10(p) 值
    observed <- -log10(sort(pvalues))
    
    # 创建数据框
    df <- data.frame(Expected = expected, Observed = observed)
    
    # 绘制QQ-plot
    ggplot(df, aes(x = Expected, y = Observed)) +
      geom_point(size = 1.5, color = "black") + # 绘制散点
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + # 45度参考线
      geom_smooth(method = "loess", color = "red", fill = "gray", se = TRUE) + # 平滑曲线和阴影
      labs(title = "Mean vigilance", 
           x = "Expected", 
           y = expression(-log[10](p))) + 
      annotate("text", x = 1, y = 8, label = "lambda = 1.0042 // red: maf < 0.05") + # 添加文本
      theme_minimal()
  }
  
  # 调用函数并传入p值列
  qq_plot(gwas_results$P)
  
}

#对FUMA的结果进行分析
{
  #annot.txt
  {
    annot.txt<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/02.FUMA/annot.txt",header = T)
    # CADD_Phred分值中，10表示score排名在前10%，20表示前1%，30表示前0.1%，因此，分值要求越低，能保留下来的位点越多。
    # 对于SNP，CADD作者建议CADD_Phred分值>15，文章中通常用10或15；InDel没有建议值。
  }
  #annov.stats.txt
  {
    annov.stats.txt<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/02.FUMA/annov.stats.txt",header = T)
  }
  #annov.txt
  {
    annov.txt<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/02.FUMA/annov.txt",header = T)
    annov.txt_filter<-annov.txt[which(annov.txt$annot=="intergenic"|annov.txt$annot=="downstream"),]
  }
  #snps.txt
  {
    snps.txt<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/02.FUMA/snps.txt",header = T)
  }
  
  {
    snps.txt<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/02.FUMA/gwascatalog.txt",header = T,sep = "\t")
  }
}
#epistasis分析
{
  #对于joint和boost，是使用fast-epistasis进行的，需要将fam文件的最后一列变为case和control
  pheno<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.pheno",header = F)
  names(pheno)<-c("FID","IID","pheno")
  fam<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/11.epistasis/01.data/siginificant_0.05_.genotype_500_50_0.2.prune.in.fam",header = F)
  names(fam)<-c("FID","IID","V3","V4","V5","V6")
  
  fam_add_pheno<-inner_join(fam,pheno)
  head(fam_add_pheno)
  fam_add_pheno_filter<-fam_add_pheno[,c(1,2,3,4,5,7)]
  table(fam_add_pheno_filter$pheno)
  #将gcta中case和control的编码修改为plink中的case和control编码
  #plink中，case用2表示，control用1表示
  fam_add_pheno_filter$pheno[fam_add_pheno_filter$pheno=='1']<-2
  fam_add_pheno_filter$pheno[fam_add_pheno_filter$pheno=='0']<-1
  write.table(fam_add_pheno_filter, file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC/11.epistasis/01.data/siginificant_0.05_.genotype_500_50_0.2.prune.in.fam", 
              quote = FALSE, row.names = FALSE,sep = "\t")
}