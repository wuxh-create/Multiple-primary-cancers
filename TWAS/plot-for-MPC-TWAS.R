


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




