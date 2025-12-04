#癌症饼图
{
  rm(list=ls())
  library(data.table)
  library(dplyr)
  library(ggplot2)
  #data<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type.txt",sep = "\t",header=F,na.strings = c("", "NA", "N/A", "NULL"))
  
  cancertype<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type.csv",sep = ",",na.strings = "")
  
  data_new<-cancertype
  
  
  {
    data_new1<-data_new[,c(1,2)] %>% as.data.frame()
    data_new1 <- na.omit(data_new1)
    
    # 对剩余的值进行计数
    count_result <- table(data_new1$V2)
    
    # 将计数结果转换为数据框
    result_df <- data.frame(
      Type = names(count_result),
      Count = as.vector(count_result)
    )
    
    result_df1<- result_df %>% group_by(Type) %>% mutate(per = scales::percent(Count / sum(result_df$Count)))
    
    cancertype<-c("breast", "prostate", "colorectal", "melanoma", "cervix", 
                  "lung", "lymphoidNeoplasms", "urinaryBladder", "kidney", 
                  "myeloidNeoplasms", "uterus", "headAndNeck", "esophagus", 
                  "ovary", "pancreas", "brain", "liver", "thyroid", "testis", 
                  "otherFemaleGenital", "tCellAndNKCellNeoplasms", "stomach", 
                  "eyeAndOrbit", "gallbladderAndBiliaryTract", "anus", 
                  "otherDigestive", "otherRespiratory", "softTissueSarcoma", 
                  "mesothelioma", "bone", "otherUrinaryOrgans", "otherMaleGenital", 
                  "otherEndocrine", "smallIntestine", "otherNervousSystem", "kaposiSarcoma")
    count<-c(21072, 13554, 9345, 6735, 6066, 4657, 3717, 2779, 2255, 1841,
             1839, 1716, 1607, 1339, 1305, 1154, 760, 642, 634, 598, 566,
             501, 490, 487, 472, 459, 434, 418, 408, 291, 226, 206, 179,
             174, 84, 19)
    
    
    df <- data.frame(x=cancertype,y=count)
    df$x <- factor(df$x,levels=c("breast", "prostate", "colorectal", "melanoma", "cervix", 
                                 "lung", "lymphoidNeoplasms", "urinaryBladder", "kidney", 
                                 "myeloidNeoplasms", "uterus", "headAndNeck", "esophagus", 
                                 "ovary", "pancreas", "brain", "liver", "thyroid", "testis", 
                                 "otherFemaleGenital", "tCellAndNKCellNeoplasms", "stomach", 
                                 "eyeAndOrbit", "gallbladderAndBiliaryTract", "anus", 
                                 "otherDigestive", "otherRespiratory", "softTissueSarcoma", 
                                 "mesothelioma", "bone", "otherUrinaryOrgans", "otherMaleGenital", 
                                 "otherEndocrine", "smallIntestine", "otherNervousSystem", "kaposiSarcoma"))
    
    #pdf("/home/wuxh/raw_data/PAS/streme_motif_zhifang.pdf",5,4)
    ggplot(data=df, aes(x=x, y=y, fill=cancertype))+
      geom_bar(stat="identity")+
      coord_flip()+ #把图像横向的参数
      guides(fill= F)+ #调节图例隐藏的参数。
      theme_bw()+
      labs(title = "Stage5")+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      xlab("count") +
      ylab("motif") 
    #dev.off()
    lab1 =as.vector(df$x)
    lab1
    lab1 =paste(lab1, "(",df$y,")",sep = "")
    lab1
    lab1 =paste(lab1, "(", round(df$y/sum(df$y)*100, 2), "%)",sep = "")
    lab1
    #生成仅有“百分比”的数据标签，并将数值小于5%的标签设置为“空”，避免占比较小的区域标签互相重叠；
    lab2<-round(df$y/sum(df$y)*100,1)
    lab2
    n<-length(lab2)
    n
    for(i in 1: n){  
      if(as.numeric(lab2[i])<2)    
        lab2[i]<-""  
      else    
        lab2[i]<-paste(lab2[i],"%",sep= "")  
    }
    lab2
    
    pdf("/home/wuxh/luohh/MPC_stage5.pdf",10,10)
    p1<-ggplot(data=df,aes(x="",y=y,fill=(x)))+geom_bar(stat="identity",width=0.4,color="white",linetype=1,size=1)
    p2<-p1+geom_text(aes(x=1,label=lab2),position= position_stack(reverse =F,vjust=0.5),size=6)
    p3<-ggplot(data=df,aes(x="",y=y,fill=x))+ geom_bar(stat="identity",width=1,color="white",linetype=1,size=1)+coord_polar(theta="y") +labs(x="",y="",title="")
    p3<-p3+geom_text(aes(x=1.25,label=lab2),position= position_stack(reverse =F,vjust=0.5),size=4)+labs(title = "Stage5")+
      theme(plot.title = element_text(hjust = 0.5))
    p4<-p3+theme_bw() + theme(axis.text=element_blank(),panel.border=element_blank(),axis.ticks=element_blank(), panel.grid=element_blank(),legend.title= element_blank(), legend.position = "right") +
      scale_fill_discrete(breaks= df$x, labels = lab1)
    p4
    dev.off()
  }
  data_new2<-data_new[,c(1,3)] %>% as.data.frame()
  
  data_new3<-data_new[,c(1,4)] %>% as.data.frame()
  
  data_new4<-data_new[,c(1,5)] %>% as.data.frame()
  
  data_new5<-data_new[,c(1,6)] %>% as.data.frame()
  
  data_new6<-data_new[,c(1,7)] %>% as.data.frame()
  
  
}
data_total<-rbind(data_new1,data_new2,data_new3,data_new4,data_new5,data_new6) %>% unique()
names(data_total)<-"Cancer_Type"
data_total1<-data_total[!is.na(data_total$Cancer_Type),] %>% as.data.frame()
names(data_total1)<-"Cancer_Type"

data<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type_gte.txt",sep = "\t",header=T,na.strings = c("", "NA", "N/A", "NULL"))

write.table(data,
            "/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type_gte.txt",
            row.names = F,col.names = F,sep = "\t",quote = FALSE,na="")





#癌症顺序分析
{
  data<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type.txt",sep = "\t",header=T,na.strings = c("", "NA", "N/A", "NULL"))
  count_dayu_4<-data[!is.na(data$V5),]
  #394个样本
  #两种癌症
  {
    # 创建包含36种癌症的数据框
    # 初始化一个空的数据框来存储组合
    new_data <- data.frame()
    
    # 使用循环生成每种癌症与其他癌症的两两组合
    for (i in 1:nrow(data_total1)) {
      for (j in 1:nrow(data_total1)) {
        if (i != j) {
          combination <- data.frame(Cancer1 = data_total1$Cancer_Type[i],
                                    Cancer2 = data_total1$Cancer_Type[j])
          new_data <- rbind(new_data, combination)
        }
      }
    }
    # 打印新的数据框
    print(new_data)
    data_new1<-data_new[!is.na(data_new$V6),] %>% as.data.frame()
    #两种癌症类型的有8540个样本信息，三种癌症的有702个样本信息
    
    data_new<-data_new1
    # 初始化一个空的列表，用于存储结果和统计数据
    result_list <- list()
    count_list <- list()
    
    # 遍历每种癌症的组合
    for (j in 1:nrow(new_data)) {
      combination <- new_data[j, ]
      
      # 初始化一个空的数据框来存储筛选后的结果
      result <- data.frame()
      
      # 对data_new中的每一行进行筛选
      for (i in 1:nrow(data_new)) {
        row <- data_new[i, ]
        
        # 检查该行是否包含当前组合的两种癌症
        if (all(combination$Cancer1 %in% row) && all(combination$Cancer2 %in% row)) {
          result <- rbind(result, row)
        }
      }
      # 去除result中的重复行
      # 仅当结果行数大于10时，才输出结果文件并计数
      if (nrow(result) > 10) {
        # 去除result中的重复行
        result <- unique(result)
        
        # 计算去重后的行数
        count <- nrow(result)
        
        # 生成文件名
        filename <- paste0(paste(combination$Cancer1, combination$Cancer2, sep="_"), "_", count, ".csv")
        
        # 保存结果到文件
        write.table(result[,-1], file.path("/home/wuxh/luohh/count2/", filename), row.names = FALSE,col.names = F,sep = "\t",quote = FALSE,na="")
        # 存储结果和统计数据到列表中
        result_list[[j]] <- result
        count_list[[j]] <- count
      }
    }
    
    # 打印统计数据
    count_list<-unlist(count_list) %>% as.data.frame()
    write.csv(count_list,"/home/wuxh/luohh/count2/00total_count", row.names = FALSE)
    
  }
  #三种癌症
  {
    # 初始化一个空的数据框来存储组合
    new_data <- data.frame()
    
    # 使用循环生成每种癌症与其他癌症的三种组合
    for (i in 1:nrow(data_total1)) {
      for (j in 1:nrow(data_total1)) {
        for (k in 1:nrow(data_total1)) {
          if (i != j && j != k && i != k) {
            combination <- data.frame(Cancer1 = data_total1$Cancer_Type[i],
                                      Cancer2 = data_total1$Cancer_Type[j],
                                      Cancer3 = data_total1$Cancer_Type[k])
            new_data <- rbind(new_data, combination)
          }
        }
      }
    }
    
    print(new_data)
    data_new1<-data_new[!is.na(data_new$V6),] %>% as.data.frame()
    #两种癌症类型的有8540个样本信息，三种癌症的有702个样本信息
    
    data_new<-data_new1
    # 初始化一个空的列表，用于存储结果和统计数据
    result_list <- list()
    count_list <- list()
    
    # 遍历每种癌症的三种组合
    for (j in 1:nrow(new_data)) {
      combination <- new_data[j, ]
      
      # 初始化一个空的数据框来存储筛选后的结果
      result <- data.frame()
      
      # 对data_new中的每一行进行筛选
      for (i in 1:nrow(data_new)) {
        row <- data_new[i, ]
        
        # 检查该行是否包含当前组合的三种癌症
        if (all(combination$Cancer1 %in% row) &&
            all(combination$Cancer2 %in% row) &&
            all(combination$Cancer3 %in% row)) {
          result <- rbind(result, row)
        }
      }
      
      # 仅当结果行数大于5时，才输出结果文件并计数
      if (nrow(result) > 2) {
        # 去除result中的重复行
        result <- unique(result)
        
        # 计算去重后的行数
        count <- nrow(result)
        
        # 生成文件名
        filename <- paste0(paste(combination$Cancer1, combination$Cancer2, combination$Cancer3, sep="_"), "_", count, ".csv")
        
        # 保存结果到文件
        write.table(result[,-1], file.path("/home/wuxh/luohh/count3/", filename), row.names = FALSE,col.names = F,sep = "\t",quote = FALSE,na="")
        
        # 存储结果和统计数据到列表中
        result_list[[j]] <- result
        count_list[[j]] <- count
      }
    }
    
    # 打印统计数据
    count_list
    count_list<-unlist(count_list) %>% as.data.frame()
    write.csv(count_list,"/home/wuxh/luohh/count3/00total_count", row.names = FALSE)
  }
  {
    
    # 初始化一个空的数据框来存储组合
    new_data <- data.frame()
    
    # 使用循环生成每种癌症与其他癌症的四种组合
    for (i in 1:5) {
      for (j in 1:5) {
        for (k in 1:5) {
          for (l in 1:5) {
            if (i != j && j != k && k != l && i != k && i != l && j != l) {
              combination <- data.frame(Cancer1 = data$Cancer_Type[i],
                                        Cancer2 = data$Cancer_Type[j],
                                        Cancer3 = data$Cancer_Type[k],
                                        Cancer4 = data$Cancer_Type[l])
              new_data <- rbind(new_data, combination)
            }
          }
        }
      }
    }
    
    
    new_data<-unite(new_data,full,c("Cancer1","Cancer2","Cancer3","Cancer4"),sep = "_")
    
    
    new_data<-separate(new_data,"full",into=c("Cancer1","Cancer2","Cancer3","Cancer4"),sep = "_")
    
    # 找到重复的行
    duplicate_rows <- duplicated(new_data_sep) | duplicated(new_data_sep, fromLast = TRUE)
    
    # 使用subset函数过滤掉重复的行
    filtered_df <- subset(new_data_sep, !duplicate_rows)
    
    # 显示结果
    print(filtered_df)
    
    
    
    
    # 打印新的数据框
    print(new_data)
    
    # 假设data_new是你的大数据框
    # 假设已经生成了combinations和new_data（前面的生成方式）
    data_new<-data_new1
    # 初始化一个空的列表，用于存储结果和统计数据
    result_list <- list()
    count_list <- list()
    
    # 遍历每种癌症的四种组合
    for (j in 1:nrow(new_data)) {
      combination <- new_data[j, ]
      
      # 初始化一个空的数据框来存储筛选后的结果
      result <- data.frame()
      
      # 对data_new中的每一行进行筛选
      for (i in 1:nrow(data_new)) {
        row <- data_new[i, ]
        
        # 检查该行是否包含当前组合的四种癌症
        if (all(combination$Cancer1 %in% row) &&
            all(combination$Cancer2 %in% row) &&
            all(combination$Cancer3 %in% row) &&
            all(combination$Cancer4 %in% row)) {
          result <- rbind(result, row)
        }
      }
      
      # 仅当结果行数大于2时，才输出结果文件并计数
      if (nrow(result) > 2) {
        # 去除result中的重复行
        result <- unique(result)
        
        # 计算去重后的行数
        count <- nrow(result)
        
        # 生成文件名
        filename <- paste0(paste(combination$Cancer1, combination$Cancer2, combination$Cancer3, combination$Cancer4, sep="_"), "_", count, ".csv")
        
        # 保存结果到文件
        write.table(result[,-1], file.path("/home/wuxh/luohh/count4/", filename), row.names = FALSE,col.names = F,sep = "\t",quote = FALSE,na="")
        # 存储结果和统计数据到列表中
        result_list[[j]] <- result
        count_list[[j]] <- count
      }
    }
    
    # 打印统计数据
    count_list<-unlist(count_list) %>% as.data.frame()
    write.csv(count_list,"/home/wuxh/luohh/count3/00total_count", row.names = FALSE)
  }
}


data_new1<-data_new[!is.na(data_new$V6),] %>% as.data.frame()

#两种癌症类型的有8540个样本信息，三种癌症的有702个样本信息

# 初始化一个新的数据框来存储匹配的行
matched_rows <- data.frame()

# 使用apply函数遍历每一行
count<-apply(data_new1, 1, function(row) {
  if ("smallIntestine" %in% row && "breast" %in% row) {
    # 如果同时包含"melanoma"和"lung"，则将整行添加到新数据框
    matched_rows <<- rbind(matched_rows, row)
  }
})

nrow(matched_rows)

matched_rows<-matched_rows[,-1]

write.table(matched_rows,
            "/home/wuxh/luohh/count/count",
            row.names = F,col.names = F,sep = "\t",quote = FALSE,na="")

