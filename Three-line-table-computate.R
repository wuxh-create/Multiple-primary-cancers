


#流行病学分析，绘制三线表

{
  rm(list=ls())
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  {
    basic<-fread("/home/luohh/UKB50wData/04.PhenoDataExtract/basic.pheno")
    #502387
    #选择两个性别一致的样本
    basic$sex_state<-ifelse(basic$sex==basic$geneticSex,"Y","N")
    basic<-basic %>% filter(sex_state=="Y")
    #487775
    basic_1001<-basic %>% filter(population=="1001")
    #430606
    #获得癌症首发年龄=诊断日期-出生年月
    clinic<-fread("/home/wuxh/SSIdb/04.UKBB/phenodata/clinic.txt",na.strings = c("", "NA", "N/A", "NULL"))
    #40005是癌症诊断日期
    #40006是ICD10编码的癌症类型
    #40008是癌症诊断年龄
    #40009是自我报告癌症发生
    #40010是自爆导致的死亡时间
    #40011是癌症肿瘤的组织学信息
    #40012癌症肿瘤的行为
    a<-names(clinic)
    clinic_select<-clinic[,2:124]
    names(clinic_select)<-a[1:123]
    #将诊断时间缺失的数据去除
    clinic_select<-clinic_select[,c(1,3,10)]
    clinic_select<-separate(clinic_select,"40005-0.0",into="nian",sep = "-")
    str(clinic_select$nian)
    str(clinic_select$`34-0.0`)
    clinic_select$nian<-as.numeric(as.character(clinic_select$nian))
    clinic_select$diagnose<-clinic_select$nian-clinic_select$`34-0.0`
    
    clinic_select<-clinic_select[!is.na(clinic_select$diagnose),]
    names(clinic_select)[1]<-"sampleID"
    #118955个数据有首发时间
    
    #统计癌症类型，筛选只发生一种多原发癌和发生大于等于三种和大于等于四种的多原发癌
    
    cancertype<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type.csv",sep = ",",na.strings = "")
    #筛选只有一种癌症类型的样本
    #将自己的统计结果和浩辉师兄的结果取差集就是单癌的样本，然后两种，三种，四种的MPC样本信息，都是用浩辉师兄所有的信息。
    
    
    multiPCancer.sample<-fread("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/multiPCancer.sample")
    
    single<-anti_join(cancertype,multiPCancer.sample)  
    names(single)[1]<-names(basic_1001)[1]
    single_1001<-inner_join(single,basic_1001) 
    #73650个单癌样本
    #添加首发时间
    single_1001_diagnose<-inner_join(single_1001,clinic_select)
    single_1001_diagnose<-single_1001_diagnose[,c(1,2,8,17)]
    names(single_1001_diagnose)<-c("sampleID","cancer","sex","diagnose")
    
    
    
    
    all_MPC<-inner_join(multiPCancer.sample,cancertype)
    names(all_MPC)[1]<-names(basic_1001)[1]
    all_MPC_1001<-inner_join(all_MPC,basic_1001) 
    #73650个单癌样本
    #添加首发时间
    all_MPC_1001_diagnose<-inner_join(all_MPC_1001,clinic_select)
    all_MPC_1001_diagnose<-all_MPC_1001_diagnose[,c(1:7,8,17)]
    names(all_MPC_1001_diagnose)<-c("sampleID","cancer1","cancer2","cancer3","cancer4","cancer5","cancer6","sex","diagnose")
    {
      anus<-all_MPC_1001_diagnose[which(all_MPC_1001_diagnose$cancer1=="anus"|
                                          all_MPC_1001_diagnose$cancer2=="anus"|
                                          all_MPC_1001_diagnose$cancer3=="anus"|
                                          all_MPC_1001_diagnose$cancer4=="anus"|
                                          all_MPC_1001_diagnose$cancer5=="anus"|
                                          all_MPC_1001_diagnose$cancer5=="anus"),]
    }
    
    
    mpc3<-fread("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.3cancers.sample")
    
    mpc4<-fread("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.4cancers.sample")
    
    mpc5<-fread("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.5cancers.sample")
    
    mpc6<-fread("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.6cancers.sample")
    
    MPC3<-rbind(mpc3,mpc4,mpc5,mpc6)
    MPC3<-inner_join(MPC3,cancertype)
    names(MPC3)[1] <- names(basic_1001)[1]
    MPC3_1001 <- inner_join(MPC3, basic_1001) 
    # 73650个单癌样本
    # 添加首发时间
    MPC3_1001_diagnose <- inner_join(MPC3_1001, clinic_select)
    MPC3_1001_diagnose <- MPC3_1001_diagnose[, c(1:7, 8, 17)]
    names(MPC3_1001_diagnose) <- c("sampleID", "cancer1", "cancer2", "cancer3", "cancer4", "cancer5", "cancer6", "sex", "diagnose")
    
    
    MPC4<-rbind(mpc4,mpc5,mpc6)
    MPC4<-inner_join(MPC4,cancertype)
    names(MPC4)[1] <- names(basic_1001)[1]
    MPC4_1001 <- inner_join(MPC4, basic_1001) 
    #73650个单癌样本
    #添加首发时间
    MPC4_1001_diagnose <- inner_join(MPC4_1001, clinic_select)
    MPC4_1001_diagnose <- MPC4_1001_diagnose[, c(1:7, 8, 17)]
    names(MPC4_1001_diagnose) <- c("sampleID", "cancer1", "cancer2", "cancer3", "cancer4", "cancer5", "cancer6", "sex", "diagnose")
    
  }
  #bmi
  {
    bmi<-fread("/home/wuxh/SSIdb/04.UKBB/phenodata/bmi_res",fill=TRUE)%>% as.data.frame() 
    names(bmi)[1]<-"sampleID"
    library(dplyr)
    # 对数据框进行处理
    new_bmi <- bmi %>%
      mutate(`23104-2.0` := ifelse(!is.na(`23104-3.0`), `23104-3.0`, `23104-2.0`),
             `23104-1.0` := ifelse(!is.na(`23104-2.0`), `23104-2.0`, `23104-1.0`),
             `23104-0.0` := ifelse(!is.na(`23104-1.0`), `23104-1.0`, `23104-0.0`)) 
    new_bmi <- bmi %>%
      mutate(`21001-2.0` := ifelse(!is.na(`21001-3.0`), `21001-3.0`, `21001-2.0`),
             `21001-1.0` := ifelse(!is.na(`21001-2.0`), `21001-2.0`, `21001-1.0`),
             `21001-0.0` := ifelse(!is.na(`21001-1.0`), `21001-1.0`, `21001-0.0`)) 
    new_bmi_select<-new_bmi[,c(1,2)]
    names(new_bmi_select)[2]<-"BMI"
    #BMI和基本信息合并
    
    single_1001_diagnose_bim<-merge(single_1001_diagnose,new_bmi_select,all.x = T)
    all_MPC_1001_diagnose_bim<-merge(all_MPC_1001_diagnose,new_bmi_select,all.x = T)
    MPC3_1001_diagnose_bim<-merge(MPC3_1001_diagnose,new_bmi_select,all.x = T)
    MPC4_1001_diagnose_bim<-merge(MPC4_1001_diagnose,new_bmi_select,all.x = T)
  }
  #smoking
  {
    smok<-fread("/home/wuxh/SSIdb/04.UKBB/phenodata/20116_smoking_res")
    new_smok <- smok %>%
      mutate(`20116-2.0` := ifelse(!is.na(`20116-3.0`), `20116-3.0`, `20116-2.0`),
             `20116-1.0` := ifelse(!is.na(`20116-2.0`), `20116-2.0`, `20116-1.0`),
             `20116-0.0` := ifelse(!is.na(`20116-1.0`), `20116-1.0`, `20116-0.0`)) 
    new_smok_select<-new_smok[,c(1,2)]
    names(new_smok_select)<-c("sampleID","smoking")
    table(new_smok_select$smoking)
    #将没有吸烟状态的样本去掉
    #bmi_basic_1001_remove_na_smok1<-bmi_basic_1001_remove_na_smok[!is.na(bmi_basic_1001_remove_na_smok$smok),]
    new_smok_select$smoking[new_smok_select$smoking==-3]<-"Prefer not to answer"
    new_smok_select$smoking[new_smok_select$smoking==0]<-"Never"
    new_smok_select$smoking[new_smok_select$smoking==1]<-"Former"
    new_smok_select$smoking[new_smok_select$smoking==2]<-"Current"
    single_1001_diagnose_bim_smok<-merge(single_1001_diagnose_bim,new_smok_select,all.x = T)
    all_MPC_1001_diagnose_bim_smok<-merge(all_MPC_1001_diagnose_bim,new_smok_select,all.x = T)
    MPC3_1001_diagnose_bim_smok<-merge(MPC3_1001_diagnose_bim,new_smok_select,all.x = T)
    MPC4_1001_diagnose_bim_smok<-merge(MPC4_1001_diagnose_bim,new_smok_select,all.x = T)
  }
  #alcohol
  {
    alcohol<-fread("/home/wuxh/SSIdb/04.UKBB/phenodata/20117_drink_alcohol_res") %>% unique() %>% as.data.frame()
    names(alcohol)<-c("sampleID","20117-0.0","20117-1.0","20117-2.0","20117-3.0")
    new_alcohol <- alcohol %>%
      mutate(`20117-2.0` := ifelse(!is.na(`20117-3.0`), `20117-3.0`, `20117-2.0`),
             `20117-1.0` := ifelse(!is.na(`20117-2.0`), `20117-2.0`, `20117-1.0`),
             `20117-0.0` := ifelse(!is.na(`20117-1.0`), `20117-1.0`, `20117-0.0`)) 
    new_alcohol_select<-new_alcohol[,c(1,2)]
    names(new_alcohol_select)<-c("sampleID","alcohol")
    
    table(new_alcohol_select$alcohol)
    
    #bmi_basic_1001_remove_na_smok_alcohol1<-bmi_basic_1001_remove_na_smok_alcohol[!is.na(bmi_basic_1001_remove_na_smok_alcohol$alcohol),]
    new_alcohol_select$alcohol[new_alcohol_select$alcohol==-3]<-"Prefer not to answer"
    new_alcohol_select$alcohol[new_alcohol_select$alcohol==0]<-"Never"
    new_alcohol_select$alcohol[new_alcohol_select$alcohol==1]<-"Former"
    new_alcohol_select$alcohol[new_alcohol_select$alcohol==2]<-"Current"
    
    single_1001_diagnose_bim_smok_alcohol<-merge(single_1001_diagnose_bim_smok,new_alcohol_select,all.x = T)
    all_MPC_1001_diagnose_bim_smok_alcohol<-merge(all_MPC_1001_diagnose_bim_smok,new_alcohol_select,all.x = T)
    MPC3_1001_diagnose_bim_smok_alcohol<-merge(MPC3_1001_diagnose_bim_smok,new_alcohol_select,all.x = T)
    MPC4_1001_diagnose_bim_smok_alcohol<-merge(MPC4_1001_diagnose_bim_smok,new_alcohol_select,all.x = T)
    
  }
  
  #样本提取
{
  library(stringr)
  library(tidyr)
  library(dplyr)
  # 设置工作目录为包含所需文件的文件夹
  setwd("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/04.cancerGroup/03.whiteData/")
  # 列出文件夹下所有以 'filtered.white.sample' 结尾的文件
  files <- list.files(pattern = "filtered\\.white\\.sample$", full.names = TRUE)
  # 读取所有文件并将它们合并成一个大的数据框
  data_list <- lapply(files, read.table, header = FALSE) # 假设文件有表头，如果没有，请设置header = FALSE
  all_cancer <- do.call(rbind, data_list) %>% unique()
  #去重后89029
  setwd("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/")
  # 列出文件夹下所有以 'filtered.white.sample' 结尾的文件
  files <- list.files(pattern = "^MPC.*cancers\\.sample$")
  # 读取所有文件并将它们合并成一个大的数据框
  data_list <- lapply(files, read.table, header = FALSE) # 假设文件有表头，如果没有，请设置header = FALSE
  multi_cancer <- do.call(rbind, data_list) %>% unique()
  #去重后14290 
  single_cancer<-anti_join(all_cancer,multi_cancer,by="V1")
  #74739
  names(multi_cancer)<-names(single_cancer)<-"sampleID"
  
  multi_cancer3<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.3cancers.sample")
  multi_cancer4<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.4cancers.sample")
  multi_cancer5<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.5cancers.sample")
  multi_cancer6<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer/MPC.6cancers.sample")
  
  multi_cancer_more3<-rbind(multi_cancer3,multi_cancer4,multi_cancer5,multi_cancer6)
  multi_cancer_more4<-rbind(multi_cancer4,multi_cancer5,multi_cancer6)
  names(multi_cancer_more3)<-names(multi_cancer_more4)<-names(single_cancer)<-"sampleID"
  
}
{
  library(tableone)
  #install.packages("table1")
  library(table1)
  #####################################################################
  {
    pvalue <- function(x, ...) {
      y <- unlist(x)
      g <- factor(rep(1:length(x), times=sapply(x, length)))
      if (is.numeric(y)) {
        p <- t.test(y ~ g)$p.value # 数值型数据用 t-test (两组比较)
      } else {
        p <- chisq.test(table(y, g))$p.value # 因子型数据用卡方
      }
      # 格式化 P 值以显示所有小数点
      p_formatted <- format(p, digits = 22)
      # 使用 HTML 实体表示小于号
      c("", sub("<", "&lt;", p_formatted))
    }
  }
  
  
  single_1001_diagnose_bim_smok_alcohol$group<-rep(1,nrow(single_1001_diagnose_bim_smok_alcohol))
  all_MPC_1001_diagnose_bim_smok_alcohol$group<-rep(0,nrow(all_MPC_1001_diagnose_bim_smok_alcohol))
  MPC3_1001_diagnose_bim_smok_alcohol$group<-rep(0,nrow(MPC3_1001_diagnose_bim_smok_alcohol))
  MPC4_1001_diagnose_bim_smok_alcohol$group<-rep(0,nrow(MPC4_1001_diagnose_bim_smok_alcohol))
  
  #分癌症画表格---单癌
  {
    anus_single_1001_diagnose_bim_smok_alcohol <- single_1001_diagnose_bim_smok_alcohol[which(single_1001_diagnose_bim_smok_alcohol$cancer == "anus"), ]
    
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    
    for (cancer in cancer_types) {
      # Generate and execute the code line
      cancer="breast"
      code_line <- sprintf("%s_single_1001_diagnose_bim_smok_alcohol <- single_1001_diagnose_bim_smok_alcohol[which(single_1001_diagnose_bim_smok_alcohol$cancer == \"%s\"), ]",
                           cancer, cancer)
      eval(parse(text=code_line))
      
      # 新增代码：创建一个新列 'cancer_sex'，将 'cancer' 和 'sex' 列的值用“-”连接起来
      combine_line <- sprintf("%s_single_1001_diagnose_bim_smok_alcohol$cancer_sex <- paste(%s_single_1001_diagnose_bim_smok_alcohol$cancer, %s_single_1001_diagnose_bim_smok_alcohol$sex, sep=\"-\")", cancer, cancer, cancer)
      eval(parse(text=combine_line))
      
      select_columns_line <- sprintf("%s_single_1001_diagnose_bim_smok_alcohol_selected <- %s_single_1001_diagnose_bim_smok_alcohol[, c(1, 8,9)]", cancer, cancer)
      eval(parse(text=select_columns_line))
    }
  }
  #分癌症画表格---all癌
  {
    anus_all_MPC_1001_diagnose_bim_smok_alcohol <- all_MPC_1001_diagnose_bim_smok_alcohol[which(all_MPC_1001_diagnose_bim_smok_alcohol$cancer1 == "anus" |
                                                                                                  all_MPC_1001_diagnose_bim_smok_alcohol$cancer2 == "anus" |
                                                                                                  all_MPC_1001_diagnose_bim_smok_alcohol$cancer3 == "anus" |
                                                                                                  all_MPC_1001_diagnose_bim_smok_alcohol$cancer4 == "anus" |
                                                                                                  all_MPC_1001_diagnose_bim_smok_alcohol$cancer5 == "anus" |
                                                                                                  all_MPC_1001_diagnose_bim_smok_alcohol$cancer6 == "anus"), ]
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    
    for (cancer in cancer_types) {
      # Generate and execute the code line
      
      cancer="breast"
      
      code_line <- sprintf("%s_all_MPC_1001_diagnose_bim_smok_alcohol <- all_MPC_1001_diagnose_bim_smok_alcohol[which(all_MPC_1001_diagnose_bim_smok_alcohol$cancer1 == \"%s\" |
                                                                                                      all_MPC_1001_diagnose_bim_smok_alcohol$cancer2 == \"%s\" |
                                                                                                      all_MPC_1001_diagnose_bim_smok_alcohol$cancer3 == \"%s\" |
                                                                                                      all_MPC_1001_diagnose_bim_smok_alcohol$cancer4 == \"%s\" |
                                                                                                      all_MPC_1001_diagnose_bim_smok_alcohol$cancer5 == \"%s\" |
                                                                                                      all_MPC_1001_diagnose_bim_smok_alcohol$cancer6 == \"%s\"), ]",
                           cancer, cancer, cancer, cancer, cancer, cancer, cancer)
      eval(parse(text=code_line))
      
      # Add the new column 'cancer' to the data frame
      new_col_line <- sprintf("%s_all_MPC_1001_diagnose_bim_smok_alcohol$cancer <- \"%s\"", cancer, cancer)
      eval(parse(text=new_col_line))
      # 新增代码：创建一个新列 'cancer_sex'，将 'cancer' 和 'sex' 列的值用“-”连接起来
      combine_line <- sprintf("%s_all_MPC_1001_diagnose_bim_smok_alcohol$cancer_sex <- paste(%s_all_MPC_1001_diagnose_bim_smok_alcohol$cancer, %s_all_MPC_1001_diagnose_bim_smok_alcohol$sex, sep=\"-\")", cancer, cancer, cancer)
      eval(parse(text=combine_line))
      select_columns_line <- sprintf("%s_all_MPC_1001_diagnose_bim_smok_alcohol_selected <- %s_all_MPC_1001_diagnose_bim_smok_alcohol[, c(1, 13, 15)]", cancer, cancer)
      eval(parse(text=select_columns_line))
    }
    
  }
  #分癌症画表格---MPC ≥ 3癌
  {
    anus_MPC3_1001_diagnose_bim_smok_alcohol <- MPC3_1001_diagnose_bim_smok_alcohol[which(MPC3_1001_diagnose_bim_smok_alcohol$cancer1 == "anus" |
                                                                                            MPC3_1001_diagnose_bim_smok_alcohol$cancer2 == "anus" |
                                                                                            MPC3_1001_diagnose_bim_smok_alcohol$cancer3 == "anus" |
                                                                                            MPC3_1001_diagnose_bim_smok_alcohol$cancer4 == "anus" |
                                                                                            MPC3_1001_diagnose_bim_smok_alcohol$cancer5 == "anus" |
                                                                                            MPC3_1001_diagnose_bim_smok_alcohol$cancer6 == "anus"), ]
    
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    
    for (cancer in cancer_types) {
      # Generate and execute the code line
      
      cancer="breast"
      
      code_line <- sprintf("%s_MPC3_1001_diagnose_bim_smok_alcohol <- MPC3_1001_diagnose_bim_smok_alcohol[which(MPC3_1001_diagnose_bim_smok_alcohol$cancer1 == \"%s\" |
                                                                                                        MPC3_1001_diagnose_bim_smok_alcohol$cancer2 == \"%s\" |
                                                                                                        MPC3_1001_diagnose_bim_smok_alcohol$cancer3 == \"%s\" |
                                                                                                        MPC3_1001_diagnose_bim_smok_alcohol$cancer4 == \"%s\" |
                                                                                                        MPC3_1001_diagnose_bim_smok_alcohol$cancer5 == \"%s\" |
                                                                                                        MPC3_1001_diagnose_bim_smok_alcohol$cancer6 == \"%s\"), ]",
                           cancer, cancer, cancer, cancer, cancer, cancer, cancer)
      eval(parse(text=code_line))
      
      # Add the new column 'cancer' to the data frame
      new_col_line <- sprintf("%s_MPC3_1001_diagnose_bim_smok_alcohol$cancer <- \"%s\"", cancer, cancer)
      eval(parse(text=new_col_line))
      # 新增代码：创建一个新列 'cancer_sex'，将 'cancer' 和 'sex' 列的值用“-”连接起来
      combine_line <- sprintf("%s_MPC3_1001_diagnose_bim_smok_alcohol$cancer_sex <- paste(%s_MPC3_1001_diagnose_bim_smok_alcohol$cancer, %s_MPC3_1001_diagnose_bim_smok_alcohol$sex, sep=\"-\")", cancer, cancer, cancer)
      eval(parse(text=combine_line))
      select_columns_line <- sprintf("%s_MPC3_1001_diagnose_bim_smok_alcohol_selected <- %s_MPC3_1001_diagnose_bim_smok_alcohol[, c(1, 13, 15)]", cancer, cancer)
      eval(parse(text=select_columns_line))
    }
    
  }
  #分癌症画表格---MPC ≥ 4癌
  {
    anus_MPC4_1001_diagnose_bim_smok_alcohol <- MPC4_1001_diagnose_bim_smok_alcohol[which(MPC4_1001_diagnose_bim_smok_alcohol$cancer1 == "anus" |
                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer2 == "anus" |
                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer3 == "anus" |
                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer4 == "anus" |
                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer5 == "anus" |
                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer6 == "anus"), ]
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    
    for (cancer in cancer_types) {
      # Generate and execute the code line
      
      cancer="breast"
      
      code_line <- sprintf("%s_MPC4_1001_diagnose_bim_smok_alcohol <- MPC4_1001_diagnose_bim_smok_alcohol[which(MPC4_1001_diagnose_bim_smok_alcohol$cancer1 == \"%s\" |
                                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer2 == \"%s\" |
                                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer3 == \"%s\" |
                                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer4 == \"%s\" |
                                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer5 == \"%s\" |
                                                                                                            MPC4_1001_diagnose_bim_smok_alcohol$cancer6 == \"%s\"), ]",
                           cancer, cancer, cancer, cancer, cancer, cancer, cancer)
      eval(parse(text=code_line))
      
      # Add the new column 'cancer' to the data frame
      new_col_line <- sprintf("%s_MPC4_1001_diagnose_bim_smok_alcohol$cancer <- \"%s\"", cancer, cancer)
      eval(parse(text=new_col_line))
      # 新增代码：创建一个新列 'cancer_sex'，将 'cancer' 和 'sex' 列的值用“-”连接起来
      combine_line <- sprintf("%s_MPC4_1001_diagnose_bim_smok_alcohol$cancer_sex <- paste(%s_MPC4_1001_diagnose_bim_smok_alcohol$cancer, %s_MPC4_1001_diagnose_bim_smok_alcohol$sex, sep=\"-\")", cancer, cancer, cancer)
      eval(parse(text=combine_line))
      select_columns_line <- sprintf("%s_MPC4_1001_diagnose_bim_smok_alcohol_selected <- %s_MPC4_1001_diagnose_bim_smok_alcohol[, c(1, 13, 15)]", cancer, cancer)
      eval(parse(text=select_columns_line))
    }
  }
  
  #分癌症画三线表:single_all
  {
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    for (cancer in cancer_types) {
      # 构造变量名称
      
      cancer<-"uterus"

      cancer_single_all <- paste(cancer, "single_all", sep = "_")
      
      # 动态构造并执行 R 代码
      eval(parse(text = paste(cancer_single_all, "<- rbind(", cancer, "_single_1001_diagnose_bim_smok_alcohol_selected, ", cancer, "_all_MPC_1001_diagnose_bim_smok_alcohol_selected)", sep = "")))
      eval(parse(text = paste(cancer_single_all, "$group <- factor(", cancer_single_all, "$group, levels = c(1, 0), labels = c('Single cancer', 'all MPC'))", sep = "")))
      
      uterus_single_all<-uterus_single_all[which(uterus_single_all$cancer_sex=="uterus-female"),]
      
      # 执行 table1 函数
      a <- eval(parse(text = paste("table1(~ cancer_sex | group, data = ", cancer_single_all, ", extra.col = list(`P-value` = pvalue), overall = F)", sep = "")))
      
      # 保存结果到 CSV 文件
      write.csv(as.data.frame(a), paste("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/", cancer_single_all, ".csv", sep = ""), row.names = F,quote=F)
    }
    # 设置工作目录为您的文件所在路径
    setwd("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/")
    
    # 获取所有以 _single_all.csv 结尾的文件名
    file_list <- list.files(pattern = "_single_all\\.csv$")
    
    # 读取并合并所有文件
    all_data <- do.call(rbind, lapply(file_list, read.csv, header = TRUE))
    
    # 如果需要，可以将合并后的数据保存为新的 CSV 文件
    write.csv(all_data, "single_all_combined_data.csv", row.names = FALSE)
    
  }
  breast 
  cervix
  #分癌症画三线表：single_mpc3
  {
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck", "kaposiSarcoma", "kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    for (cancer in cancer_types) {
      # 构造变量名称
      
      cancer<-"uterus"
      
      cancer_single_MPC3 <- paste(cancer, "single_MPC3", sep = "_")
      
      # 动态构造并执行 R 代码
      eval(parse(text = paste(cancer_single_MPC3, "<- rbind(", cancer, "_single_1001_diagnose_bim_smok_alcohol_selected, ", cancer, "_MPC3_1001_diagnose_bim_smok_alcohol_selected)", sep = "")))
      eval(parse(text = paste(cancer_single_MPC3, "$group <- factor(", cancer_single_MPC3, "$group, levels = c(1, 0), labels = c('Single cancer', 'MPC ≥ 3'))", sep = "")))
      
      uterus_single_MPC3<-uterus_single_MPC3[which(uterus_single_MPC3$cancer_sex=="uterus-female"),]
      
      # 执行 table1 函数
      a <- eval(parse(text = paste("table1(~ cancer_sex | group, data = ", cancer_single_MPC3, ", extra.col = list(`P-value` = pvalue), overall = F)", sep = "")))
      
      # 保存结果到 CSV 文件
      write.csv(as.data.frame(a), paste("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/", cancer_single_MPC3, ".csv", sep = ""), row.names = F,quote=F)
      
    }
    
    # 设置工作目录为您的文件所在路径
    setwd("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/")
    
    # 获取所有以 _single_all.csv 结尾的文件名
    file_list <- list.files(pattern = "_single_MPC3\\.csv$")
    
    # 读取并合并所有文件
    all_data <- do.call(rbind, lapply(file_list, read.csv, header = TRUE))
    
    # 如果需要，可以将合并后的数据保存为新的 CSV 文件
    write.csv(all_data, "single_MPC3_combined_data.csv", row.names = FALSE)
  }
  #分癌症画三线表：single_mpc4
  {
    cancer_types <- c("anus", "bone", "brain", "breast", "cervix",
                      "colorectal", "esophagus", "eyeAndOrbit", "gallbladderAndBiliaryTract", 
                      "headAndNeck","kidney", "liver", 
                      "lung", "lymphoidNeoplasms", "melanoma", "mesothelioma", 
                      "myeloidNeoplasms", "otherDigestive", "otherEndocrine", 
                      "otherFemaleGenital", "otherMaleGenital", "otherNervousSystem", 
                      "otherRespiratory", "otherUrinaryOrgans", "ovary", 
                      "pancreas", "prostate", "smallIntestine", "softTissueSarcoma", 
                      "stomach", "tCellAndNKCellNeoplasms", "testis", "thyroid", 
                      "urinaryBladder", "uterus")
    # "kaposiSarcoma"
    for (cancer in cancer_types) {
      # 构造变量名称
      
      cancer<-"uterus"
      
      cancer_single_MPC4 <- paste(cancer, "single_MPC4", sep = "_")
      
      # 动态构造并执行 R 代码
      eval(parse(text = paste(cancer_single_MPC4, "<- rbind(", cancer, "_single_1001_diagnose_bim_smok_alcohol_selected, ", cancer, "_MPC4_1001_diagnose_bim_smok_alcohol_selected)", sep = "")))
      eval(parse(text = paste(cancer_single_MPC4, "$group <- factor(", cancer_single_MPC4, "$group, levels = c(1, 0), labels = c('Single cancer', 'MPC ≥ 4'))", sep = "")))
      
      uterus_single_MPC4<-uterus_single_MPC4[which(uterus_single_MPC4$cancer_sex=="uterus-female"),]
      
      # 执行 table1 函数
      a <- eval(parse(text = paste("table1(~ cancer_sex | group, data = ", cancer_single_MPC4, ", extra.col = list(`P-value` = pvalue), overall = F)", sep = "")))
      
      # 保存结果到 CSV 文件
      write.csv(as.data.frame(a), paste("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/", cancer_single_MPC4, ".csv", sep = ""), row.names = F,quote=F)
    }
    # 设置工作目录为您的文件所在路径
    setwd("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/")
    
    # 获取所有以 _single_all.csv 结尾的文件名
    file_list <- list.files(pattern = "_single_MPC4\\.csv$")
    
    # 读取并合并所有文件
    all_data <- do.call(rbind, lapply(file_list, read.csv, header = TRUE))
    
    # 如果需要，可以将合并后的数据保存为新的 CSV 文件
    write.csv(all_data, "single_MPC4_combined_data.csv", row.names = FALSE)
  }
  
  {
    single_1001_diagnose_bim_smok_alcohol_select<-single_1001_diagnose_bim_smok_alcohol[,c(1,4,5,6,7,8)]
    all_MPC_1001_diagnose_bim_smok_alcohol_select<-all_MPC_1001_diagnose_bim_smok_alcohol[,c(1,9,10,11,12,13)]
    MPC3_1001_diagnose_bim_smok_alcohol_select<-MPC3_1001_diagnose_bim_smok_alcohol[,c(1,9,10,11,12,13)]
    MPC4_1001_diagnose_bim_smok_alcohol_select<-MPC4_1001_diagnose_bim_smok_alcohol[,c(1,9,10,11,12,13)]
    
    
    
    {
      single_and_all_basic_all<-rbind(single_1001_diagnose_bim_smok_alcohol_select,all_MPC_1001_diagnose_bim_smok_alcohol_select) %>% unique()
      #332548
      single_and_all_basic_all$group<- factor(single_and_all_basic_all$group,
                                              levels=c(1,0), labels=c("Single cancer", "all MPC"))
      
      
      single_and_all_basic_all1<-single_and_all_basic_all[!is.na(single_and_all_basic_all$BMI),]
      
      single_and_all_basic_all1<-single_and_all_basic_all1[!is.na(single_and_all_basic_all1$alcohol),]
      
      
      # 假设你的数据框叫做 df，并且有一个列名叫做 'age'
      # df <- ... # 你的数据框
      # 定义年龄区间
      breaks <- c(0, 30, 50,70,Inf)
      labels <- c("0-30", "30-50", "50-70","70+")
      # 使用cut函数进行年龄区间划分
      single_and_all_basic_all1$diagnose <- cut(single_and_all_basic_all1$diagnose, breaks = breaks, labels = labels, right = FALSE)
  
      library(dplyr)
      
      # 假设您的数据框是df，需要修改的列名是column_name
      single_and_all_basic_all1 <- single_and_all_basic_all1 %>%
        mutate(BMI = case_when(
          BMI < 18.5 ~ "Underweight(<18.5)",
          BMI >= 18.5 & BMI < 23.0 ~ "Normal(18.5-23.0)",
          BMI >= 23.0 & BMI < 25.0 ~ "Overweight(23.0-<25.0)",
          BMI >= 25 ~ "Obesity(≥25)",
          TRUE ~ as.character(BMI)  # 处理其他情况，例如NA值
        ))
      
      names(single_and_all_basic_all1)
      ######################################################################################
      table1(~ diagnose + BMI + smoking + alcohol|group,
             data=single_and_all_basic_all1,
             extra.col=list(`P-value`=pvalue),overall=F)
    }
    {
      single_and_MPC3_basic_all<-rbind(single_1001_diagnose_bim_smok_alcohol_select,MPC3_1001_diagnose_bim_smok_alcohol_select) %>% unique()
      #332548
      single_and_MPC3_basic_all$group<- factor(single_and_MPC3_basic_all$group,
                                              levels=c(1,0), labels=c("Single cancer", "MPC ≥ 3"))
      
      
      single_and_MPC3_basic_all1<-single_and_MPC3_basic_all[!is.na(single_and_MPC3_basic_all$BMI),]
      
      single_and_MPC3_basic_all1<-single_and_MPC3_basic_all1[!is.na(single_and_MPC3_basic_all1$alcohol),]
      
      
      # 假设你的数据框叫做 df，并且有一个列名叫做 'age'
      # df <- ... # 你的数据框
      # 定义年龄区间
      breaks <- c(0, 30, 50,70,Inf)
      labels <- c("0-30", "30-50", "50-70","70+")
      # 使用cut函数进行年龄区间划分
      single_and_MPC3_basic_all1$diagnose <- cut(single_and_MPC3_basic_all1$diagnose, breaks = breaks, labels = labels, right = FALSE)
      
      single_and_MPC3_basic_all1 <- single_and_MPC3_basic_all1 %>%
        mutate(BMI = case_when(
          BMI < 18.5 ~ "Underweight(<18.5)",
          BMI >= 18.5 & BMI < 23.0 ~ "Normal(18.5-23.0)",
          BMI >= 23.0 & BMI < 25.0 ~ "Overweight(23.0-<25.0)",
          BMI >= 25 ~ "Obesity(≥25)",
          TRUE ~ as.character(BMI)  # 处理其他情况，例如NA值
        ))
      
      
      names(single_and_MPC3_basic_all1)
      ######################################################################################
      table1(~ diagnose + BMI + smoking + alcohol|group,
             data=single_and_MPC3_basic_all1,
             extra.col=list(`P-value`=pvalue),overall=F)
    }
    {
      single_and_MPC4_basic_all<-rbind(single_1001_diagnose_bim_smok_alcohol_select,MPC4_1001_diagnose_bim_smok_alcohol_select) %>% unique()
      #332548
      single_and_MPC4_basic_all$group<- factor(single_and_MPC4_basic_all$group,
                                               levels=c(1,0), labels=c("Single cancer", "MPC ≥ 4"))
      
      
      single_and_MPC4_basic_all1<-single_and_MPC4_basic_all[!is.na(single_and_MPC4_basic_all$BMI),]
      
      single_and_MPC4_basic_all1<-single_and_MPC4_basic_all1[!is.na(single_and_MPC4_basic_all1$alcohol),]
      
      
      # 假设你的数据框叫做 df，并且有一个列名叫做 'age'
      # df <- ... # 你的数据框
      # 定义年龄区间
      breaks <- c(0, 30, 50,70,Inf)
      labels <- c("0-30", "30-50", "50-70","70+")
      # 使用cut函数进行年龄区间划分
      single_and_MPC4_basic_all1$diagnose <- cut(single_and_MPC4_basic_all1$diagnose, breaks = breaks, labels = labels, right = FALSE)
      
      single_and_MPC4_basic_all1 <- single_and_MPC4_basic_all1 %>%
        mutate(BMI = case_when(
          BMI < 18.5 ~ "Underweight(<18.5)",
          BMI >= 18.5 & BMI < 23.0 ~ "Normal(18.5-23.0)",
          BMI >= 23.0 ~ "Overweight(23.0-<25.0)",
          BMI >= 25 ~ "Obesity(≥25)",
          TRUE ~ as.character(BMI)  # 处理其他情况，例如NA值
        ))
      
      names(single_and_MPC4_basic_all1)
      ######################################################################################
      table1(~ diagnose + BMI + smoking + alcohol|group,
             data=single_and_MPC4_basic_all1,
             extra.col=list(`P-value`=pvalue),overall=F)
    }
  }
  
  {
    single_and_three_basic_all<-rbind(cancertype_single_diagnose_white_bim_smok_alcohol,cancertype_three_diagnose_white_bim_smok_alcohol)
    #332548
    single_and_three_basic_all$group<- factor(single_and_three_basic_all$group,
                                              levels=c(1,0), labels=c("Single cancer", "MPC ≥ 3"))
    # 假设你的数据框叫做 df，并且有一个列名叫做 'age'
    # df <- ... # 你的数据框
    # 定义年龄区间
    breaks <- c(0, 30, 40, 50, 60,70,80,Inf)
    labels <- c("0-30", "30-40", "40-50", "50-60","60-70","70-80","80+")
    # 使用cut函数进行年龄区间划分
    single_and_three_basic_all$diagnose <- cut(single_and_three_basic_all$diagnose, breaks = breaks, labels = labels, right = FALSE)
    names(single_and_three_basic_all)
    ######################################################################################
    table1(~ cancer_sex + diagnose + BMI + smoking + alcohol|group,
           data=single_and_three_basic_all,
           extra.col=list(`P-value`=pvalue),overall=F)
  }
  
  
  {
    single_and_four_basic_all <- rbind(cancertype_single_diagnose_white_bim_smok_alcohol, cancertype_four_diagnose_white_bim_smok_alcohol)
    #332548
    single_and_four_basic_all$group <- factor(single_and_four_basic_all$group,
                                              levels = c(1, 0), labels = c("Single cancer", "MPC ≥ 4"))
    # 假设你的数据框叫做 df，并且有一个列名叫做 'age'
    # df <- ... # 你的数据框
    # 定义年龄区间
    breaks <- c(0, 30, 40, 50, 60, 70, 80, Inf)
    labels <- c("0-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80+")
    # 使用cut函数进行年龄区间划分
    single_and_four_basic_all$diagnose <- cut(single_and_four_basic_all$diagnose, breaks = breaks, labels = labels, right = FALSE)
    names(single_and_four_basic_all)
    ######################################################################################
    table1(~ cancer_sex + diagnose + BMI + smoking + alcohol | group,
           data = single_and_four_basic_all,
           extra.col = list(`P-value` = pvalue), overall = F)
  }
}

}

save.image("/home/wuxh/luohh/3.cancertype_no_time/sanxian_res/sanxian.Rdata")
