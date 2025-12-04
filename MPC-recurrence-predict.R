


#浩辉师兄数据
{
  #两种，三种，以上多原发癌样本：/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/03.multiPcancer
  #癌症诊断时间数据：/home/luohh/UKB50wMultiPcancer2/05.otherAnalysis/01.cancerPair/02.result/cancerDiagnose.allPcancer.uniqueCancer.date
}

##part1----------------------------------------------------------------------------------------------------------------------------------
##载入包
{
  rm(list=ls())
  library(data.table)
  library(dplyr)
  library(tidyr)
}
##part2----------------------------------------------------------------------------------------------------------------------------------
#提取case和control，筛选白人样本，添加首发时间和诊断时间数据，提取不同类型PC之间的顺序情况，提取不同PC顺序对应的样本
{
  #现在一种癌症复发的情况，别人已经做过了，我们现在要做一种癌症转移为另一种癌症的情况
  #思路是根据已有的研究复发的论文学习如何研究一种癌症转移为另一种癌症的情况。
  #在研究的过程中使用机器学习的算法
  #研究过程的思路是：已经复发了的作为case，超过5年/10年未复发的作为control。
  {
    ##part2.1----------------------------------------------------------------------------------------------------------------------------
    #提取case和control，筛选白人。
    #提取不同癌症复发的作为case，multiple_cancer_recurrence
    {
      # all_recurrence<-read.table("/home/wuxh/luohh/cancerDiagnose.allPcancer.uniqueCancer.date.csv",sep = ",")
      # all_recurrence$V3<-as.Date(all_recurrence$V3, format="%Y/%m/%d")
      # all_recurrence$V5<-as.Date(all_recurrence$V5, format="%Y/%m/%d")
      # all_recurrence$first_diagnose_cha <- as.numeric(all_recurrence$V5 - all_recurrence$V3)
      # all_recurrence$first_diagnose_cha <- all_recurrence$first_diagnose_cha / 365
      # #no_recurrence<-all_recurrence[all_recurrence$first_diagnose_cha > 10,]
      # multiple_cancer_recurrence<-all_recurrence[all_recurrence$V2 != all_recurrence$V4,]
      # single_cancer_recurrence<-all_recurrence[all_recurrence$V2 == all_recurrence$V4,]
      # write.table(multiple_cancer_recurrence,
      #             "/home/wuxh/luohh/5.predict.MPC.recurrence/multiple_cancer_recurrence",
      #             row.names = F,col.names = T,sep = "\t",quote = FALSE)
      # write.table(single_cancer_recurrence,
      #             "/home/wuxh/luohh/5.predict.MPC.recurrence/single_cancer_recurrence",
      #             row.names = F,col.names = T,sep = "\t",quote = FALSE)
      multiple_cancer_recurrence<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/multiple_cancer_recurrence",sep = "\t",header = T)
      multiple_cancer_recurrence$V1<-as.numeric(multiple_cancer_recurrence$V1)
      
      single_cancer_recurrence<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/single_cancer_recurrence",sep = "\t",header = T)
      
      all_recurrence<-rbind(multiple_cancer_recurrence,single_cancer_recurrence)
      
      #共计存在复发的诊断日期的样本8540个，其中相同癌症的有1122个，不同癌症的有7418个
      #不同癌症复发的样本作为case
    }
    #筛选单癌样本作control【后面再筛选5年内没有复发的】
    {
      all_sample<-read.table("/home/wuxh/luohh/3.cancertype_no_time/count_cancer_type.csv",sep = ",",na.strings = "")
      single_sample<-all_sample[is.na(all_sample$V3),]
      #筛选乳腺癌的单癌样本
      breast<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/04.cancerGroup/03.whiteData/breast.filtered.white.sample")
      all_single_cancer<-read.table("/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/08.singleCancer/singleCacner.filtered.white.sample")
      
      breast_filter<-inner_join(breast,all_single_cancer)
      
      single_sample<-inner_join(breast_filter,single_sample)
      
      # all_MPC<-all_sample[!is.na(all_sample$V3),]
      
      all_MPC<-all_recurrence
      
    }
    #将single和all_MPC的样本结合白人样本进行筛选
    {
      #basic.info信息初步过滤
      {
        basic<-read.table("/home/luohh/UKB50wData/04.PhenoDataExtract/basic.pheno",sep="\t",header = T)
        #502387
        #选择两个性别一致的样本
        basic$sex_state<-ifelse(basic$sex==basic$geneticSex,"Y","N")
        basic<-basic %>% filter(sex_state=="Y")
        #487775
        basic_white<-basic %>% filter(population=="1001"|population=="1002"|population=="1003"|population=="1")
        #459680
        basic_white$sampleID<-as.integer(basic_white$sampleID)
      }
      #不同
      #multiple_cancer_recurrence<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/predict_MPC/cancerDiagnose.diffCancer.date.csv",sep = ",",fill = T)
      names(multiple_cancer_recurrence)[1]<-names(basic_white)[1]
      multi_cancer_white<-inner_join(multiple_cancer_recurrence,basic_white)
      
      #single
      names(single_sample)[1]<-names(basic_white)[1]
      single_sample_white<-inner_join(single_sample,basic_white)
      #all_MPC
      names(all_MPC)[1]<-names(basic_white)[1]
      all_MPC_white<-inner_join(all_MPC,basic_white)
    }
    ##part2.2----------------------------------------------------------------------------------------------------------------------------
    #为single添加首发时间,为MPC添加诊断时间数据
    {
      #诊断时间数据预处理
      {
        #获得癌症首发年龄=诊断日期-出生年月
        clinic<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/clinic.txt",sep="\t",quote="",na.strings = c("", "NA", "N/A", "NULL"), fill=TRUE)
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
        clinic_select<-clinic[,2:124]
        names(clinic_select)<-a[1:123]
        #将诊断时间缺失的数据去除
        clinic_select<-clinic_select[,c(1,3,10)]
        clinic_select<-separate(clinic_select,"40005-0.0",into="Date_of_cancer_diagnosis",sep = "-")
        names(clinic_select)<-c("sampleID","Year_of_birth","Date_of_cancer_diagnosis")
        str(clinic_select$Date_of_cancer_diagnosis)
        str(clinic_select$Year_of_birth)
        clinic_select$Date_of_cancer_diagnosis<-as.numeric(as.character(clinic_select$Date_of_cancer_diagnosis))
        clinic_select$sampleID<-as.numeric(as.character(clinic_select$sampleID))
        clinic_select$Year_of_birth<-as.numeric(as.character(clinic_select$Year_of_birth))
      }
      #筛选5年内未复发的单癌样本作为control
      {
        #单癌添加首发时间
        single_sample_white_clinic<-inner_join(single_sample_white,clinic_select)
        #去除5年内未复发的多原发癌
        #找首次发病到现在时间大于5年的
        #首次发病的时间就是诊断时间
        #现在的时间就是出生日期加现在的年龄
        #首发年龄就是癌症的诊断日期减去出生日期
        single_sample_white_clinic$diagnose_first<-single_sample_white_clinic$Date_of_cancer_diagnosis-single_sample_white_clinic$Year_of_birth
        #现在的时间就是年龄+出生日期
        single_sample_white_clinic$Now_time<-single_sample_white_clinic$age+single_sample_white_clinic$Year_of_birth
        #现在的时间-癌症诊断时间=患原发癌之后的观察时间
        single_sample_white_clinic$Now_jian_first_diagnose<-single_sample_white_clinic$Now_time-single_sample_white_clinic$Date_of_cancer_diagnosis
        single_sample_white_clinic1<-single_sample_white_clinic[!is.na(single_sample_white_clinic$Now_jian_first_diagnose),]
        single_sample_white_clinic_remove <- single_sample_white_clinic1[abs(single_sample_white_clinic1$Now_jian_first_diagnose) <= 5, ] %>% unique()
        single_sample_white_clinic_keep <- single_sample_white_clinic1[abs(single_sample_white_clinic1$Now_jian_first_diagnose) > 5, ] %>% unique()
        single_sample_white_clinic_keep<-single_sample_white_clinic_keep[which(single_sample_white_clinic_keep$geneticSex=='female'),]
        write.table(single_sample_white_clinic_keep,
                    "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/single_sample_1001_clinic_keep",
                    row.names = F,col.names = T,sep = "\t",quote = FALSE)
      }
    }
    #为不同类型MPC筛选白人样本，并添加诊断时间数据
    {
      names(all_recurrence)[1]<-"sampleID"
      all_MPC_white_select<-all_MPC_white[,c(1,20,21)] %>% as.data.frame()
      names(all_MPC_white_select)<-c("sampleID","sex","age")
      all_MPC_white_diagnose_filter_timeorder<-inner_join(all_recurrence,all_MPC_white_select)
      all_MPC_white_diagnose_filter_timeorder2<-inner_join(clinic_select,all_MPC_white_diagnose_filter_timeorder)
      all_MPC_white_diagnose_filter_timeorder2$diagnose_first<-all_MPC_white_diagnose_filter_timeorder2$Date_of_cancer_diagnosis-all_MPC_white_diagnose_filter_timeorder2$Year_of_birth
      write.table(all_MPC_white_diagnose_filter_timeorder2,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/all_MPC_white_diagnose_filter_timeorder2",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
    }
    ##part2.3----------------------------------------------------------------------------------------------------------------------------
    #统计癌症发生顺序
    #具体包括12，23，34，45是乳腺癌到其他癌症的，进行统计
    {
      #12
      {
        all_recurrence1<-inner_join(all_recurrence,basic)
        all_recurrence1<-all_recurrence1[which(all_recurrence1$sex=='female'),]
        all_recurrence<-all_recurrence1[,1:18]
        multiple_cancer_recurrence_unite<-unite(all_recurrence,"cancer12",c("V2","V4"),sep = "_",remove=F)
        multiple_cancer_recurrence_unite_white<-inner_join(multiple_cancer_recurrence_unite,basic_white)
        COUNT<- table(multiple_cancer_recurrence_unite_white$cancer12) %>% as.data.frame()
        
        write.table(COUNT,
                    file = "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/12/zz_cancertype_count.txt",
                    row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        
        a <- unique(multiple_cancer_recurrence_unite_white$cancer) 
        #找到癌症发生顺序对应的样本
        # 遍历每个癌症类型
        for(cancer_type in a) {
          # 筛选出对应癌症类型的样本信息
          subset_df <- multiple_cancer_recurrence_unite_white[multiple_cancer_recurrence_unite_white$cancer == cancer_type, ]
          # 将筛选出来的样本信息保存到文本文件中
          write.table(subset_df, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/12/',cancer_type, '.txt'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
          sample<-subset_df[,1]
          write.table(sample, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/12/',cancer_type, '.sample'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        }
      }
      #23
      {
        multiple_cancer_recurrence_unite<-unite(multiple_cancer_recurrence,"cancer23",c("V4","V6"),sep = "_",remove=F)
        multiple_cancer_recurrence_unite_white<-inner_join(multiple_cancer_recurrence_unite,basic_white)
        COUNT<- table(multiple_cancer_recurrence_unite_white$cancer23) %>% as.data.frame()
        
        write.table(COUNT,
                    file = "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/23/zz_cancertype_count.txt",
                    row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        
        a <- unique(multiple_cancer_recurrence_unite_white$cancer) 
        #找到癌症发生顺序对应的样本
        # 遍历每个癌症类型
        for(cancer_type in a) {
          # 筛选出对应癌症类型的样本信息
          subset_df <- multiple_cancer_recurrence_unite_white[multiple_cancer_recurrence_unite_white$cancer == cancer_type, ]
          # 将筛选出来的样本信息保存到文本文件中
          write.table(subset_df, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/23/',cancer_type, '.txt'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
          sample<-subset_df[,1]
          write.table(sample, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/23/',cancer_type, '.sample'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        }
      }
      #34
      {
        multiple_cancer_recurrence_unite<-unite(multiple_cancer_recurrence,"cancer34",c("V6","V8"),sep = "_",remove=F)
        multiple_cancer_recurrence_unite_white<-inner_join(multiple_cancer_recurrence_unite,basic_white)
        COUNT<- table(multiple_cancer_recurrence_unite_white$cancer34) %>% as.data.frame()
        
        write.table(COUNT,
                    file = "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/34/zz_cancertype_count.txt",
                    row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        
        a <- unique(multiple_cancer_recurrence_unite_white$cancer) 
        #找到癌症发生顺序对应的样本
        # 遍历每个癌症类型
        for(cancer_type in a) {
          # 筛选出对应癌症类型的样本信息
          subset_df <- multiple_cancer_recurrence_unite_white[multiple_cancer_recurrence_unite_white$cancer == cancer_type, ]
          # 将筛选出来的样本信息保存到文本文件中
          write.table(subset_df, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/34/',cancer_type, '.txt'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
          sample<-subset_df[,1]
          write.table(sample, 
                      file = paste0('/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/34/',cancer_type, '.sample'), 
                      row.names = FALSE, col.names = TRUE, sep = "\t",quote=F)
        }
      }
      #12，23，34数量统计直方图
      {
        yier<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/12/zz_cancertype_count.txt",header = T)
        ersan<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/23/zz_cancertype_count.txt",header = T)
        sansi<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/34/zz_cancertype_count.txt",header = T)
        total<-rbind(yier,ersan,sansi)
        library(dplyr)
        total_count<-total %>% 
          group_by(Var1) %>% 
          summarise(Freq=sum(Freq))
        total_count_sep<-separate(total_count,"Var1",into = c("cancer1","cancer2"),sep = "_")
        total_count_sep_filter <- total_count_sep %>%
          filter(cancer1 == "breast", !is.na(cancer2) & cancer2 != "") %>% 
          unite("type",c("cancer1","cancer2"),sep = "->")
        
        #画图
        {
          library(ggplot2)
          # 假设您的数据框名为df
          # 使用ggplot绘制柱形图
          total_count_sep_filter$color_group <- ifelse(total_count_sep_filter$Freq > 100, "#FFD180",
                                   ifelse(total_count_sep_filter$Freq > 50, "#99CCFF", "#F7A8A8"))
          
          # 使用ggplot绘制柱形图，并使用新的颜色变量
          p<-ggplot(total_count_sep_filter, aes(x=reorder(type, Freq), y=Freq, fill=color_group)) + 
            geom_bar(stat="identity", show.legend=FALSE) +  # 不显示图例
            coord_flip() +  # 翻转坐标轴，使条形图水平显示
            labs(x="Cancer Type", y="Frequency", title="Frequency of Cancer Types") +
            theme_minimal() +  # 使用简洁的主题
            theme(panel.grid.major = element_blank(),  # 去除主要的背景线
                  panel.grid.minor = element_blank(),  # 去除次要的背景线
                  panel.border = element_rect(colour = "black", fill=NA, size=1),  # 添加边框
                  axis.text.x = element_text(size=12),  # 放大x轴标签字体
                  axis.text.y = element_text(size=12)) +  # 放大y轴标签字体
            geom_text(aes(label=Freq), hjust=-0.1, size=3.5) +  # 添加每种type对应的数量
            scale_fill_manual(values = c("#FFD180", "#99CCFF", "#F7A8A8"))  # 手动设置颜色
          
          pdf_file <- "cancer_types_frequency.pdf"
          
          # 保存ggplot图形为PDF
          ggsave(pdf_file, plot = p, device = "pdf", width = 9, height = 8, path = "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/")
          
        }
        
      }
    }
  }
}
##part3---------------------------------------------------------------------------------------------------------------------------------- 
#将不同PC顺序的多原发癌样本和单癌样本进行PSM匹配，比例为1：1
{
  ##part3.1---------------------------------------------------------------------------------------------------------------------------------- 
  #数据准备
  {
    {
      args=commandArgs(TRUE)
      Cond1="12"
      Cond2="23"
      Cond3="34"
      
      Cond4="breast_colorectal"
      # Cond5="breast_headAndNeck"
      # Cond6="breast_kidney"
      Cond7="breast_lung"
      # Cond8="breast_lymphoidNeoplasms"
      Cond9="breast_melanoma"
      # Cond10="breast_myeloidNeoplasms"
      # Cond11="breast_ovary"
      # Cond12="breast_pancreas"
      # Cond13="breast_urinaryBladder"
      Cond14="breast_uterus"
      Cond15="breast_breast"
    }
    #不同顺序的PC样本
    {
      #breast_colorectal
      {
        yi_er<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond1,"/",Cond4,".sample"),header = T)
        er_san<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond2,"/",Cond4,".sample"),header = T)
        san_si<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond3,"/",Cond4,".sample"),header = T)
        breast_colorectal_sample<-rbind(yi_er,er_san,san_si)%>% unique() %>% as.data.frame()
        names(breast_colorectal_sample)[1]<-names(all_MPC_white_diagnose_filter_timeorder2)[1]
        all_MPC_white_diagnose_filter_timeorder2$sampleID<-as.integer(all_MPC_white_diagnose_filter_timeorder2$sampleID)
        breast_colorectal_sample1<-inner_join(breast_colorectal_sample,all_MPC_white_diagnose_filter_timeorder2)
      }
      #Cond7="breast_lung"
      {
        yi_er<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond1,"/",Cond7,".sample"),header = T)
        er_san<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond2,"/",Cond7,".sample"),header = T)
        # san_si<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond3,"/",Cond7,".sample"),header = T)
        breast_lung_sample<-rbind(yi_er,er_san)%>% unique() %>% as.data.frame()
        names(breast_lung_sample)[1]<-names(all_MPC_white_diagnose_filter_timeorder2)[1]
        all_MPC_white_diagnose_filter_timeorder2$sampleID<-as.integer(all_MPC_white_diagnose_filter_timeorder2$sampleID)
        breast_lung_sample1<-inner_join(breast_lung_sample,all_MPC_white_diagnose_filter_timeorder2)
      }
      #Cond9="breast_melanoma"
      {
        yi_er<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond1,"/",Cond9,".sample"),header = T)
        er_san<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond2,"/",Cond9,".sample"),header = T)
        # san_si<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond3,"/",Cond9,".sample"),header = T)
        breast_melanoma_sample<-rbind(yi_er,er_san)%>% unique() %>% as.data.frame()
        names(breast_melanoma_sample)[1]<-names(all_MPC_white_diagnose_filter_timeorder2)[1]
        all_MPC_white_diagnose_filter_timeorder2$sampleID<-as.integer(all_MPC_white_diagnose_filter_timeorder2$sampleID)
        breast_melanoma_sample1<-inner_join(breast_melanoma_sample,all_MPC_white_diagnose_filter_timeorder2)
      }
      #Cond14="breast_uterus"
      {
        yi_er<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond1,"/",Cond14,".sample"),header = T)
        er_san<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond2,"/",Cond14,".sample"),header = T)
        # san_si<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond3,"/",Cond14,".sample"),header = T)
        breast_uterus_sample<-rbind(yi_er,er_san)%>% unique() %>% as.data.frame()
        names(breast_uterus_sample)[1]<-names(all_MPC_white_diagnose_filter_timeorder2)[1]
        all_MPC_white_diagnose_filter_timeorder2$sampleID<-as.integer(all_MPC_white_diagnose_filter_timeorder2$sampleID)
        breast_uterus_sample1<-inner_join(breast_uterus_sample,all_MPC_white_diagnose_filter_timeorder2)
      }
      #Cond15="breast_breast"
      {
        yi_er<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond1,"/",Cond15,".sample"),header = T)
        # er_san<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond2,"/",Cond15,".sample"),header = T)
        # san_si<-read.table(paste0("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/",Cond3,"/",Cond14,".sample"),header = T)
        breast_breast_sample<-yi_er%>% unique() %>% as.data.frame()
        names(breast_breast_sample)[1]<-names(all_MPC_white_diagnose_filter_timeorder2)[1]
        all_MPC_white_diagnose_filter_timeorder2$sampleID<-as.integer(all_MPC_white_diagnose_filter_timeorder2$sampleID)
        breast_breast_sample1<-inner_join(breast_breast_sample,all_MPC_white_diagnose_filter_timeorder2)
      }
      
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_headAndNeck_sample1,breast_kidney_sample1,breast_lung_sample1,breast_lymphoidNeoplasms_sample1,
      #            breast_melanoma_sample1,breast_myeloidNeoplasms_sample1,breast_ovary_sample1,breast_pancreas_sample1,breast_urinaryBladder_sample1,
      #            breast_uterus_sample1,breast_breast_sample1)
      MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
                         breast_melanoma_sample1,
                         breast_uterus_sample1,breast_breast_sample1)
    }
    #single样本
    {
      single<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/single_sample_1001_clinic_keep",sep = "\t",header = T)
      sample<-single[,1] %>% unique()
      single_sample<-single
    }
    #bmi
    {
      bmi<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/bmi_res",fill=TRUE)
      names(bmi)<-bmi[1,]
      bmi<-bmi[-1,]
      names(bmi)[1]<-"sampleID"
      names(bmi)
      str(bmi)
      bmi <- data.frame(lapply(bmi, function(x) as.integer(as.character(x))))
      str(bmi)
      names(bmi)<-c("sampleID","21001-0.0","21001-1.0","21001-2.0","21001-3.0","23104-0.0","23104-1.0","23104-2.0","23104-3.0")
      
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
      
      single_sample_keep_bim<-merge(single_sample,new_bmi_select,all.x = T)
      MPC_sample_keep_bim<-merge(MPC_sample1,new_bmi_select,all.x = T)
    }
    #smoking
    {
      smok<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/20116_smoking_res",fill=TRUE)
      names(smok)<-smok[1,]
      smok<-smok[-1,]
      names(smok)[1]<-"sampleID"
      names(smok)
      str(smok)
      smok <- data.frame(lapply(smok, function(x) as.integer(as.character(x))))
      str(smok)
      names(smok)<-c("sampleID","20116-0.0","20116-1.0","20116-2.0","20116-3.0")
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
      
      
      single_sample_keep_bim_smok<-merge(single_sample_keep_bim,new_smok_select,all.x = T)
      MPC_sample_keep_bim_smok<-merge(MPC_sample_keep_bim,new_smok_select,all.x = T)
    }
    #alcohol
    {
      alcohol<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/20117_drink_alcohol_res",fill=TRUE) 
      names(alcohol)<-alcohol[1,]
      alcohol<-alcohol[-1,]
      str(alcohol)
      names(alcohol)
      alcohol <- data.frame(lapply(alcohol, function(x) as.integer(as.character(x))))
      str(alcohol)
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
      
      
      single_sample_keep_bim_smok_alcohol<-merge(single_sample_keep_bim_smok,new_alcohol_select,all.x = T)
      MPC_sample_keep_bim_smok_alcohol<-merge(MPC_sample_keep_bim_smok,new_alcohol_select,all.x = T)
      
  }
  ##part3.2---------------------------------------------------------------------------------------------------------------------------------- 
  #从大数据框中挑选出感兴趣的列，进行样本的倾向性评分匹配
  {
    single_sample_psm<-single_sample_keep_bim_smok[,c(1,2,8,20,21,17)]
    names(single_sample_psm)<-c("sampleID","cancer","sex","BMI","smoking","diagnose_age")
    MPC_sample_psm<-MPC_sample_keep_bim_smok[,c(1,4,21,24,25,23)]
    names(MPC_sample_psm)<-c("sampleID","cancer","sex","BMI","smoking","diagnose_age")
    single_sample_psm$group<-rep("single",nrow(single_sample_psm))
    MPC_sample_psm$group<-rep("MPC",nrow(MPC_sample_psm))
    MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,breast_melanoma_sample1,breast_uterus_sample1,breast_breast_sample1)
    {
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
      #                    breast_melanoma_sample1,
      #                    breast_uterus_sample1,breast_breast_sample1)
      breast_colorectal_sample2<-breast_colorectal_sample1[,1] %>% as.data.frame()
      names(breast_colorectal_sample2)<-'sampleID'
      breast_colorectal_case<-inner_join(MPC_sample_psm,breast_colorectal_sample2)
      nrow(breast_colorectal_case)
      breast_colorectal_control<-single_sample_psm[1:1910,]
      breast_colorectal_all<-rbind(breast_colorectal_case,breast_colorectal_control)
      nrow(breast_colorectal_all)
      a<-table1(~ BMI+smoking+diagnose_age|group,data=breast_colorectal_all,
             extra.col=list(`P-value`=pvalue),overall=F)
      a
      
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_colorectal_191/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_colorectal_all,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_colorectal_191/breast_colorectal_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_colorectal_all[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_colorectal_191/breast_colorectal_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
    }
    {
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
      #                    breast_melanoma_sample1,
      #                    breast_uterus_sample1,breast_breast_sample1)
      breast_lung_sample2<-breast_lung_sample1[,1] %>% as.data.frame()
      names(breast_lung_sample2)<-'sampleID'
      nrow(breast_lung_sample2)
      breast_lung_case<-inner_join(MPC_sample_psm,breast_lung_sample2)
      breast_lung_control<-single_sample_psm[1911:3710,]
      nrow(breast_lung_control)
      breast_lung_all<-rbind(breast_lung_case,breast_lung_control)
      nrow(breast_lung_all)
      a<-table1(~ BMI+smoking+diagnose_age|group,data=breast_lung_all,
             extra.col=list(`P-value`=pvalue),overall=F)
      a
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_lung_180/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_lung_all,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_lung_180/breast_lung_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_lung_all[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_lung_180/breast_lung_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
    }
    {
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
      #                    breast_melanoma_sample1,
      #                    breast_uterus_sample1,breast_breast_sample1)
      breast_melanoma_sample2<-breast_melanoma_sample1[,1] %>% as.data.frame()
      names(breast_melanoma_sample2)<-'sampleID'
      nrow(breast_melanoma_sample2)
      breast_melanoma_case<-inner_join(MPC_sample_psm,breast_melanoma_sample2)
      breast_melanoma_control<-single_sample_psm[3711:5320,]
      nrow(breast_melanoma_control)
      breast_melanoma_all<-rbind(breast_melanoma_case,breast_melanoma_control)
      nrow(breast_melanoma_all)
      a<-table1(~ BMI+smoking+diagnose_age|group,data=breast_melanoma_all,
             extra.col=list(`P-value`=pvalue),overall=F)
      a
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_melanoma_161/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_melanoma_all,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_melanoma_161/breast_melanoma_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_melanoma_all[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_melanoma_161/breast_melanoma_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
    }
    {
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
      #                    breast_melanoma_sample1,
      #                    breast_uterus_sample1,breast_breast_sample1)
      breast_uterus_sample2<-breast_uterus_sample1[,1] %>% as.data.frame()
      names(breast_uterus_sample2)<-'sampleID'
      nrow(breast_uterus_sample2)
      breast_uterus_case<-inner_join(MPC_sample_psm,breast_uterus_sample2)
      breast_uterus_control<-single_sample_psm[5321:7010,]
      nrow(breast_uterus_control)
      breast_uterus_all<-rbind(breast_uterus_case,breast_uterus_control)
      nrow(breast_uterus_all)
      a<-table1(~ BMI+smoking+diagnose_age|group,data=breast_uterus_all,
             extra.col=list(`P-value`=pvalue),overall=F)
      a
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_uterus_all,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/breast_uterus_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_uterus_all[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/breast_uterus_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
    }
    {
      # MPC_sample1<-rbind(breast_colorectal_sample1,breast_lung_sample1,
      #                    breast_melanoma_sample1,
      #                    breast_uterus_sample1,breast_breast_sample1)
      breast_breast_sample2<-breast_breast_sample1[,1] %>% as.data.frame()
      names(breast_breast_sample2)<-'sampleID'
      nrow(breast_breast_sample2)
      breast_breast_case<-inner_join(MPC_sample_psm,breast_breast_sample2)
      breast_breast_control<-single_sample_psm[4096:8505,]
      nrow(breast_breast_control)
      breast_breast_all<-rbind(breast_breast_case,breast_breast_control)
      nrow(breast_breast_all)
      a<-table1(~ BMI+smoking+diagnose_age|group,data=breast_breast_all,
             extra.col=list(`P-value`=pvalue),overall=F)
      a
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_breast_441/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_breast_all,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_breast_441/breast_breast_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(breast_breast_all[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_breast_441/breast_breast_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
      
    }
  }    
  ##part3.3---------------------------------------------------------------------------------------------------------------------------------- 
  #从casecaseGWAS结果中提取P小于等于10-5的SNP
  {
     
  }
  ##part3.4---------------------------------------------------------------------------------------------------------------------------------- 
  #提取特定样本特定SNP的基因型数据
  {
    # genotype path:/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/01.genoData
    /home/wuxh/software/plink2 --bfile /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/01.genoData/merge \
    --keep /home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/breast_uterus_sample \
    --extract /home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_fastGWA_p_1e-5_SNP_uniq \
    --make-bed --out /home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/genotype 
    #去除LD
    cd /home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169
    
    /home/wuxh/software/plink2 --bfile genotype --indep-pairwise 50 5 0.2 --out pruned_data
    /home/wuxh/software/plink2 --bfile genotype --extract pruned_data.prune.in --make-bed --out genotype_LD_0.2
    /home/wuxh/software/plink2 --bfile genotype --extract pruned_data.prune.in --export A --out genotype_LD_0.2
    
    
    
    
      
    cut -f1,7- genotype_LD_0.2.raw > genotype_LD_0.2_only
    
    sed 's/\bNA\b/-1/g' genotype_LD_0.2_only > genotype_LD_0.2_only_trans_NA_-1
    cut -f1 genotype_LD_0.2_only_trans_NA_-1 > genotype_LD_0.2_only_trans_NA_-1_sample
  }  
    
    
    
    
    
    {
      
      
      single_all<-rbind(single_sample_psm,MPC_sample_psm)
      single_all_remove<-single_all[which(single_all$sex=="male"),]
      single_all_na<-single_all[!is.na(single_all$BMI) & 
                                  !is.na(single_all$smoking) &
                                  !is.na(single_all$diagnose_age) ,]
      single_all_na<-single_all[!is.na(single_all$sampleID) & 
                                  !is.na(single_all$cancer) & 
                                  !is.na(single_all$sex) & 
                                  !is.na(single_all$BMI) & 
                                  !is.na(single_all$smoking) &
                                  !is.na(single_all$diagnose_age) ,]
      head(single_all_na)
      set.seed(123)
      library(tableone)
      #install.packages("table1")
      library(table1)
      #####################################################################
      {
        pvalue <- function(x, ...) {
          y <- unlist(x)
          g <- factor(rep(1:length(x), times=sapply(x, length)))
          if (is.numeric(y)) {
            p <- t.test(y ~ g)$p.value #数值型数据用t-test(两组比较)
          } else {
            p <- chisq.test(table(y, g))$p.value #因子型数据用卡方
          }
          c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
        }
      }
      ######################################################################################
      table1(~ sex+BMI+smoking+diagnose_age|group,data=single_all_na,
             extra.col=list(`P-value`=pvalue),overall=F)
      table1(~ sex+BMI+smoking+diagnose_age|group,data=single_all_na,
             overall=F)
      library(MatchIt)
      #PSM
      #####################################################################
      single_all_na$group<- factor(single_all_na$group,levels=c("single","MPC"), labels=c(1,0))
      breaks <- c(0, 30, 50,70,Inf)
      labels <- c("0-30", "30-50", "50-70","70+")
      # 使用cut函数进行年龄区间划分
      single_all_na$diagnose_age <- cut(single_all_na$diagnose_age, breaks = breaks, labels = labels, right = FALSE)
      single_all_na1 <- single_all_na %>%
        mutate(BMI = case_when(
          BMI < 18.5 ~ "Underweight(<18.5)",
          BMI >= 18.5 & BMI < 23.0 ~ "Normal(18.5-23.0)",
          BMI >= 23.0 & BMI < 25.0 ~ "Overweight(23.0-<25.0)",
          BMI >= 25 ~ "Obesity(≥25)",
          TRUE ~ as.character(BMI)  # 处理其他情况，例如NA值
        ))
      names(single_all_na1)
      ######################################################################################
      table1(~ BMI+smoking+diagnose_age|group,
             data=single_all_na,
             extra.col=list(`P-value`=pvalue),overall=F)
      single_all_na1$diagnose_age <- as.character(single_all_na1$diagnose_age)
      # match.it <- matchit(group ~ sex+BMI+smoking+alcohol+diagnose_age, data = single_all_na1, 
      #                     method = "nearest", distance = "logit", replace = FALSE, ratio = 1, caliper = 0.2) 
      #
      match.it <- matchit(group ~ BMI+smoking+diagnose_age, data = single_all_na1, 
                          method="nearest", distance = "logit", replace = TRUE, ratio=5)
      a <- summary(match.it)
      #plot(match.it, type="hist")
      lalonde_matched<-match.data(match.it)
      #####################################################
      lalonde_matched$group<- factor(lalonde_matched$group,levels=c(1,0), labels=c("single","MPC"))
      table1(~ BMI+smoking+diagnose_age|group,data=lalonde_matched,extra.col=list(`P-value`=pvalue),overall=F)
      table1(~ BMI+smoking+diagnose_age|group,data=lalonde_matched,overall=F)
      a<-table1(~ sex+BMI+smoking+diagnose_age|group,data=lalonde_matched,extra.col=list(`P-value`=pvalue),overall=F)
      write.table(a,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/total/PSM_table",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      
      write.table(lalonde_matched,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/total/PSM_match_all",
                  row.names = F,col.names = T,sep = "\t",quote = FALSE)
      write.table(lalonde_matched[,c(1,1)],
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/total/PSM_match_all_sample",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
      write.table(single_all_na1,
                  "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/total/single_MPC_metainfo",
                  row.names = F,col.names = F,sep = "\t",quote = FALSE)
    }
  }
}
##part4---------------------------------------------------------------------------------------------------------------------------------- 
#提取年龄(首次患癌的年龄)【84】，性别【22001】，身高【12144】，体重【21002】，吸烟【22506】，饮酒【20117】
{
  lalonde_matched<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/breast_uterus_all",
                              sep = "\t",header = T)
  sample<-lalonde_matched[,1] %>% as.data.frame()
  #Weight
  {
    Weight<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/21002_Weight.txt",sep = "\t",header = T)
    names(Weight)<-c("eid","21002-0.0","21002-1.0","21002-2.0","21002-3.0")
    new_Weight <- Weight %>%
      mutate(`21002-2.0` := ifelse(!is.na(`21002-3.0`), `21002-3.0`, `21002-2.0`),
             `21002-1.0` := ifelse(!is.na(`21002-2.0`), `21002-2.0`, `21002-1.0`),
             `21002-0.0` := ifelse(!is.na(`21002-1.0`), `21002-1.0`, `21002-0.0`)) 
    new_Weight_select<-new_Weight[,c(1,3)]
    names(new_Weight_select)[2]<-"Weight"
    names(sample)<-names(new_Weight_select)[1]
    sample_Weight<-merge(sample,new_Weight_select)
  }
  #Alcohol_drinker
  {
    Alcohol_drinker<-read.table("/home/wuxh/SSIdb/04.UKBB/phenodata/20117_Alcohol_drinker_status.txt",sep = "\t",header = T)
    names(Alcohol_drinker)<-c("eid","20117-0.0","20117-1.0","20117-2.0","20117-3.0")
    new_Alcohol_drinker <- Alcohol_drinker %>%
      mutate(`20117-2.0` := ifelse(!is.na(`20117-3.0`), `20117-3.0`, `20117-2.0`),
             `20117-1.0` := ifelse(!is.na(`20117-2.0`), `20117-2.0`, `20117-1.0`),
             `20117-0.0` := ifelse(!is.na(`20117-1.0`), `20117-1.0`, `20117-0.0`)) 
    new_Alcohol_drinker_select<-new_Alcohol_drinker[,c(1,2)]
    names(new_Alcohol_drinker_select)[2]<-"Alcohol"
    sample_Weight_Alcohol<-merge(sample_Weight,new_Alcohol_drinker_select)
    
    single_all_sex_smoking_diagnose_group<-lalonde_matched[,c(1,3,5,6,7)]
    
    names(single_all_sex_smoking_diagnose_group)[1]<-names(sample_Weight_Alcohol)[1]
    sample_Weight_Alcohol_sex_smoking_diagnose_group<-merge(sample_Weight_Alcohol,single_all_sex_smoking_diagnose_group)
    sample_Weight_Alcohol_sex_smoking_diagnose_group_order<-sample_Weight_Alcohol_sex_smoking_diagnose_group[,c(1,2,4,6,3,5,7)]
    
  }
  #将female替换为0，male替换为1
  #smoking中,将Prefer not to answer替换为-3，Never替换为0，Former替换为1，Current替换为2
  {
    sample_Weight_Alcohol_sex_smoking_diagnose_group_order <- sample_Weight_Alcohol_sex_smoking_diagnose_group_order %>%
      mutate(
        sex = case_when(
          sex == "female" ~ "0",
          sex == "male" ~ "1",
          TRUE ~ sex # 默认情况，保留原始值
        ),
        smoking = case_when(
          smoking == "Prefer not to answer" ~ "-3",
          smoking == "Never" ~ "0",
          smoking == "Former" ~ "1",
          smoking == "Current" ~ "2",
          TRUE ~ smoking # 默认情况，保留原始值
        ),
        group = case_when(
          group == "single" ~ "0",
          group =="MPC" ~ "1",
          TRUE ~ group #默认情况，保留原始值
        )
      )
  }
  #添加pc40
  {
    pc40<-read.table("/home/luohh/UKB50wData/04.PhenoDataExtract/pca40.snv.eigenvec",header = F,sep="\t")
    pc40<-pc40[,-2]
    names(pc40)[1]<-names(sample_Weight_Alcohol_sex_smoking_diagnose_group_order)[1]
    sample_Weight_Alcohol_sex_smoking_diagnose_group_order_pc40<-inner_join(sample_Weight_Alcohol_sex_smoking_diagnose_group_order,pc40)
    names(sample_Weight_Alcohol_sex_smoking_diagnose_group_order_pc40)<-c("sampleID","Weight","sex","diagnose_age","Alcohol","smoking","group",
                                                                          "pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","pc11",
                                                                          "pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20","pc21",
                                                                          "pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30","pc31",
                                                                          "pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")
  }
  #保证顺序
  #thousand
  {
    sample<-read.table("/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/genotype_LD_0.2_only_trans_NA_-1_sample",
                       sep="\t",header=T) %>% as.data.frame()
    names(sample)[1]<-names(sample_Weight_Alcohol_sex_smoking_diagnose_group_order_pc40)[1]
    sample_order<-inner_join(sample,sample_Weight_Alcohol_sex_smoking_diagnose_group_order_pc40)
    write.table(sample_order,
                "/home/wuxh/luohh/5.predict.MPC.recurrence/2predict_MPC/psm/breast_uterus_169/genotype_LD_0.2_only_trans_NA_-1_phenotype",
                row.names = F,col.names = T,sep = "\t",quote = FALSE) 
  }
}





