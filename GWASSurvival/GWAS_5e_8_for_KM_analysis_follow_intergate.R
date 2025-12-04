


rm(list=ls())
gc()
GATE<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/03.result/step2Results.GWAS_5e-8.txt",
                 header = F)
names(GATE)<-c("CHR","POS","SNP","Allele1","Allele2","AC_Allele2","AF_Allele2","imputationInfo","N", 
               "BETA","SE","Tstat","p.value","p.value.NA","Is.SPA.converge","varT","varTstar","AF.Events","AF.Censored","N.Events","N.Censored")
KM<-read.table("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/06.GWASSurvival_5e_8_SNP_Survival/z_snp_p_values",
               header = T)
names(KM)
GATE_KM<-left_join(GATE,KM)
GATE_KM$OR<-exp(GATE_KM$BETA)
GATE_KM_select<-GATE_KM[,c(3,4,5,7,9,10,11,13,14,23,20,21,22)]
names(GATE_KM_select)[13]<-"KM P_value"
write.table(GATE_KM_select, 
            file = "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/GATE_5e_8_KM_intergate_res.txt", 
            quote = FALSE, row.names = FALSE,sep = "\t")
