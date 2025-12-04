#!/bin/bash
# 定义基础路径
base_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/01.genoData"
# 定义输出路径
output_base="/home/luohh/Wuxh/MPC/gwasSurvival"
# 遍历所有子文件夹
# 01.breast 
for folder in 03.prostate 04.lung 05.melanoma 06.urinaryBladder 07.lymphoidNeoplasms 08.cervix 09.kidney 10.ovary 11.otherDigestive 12.uterus 13.esophagus 14.stomach 15.myeloidNeoplasms 16.headAndNeck 17.pancreas 18.liver 19.otherFemaleGenital 20.tCellAndNKCellNeoplasms 21.otherRespiratory 22.gallbladderAndBiliaryTract 23.softTissueSarcoma 24.otherUrinaryOrgans 25.thyroid 26.smallIntestine 27.brain 28.anus 29.eyeAndOrbit 30.bone 31.testis 32.mesothelioma 33.otherEndocrine; do
    # 构建完整的路径
    tissue=$(echo$folder | cut -d '.' -f 2)
    plink_file="${output_base}/${folder}/genotype"
    # 构建输出前缀
    output_prefix="${output_base}/${folder}/${folder}"
    mkdir -p "$output_prefix"
    # 运行第一步的Rscript命令
    Rscript /home/luohh/Biosoft/GATE/extdata/step1_fitNULLGLMM.R \
        --plinkFile="$plink_file" \
        --phenoFile="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData/${folder}/${tissue}.mpc.sampleID.pheno" \
        --phenoCol=caseControl \
        --covarColList=sex,age,ageAtDiagnosis,stage,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 \
        --eventTimeCol=AgeOfEventFinal \
        --eventTimeBinSize=1 \
        --sampleIDColinphenoFile=sampleID \
        --traitType=survival \
        --outputPrefix="$output_prefix" \
        --nThreads=100 \
        --LOCO=FALSE \
        --minMAFforGRM=0.1 \
        --skipModelFitting=FALSE \
        --tauInit=1,0 \
        --pcgforUhatforSurvAnalysis=FALSE \
        --numRandomMarkerforVarianceRatio=15000
done
