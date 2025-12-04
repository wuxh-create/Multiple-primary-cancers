#!/bin/bash
# 定义基础路径
base_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data"
# 定义输出路径
output_base="/home/luohh/Wuxh/MPC/gwasSurvival"

# 01.breast

# 遍历所有子文件夹
for folder in 06.urinaryBladder; do 
# 02.colorectal 03.prostate 04.lung 05.melanoma 06.urinaryBladder 07.lymphoidNeoplasms 08.cervix 09.kidney 10.ovary 11.otherDigestive 12.uterus 13.esophagus 
# 14.stomach 15.myeloidNeoplasms 16.headAndNeck 17.pancreas 18.liver 19.otherFemaleGenital 20.tCellAndNKCellNeoplasms 21.otherRespiratory 22.gallbladderAndBiliaryTract 
# 23.softTissueSarcoma 24.otherUrinaryOrgans 25.thyroid 26.smallIntestine 27.brain 28.anus 29.eyeAndOrbit 30.bone 31.testis 32.mesothelioma 33.otherEndocrine
    # 构建完整的路径
    folder_name=$(echo "$folder" | sed 's/^[0-9]*\.//')
    echo "${folder_name}"

    plink_file="${output_base}/${folder}/genotype.bgen"
    plink_bgi_file="${output_base}/${folder}/genotype.bgen.bgi"
    plink_sampleID_file="${base_path}/02.phenoData/${folder}/${folder_name}.mpc.sampleID.txt"

    # 构建输出目录
    output_dir="${output_base}/${folder}"
    mkdir -p "$output_dir"

    # 运行第二步的Rscript命令
    Rscript /home/luohh/Biosoft/GATE/extdata/step2_SPAtests.R \
        --bgenFile="$plink_file" \
        --bgenFileIndex="$plink_bgi_file" \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile="$plink_sampleID_file" \
        --GMMATmodelFile="${output_dir}/${folder}.rda" \
        --varianceRatioFile="${output_dir}/${folder}.varianceRatio.txt" \
        --SAIGEOutputFile="${output_dir}/${folder}.step2Results.GWAS.txt" \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsSPAfast=TRUE
done
