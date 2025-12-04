#!/bin/bash

# 定义根路径
base_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/03.result"

# 定义所有子文件夹名称
folders=(
    "01.breast" "02.colorectal" "03.prostate" "04.lung" "05.melanoma" "06.urinaryBladder" 
    "07.lymphoidNeoplasms" "08.cervix" "09.kidney" "10.ovary" "11.otherDigestive" "12.uterus" 
    "13.esophagus" "14.stomach" "15.myeloidNeoplasms" "16.headAndNeck" "17.pancreas" "18.liver" 
    "19.otherFemaleGenital" "20.tCellAndNKCellNeoplasms" "21.otherRespiratory" "22.gallbladderAndBiliaryTract" 
    "23.softTissueSarcoma" "24.otherUrinaryOrgans" "25.thyroid" "26.smallIntestine" "27.brain" 
    "28.anus" "29.eyeAndOrbit" "30.bone" "31.testis" "32.mesothelioma" "33.otherEndocrine"
)

# 遍历每个子文件夹
for folder in "${folders[@]}"; do
    
    mkdir -p "/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/04.geneBasedGWAS/${folder}"
    
    # 构建输入文件路径
    input_file="$base_path/$folder/fastGWA_GLMM_final.fastGWA"
    
    # 构建输出文件路径
    output_file="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/04.geneBasedGWAS/${folder}/assoc.geneBased.txt"
    
    # 检查输入文件是否存在
    if [[ -f "$input_file" ]]; then
        # 执行awk命令
        awk '{print $2"\t"$4"\t"$5"\t"$7"\t"$11"\t"$12"\t"$13"\t"$6}' "$input_file" > "$output_file"
        echo "Processed: $folder"
    else
        echo "File not found: $input_file"
    fi
done