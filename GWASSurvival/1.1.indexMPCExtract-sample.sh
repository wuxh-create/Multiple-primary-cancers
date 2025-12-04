#!/bin/bash

# 定义源文件夹和目标文件夹
source_dir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/03.SpecificMPC"
target_dir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData"

# 定义子文件夹列表
subfolders=("01.breast" "07.lymphoidNeoplasms" "13.esophagus" "19.otherFemaleGenital" "25.thyroid" "31.testis"
            "02.colorectal" "08.cervix" "14.stomach" "20.tCellAndNKCellNeoplasms" "26.smallIntestine" "32.mesothelioma"
            "03.prostate" "09.kidney" "15.myeloidNeoplasms" "21.otherRespiratory" "27.brain" "33.otherEndocrine"
            "04.lung" "10.ovary" "16.headAndNeck" "22.gallbladderAndBiliaryTract" "28.anus"
            "05.melanoma" "11.otherDigestive" "17.pancreas" "23.softTissueSarcoma" "29.eyeAndOrbit"
            "06.urinaryBladder" "12.uterus" "18.liver" "24.otherUrinaryOrgans" "30.bone")

# 遍历子文件夹
for subfolder in "${subfolders[@]}"; do
    # 提取子文件夹名称（去掉数字和点）
    folder_name=$(echo "$subfolder" | sed 's/^[0-9]*\.//')

    echo "${folder_name}"

    # 定义源文件和目标文件路径
    source_file="$source_dir/$subfolder/${folder_name}.mpc.pheno"

    echo "${source_file}"

    target_file="$target_dir/$subfolder/${folder_name}.mpc.sample"
    target_fileID="$target_dir/$subfolder/${folder_name}.mpc.sampleID.txt"

    # 创建目标文件夹（如果不存在）
    mkdir -p "$(dirname "$target_file")"
    
    # 处理文件
    awk '{print $1 "\t"$2}' "$source_file" > "$target_file"

    awk '{print $1}' "$source_file" > "$target_fileID" 

done

echo "处理完成。"
