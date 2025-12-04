# #!/bin/bash
# # PLINK 2 executable path
# PLINK2="/home/wuxh/software/plink2"  # 请替换为您的PLINK 2可执行文件路径
# # 基因型文件路径
# genofile="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/merge"
# # 子文件夹列表
# subfolders=(
#     "01.breast" "07.lymphoidNeoplasms" "13.esophagus" "19.otherFemaleGenital" "25.thyroid" "31.testis"
#     "02.colorectal" "08.cervix" "14.stomach" "20.tCellAndNKCellNeoplasms" "26.smallIntestine" "32.mesothelioma"
#     "03.prostate" "09.kidney" "15.myeloidNeoplasms" "21.otherRespiratory" "27.brain" "33.otherEndocrine"
#     "04.lung" "10.ovary" "16.headAndNeck" "22.gallbladderAndBiliaryTract" "28.anus"
#     "05.melanoma" "11.otherDigestive" "17.pancreas" "23.softTissueSarcoma" "29.eyeAndOrbit"
#     "06.urinaryBladder" "12.uterus" "18.liver" "24.otherUrinaryOrgans" "30.bone"
# )
# #phenodir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData"
# phenodir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData"
# gendir="/home/luohh/Wuxh/MPC/gwasSurvival"
# # 遍历子文件夹
# for subfolder in "${subfolders[@]}"; do
#     # 读取子文件夹中的样本ID
#     folder_name=$(echo "$subfolder" | sed 's/^[0-9]*\.//')
#     echo "${folder_name}"
#     # 定义源文件和目标文件路径
#     samplelist="$phenodir/$subfolder/${folder_name}.mpc.sample"
#     echo "${samplelist}"
#     # 提取基因型
#     outprefix="${gendir}/${subfolder}"
#     mkdir -p "${gendir}/${subfolder}"
    
#     rm "$outprefix/genotype.bgen.bgi"

#     plink2 --bfile "$genofile" --keep "$samplelist" --maf 0.001 --max-maf 0.9999 --geno 0.05 --export bgen-1.2 bits=8 --out "$outprefix/genotype"
#     /home/luohh/Biosoft/miniconda3/envs/GCTA/bin/bgenix -g $outprefix/genotype.bgen -index

#     plink2 --bfile "$genofile" --keep "$samplelist" --maf 0.001 --max-maf 0.9999 --geno 0.05 --make-bed --out "$outprefix/genotype"

#     # 删除临时样本ID文件
#     # rm "$samplelist"
# done

# echo "基因型提取完成。"



#!/bin/bash
# PLINK 2 executable path
PLINK2="/home/wuxh/software/plink2"  # 请替换为您的PLINK 2可执行文件路径
# 基因型文件路径
genofile="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/merge"
# 子文件夹列表

# "01.breast" "02.colorectal" "03.prostate" "04.lung" "05.melanoma" "06.urinaryBladder" "07.lymphoidNeoplasms" "08.cervix"
# "09.kidney" "10.ovary" "11.otherDigestive" "12.uterus" "13.esophagus" "14.stomach" "15.myeloidNeoplasms" "16.headAndNeck"
# "17.pancreas" "18.liver" "19.otherFemaleGenital" "20.tCellAndNKCellNeoplasms" "21.otherRespiratory" "22.gallbladderAndBiliaryTract" 
# "23.softTissueSarcoma" "24.otherUrinaryOrgans" "25.thyroid" "26.smallIntestine" "27.brain" "28.anus" "29.eyeAndOrbit" "30.bone" 
# "31.testis" 

subfolders=(
    "32.mesothelioma" "33.otherEndocrine"
)
#phenodir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData"
phenodir="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/02.phenoData"
gendir="/home/luohh/Wuxh/MPC/gwasSurvival"

# 最大并行任务数
MAX_JOBS=2

# 遍历子文件夹
for subfolder in "${subfolders[@]}"; do
    # 检查当前运行的作业数量
    while (( $(jobs -r | wc -l) >= MAX_JOBS )); do
        sleep 10  # 如果达到最大并行任务数，等待10秒后再次检查
    done

    # 启动后台任务
    (
        # 读取子文件夹中的样本ID
        folder_name=$(echo "$subfolder" | sed 's/^[0-9]*\.//')
        echo "正在处理: ${folder_name}"
        echo "样本列表: ${phenodir}/${subfolder}/${folder_name}.mpc.sample"

        # 定义源文件和目标文件路径
        outprefix="${gendir}/${subfolder}"
        mkdir -p "${gendir}/${subfolder}"

        # 删除旧的输出文件
        rm -f "$outprefix/genotype.bgen.bgi"

        # 提取基因型并导出为bgen格式
        plink2 --bfile "$genofile" --keep "${phenodir}/${subfolder}/${folder_name}.mpc.sample" \
            --maf 0.001 --max-maf 0.9999 --geno 0.05 \
            --export bgen-1.2 bits=8 --out "$outprefix/genotype"
        /home/luohh/Biosoft/miniconda3/envs/GCTA/bin/bgenix -g "$outprefix/genotype.bgen" -index

        # 提取基因型并导出为bed格式
        plink2 --bfile "$genofile" --keep "${phenodir}/${subfolder}/${folder_name}.mpc.sample" \
            --maf 0.001 --max-maf 0.9999 --geno 0.05 \
            --make-bed --out "$outprefix/genotype_bed"

        # 删除临时样本ID文件（如果需要）
        # rm "${phenodir}/${subfolder}/${folder_name}.mpc.sample"
    ) &

    echo "已启动任务: ${subfolder}"
done

# 等待所有后台任务完成
wait

echo "所有基因型提取任务已完成。"