# Generate 200 GRM segments using GCTA. Later, these GRM-segments will be merged into a GRM. 
# See details at https://yanglab.westlake.edu.cn/software/gcta/index.html#MakingaGRM

gcta="/home/luohh/Biosoft/miniconda3/envs/myenv/bin/gcta64"

##  Step 1
# For step 1, please use job-arrays on HPC. The script is designed for dataset like the UK Biobank (N=~450k).
# User may reduce the number of segments accordingly if a smaller dataset is being used. 

#  a list of paths to the BGEN files (CHR 1 to 22):
ukb_bfile="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/merge"
#  the .sample file for the BGEN files (only one .sample file is needed)
#ukb_sample_file="ukb22828_c1_b0_v3_s487281.sample"

#  the list for QCed high-quality model-SNPs to construct the GRM 
#model_snp="/home/luohh/UKB50w_GWAS/01.Data/01.rawData/QC.01.rsid"

#  The inferred EUR individual list
#EUR_list="/home/luohh/UKB50w_GWAS/01.Data/02.phenoData/01.multiCancer/GWAS.gcta.sample"

output_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/03.grmData/01.MPC/"
# job_id=${SGE_TASK_ID}
# echo "Processing Part No. ${job_id} of the GRM...."

# # 设置多线程
# trap "exec 1000>&-;exec 1000<&-;exit 0" 2

# mkfifo genodata
# exec 1000<>genodata
# rm -rf genodata

# for((n=1;n<=$1;n++))
# do
# 	echo >&1000
# done

# date

# for i in $(seq 49 69)
# do
# 	read -u1000
# 	{
# 		${gcta} --bfile ${ukb_bfile} \
# 				--make-grm-part 100 $i \
# 				--thread-num 5 \
# 				--out ${output_path}UKB_GRM_parts
# 		echo >&1000
# 	}&
# done

# wait
# date
# exec 1000>&-
# exec 1000<&-

#if [[ $? -eq 0 ]]
#then   
#  echo "Completed Part No. ${job_id} of the GRM! \n"
#fi


# ##  Step 2
#  Merge the 200 parts into a complete GRM
# cat ${output_path}UKB_GRM_parts.part_100_*.grm.id > ${output_path}UKB_GRM.grm.id
# cat ${output_path}UKB_GRM_parts.part_100_*.grm.bin > ${output_path}UKB_GRM.grm.bin
# cat ${output_path}UKB_GRM_parts.part_100_*.grm.N.bin > ${output_path}UKB_GRM.grm.N.bin


# # ##  Step 3
# # # Get the sparse GRM
${gcta} --grm ${output_path}UKB_GRM --make-bK-sparse 0.05 --out ${output_path}Sparse_GRM_UKB

# # Get the unrelated individual list
# ${gcta} --grm ${output_path}UKB_GRM --grm-cutoff 0.05 --make-grm --out ${output_path}unrelated_0.05_UKB

# Get covar
#${gcta} --grm ${output_path}UKB_GRM --pca 10 --out ${output_path}pca10 --thread-num 40
