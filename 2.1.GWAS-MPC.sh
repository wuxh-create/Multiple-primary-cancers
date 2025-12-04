# Perform GWAS with fastGWA-GLMM

gcta="/home/luohh/Biosoft/miniconda3/envs/myenv/bin/gcta64"
geno_path="/home/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/merge"
path_to_sparse_GRM="/home/luohh/UKB50wMultiPcancer/01.data/03.grmData/01.MPC/Sparse_GRM_UKB"
qcovar_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.qcovar"
covar_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.covar"
pheno_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.pheno"
output_path="/home/luohh/UKB50wMultiPcancer/03.result/01.MPC/"

# run GCTA-fastGWA-GLMM
${gcta} --bfile ${geno_path} \
 --grm-sparse ${path_to_sparse_GRM} \
 --fastGWA-mlm-binary \
 --qcovar ${qcovar_file} \
 --covar ${covar_file} \
 --pheno ${pheno_file} \
 --threads 80 \
 --out ${output_path}fastGWA_GLMM_final

