# Perform GWAS with fastGWA-GLMM

gcta="gcta64"
geno_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/06.5cancersMPCAllSite/01.data/5cancerMPCAll-merge"
path_to_sparse_GRM="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/03.grmData/01.MPC/Sparse_GRM_UKB"
qcovar_file="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/06.5cancerMPC/GWAS.qcovar"
covar_file="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/06.5cancerMPC/GWAS.covar"
pheno_file="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/06.5cancerMPC/GWAS.pheno"
output_path="/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/06.5cancersMPCAllSite/03.result/"

# run GCTA-fastGWA-GLMM
${gcta} --bfile ${geno_path} \
 --grm-sparse ${path_to_sparse_GRM} \
 --fastGWA-mlm-binary \
 --qcovar ${qcovar_file} \
 --covar ${covar_file} \
 --pheno ${pheno_file} \
 --threads 80 \
 --out ${output_path}fastGWA_GLMM_final

