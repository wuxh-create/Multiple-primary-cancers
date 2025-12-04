# Perform GWAS with fastGWA-GLMM

gcta="/home/luohh/Biosoft/miniconda3/envs/myenv/bin/gcta64"
geno_path="/home/luohh/UKB50wMultiPcancer/01.data/01.genoData/03.3cancersMPC/allbgen.list"
genoSample_path="/home/luohh/UKB50wMultiPcancer/01.data/01.genoData/03.3cancersMPC/chr1.sample"
sample_path="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/05.4cancersMPC/GWAS.gcta.sample"
path_to_sparse_GRM="/home/luohh/UKB50wMultiPcancer/01.data/03.grmData/01.MPC/Sparse_GRM_UKB"
qcovar_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/05.4cancersMPC/GWAS.qcovar"
covar_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/05.4cancersMPC/GWAS.covar"
pheno_file="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/05.4cancersMPC/GWAS.pheno"
output_path="/home/luohh/UKB50wMultiPcancer/03.result/05.4cancersMPCAllSite/"

# run GCTA-fastGWA-GLMM
${gcta} --mbgen ${geno_path} \
 --grm-sparse ${path_to_sparse_GRM} \
 --sample ${genoSample_path} \
 --maf 0.000001 \
 --keep ${sample_path} \
 --fastGWA-mlm-binary \
 --qcovar ${qcovar_file} \
 --covar ${covar_file} \
 --pheno ${pheno_file} \
 --threads 30 \
 --out ${output_path}fastGWA_GLMM_final

