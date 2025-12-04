import subprocess

variable_names = [
    "PrediXcan_Adipose_Subcutaneous","PrediXcan_Brain_Spinal_cord_cervical_c-1",
    "PrediXcan_Nerve_Tibial","PrediXcan_Adipose_Visceral_Omentum",
    "PrediXcan_Brain_Substantia_nigra","PrediXcan_Ovary",
    "PrediXcan_Adrenal_Gland","PrediXcan_Breast_Mammary_Tissue",
    "PrediXcan_Pancreas","PrediXcan_Artery_Aorta",
    "PrediXcan_Cells_Cultured_fibroblasts","PrediXcan_Pituitary",
    "PrediXcan_Artery_Coronary","PrediXcan_Cells_EBV-transformed_lymphocytes",
    "PrediXcan_Prostate","PrediXcan_Artery_Tibial",
    "PrediXcan_Colon_Sigmoid","PrediXcan_Skin_Not_Sun_Exposed_Suprapubic",
    "PrediXcan_Brain_Amygdala","PrediXcan_Colon_Transverse",
    "PrediXcan_Skin_Sun_Exposed_Lower_leg","PrediXcan_Brain_Anterior_cingulate_cortex_BA24",
    "PrediXcan_Esophagus_Gastroesophageal_Junction","PrediXcan_Small_Intestine_Terminal_Ileum",
    "PrediXcan_Brain_Caudate_basal_ganglia","PrediXcan_Esophagus_Mucosa",
    "PrediXcan_Spleen","PrediXcan_Brain_Cerebellar_Hemisphere",
    "PrediXcan_Esophagus_Muscularis","PrediXcan_Stomach",
    "PrediXcan_Brain_Cerebellum","PrediXcan_Heart_Atrial_Appendage",
    "PrediXcan_Testis","PrediXcan_Brain_Cortex",
    "PrediXcan_Heart_Left_Ventricle","PrediXcan_Thyroid",
    "PrediXcan_Brain_Frontal_Cortex_BA9","PrediXcan_Kidney_Cortex",
    "PrediXcan_Uterus","PrediXcan_Brain_Hippocampus",
    "PrediXcan_Liver","PrediXcan_Vagina",
    "PrediXcan_Brain_Hypothalamus","PrediXcan_Lung",
    "PrediXcan_Whole_Blood","PrediXcan_Brain_Nucleus_accumbens_basal_ganglia",
    "PrediXcan_Minor_Salivary_Gland","PrediXcan_Brain_Putamen_basal_ganglia",
    "PrediXcan_Muscle_Skeletal"
]
base_command = (
    "/home/wuxh/luohh/07.TWAS/01.install.software/MetaXcan/software/SPrediXcan.py "
    "--model_db_path /home/wuxh/03.MultiPCancer/05.TWAS/01.install.software/MR-JTI/pre-training/{db_path} "
    "--covariance /home/wuxh/03.MultiPCancer/05.TWAS/01.install.software/MR-JTI/pre-training/{cov_path}.txt.gz "
    "--gwas_folder /home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/03.result/01.MPC "
    "--gwas_file_pattern '.*5e-8.header.results' "
    "--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P "
    "--output_file /home/wuxh/luohh/07.TWAS/03.result/{output_file}"
)
# 遍历变量名称列表，构建并执行命令
for var_name in variable_names:
    db_path = f"{var_name}.db"
    cov_path = f"{var_name}"
    output_file = f"{var_name}.csv"
    command = base_command.format(db_path=db_path, cov_path=var_name, output_file=output_file)
    print(command)  # 打印命令，实际使用时可以替换为os.system(command)来执行命令
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command executed successfully for {var_name}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing command for {var_name}: {e}")

