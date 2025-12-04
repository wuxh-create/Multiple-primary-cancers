gcta="/home/luohh/Biosoft/miniconda3/envs/myenv/bin/gcta64"

allPcancerGenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/01.genoData/merge"
phenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData"
grmPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/04.grmData"
outputPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/03.result"

cancerClass=("breast" "colorectal" "prostate" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "cervix" "kidney" "ovary" "otherDigestive" "uterus" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "otherFemaleGenital" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "testis" "mesothelioma" "otherEndocrine")
allSexCancerClass=("colorectal" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "kidney" "otherDigestive" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "mesothelioma" "otherEndocrine")
singleSexCancerClass=("breast" "prostate" "cervix" "ovary" "uterus" "otherFemaleGenital" "prostate" "testis")

trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo genodata
exec 1000<>genodata
rm -rf genodata

for((n=1;n<=$1;n++))
do
	echo >&1000
done

n=0
for cancer in ${cancerClass[@]}
do
	n=$(($n+1))
	dirID=$(printf "%02d\n" "$n").$cancer

	mkdir -p ${outputPath}/$dirID

	for allSexCancer in ${allSexCancerClass[@]}
	do
		if [ $cancer = $allSexCancer ]
		then
			read -u1000
			{
				# run GCTA-fastGWA-GLMM
				${gcta} --bfile ${allPcancerGenoPath} \
				  --grm-sparse ${grmPath}/Sparse_GRM_UKB \
				  --keep ${phenoPath}/$dirID/GWAS.last.gcta.sample \
				  --fastGWA-mlm-binary \
				  --qcovar ${phenoPath}/$dirID/GWAS.last.qcovar \
				  --covar ${phenoPath}/$dirID/GWAS.last.covar \
				  --pheno ${phenoPath}/$dirID/GWAS.last.pheno \
				  --threads 80 \
				  --out ${outputPath}/$dirID/fastGWA_GLMM_final
				echo >&1000
			}&
		fi
	done	
	
	for singleSexCancer in ${singleSexCancerClass[@]}
	do
		if [ $cancer = $singleSexCancer ]
		then
			read -u1000
			{
				# run GCTA-fastGWA-GLMM
				${gcta} --bfile ${allPcancerGenoPath} \
				  --grm-sparse ${grmPath}/Sparse_GRM_UKB \
				  --keep ${phenoPath}/$dirID/GWAS.last.gcta.sample \
				  --fastGWA-mlm-binary \
				  --qcovar ${phenoPath}/$dirID/GWAS.last.qcovar \
				  --pheno ${phenoPath}/$dirID/GWAS.last.pheno \
				  --threads 80 \
				  --out ${outputPath}/$dirID/fastGWA_GLMM_final
				echo >&1000
			}&
		fi
	done
done
wait

exec 1000>&-
exec 1000<&-

