gcta="/home/luohh/Biosoft/miniconda3/envs/myenv/bin/gcta64"

allPcancerGenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/01.genoData/allbgen.list"
samplePath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/01.genoData/chr1.sample"
phenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData"
grmPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/03.grmData"
outputPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/03.result"

cancerClass=("breast" "prostate" "colorectal" "melanoma" "cervix" "lung" "lymphoidNeoplasms" "urinaryBladder" "uterus" "myeloidNeoplasms" "kidney" "headAndNeck" "ovary" "pancreas" "brain" "esophagus" "testis" "thyroid" "liver" "tCellAndNKCellNeoplasms" "stomach" "softTissueSarcoma" "otherFemaleGenital" "eyeAndOrbit" "otherRespiratory" "mesothelioma" "anus" "otherMaleGenital" "smallIntestine" "otherEndocrine" "otherDigestive" "gallbladderAndBiliaryTract")

allSexCancerClass=("colorectal" "melanoma" "lung" "lymphoidNeoplasms" "urinaryBladder" "myeloidNeoplasms" "kidney" "headAndNeck" "pancreas" "brain" "esophagus" "thyroid" "liver" "tCellAndNKCellNeoplasms" "stomach" "softTissueSarcoma" "eyeAndOrbit" "otherRespiratory" "mesothelioma" "anus" "smallIntestine" "otherEndocrine" "otherDigestive" "gallbladderAndBiliaryTract")
singleSexCancerClass=("breast" "cervix" "ovary" "uterus" "otherFemaleGenital" "prostate" "testis" "otherMaleGenital")

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
				${gcta} --mbgen ${allPcancerGenoPath} \
                  --sample ${samplePath} \
				  --grm-sparse ${grmPath}/Sparse_GRM_UKB \
				  --keep ${phenoPath}/$dirID/GWAS.last.gcta.sample \
				  --fastGWA-mlm-binary \
				  --qcovar ${phenoPath}/$dirID/GWAS.last.qcovar \
				  --covar ${phenoPath}/$dirID/GWAS.last.covar \
				  --pheno ${phenoPath}/$dirID/GWAS.last.pheno \
				  --threads 10 \
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
				${gcta} --mbgen ${allPcancerGenoPath} \
                  --sample ${samplePath} \
				  --grm-sparse ${grmPath}/Sparse_GRM_UKB \
				  --keep ${phenoPath}/$dirID/GWAS.last.gcta.sample \
				  --fastGWA-mlm-binary \
				  --qcovar ${phenoPath}/$dirID/GWAS.last.qcovar \
				  --pheno ${phenoPath}/$dirID/GWAS.last.pheno \
				  --threads 10 \
				  --out ${outputPath}/$dirID/fastGWA_GLMM_final
				echo >&1000
			}&
		fi
	done
done
wait

exec 1000>&-
exec 1000<&-

