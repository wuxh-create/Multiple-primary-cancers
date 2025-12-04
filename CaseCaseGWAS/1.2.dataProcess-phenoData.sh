phenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData"
samplePath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/03.sampleData"
gwasPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/03.result"

cancerClass=("breast" "colorectal" "prostate" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "cervix" "kidney" "ovary" "otherDigestive" "uterus" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "otherFemaleGenital" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "testis" "mesothelioma" "otherEndocrine")
allSexCancerClass=("colorectal" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "kidney" "otherDigestive" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "mesothelioma" "otherEndocrine")
femaleCancerClass=("breast" "cervix" "ovary" "uterus" "otherFemaleGenital")
maleCancerClass=("prostate" "testis")


for cancer in ${cancerClass[@]}
do
	n=$(($n+1))
	dirID=$(printf "%02d\n" "$n").$cancer

    mkdir -p $phenoPath/$dirID
    mkdir -p $gwasPath/$dirID

    cp $samplePath/01.mpcSample/${cancer}.mpc.sample $phenoPath/$dirID/case.sample
    cp $samplePath/02.singleCancerSample/${cancer}.singleCancer.sample $phenoPath/$dirID/control.sample

    cat $phenoPath/$dirID/case.sample $phenoPath/$dirID/control.sample > $phenoPath/$dirID/GWAS.sample
    awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.sample > $phenoPath/$dirID/GWAS.gcta.sample

    awk '{print $1"\t"$1"\t1"}' $phenoPath/$dirID/case.sample > $phenoPath/$dirID/case.pheno
    awk '{print $1"\t"$1"\t0"}' $phenoPath/$dirID/control.sample > $phenoPath/$dirID/control.pheno

    cat $phenoPath/$dirID/case.pheno $phenoPath/$dirID/control.pheno > $phenoPath/$dirID/GWAS.pheno

    ### covar qcovar pheno提取
	for allSexCancer in ${allSexCancerClass[@]}
	do
		if [ $cancer = $allSexCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/selectSexSample.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getCovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getQcovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done
	
	for femaleCancer in ${femaleCancerClass[@]}
	do
		if [ $cancer = $femaleCancer ]	
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/selectSexSample.py $dirID female
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done

	for maleCancer in ${maleCancerClass[@]}
	do
		if [ $cancer = $maleCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/selectSexSample.py $dirID male
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done

    awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.last.sample > $phenoPath/$dirID/GWAS.last.gcta.sample

done