#!/usr/bin/perl

samplePath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/04.sampleData"
phenoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData"
cancerClass=("breast" "prostate" "colorectal" "melanoma" "cervix" "lung" "lymphoidNeoplasms" "urinaryBladder" "uterus" "myeloidNeoplasms" "kidney" "headAndNeck" "ovary" "pancreas" "brain" "esophagus" "testis" "thyroid" "liver" "tCellAndNKCellNeoplasms" "stomach" "softTissueSarcoma" "otherFemaleGenital" "eyeAndOrbit" "otherRespiratory" "mesothelioma" "anus" "otherMaleGenital" "smallIntestine" "otherEndocrine" "otherDigestive" "gallbladderAndBiliaryTract")

allSexCancerClass=("colorectal" "melanoma" "lung" "lymphoidNeoplasms" "urinaryBladder" "myeloidNeoplasms" "kidney" "headAndNeck" "pancreas" "brain" "esophagus" "thyroid" "liver" "tCellAndNKCellNeoplasms" "stomach" "softTissueSarcoma" "eyeAndOrbit" "otherRespiratory" "mesothelioma" "anus" "smallIntestine" "otherEndocrine" "otherDigestive" "gallbladderAndBiliaryTract")
femaleCancerClass=("breast" "cervix" "ovary" "uterus" "otherFemaleGenital")
maleCancerClass=("prostate" "testis" "otherMaleGenital")

n=0
for cancer in ${cancerClass[@]}
do
	n=$(($n+1))
	dirID=$(printf "%02d\n" "$n").$cancer
	mkdir -p $phenoPath/$dirID

	cp $samplePath/$cancer.singleCancer.sample $phenoPath/$dirID
	awk '{print $1"\t"$1"\t1"}' $phenoPath/$dirID/$cancer.singleCancer.sample > $phenoPath/$dirID/$cancer.singleCancer.pheno
	cat $phenoPath/$dirID/$cancer.singleCancer.sample $phenoPath/nonCancer.filtered.white.sample > $phenoPath/$dirID/GWAS.sample
	cat $phenoPath/$dirID/$cancer.singleCancer.pheno $phenoPath/nonCancer.filtered.white.pheno > $phenoPath/$dirID/GWAS.pheno
	awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.sample > $phenoPath/$dirID/GWAS.gcta.sample
	for allSexCancer in ${allSexCancerClass[@]}
	do
		if [ $cancer = $allSexCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/selectSexSample.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getCovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getQcovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done
	
	for femaleCancer in ${femaleCancerClass[@]}
	do
		if [ $cancer = $femaleCancer ]	
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/selectSexSample.py $dirID female
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done

	for maleCancer in ${maleCancerClass[@]}
	do
		if [ $cancer = $maleCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/selectSexSample.py $dirID male
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/02.script/OtherScript/getPheno.py $dirID
		fi
	done
	
	awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.last.sample > $phenoPath/$dirID/GWAS.last.gcta.sample			
done
