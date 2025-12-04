#!/usr/bin/perl

samplePath="/home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC"
phenoPath="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/03.SpecificMPC"
cancerClass=("breast" "colorectal" "prostate" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "cervix" "kidney" "ovary" "otherDigestive" "uterus" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "otherFemaleGenital" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "testis" "mesothelioma" "otherEndocrine")

allSexCancerClass=("colorectal" "lung" "melanoma" "urinaryBladder" "lymphoidNeoplasms" "kidney" "otherDigestive" "esophagus" "stomach" "myeloidNeoplasms" "headAndNeck" "pancreas" "liver" "tCellAndNKCellNeoplasms" "otherRespiratory" "gallbladderAndBiliaryTract" "softTissueSarcoma" "otherUrinaryOrgans" "thyroid" "smallIntestine" "brain" "anus" "eyeAndOrbit" "bone" "mesothelioma" "otherEndocrine")
femaleCancerClass=("breast" "cervix" "ovary" "uterus" "otherFemaleGenital")
maleCancerClass=("prostate" "testis")

n=0
for cancer in ${cancerClass[@]}
do
	n=$(($n+1))
	dirID=$(printf "%02d\n" "$n").$cancer

	mkdir -p $phenoPath/$dirID

	cp $samplePath/$cancer.mpc.sample $phenoPath/$dirID
	awk '{print $1"\t"$1"\t1"}' $phenoPath/$dirID/$cancer.mpc.sample > $phenoPath/$dirID/$cancer.mpc.pheno
	cat $phenoPath/$dirID/$cancer.mpc.sample $phenoPath/../nonCancer.filtered.white.sample > $phenoPath/$dirID/GWAS.sample
	cat $phenoPath/$dirID/$cancer.mpc.pheno $phenoPath/../nonCancer.filtered.white.pheno > $phenoPath/$dirID/GWAS.pheno
	awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.sample > $phenoPath/$dirID/GWAS.gcta.sample
	for allSexCancer in ${allSexCancerClass[@]}
	do
		if [ $cancer = $allSexCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/selectSexSample.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getCovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getQcovar.py $dirID allSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getPheno.py $dirID
		fi
	done
	
	for femaleCancer in ${femaleCancerClass[@]}
	do
		if [ $cancer = $femaleCancer ]	
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/selectSexSample.py $dirID female
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getPheno.py $dirID
		fi
	done

	for maleCancer in ${maleCancerClass[@]}
	do
		if [ $cancer = $maleCancer ]
		then
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/selectSexSample.py $dirID male
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getCovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getQcovar.py $dirID singleSex
			/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/python /home/luohh/UKB50wMultiPcancer/02.script/OtherScript/getPheno.py $dirID
		fi
	done
	
	awk '{print $1"\t"$1}' $phenoPath/$dirID/GWAS.last.sample > $phenoPath/$dirID/GWAS.last.gcta.sample			
done
