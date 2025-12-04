#!/bin/bash
# 从UKB数据中提取GWAS数据

UKBGenoData='/home/luohh/UKB50wData/01.Download/02.ImputationData'
keepSample='/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.gcta.sample'
excludeVar='/home/luohh/UKB50wMultiPcancer/01.data/05.genotypeQC/hardy6.impuScore4.MHC.multiAllele.exclude.var'
outputPath='/home/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC'

trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo genodata
exec 1000<>genodata
rm -rf genodata

for((n=1;n<=$1;n++))
do
	echo >&1000
done

date
for chr in $(seq 1 22)
do
read -u1000
{
	plink2 --bgen $UKBGenoData/ukb_imp_chr${chr}_v3.bgen 'ref-first' --sample $UKBGenoData/ukb_imp_chr${chr}_v3.sample \
	--exclude ${excludeVar} --keep ${keepSample}  \
	--maf 0.0001 --max-maf 0.9999 \
	--geno 0.05 --make-bed --out ${outputPath}/chr${chr}
	echo >&1000
}&
done
wait

exec 1000>&-
exec 1000<&-

date

### 合并bfile
for chr in $(seq 1 22)
do
	echo $outputPath/chr$chr >> $outputPath/allBfile.list
done

plink2 --pmerge-list $outputPath/allBfile.list bfile --make-bed --out $outputPath/merge
#rm $outputPath/chr* $outputPath/allBfile.list $outputPath/merge.pgen $outputPath/merge.psam $outputPath/merge.pvar
