samplePath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/03.sampleData"
genoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/01.genoData"
UKBGenoData='/home/luohh/UKB50wData/01.Download/02.ImputationData'
excludeVar='/home/luohh/UKB50wMultiPcancer/01.data/05.genotypeQC/hardy6.impuScore4.MHC.multiAllele.exclude.var'

### 数据准备
mkdir -p $samplePath/01.mpcSample
mkdir -p $samplePath/02.singleCancerSample

cp /home/luohh/UKB50wMultiPcancer/01.data/04.sampleData/07.SpecificMPC/* $samplePath/01.mpcSample
cp /home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/04.sampleData/* $samplePath/02.singleCancerSample

### 获取所有case case study 样本
cat $samplePath/01.mpcSample/* $samplePath/02.singleCancerSample/* | sort | uniq > $samplePath/caseCase.all.sample
awk '{print $1"\t"$1}' $samplePath/caseCase.all.sample > $samplePath/caseCase.all.gcta.sample

### 基因型数据准备
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
	plink2 --bgen $UKBGenoData/ukb_imp_chr${chr}_v3.bgen 'ref-first' --sample $UKBGenoData/ukb_imp_chr${chr}_v3.sample --exclude ${excludeVar} --keep ${samplePath}/caseCase.all.gcta.sample  --maf 0.0001 --max-maf 0.9999 --geno 0.05 --make-bed --out ${genoPath}/chr${chr}
	echo >&1000
}&
done
wait

exec 1000>&-
exec 1000<&-

### 合并bfile
for chr in $(seq 1 22)
do
	echo $genoPath/chr$chr >> $genoPath/allBfile.list
done

plink2 --pmerge-list $genoPath/allBfile.list bfile --make-bed --out $genoPath/merge
rm $genoPath/chr* $genoPath/allBfile.list $genoPath/merge.pgen $genoPath/merge.psam $genoPath/merge.pvar
