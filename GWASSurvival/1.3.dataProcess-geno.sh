genoPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/01.genoData/merge"
phenoPath="/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/multiPCancer.gcta.sample"
outputPath="/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/04.gwasSurvival/01.data/01.genoData/15379Sample"

plink2 --bfile $genoPath --keep $phenoPath --maf 0.0001 --max-maf 0.9999 --geno 0.05 --export bgen-1.2 bits=8 --out $outputPath/mpcOnly
/home/luohh/Biosoft/miniconda3/envs/GCTA/bin/bgenix -g ${outputPath}/mpcOnly.bgen -index

plink2 --bgen $outputPath/mpcOnly.bgen "ref-first" --sample $outputPath/mpcOnly.sample --make-bed --out $outputPath/mpcOnly
