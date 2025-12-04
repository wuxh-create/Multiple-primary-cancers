gcta="gcta64"

genoPath="/home/wuxh/MPC/04.3cancersMPCAllSite/02.GeneBaseGWAS/02.GeneBaseGWAS/01.data"
genePath="/home/wuxh/MPC/04.3cancersMPCAllSite/02.GeneBaseGWAS/02.GeneBaseGWAS/01.data/glist-hg19.autosome.txt"
summaryPath="/home/wuxh/MPC/04.3cancersMPCAllSite/02.GeneBaseGWAS/02.GeneBaseGWAS/01.data"
outputPath="/home/wuxh/MPC/04.3cancersMPCAllSite/02.GeneBaseGWAS/02.GeneBaseGWAS/03.result"


${gcta} --bfile ${genoPath}/merge \
    --fastBAT ${summaryPath}/assoc.geneBased.txt \
    --threads 20 \
    --fastBAT-gene-list ${genePath} \
    --out ${outputPath}/genebase
