import sys

keepSample = {}
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.sample","r") as fi:
    for line in fi:
        data = line.strip()
        keepSample[data] = 0
fi.close()

fo = open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.pheno","w")
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.pheno","r") as fi:
    for line in fi:
        data = line.strip().split()
        if data[0] in keepSample:
            fo.write(line)
fi.close()
fo.close()
