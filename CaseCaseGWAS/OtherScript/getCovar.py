import sys

### 性别
sexDict = {}
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/basic.pheno","r") as fi:
    for line in fi:
        data = line.strip().split("\t")
        if data[1] == "male":
            sexDict[data[0]] = "M"
        elif data[1] == "female":
            sexDict[data[0]] = "F"
fi.close()

fo = open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.covar","w")
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.sample","r") as fi:
    for line in fi:
        sample = line.strip()
        if sys.argv[2] != "singleSex":
            fo.write("\t".join((sample,sample,sexDict[sample]))+"\n")
        else:
            fo.write("\t".join((sample,sample))+"\n")
fi.close()
fo.close()
