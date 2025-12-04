import sys

### 取特有性别样本
sexDict = {}
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/basic.pheno","r") as fi:
    for line in fi:
        data = line.strip().split("\t")
        if data[1] == "male":
            sexDict[data[0]] = "M"
        elif data[1] == "female":
            sexDict[data[0]] = "F"
fi.close()

fo = open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.sample","w")
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/05.caseCaseGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.sample","r") as fi:
    for line in fi:
        sampleid = line.strip()
        if sys.argv[2] == "female" and sexDict[sampleid] == "F":
            fo.write(sampleid+"\n")
        elif sys.argv[2] == "male" and sexDict[sampleid] == "M":
            fo.write(sampleid+"\n")
        elif sys.argv[2] == "allSex":
            fo.write(sampleid+"\n")
fi.close()
fo.close()
