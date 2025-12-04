import sys

### 年龄 SNV
ageDict = {}
with open("/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/basic.pheno","r") as fi:
    fi.readline()
    for line in fi:
        data = line.strip().split("\t")
        ageDict[data[0]] = eval(data[3])
fi.close()

### 性别 SNV
sexDict = {}
with open("/home/luohh/UKB50wMultiPcancer/01.data/02.phenoData/basic.pheno","r") as fi:
    fi.readline()
    for line in fi:
        data = line.strip().split("\t")
        if data[1] == "female":
            sex = 2
        elif data[1] == "male":
            sex = 1
        sexDict[data[0]] = sex
fi.close()

### 基因型协变量 SNV
snvPcaInfoDict = {}
with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData/pca40.snv.eigenvec","r") as fi:
    for line in fi:
        data = line.strip().split("\t")
        snvPcaInfoDict[data[0]] = data[:22]
fi.close()

### 生成数量性状协变量文件
if sys.argv[2] == "singleSex":
    fo = open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.qcovar","w")
    with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.sample","r") as fi:
        for line in fi:
            sample = line.strip()
            mylist = snvPcaInfoDict[sample][:]
            mylist.insert(2,str(ageDict[sample]))
            mylist.insert(3,str(ageDict[sample]*ageDict[sample]))
            fo.write("\t".join(mylist)+"\n")
    fi.close()
    fo.close()
else:
    fo = open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.qcovar","w")
    with open("/home/luohh/UKB50wMultiPcancer/05.otherAnalysis/02.singleCancerGWAS/01.data/02.phenoData/"+sys.argv[1]+"/GWAS.last.sample","r") as fi:
        for line in fi:
            sample = line.strip()
            mylist = snvPcaInfoDict[sample][:]
            mylist.insert(2,str(ageDict[sample]))
            mylist.insert(3,str(ageDict[sample]*ageDict[sample]))
            mylist.insert(4,str(ageDict[sample]*sexDict[sample]))
            mylist.insert(5,str(sexDict[sample]*ageDict[sample]*ageDict[sample]))
            fo.write("\t".join(mylist)+"\n")
    fi.close()
    fo.close()
