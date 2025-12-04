import datetime
import sys
import os

def dateSubtraction(date1,date2):
    d0 = datetime.datetime.strptime(date1, "%Y-%m-%d")
    df = datetime.datetime.strptime(date2, "%Y-%m-%d")
    return (df-d0).days

ICD10CancerDict = {}
for i in range(98):
    ICD10CancerDict["C{:0>2d}".format(i)] = []
for i in ["D00","D01","D02","D03","D04","D05","D06","D07","D08","D09"]:
    ICD10CancerDict[i] = []
for i in range(98):
    for j in ["","0","1","2","3","4","5","6","7","8","9"]:
        ICD10CancerDict["C{:0>2d}".format(i)].append("C{:0>2d}".format(i)+j)
for i in ["D00","D01","D02","D03","D04","D05","D06","D07","D08","D09"]:
    for j in ["","0","1","2","3","4","5","6","7","8","9"]:
        ICD10CancerDict[i].append(i+j)

### 各个分类癌症编码
clinicCancerDict = {}
clinicCancerDict["headAndNeck"] = ICD10CancerDict["C00"] + ICD10CancerDict["C01"] + ICD10CancerDict["C02"] + ICD10CancerDict["C03"] + ICD10CancerDict["C04"] + ICD10CancerDict["C05"] + ICD10CancerDict["C06"] + ICD10CancerDict["C07"] + ICD10CancerDict["C08"] + ICD10CancerDict["C09"] + ICD10CancerDict["C10"] + ICD10CancerDict["C11"] + ICD10CancerDict["C12"] + ICD10CancerDict["C13"] + ICD10CancerDict["C14"] + ["D000"]
clinicCancerDict["esophagus"] = ICD10CancerDict["C15"] + ["D001"]
clinicCancerDict["stomach"] = ICD10CancerDict["C16"] + ["D002"]
clinicCancerDict["smallIntestine"] = ICD10CancerDict["C17"]
clinicCancerDict["colorectal"] = ICD10CancerDict["C18"] + ICD10CancerDict["C19"] + ICD10CancerDict["C20"] + ["D010","D011","D012"]
clinicCancerDict["anus"] = ICD10CancerDict["C21"] + ["D013"]
clinicCancerDict["liver"] = ICD10CancerDict["C22"]
clinicCancerDict["gallbladderAndBiliaryTract"] = ICD10CancerDict["C23"] + ICD10CancerDict["C24"]
clinicCancerDict["pancreas"] = ICD10CancerDict["C25"]
clinicCancerDict["otherDigestive"] = ICD10CancerDict["C26"] + ICD10CancerDict["C48"] + ["D014","D017","D019"]
clinicCancerDict["otherRespiratory"] = ICD10CancerDict["C30"] + ICD10CancerDict["C31"] + ICD10CancerDict["C32"] + ICD10CancerDict["C33"] + ICD10CancerDict["C39"] + ["C381","C382","C383","C384","D020","D021","D023"]
clinicCancerDict["lung"] = ICD10CancerDict["C34"] + ["D022"]
clinicCancerDict["otherEndocrine"] = ICD10CancerDict["C34"] + ["D022"]
clinicCancerDict["bone"] = ICD10CancerDict["C40"] + ICD10CancerDict["C41"]
clinicCancerDict["melanoma"] = ICD10CancerDict["C43"] + ICD10CancerDict["D03"]
clinicCancerDict["mesothelioma"] = ICD10CancerDict["C45"]
clinicCancerDict["kaposiSarcoma"] = ICD10CancerDict["C46"]
clinicCancerDict["softTissueSarcoma"] = ICD10CancerDict["C47"] + ICD10CancerDict["C49"]
clinicCancerDict["breast"] = ICD10CancerDict["C50"] + ICD10CancerDict["D05"]
clinicCancerDict["otherFemaleGenital"] = ICD10CancerDict["C51"] + ICD10CancerDict["C52"] + ICD10CancerDict["C57"] + ICD10CancerDict["C58"] + ["D071","D072","D073"]
clinicCancerDict["cervix"] = ICD10CancerDict["C53"] + ICD10CancerDict["D06"]
clinicCancerDict["uterus"] = ICD10CancerDict["C54"] + ICD10CancerDict["C55"] + ["D070"]
clinicCancerDict["ovary"] = ICD10CancerDict["C56"]
clinicCancerDict["otherMaleGenital"] = ICD10CancerDict["C60"] + ICD10CancerDict["C63"] + ["D074","D076"]
clinicCancerDict["prostate"] = ICD10CancerDict["C61"] + ["D075"]
clinicCancerDict["testis"] = ICD10CancerDict["C62"]
clinicCancerDict["kidney"] = ICD10CancerDict["C64"] + ICD10CancerDict["C65"]
clinicCancerDict["otherUrinaryOrgans"] = ICD10CancerDict["C66"] + ICD10CancerDict["C68"] + ["D091"]
clinicCancerDict["urinaryBladder"] = ICD10CancerDict["C67"] + ["D090"]
clinicCancerDict["eyeAndOrbit"] = ICD10CancerDict["C69"] + ["D092"]
clinicCancerDict["otherNervousSystem"] = ICD10CancerDict["C70"] + ICD10CancerDict["C72"]
clinicCancerDict["brain"] = ICD10CancerDict["C71"]
clinicCancerDict["thyroid"] = ICD10CancerDict["C73"]
clinicCancerDict["tCellAndNKCellNeoplasms"] = ICD10CancerDict["C81"] + ICD10CancerDict["C84"] + ICD10CancerDict["C86"]
clinicCancerDict["lymphoidNeoplasms"] = ICD10CancerDict["C82"] + ICD10CancerDict["C83"] + ICD10CancerDict["C91"] + ["C851","C884"]
clinicCancerDict["myeloidNeoplasms"] = ICD10CancerDict["C90"] + ICD10CancerDict["C92"] + ["C930","C931"]

clinicCancerList = []
for key,value in clinicCancerDict.items():
    clinicCancerList += value

cancerSampleList = []
with open("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/02.phenoData/01.MPC/GWAS.newsample.pheno","r") as fi:
    for line in fi:
        data = line.strip().split("\t")
        if data[2] == "1":
            cancerSampleList.append(data[0])
fi.close()

fo = open("/home/wuxh/03.MultiPCancer/luohh/UKB50wMultiPcancer/01.data/01.genoData/01.MPC/newGenotype/GATE/clinic.tsv","w")
fo.write("\t".join(["sampleID","sex","AgeOfEventFinal","birthYear","population","osDays","status","ageAtDiagnosis","stage"])+"\n")

with open('/home/wuxh/UKB50wData/01.Download/04.PhenotypeData/clinic.txt',"r",encoding='gbk') as fi:
    header = fi.readline().strip("\n").lstrip().split("\t")
    
    # 根据header提取各个目标表型所在的index，方便后续提取对应表型数据
    # 31：性别、34：出生年份、52：出生月份
    # 21000：人群
    # 40000：死亡日期、40005：诊断癌症日期、40006：患病记录

    for index in range(len(header)):
        if "31" == header[index].split("-")[0]:
            sexIndex = index 
        if "34" == header[index].split("-")[0]:
            birthYearIndex = index
        if "52" == header[index].split("-")[0]:
            birthMouthIndex = index
    
    popIndexList = []
    for index in range(len(header)):
        if "21000" == header[index].split("-")[0]:
            popIndexList.append(index)
    
    deathIndexList = []
    for index in range(len(header)):
        if "40000" == header[index].split("-")[0]:
            deathIndexList.append(index)
    
    diagIndexList = []
    for index in range(len(header)):
        if "40005" == header[index].split("-")[0]:
            diagIndexList.append(index)

    cancerIndexList = []
    for index in range(len(header)):
        if "40006" == header[index].split("-")[0]:
            cancerIndexList.append(index)

    cancerStageList = []
    for index in range(len(header)):
        if "40012" == header[index].split("-")[0]:
            cancerStageList.append(index)

    # 读取表型文件内容
    for line in fi:
        data = line.strip("\n").lstrip().split("\t")
        #try:

        # 判断这行数据的样本是否患有目标疾病，没有就跳过，提高运行速率
        if data[0] in cancerSampleList:

            # 获取每个样本的出生日期，由于没有具体日期记录，所以我们选取每个月第一天作为他们的出生日期 
            birthDate = "-".join((data[birthYearIndex],data[birthMouthIndex],"01"))

            # 该样本的具体性别
            sexData = data[sexIndex]

            # 该样本得人群信息
            popData = [data[index] for index in popIndexList]

            # 该样本的死亡日期、诊断该癌症的日期以及对应第几次癌症记录
            targetData1 = [data[index] for index in deathIndexList]
            targetData2 = [data[index] for index in diagIndexList]
            targetData3 = [data[index] for index in cancerIndexList]
            targetData4 = [data[index] for index in cancerStageList]

            # 如果同一疾病有多个诊断日期，取最近一次日期作为诊断日期
            diagDateList = []
            stageList = []
            for i in range(len(targetData3)):
                if targetData3[i] in clinicCancerList:
                    diagDateList.append(targetData2[i])
                    stageList.append(targetData4[i])
            
            # 判断该癌症样本是否有预后信息
            if diagDateList == []:
                continue
            else:
                indexList = [i for i in sorted(enumerate(diagDateList),key=lambda x:x[1],reverse=True)]
                stageSortedList = []
                for i in indexList:
                    stageSortedList.append(stageList[i[0]])
                diagDateList.sort(reverse=True)

            # 获取诊断时的年龄，返回结果以天为单位
            ageAtDiagnosis = dateSubtraction(birthDate,diagDateList[0])/365

            birthYear = data[birthYearIndex]

            # 第一个癌症的分期作为MPC分析
            cancerStage = stageSortedList[0]
            if cancerStage == "":
                continue

            # 将性别编码转为性别
            if sexData == "0":
                sex = "female"
            else:
                sex = "male"

            popDict = {"1":"White", "1001":"White", "1002":"White", "1003":"White", "2":"Mixed", "2001":"Mixed", "2002":"Mixed", "2003":"Mixed", "2004":"Mixed",
                        "3":"Asian", "3001":"Asian", "3002":"Asian", "3003":"Asian", "3004":"Asian", "4":"Black", "4001":"Black", "4002":"Black", "4003":"Black",
                        "5":"Chinese", "6":"Other", "-1":"Unknown", "-3":"Unknown"}
            # 获取人群信息
            if popData[0] == "" and popData[1] == "" and popData[2] == "":
                popData = popDict["-1"]
            else:
                if popData[0] != "":
                    popData = popDict[popData[0]]
                elif popData[1] != "":
                    popData = popDict[popData[1]]
                else:
                    popData = popDict[popData[2]]

            # 以 TCGA临床文件的格式，输出临床文件
            # 如果该人有死亡时间记录，我们则认为他以及死了，否则则活着。
            if targetData1[0] == "":

                # 由于没有我们不知道我们拿到的表型文件最后一次随访日期是多少，
                # 我们用最近一个病人死亡时间近似估计最后一次随访时间，即，2021-11-12
                osDays = dateSubtraction(diagDateList[0],"2021-11-12")
                # 获取年龄，返回结果以年为单位
                censorAge = str(round(dateSubtraction(birthDate,"2021-11-12")/365))
                fo.write("\t".join((data[0],sex,censorAge,birthYear,popData,str(osDays),"alive",str(ageAtDiagnosis),cancerStage))+"\n")
            else:
                osDays = dateSubtraction(diagDateList[0],targetData1[0])
                censorAge = str(round(dateSubtraction(birthDate,targetData1[0])/365))
                fo.write("\t".join((data[0],sex,censorAge,birthYear,popData,str(osDays),"dead",str(ageAtDiagnosis),cancerStage))+"\n")
        # except:
        #     print(data)
fi.close()
fo.close()
