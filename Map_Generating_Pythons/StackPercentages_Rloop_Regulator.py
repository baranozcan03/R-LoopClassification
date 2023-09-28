import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

rloop_intersections= sys.argv[1]
reg_intersection_list=[sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13],sys.argv[14],sys.argv[15]]
numberOfregulators = int(sys.argv[16])
stackPercentages_path_csv= sys.argv[17]
reg_list=["ATR","BRCA1","DDX21","FIP1L1","NFAT5","PRMT5","RAD51","RPA1","SETX","SIRT7","SMC3","SRPK2","SRSF1","XRN2"]
segment_file_path=sys.argv[18]
pngFilePath = sys.argv[19]

def dicBuilder(segment_file_path):
    MyDic={}
    with open( segment_file_path,"r") as inputs:
        for line in inputs:
            line = line.strip().split("\t")
            MyDic[(line[1],line[2],line[3],line[0])] = []
    return MyDic
    

BinDic= dicBuilder(segment_file_path)


def dicMutator(intersectionPath):
    with open(intersectionPath, "r") as inputs:
        for line in inputs:
            line = line.strip().split("\t")
            oldValue = (BinDic.get((line[1],line[2],line[3],line[0]))).copy()
            oldValue.append(int(line[4]))
            BinDic.update({(line[1],line[2],line[3],line[0]):oldValue})


dicMutator(rloop_intersections)


for index in range (0,numberOfregulators): 
    dicMutator(reg_intersection_list[index])


BinDicSeg1 = dict()
BinDicSeg2 = dict()
BinDicSeg3 = dict()
BinDicSeg4 = dict()
BinDicSeg5 = dict()
BinDicSeg6 = dict()
BinDicSeg7 = dict()
BinDicSeg8 = dict()
BinDicSeg9 = dict()
BinDicSeg10 = dict()



for k in BinDic:
    if k[2] == "E1":
        value = BinDic[k]
        BinDicSeg1[k] = value
    elif k[2] == "E2":
        value = BinDic[k]
        BinDicSeg2[k] = value
    elif k[2] == "E3":
        value = BinDic[k]
        BinDicSeg3[k] = value
    elif k[2] == "E4":
        value = BinDic[k]
        BinDicSeg4[k] = value
    elif k[2] == "E5":
        value = BinDic[k]
        BinDicSeg5[k] = value
    elif k[2] == "E6":
        value = BinDic[k]
        BinDicSeg6[k] = value
    elif k[2] == "E7":
        value = BinDic[k]
        BinDicSeg7[k] = value
    elif k[2] == "E8":
        value = BinDic[k]
        BinDicSeg8[k] = value
    elif k[2] == "E9":
        value = BinDic[k]
        BinDicSeg9[k] = value
    elif k[2] == "E10":
        value = BinDic[k]
        BinDicSeg10[k] = value
        

segmentPercentages=[0,0,0,0,0,0,0,0,0,0]
totalBins=0
rBins=0
for key in BinDicSeg1:
    totalBins+=1
    value = BinDicSeg1[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[0] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg2:
    totalBins+=1
    value = BinDicSeg2[key] 
    if value[0] != 0:
        rBins+=1
segmentPercentages[1] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg3:
    totalBins+=1
    value = BinDicSeg3[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[2] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg4:
    totalBins+=1
    value = BinDicSeg4[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[3] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg5:
    totalBins+=1
    value = BinDicSeg5[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[4] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg6:
    totalBins+=1
    value = BinDicSeg6[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[5] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg7:
    totalBins+=1
    value = BinDicSeg7[key] 
    if value[0] != 0:
        rBins+=1
segmentPercentages[6] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg8:
    totalBins+=1
    value = BinDicSeg8[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[7] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg9:
    totalBins+=1
    value = BinDicSeg9[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[8] = rBins*100/totalBins
totalBins=0
rBins=0
for key in BinDicSeg10:
    totalBins+=1
    value = BinDicSeg10[key] 
    if value[0] != 0: 
        rBins+=1
segmentPercentages[9] = rBins*100/totalBins




#number of bins, presence of r loops and regulators in binary. 
#"00" if none, "01" only r-loop, "10" only reg, "11" both

def SegmentDicAssessor(dic,regList,segmentName):
    Output = open(stackPercentages_path_csv+"/"+segmentName+".csv" , "w")
    Output.write("RegulatorName,None,Only_Regulator,Only_RLoop,Both" + "\n")
    for regName in regList:
        NumsBin=0
        Nums00=0
        Nums01=0
        Nums10=0
        Nums11=0
        rloop=0
        for key in dic:
            value = dic[key]
            NumsBin += 1
            if (value[0] == 0) and (value[(regList.index(regName)+1)] == 0):
                Nums00+=1
            elif (value[0] == 0) and (value[(regList.index(regName)+1)] != 0):
                Nums01+=1
            elif (value[0] != 0) and (value[(regList.index(regName)+1)] == 0):
                Nums10+=1
            elif (value[0] != 0) and (value[(regList.index(regName)+1)] != 0):
                Nums11+=1
        percent00=Nums00*100/NumsBin
        percent01=Nums01*100/NumsBin
        percent10=Nums10*100/NumsBin
        percent11=Nums11*100/NumsBin
        Output.write(str(regName)+","+str(percent00)+","+str(percent01)+","+str(percent10)+","+str(percent11)+"\n")
    Output.close()

SegmentDicAssessor(BinDicSeg1,reg_list,"E1")
SegmentDicAssessor(BinDicSeg2,reg_list,"E2")
SegmentDicAssessor(BinDicSeg3,reg_list,"E3")
SegmentDicAssessor(BinDicSeg4,reg_list,"E4")
SegmentDicAssessor(BinDicSeg5,reg_list,"E5")
SegmentDicAssessor(BinDicSeg6,reg_list,"E6")
SegmentDicAssessor(BinDicSeg7,reg_list,"E7")
SegmentDicAssessor(BinDicSeg8,reg_list,"E8")
SegmentDicAssessor(BinDicSeg9,reg_list,"E9")
SegmentDicAssessor(BinDicSeg10,reg_list,"E10")




segments=["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"]
with sns.axes_style("whitegrid"):
    plot =sns.barplot(x=segments, y=segmentPercentages)
figure = plot.get_figure()
figure.savefig(pngFilePath+"/_Rloop_percentages.png")
plt.clf()


for segment in segments:
    df = pd.read_csv(stackPercentages_path_csv+"/"+segment+".csv")
    df.set_index('RegulatorName',inplace=True)
    print(df)
    plot =df.plot(kind='barh',stacked=True,color=['darkred','royalblue','plum','forestgreen'])
    figure = plot.get_figure()
    figure.savefig(pngFilePath+"/"+segment+"_Rloop-Reg_percentages.png")
    plt.clf()
    df_reg_r_loop= df[['Only_Regulator','Only_RLoop','Both']].copy()
    print(df_reg_r_loop)
    plot =df_reg_r_loop.plot(kind='barh',stacked=True,color=['royalblue','plum','forestgreen'])
    figure = plot.get_figure()
    figure.savefig(pngFilePath+"/"+segment+"_Rloop-Reg_percentages_ecpetNone.png")
    plt.clf()