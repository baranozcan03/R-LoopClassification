import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path = sys.argv[1]
histogram_path_png = sys.argv[2]
histogram_path_csv = sys.argv[3]
segment_file_path = sys.argv[4]

dir_list = os.listdir(path)
for i in range(0, len(dir_list)):
    dir_list[i] = [dir_list[i],path+"/" +dir_list[i]]
print(dir_list)

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

for index in range (0,len(dir_list)): 
    dicMutator(dir_list[index][1])

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
        

###Here we have a dic for every segmetn with a key as the (start,end,segment,chr)= 


reg_list = ["RAD51","R-Loops","ATR","DS-seq","DDX21","XRN2","SRSF1","SIRT7","DS_SIM-seq","XR-seq","PRMT5","SETX","ATAC-seq","RPA1","XR_SIM-seq","FIP1L1","NFAT-5","RNA-seq","SRPK2","SMC3","BRCA1",]
first_line=""
for i in reg_list[:-1]:
    first_line = first_line + i + "," 
first_line = first_line + reg_list[-1]+ "\n"
print(first_line)


def SegmentDicAssessor(dic,reg_interseciton_List,segmentName):
    Output = open(histogram_path_csv+"/"+segmentName+".csv" , "w")
    Output.write(first_line)
    for key in dic:
       for i in range(0,len(reg_interseciton_List)):
       with open(reg_interseciton_List[i][1], "r") as inputs:
         
    
    
    for i in range(0,len(reg_interseciton_List)):
        for key in dic:
            if key[2] == segmentName:
                value = dic[key]
                with open(reg_interseciton_List[i][1], "r") as inputs:
                    for line in inputs:
                        line = line.strip().split("\t")
                        if line[1] == key[0]:
                            value[i] = line[4]  
                newline= ""
                for l in range(0,len(reg_interseciton_List)-1):
                    newline= newline+str(value[l])+","
                newline = newline + str(value[-1]) + "\n"
                Output.write(newline)
    Output.close()

SegmentDicAssessor(BinDicSeg1,dir_list,"E1")
SegmentDicAssessor(BinDicSeg2,dir_list,"E2")
SegmentDicAssessor(BinDicSeg3,dir_list,"E3")
SegmentDicAssessor(BinDicSeg4,dir_list,"E4")
SegmentDicAssessor(BinDicSeg5,dir_list,"E5")
SegmentDicAssessor(BinDicSeg6,dir_list,"E6")
SegmentDicAssessor(BinDicSeg7,dir_list,"E7")
SegmentDicAssessor(BinDicSeg8,dir_list,"E8")
SegmentDicAssessor(BinDicSeg9,dir_list,"E9")
SegmentDicAssessor(BinDicSeg10,dir_list,"E10")

segments=["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"]
for segment in segments:
    df = pd.read_csv(histogram_path_csv+"/"+segment+".csv")
    print(df)
    for i in range (0,len(reg_list)):
        mean = df[reg_list[i]].mean()
        standart_deviation=df[reg_list[i]].std()
        one_STD_aboveTheMean = mean + standart_deviation
        one_STD_belowTheMean = mean - standart_deviation
        value_counts = df[reg_list[i]].value_counts()
        value_counts_df = pd.DataFrame({'value': value_counts.index, 'count': value_counts.values})
        plot_counts = sns.distplot(data=value_counts_df,x='value',y="count")
        plt.axvline(x=one_STD_aboveTheMean, ymin=0, ymax=1000,linestyle = "--", color="green")
        plt.axvline(x=one_STD_belowTheMean, ymin=0, ymax=1000,linestyle = "--", color="green")
        Figure= plot_counts.get_figure()
        Figure.savefig(histogram_path_png+"/"+reg_list[i]+"_"+ segment  +"_histogram.png")
        plt.clf()



'''
def histogrammer(intersectionPath,regCloumnName):
    for i in range (0,segment_num):
        reads=[]
        with open( intersectionPath[1],"r") as inputs:
            for line in inputs:
                line = line.strip().split("\t")
                if (int(line[3][1:]) == (i +1) ):
                    times= int((int(line[2])-int(line[1])) /5000)
                    for m in range(0,times):
                        if times != 0:
                            reads.append(float(int(line[4])/times))
        Output = open(histogram_path_csv+"/"+intersectionPath[0]+"_E"+str(i+1)+".csv" , "w")
        Output.write("Bin_ID,Counts"+"\n")
        for k in range(0,len(reads)):
            Output.write("Bin_ID_"+str(k+1)+","+str(reads[k])+"\n") 
        Output.close()
        
        df = pd.read_csv(histogram_path_csv+"/"+intersectionPath[0]+"_E"+str(i+1)+".csv")
        #df_sorted = df.sort_values(['Counts'],ascending=False)
        standart_deviation= df['Counts'].std()
        mean= df['Count'].mean()
        one_STD_aboveTheMean=mean+standart_deviation
        one_STD_belowTheMean=mean-standart_deviation
        
        #Plot = sns.histplot(data=df,x='Bin_ID',y="Counts")
        #plt.axhline(y=one_STD_aboveTheMean, xmin=df['Counts'].min(), xmac=df['Counts'].max(),linestyle = "--", color="green")
        #plt.axhline(y=one_STD_belowTheMean, xmin=df['Counts'].min(), xmac=df['Counts'].max(),linestyle = "--", color="green")
        #Figure= Plot.get_figure()
        #Figure.savefig(histogram_path_png+"/"+intersectionPath[0]+"_E"+str(i+1)+".png")
        #plt.clf()
        
        #Plot_sorted = sns.histplot(data=df_sorted,x='Bin_ID',y="Counts")
        #plt.axhline(y=one_STD_aboveTheMean, xmin=df['Counts'].min(), xmac=df['Counts'].max(),linestyle = "--", color="green")
        #plt.axhline(y=one_STD_belowTheMean, xmin=df['Counts'].min(), xmac=df['Counts'].max(),linestyle = "--", color="green")
        #Figure_sorted= Plot_sorted.get_figure()
        #Figure_sorted.savefig(histogram_path_png+"/"+intersectionPath[0]+"_E"+str(i+1)+"_sorted.png")
        #plt.clf()

        value_counts = df['Counts'].value_counts()
        value_counts_df = pd.DataFrame({'value': value_counts.index, 'count': value_counts.values})
        plot_counts = sns.histplot(data=value_counts_df,x='value',y="count")
        plt.axhline(y=one_STD_aboveTheMean, ymin=df['Counts'].min(), ymax=df['Counts'].max(),linestyle = "--", color="green")
        plt.axhline(y=one_STD_belowTheMean, ymin=df['Counts'].min(), ymax=df['Counts'].max(),linestyle = "--", color="green")
        Figure= plot_counts.get_figure()
        Figure.savefig(histogram_path_png+"/"+intersectionPath[0]+"_E"+str(i+1)+"_histogram.png")
        plt.clf()

        axes = df.hist(['Counts'], by='Bin_ID',bins=15, layout=(2,2), legend=True, yrot=90,sharex=True,sharey=True, log=True, figsize=(6,6))
        for ax in axes.flatten():
            ax.set_xlabel('N')
            ax.set_ylabel('Count')
            ax.set_ylim(bottom=1,top=100) 
        Figure= axes.get_figure()
        Figure.savefig(histogram_path_png+"/"+intersectionPath[0]+"_E"+str(i+1)+"_idk.png")
        plt.clf()

        



for i in range(0,len(dir_list)) :
    histogrammer(dir_list[i],10)
'''