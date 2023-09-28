import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

geneInputFilname = sys.argv[1]
ATRintersection = sys.argv[2]
BRCA1intersection = sys.argv[3]
DDX21intersection = sys.argv[4]
FIP1L1intersection = sys.argv[5]
NFAT_5intersection = sys.argv[6]
PRMT5intersection = sys.argv[7]
RAD51intersection = sys.argv[8]
RPA1intersection = sys.argv[9]
SETXintersection = sys.argv[10]
SIRT7intersection = sys.argv[11]
SMC3intersection = sys.argv[12]
SRPK2intersection = sys.argv[13]
SRSF1intersection = sys.argv[14]
XRN2intersection = sys.argv[15]
ATRcsv = sys.argv[16]
BRCA1csv = sys.argv[17]
DDX21csv = sys.argv[18]
FIP1L1csv = sys.argv[19]
NFAT_5csv = sys.argv[20]
PRMT5csv = sys.argv[21]
RAD51csv = sys.argv[22]
RPA1csv = sys.argv[23]
SETXcsv = sys.argv[24]
SIRT7csv = sys.argv[25]
SMC3csv = sys.argv[26]
SRPK2csv = sys.argv[27]
SRSF1csv = sys.argv[28]
XRN2csv = sys.argv[29]
ATRpng = sys.argv[30]
BRCA1png = sys.argv[31]
DDX21png = sys.argv[32]
FIP1L1png = sys.argv[33]
NFAT_5png = sys.argv[34]
PRMT5png = sys.argv[35]
RAD51png = sys.argv[36]
RPA1png = sys.argv[37]
SETXpng = sys.argv[38]
SIRT7png = sys.argv[39]
SMC3png = sys.argv[40]
SRPK2png = sys.argv[41]
SRSF1png = sys.argv[42]
XRN2png = sys.argv[43]
line1=int(sys.argv[44])
line2=int(sys.argv[45])
line3=int(sys.argv[46])
line4=int(sys.argv[47])
line5=int(sys.argv[48])
line6=int(sys.argv[49])
line7=int(sys.argv[50])
line8=int(sys.argv[51])
line9=int(sys.argv[52])
line10=int(sys.argv[53])
line11=int(sys.argv[54])
line12=int(sys.argv[55])
line13=int(sys.argv[56])
line14=int(sys.argv[57])
totalcsvPath=sys.argv[58]
totalpngPat=sys.argv[59]
ATACintersection=sys.argv[60]
ATACcsv=sys.argv[61]
totalnormalizedcsvPat_atac=sys.argv[62]
totalnormalizedPNG_atac=sys.argv[63]
line15=int(sys.argv[64])
RNAintersection=sys.argv[65]
RNAcsv=sys.argv[66]
totalnormalizedcsvPat_rna=sys.argv[67]
totalnormalizedPNG_rna=sys.argv[68]
line16=int(sys.argv[69])



segmentLengths= [0,0,0,0,0,0,0,0,0,0]
ATRnums= [0,0,0,0,0,0,0,0,0,0]
BRCA1nums= [0,0,0,0,0,0,0,0,0,0]
DDX21nums= [0,0,0,0,0,0,0,0,0,0]
FIP1L1nums= [0,0,0,0,0,0,0,0,0,0]
NFAT_5nums= [0,0,0,0,0,0,0,0,0,0]
PRMT5nums= [0,0,0,0,0,0,0,0,0,0]
RAD51nums= [0,0,0,0,0,0,0,0,0,0]
RPA1nums= [0,0,0,0,0,0,0,0,0,0]
SETXnums= [0,0,0,0,0,0,0,0,0,0]
SIRT7nums= [0,0,0,0,0,0,0,0,0,0]
SMC3nums= [0,0,0,0,0,0,0,0,0,0]
SRPK2nums= [0,0,0,0,0,0,0,0,0,0]
SRSF1nums= [0,0,0,0,0,0,0,0,0,0]
XRN2nums= [0,0,0,0,0,0,0,0,0,0]
ATACnums= [0,0,0,0,0,0,0,0,0,0]
RNAnums= [0,0,0,0,0,0,0,0,0,0]

def calculateNums(intersectionPath , numList ):
    with open( intersectionPath , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        numList[index] += int(line[4])

calculateNums(ATRintersection,ATRnums)
calculateNums(BRCA1intersection,BRCA1nums)
calculateNums(DDX21intersection,DDX21nums)
calculateNums(FIP1L1intersection,FIP1L1nums)
calculateNums(NFAT_5intersection,NFAT_5nums)
calculateNums(PRMT5intersection,PRMT5nums)
calculateNums(RAD51intersection,RAD51nums)
calculateNums(RPA1intersection,RPA1nums)
calculateNums(SETXintersection,SETXnums)
calculateNums(SIRT7intersection,SIRT7nums)
calculateNums(SMC3intersection,SMC3nums)
calculateNums(SRPK2intersection,SRPK2nums)
calculateNums(SRSF1intersection,SRSF1nums)
calculateNums(XRN2intersection,XRN2nums)
calculateNums(ATACintersection,ATACnums)
calculateNums(RNAintersection,RNAnums)


with open( geneInputFilname , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        segmentLengths[index] += 5000


def csvCalculator(csvPath,numList,SegLen,readNum,regName):
    Output = open(csvPath , "w")
    for i in range(0,10):
        rpkm = numList[i] / ( ( SegLen[i] / 1000) * (readNum /1000000))
        numList[i] = rpkm
        rpkmF = "{:.2f}".format(rpkm)
        newline = regName + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
        Output.write(newline)
    Output.close()

csvCalculator(ATRcsv,ATRnums,segmentLengths,line1,"ATR")
csvCalculator(BRCA1csv,BRCA1nums,segmentLengths,line2,"BRCA1")
csvCalculator(DDX21csv,DDX21nums,segmentLengths,line3,"DDX21")
csvCalculator(FIP1L1csv,FIP1L1nums,segmentLengths,line4,"FIP1L1")
csvCalculator(NFAT_5csv,NFAT_5nums,segmentLengths,line5,"NFAT_5")
csvCalculator(PRMT5csv,PRMT5nums,segmentLengths,line6,"PRMT5")
csvCalculator(RAD51csv,RAD51nums,segmentLengths,line7,"RAD51")
csvCalculator(RPA1csv,RPA1nums,segmentLengths,line8,"RPA1")
csvCalculator(SETXcsv,SETXnums,segmentLengths,line9,"SETX")
csvCalculator(SIRT7csv,SIRT7nums,segmentLengths,line10,"SIRT7")
csvCalculator(SMC3csv,SMC3nums,segmentLengths,line11,"SMC3")
csvCalculator(SRPK2csv,SRPK2nums,segmentLengths,line12,"SRPK2")
csvCalculator(SRSF1csv,SRSF1nums,segmentLengths,line13,"SRSF1")
csvCalculator(XRN2csv,XRN2nums,segmentLengths,line14,"XRN2")
csvCalculator(ATACcsv,ATACnums,segmentLengths,line15,"ATAC-seq")
csvCalculator(RNAcsv,RNAnums,segmentLengths,line16,"RNA-seq")

def writeOnTop(fromcsv):
    with open( fromcsv , "r" ) as inputs:
        for line in inputs:
            line = line.strip() + "\n"
            Output.write(line)

Output = open(totalcsvPath, "w")
Output.write("Data-Type,Segment,RPKM"+"\n")
writeOnTop(ATRcsv)
writeOnTop(BRCA1csv)
writeOnTop(DDX21csv)
writeOnTop(FIP1L1csv)
writeOnTop(NFAT_5csv)
writeOnTop(PRMT5csv)
writeOnTop(RAD51csv)
writeOnTop(RPA1csv)
writeOnTop(SETXcsv)
writeOnTop(SIRT7csv)
writeOnTop(SMC3csv)
writeOnTop(SRPK2csv)
writeOnTop(SRSF1csv)
writeOnTop(XRN2csv)
Output.close()


df = pd.read_csv(totalcsvPath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalpngPat)
plt.clf()

df = pd.read_csv(totalcsvPath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalpngPat[:-4]+"_sorted.png")
plt.clf()

'''
ATRnums_normalised= [0,0,0,0,0,0,0,0,0,0]
BRCA1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
DDX21nums_normalised= [0,0,0,0,0,0,0,0,0,0]
FIP1L1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
NFAT_5num_normaliseds= [0,0,0,0,0,0,0,0,0,0]
PRMT5nums_normalised= [0,0,0,0,0,0,0,0,0,0]
RAD51nums_normalised= [0,0,0,0,0,0,0,0,0,0]
RPA1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SETXnums_normalised= [0,0,0,0,0,0,0,0,0,0]
SIRT7nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SMC3nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SRPK2nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SRSF1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
XRN2nums_normalised= [0,0,0,0,0,0,0,0,0,0]


def normalise(raw,reference,new):
    for i in range(0,10):
        new[i] = raw[i] / reference[i]



normalise(ATRnums,ATACnums,ATRnums_normalised)
normalise(BRCA1nums,ATACnums,BRCA1nums_normalised)
normalise(DDX21nums,ATACnums,DDX21nums_normalised)
normalise(FIP1L1nums,ATACnums,FIP1L1nums_normalised)
normalise(NFAT_5nums,ATACnums,NFAT_5num_normaliseds)
normalise(PRMT5nums,ATACnums,PRMT5nums_normalised)
normalise(RAD51nums,ATACnums,RAD51nums_normalised)
normalise(RPA1nums,ATACnums,RPA1nums_normalised)
normalise(SETXnums,ATACnums,SETXnums_normalised)
normalise(SIRT7nums,ATACnums,SIRT7nums_normalised)
normalise(SMC3nums,ATACnums,SMC3nums_normalised)
normalise(SRPK2nums,ATACnums,SRPK2nums_normalised)
normalise(SRSF1nums,ATACnums,SRSF1nums_normalised)
normalise(XRN2nums,ATACnums,XRN2nums_normalised)


def writerList(list,protein):
    for i in range (0,10):
        newline= protein + "," + "E" + str(i+1)+ "," + str(list[i]) + "\n"
        Output.write(newline)


Output = open(totalnormalizedcsvPat_atac, "w")
Output.write("Data-Type,Segment,Relative_ATAC"+"\n")
writerList(ATRnums_normalised,"ATR")
writerList(BRCA1nums_normalised,"BRCA1")
writerList(DDX21nums_normalised,"DDX21")
writerList(FIP1L1nums_normalised,"FIP1L1")
writerList(NFAT_5num_normaliseds,"NFAT-5")
writerList(PRMT5nums_normalised,"PRMT5")
writerList(RAD51nums_normalised,"RAD51")
writerList(RPA1nums_normalised,"RPA1")
writerList(SETXnums_normalised,"SETX")
writerList(SIRT7nums_normalised,"SIRT7")
writerList(SMC3nums_normalised,"SMC3")
writerList(SRPK2nums_normalised,"SRPK2")
writerList(SRSF1nums_normalised,"SRSF1")
writerList(XRN2nums_normalised,"XRN2")
Output.close()



df = pd.read_csv(totalnormalizedcsvPat_atac)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative_ATAC'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalnormalizedPNG_atac)
plt.clf()

df = pd.read_csv(totalnormalizedcsvPat_atac)
df = df.sort_values(['Relative_ATAC'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative_ATAC'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalnormalizedPNG_atac[:-4]+"_sorted.png")
plt.clf()




ATRnums_normalised= [0,0,0,0,0,0,0,0,0,0]
BRCA1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
DDX21nums_normalised= [0,0,0,0,0,0,0,0,0,0]
FIP1L1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
NFAT_5num_normaliseds= [0,0,0,0,0,0,0,0,0,0]
PRMT5nums_normalised= [0,0,0,0,0,0,0,0,0,0]
RAD51nums_normalised= [0,0,0,0,0,0,0,0,0,0]
RPA1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SETXnums_normalised= [0,0,0,0,0,0,0,0,0,0]
SIRT7nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SMC3nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SRPK2nums_normalised= [0,0,0,0,0,0,0,0,0,0]
SRSF1nums_normalised= [0,0,0,0,0,0,0,0,0,0]
XRN2nums_normalised= [0,0,0,0,0,0,0,0,0,0]
normalise(ATRnums,RNAnums,ATRnums_normalised)
normalise(BRCA1nums,RNAnums,BRCA1nums_normalised)
normalise(DDX21nums,RNAnums,DDX21nums_normalised)
normalise(FIP1L1nums,RNAnums,FIP1L1nums_normalised)
normalise(NFAT_5nums,RNAnums,NFAT_5num_normaliseds)
normalise(PRMT5nums,RNAnums,PRMT5nums_normalised)
normalise(RAD51nums,RNAnums,RAD51nums_normalised)
normalise(RPA1nums,RNAnums,RPA1nums_normalised)
normalise(SETXnums,RNAnums,SETXnums_normalised)
normalise(SIRT7nums,RNAnums,SIRT7nums_normalised)
normalise(SMC3nums,RNAnums,SMC3nums_normalised)
normalise(SRPK2nums,RNAnums,SRPK2nums_normalised)
normalise(SRSF1nums,RNAnums,SRSF1nums_normalised)
normalise(XRN2nums,RNAnums,XRN2nums_normalised)
Output = open(totalnormalizedcsvPat_rna, "w")
Output.write("Data-Type,Segment,Relative_RNA"+"\n")
writerList(ATRnums_normalised,"ATR")
writerList(BRCA1nums_normalised,"BRCA1")
writerList(DDX21nums_normalised,"DDX21")
writerList(FIP1L1nums_normalised,"FIP1L1")
writerList(NFAT_5num_normaliseds,"NFAT-5")
writerList(PRMT5nums_normalised,"PRMT5")
writerList(RAD51nums_normalised,"RAD51")
writerList(RPA1nums_normalised,"RPA1")
writerList(SETXnums_normalised,"SETX")
writerList(SIRT7nums_normalised,"SIRT7")
writerList(SMC3nums_normalised,"SMC3")
writerList(SRPK2nums_normalised,"SRPK2")
writerList(SRSF1nums_normalised,"SRSF1")
writerList(XRN2nums_normalised,"XRN2")
Output.close()

df = pd.read_csv(totalnormalizedcsvPat_rna)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative_RNA'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalnormalizedPNG_rna)
plt.clf()

df = pd.read_csv(totalnormalizedcsvPat_rna)
df = df.sort_values(['Relative_RNA'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative_RNA'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalnormalizedPNG_rna[:-4]+"_sorted.png")
plt.clf()
'''

Output = open(totalcsvPath[:-4] +"_median.csv" , "w")
Output.write("Data-Type,Regulator,CHR,Start,End,Segment,RPKM"+"\n")
with open( ATRintersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line1 /1000000))
        newline = "RPKM_median"+","+"ATR"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( BRCA1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line2 /1000000))
        newline = "RPKM_median"+","+"BRCA1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( DDX21intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line3 /1000000))
        newline = "RPKM_median"+","+"DDX21"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( FIP1L1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line4 /1000000))
        newline = "RPKM_median"+","+"FIP1L1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( NFAT_5intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line5 /1000000))
        newline = "RPKM_median"+","+"NFAT_5"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( PRMT5intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line6 /1000000))
        newline = "RPKM_median"+","+"PRMT5"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( RAD51intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line7 /1000000))
        newline = "RPKM_median"+","+"RAD51"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( RPA1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line8 /1000000))
        newline = "RPKM_median"+","+"RPA1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SETXintersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line9 /1000000))
        newline = "RPKM_median"+","+"SETX"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SIRT7intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line10 /1000000))
        newline = "RPKM_median"+","+"SIRT7"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SMC3intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line11 /1000000))
        newline = "RPKM_median"+","+"SMC3"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SRPK2intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line12 /1000000))
        newline = "RPKM_median"+","+"SRPK2"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SRSF1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line13 /1000000))
        newline = "RPKM_median"+","+"SRSF1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( XRN2intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line15 /1000000))
        newline = "RPKM_median"+","+"XRN2"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()

df = pd.read_csv(totalcsvPath[:-4] +"_median.csv")
df = df.groupby(['Segment','Regulator'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)
Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalpngPat[:-4]+"_median.png")
plt.clf()












Output = open(totalcsvPath[:-4] +"_mean.csv" , "w")
Output.write("Data-Type,Regulator,CHR,Start,End,Segment,RPKM"+"\n")
with open( ATRintersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line1 /1000000))
        newline = "RPKM_mean"+","+"ATR"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( BRCA1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line2 /1000000))
        newline = "RPKM_mean"+","+"BRCA1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( DDX21intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line3 /1000000))
        newline = "RPKM_mean"+","+"DDX21"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( FIP1L1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line4 /1000000))
        newline = "RPKM_mean"+","+"FIP1L1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( NFAT_5intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line5 /1000000))
        newline = "RPKM_mean"+","+"NFAT_5"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( PRMT5intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line6 /1000000))
        newline = "RPKM_mean"+","+"PRMT5"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( RAD51intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line7 /1000000))
        newline = "RPKM_mean"+","+"RAD51"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( RPA1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line8 /1000000))
        newline = "RPKM_mean"+","+"RPA1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SETXintersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line9 /1000000))
        newline = "RPKM_mean"+","+"SETX"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SIRT7intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line10 /1000000))
        newline = "RPKM_mean"+","+"SIRT7"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SMC3intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line11 /1000000))
        newline = "RPKM_mean"+","+"SMC3"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SRPK2intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line12 /1000000))
        newline = "RPKM_mean"+","+"SRPK2"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( SRSF1intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line13 /1000000))
        newline = "RPKM_mean"+","+"SRSF1"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
with open( XRN2intersection , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (line15 /1000000))
        newline = "RPKM_mean"+","+"XRN2"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()

df = pd.read_csv(totalcsvPath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Regulator'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)
Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.4,'pad':0.04,'aspect': 7},annot_kws={"size":7})
Figure= Plot.get_figure()
Figure.savefig(totalpngPat[:-4]+"_mean.png")
plt.clf()
