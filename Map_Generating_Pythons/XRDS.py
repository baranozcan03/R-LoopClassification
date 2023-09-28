import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

geneInputFilname = sys.argv[1]
xrIntersectionFilename = sys.argv[2]
xrSimIntersectionFilename = sys.argv[3]
dsIntersectionFilename = sys.argv[4]
dsSimIntersectionFilename = sys.argv[5]
xrcsvFilePath = sys.argv[6]
xrSimcsvFilePath = sys.argv[7]
dscsvFilePath = sys.argv[8]
dsSimcsvFilePath = sys.argv[9]
xrpngFilePath = sys.argv[10]
xrSimpngFilePath = sys.argv[11]
dspngFilePath= sys.argv[12]
dsSimpngFilePath= sys.argv[13]
xrreadNum= int(sys.argv[14])
xrSimreadNum= int(sys.argv[15])
dsreadNum= int(sys.argv[16])
dsSimreadNum= int(sys.argv[17])
relativeCSV= sys.argv[18]
ratioCSV= sys.argv[19]
relativePNG = sys.argv[20]
ratioPNG =sys.argv[21]

segmentLengths= [0,0,0,0,0,0,0,0,0,0]
xrNums= [0,0,0,0,0,0,0,0,0,0]
xrSimNums = [0,0,0,0,0,0,0,0,0,0]
dsNums= [0,0,0,0,0,0,0,0,0,0]
dsSimNums= [0,0,0,0,0,0,0,0,0,0]


with open( xrIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        xrNums[index] = xrNums[index] + int(line[4])

with open( xrSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        xrSimNums[index] = xrSimNums[index] + int(line[4])

with open( dsIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        dsNums[index] = dsNums[index] + int(line[4])

with open( dsSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        dsSimNums[index] = dsSimNums[index] + int(line[4])

with open( geneInputFilname , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        segmentLengths[index] += 5000


Output = open(xrcsvFilePath , "w")

Output.write("Data-Type,Segment,RPKM"+"\n")
for i in range(0,10):
    rpkm = xrNums[i] / ( ( segmentLengths[i] / 1000) * (xrreadNum /1000000))
    rpkmF = "{:.2f}".format(rpkm)
    newline = "XR-seq" + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
    Output.write(newline)

Output.close()


Output = open(xrSimcsvFilePath , "w")

Output.write("Data-Type,Segment,RPKM"+"\n")
for i in range(0,10):
    rpkm = xrSimNums[i] / ( ( segmentLengths[i] / 1000) * (xrSimreadNum /1000000))
    rpkmF = "{:.2f}".format(rpkm)
    newline = "XR-sim-seq" + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
    Output.write(newline)

Output.close()


Output = open(dscsvFilePath , "w")

Output.write("Data-Type,Segment,RPKM"+"\n")
for i in range(0,10):
    rpkm = dsNums[i] / ( ( segmentLengths[i] / 1000) * (dsreadNum /1000000))
    rpkmF = "{:.2f}".format(rpkm)
    newline = "DS-seq" + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
    Output.write(newline)

Output.close()


Output = open(dsSimcsvFilePath , "w")

Output.write("Data-Type,Segment,RPKM"+"\n")
for i in range(0,10):
    rpkm = dsSimNums[i] / ( ( segmentLengths[i] / 1000) * (dsSimreadNum  /1000000))
    rpkmF = "{:.2f}".format(rpkm)
    newline = "DS-sim-seq" + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
    Output.write(newline)

Output.close()



df = pd.read_csv(xrcsvFilePath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrpngFilePath)
plt.clf()

df = pd.read_csv(xrcsvFilePath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrpngFilePath[:-4]+"_sorted.png")
plt.clf()


df = pd.read_csv(xrSimcsvFilePath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrSimpngFilePath)
plt.clf()

df = pd.read_csv(xrSimcsvFilePath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrSimpngFilePath[:-4]+"_sorted.png")
plt.clf()


df = pd.read_csv(dscsvFilePath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dspngFilePath)
plt.clf()

df = pd.read_csv(dscsvFilePath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dspngFilePath[:-4]+"_sorted.png")
plt.clf()


df = pd.read_csv(dsSimcsvFilePath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dsSimpngFilePath)
plt.clf()

df = pd.read_csv(dsSimcsvFilePath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dsSimpngFilePath[:-4]+"_sorted.png")
plt.clf()


Output = open(relativeCSV , "w")

Output.write("Data-Type,Segment,Relative-Repair"+"\n")
for i in range(0,10):
    rpkmX = xrNums[i] / ( ( segmentLengths[i] / 1000) * (xrreadNum /1000000))
    rpkmSX = xrSimNums[i] / ( ( segmentLengths[i] / 1000) * (xrSimreadNum /1000000))
    rpkmD = dsNums[i] / ( ( segmentLengths[i] / 1000) * (dsreadNum /1000000))
    rpkmDX = dsSimNums[i] / ( ( segmentLengths[i] / 1000) * (dsSimreadNum /1000000))
    relative = ( rpkmX / rpkmSX ) / (rpkmD / rpkmDX)
    relative = "{:.2f}".format(relative)
    newline = "Relative_Repair" + "," + "E" +str(i+1)  + "," + relative + "\n"
    Output.write(newline)

Output.close()



Output = open(ratioCSV , "w")

Output.write("Data-Type,Segment,Normalised-Damage"+"\n")
for i in range(0,10):
    rpkmD = dsNums[i] / ( ( segmentLengths[i] / 1000) * (dsreadNum /1000000))
    rpkmDX = dsSimNums[i] / ( ( segmentLengths[i] / 1000) * (dsSimreadNum /1000000))
    ratio = rpkmD / rpkmDX
    ratio = "{:.2f}".format(ratio)
    newline = "Normalised_Damage" + "," + "E" +str(i+1)  + "," + ratio + "\n"
    Output.write(newline)

Output.close()



df = pd.read_csv(relativeCSV)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative-Repair'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(relativePNG)
plt.clf()



df = pd.read_csv(relativeCSV)
df = df.sort_values(['Relative-Repair'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Relative-Repair'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(relativePNG[:-4]+"_sorted.png")
plt.clf()


df = pd.read_csv(ratioCSV)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Normalised-Damage'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(ratioPNG)
plt.clf()


df = pd.read_csv(ratioCSV)
df = df.sort_values(['Normalised-Damage'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['Normalised-Damage'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(ratioPNG[:-4]+"_sorted.png")
plt.clf()


Output = open(xrcsvFilePath[:-4] +"_mean.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( xrIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (xrreadNum /1000000))
        newline = "RPKM_mean"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(xrcsvFilePath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Data'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrpngFilePath[:-4]+"_mean.png")
plt.clf()

Output = open(xrSimcsvFilePath[:-4] +"_mean.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( xrSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (xrSimreadNum /1000000))
        newline = "RPKM_mean"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(xrSimcsvFilePath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Data'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrSimpngFilePath[:-4]+"_mean.png")
plt.clf()

Output = open(dscsvFilePath[:-4] +"_mean.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( dsIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (dsreadNum /1000000))
        newline = "RPKM_mean"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(dscsvFilePath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Data'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dspngFilePath[:-4]+"_mean.png")
plt.clf()

Output = open(dsSimcsvFilePath[:-4] +"_mean.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( dsSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (dsSimreadNum /1000000))
        newline = "RPKM_mean"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(dsSimcsvFilePath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Data'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dsSimpngFilePath[:-4]+"_mean.png")
plt.clf()





Output = open(xrcsvFilePath[:-4] +"_median.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( xrIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (xrreadNum /1000000))
        newline = "RPKM_median"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(xrcsvFilePath[:-4] +"_median.csv")
df = df.groupby(['Segment','Data'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrpngFilePath[:-4]+"_median.png")
plt.clf()

Output = open(xrSimcsvFilePath[:-4] +"_median.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( xrSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (xrSimreadNum /1000000))
        newline = "RPKM_median"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(xrSimcsvFilePath[:-4] +"_median.csv")
df = df.groupby(['Segment','Data'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(xrSimpngFilePath[:-4]+"_median.png")
plt.clf()

Output = open(dscsvFilePath[:-4] +"_median.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( dsIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (dsreadNum /1000000))
        newline = "RPKM_median"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(dsIntersectionFilename[:-4] +"_median.csv")
df = df.groupby(['Segment','Data'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dspngFilePath[:-4]+"_median.png")
plt.clf()

Output = open(dsSimcsvFilePath[:-4] +"_median.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( dsSimIntersectionFilename , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * (dsSimreadNum /1000000))
        newline = "RPKM_median"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(dsSimcsvFilePath[:-4] +"_median.csv")
df = df.groupby(['Segment','Data'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(dsSimpngFilePath[:-4]+"_median.png")
plt.clf()