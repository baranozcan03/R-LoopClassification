import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

segmentsFilePath= sys.argv[1]
intersectionFilePath= sys.argv[2]
csvFilePath= sys.argv[3]
pngFilePath= sys.argv[4]
readNum= int(sys.argv[5])

stateLength=[0,0,0,0,0,0,0,0,0,0]
StateAppear=[0,0,0,0,0,0,0,0,0,0]

with open( intersectionFilePath , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        index=int(line[3][1:])-1
        StateAppear[index] += int(line[4])
        stateLength[index] += 5000


Output = open(csvFilePath , "w")

Output.write("Data-Type,Segment,RPKM"+"\n")
for i in range(0,10):
    rpkm = StateAppear[i] / ( ( stateLength[i] / 1000) * 100000)
    rpkmF = "{:.2f}".format(rpkm)
    newline = "R-loop" + "," + "E" +str(i+1)  + "," + rpkmF + "\n"
    Output.write(newline)

Output.close()

df = pd.read_csv(csvFilePath)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngFilePath)
plt.clf()

df = pd.read_csv(csvFilePath)
df = df.sort_values(['RPKM'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RPKM'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngFilePath[:-4]+"_sorted.png")
plt.clf()

Output = open(csvFilePath[:-4] +"_median.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( intersectionFilePath , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * 1000000)
        newline = "RPKM_median"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(csvFilePath[:-4] +"_median.csv")
df = df.groupby(['Segment','Data'])['RPKM'].median().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngFilePath[:-4]+"_median.png")
plt.clf()


Output = open(csvFilePath[:-4] +"_mean.csv" , "w")
Output.write("Data,Chr,Start,End,Segment,RPKM"+"\n")
with open( intersectionFilePath , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        rpkm = int(line[4]) / ( ( 5000 / 1000) * 100000)
        newline = "RPKM_mean"+","+str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + "," +str(line[3][1:]) +","+ str(rpkm) + "\n"
        Output.write(newline)
Output.close()
df = pd.read_csv(csvFilePath[:-4] +"_mean.csv")
df = df.groupby(['Segment','Data'])['RPKM'].mean().unstack('Segment')
rpkm_df = pd.DataFrame(df)
rpkm_df = rpkm_df.sort_index(ascending=True)
print(rpkm_df)

Plot = sns.heatmap(rpkm_df,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngFilePath[:-4]+"_mean.png")
plt.clf()