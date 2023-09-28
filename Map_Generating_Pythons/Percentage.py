import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

segmentsFilePath= sys.argv[1]
outfilePath= sys.argv[2]
pngOutPath= sys.argv[3]
total=0
segments=[0,0,0,0,0,0,0,0,0,0]

with open( segmentsFilePath , "r" ) as inputs:
     for line in inputs:
        line=line.strip().split("\t")
        index=int(line[3][1:])-1
        segments[index]+=5000
        total+=5000


Output = open(outfilePath , "w")

Output.write("Data-Type,Segment,RoundedPercentages"+"\n")
for i in range (0,10):
    newline = "Percentages" + "," +"E" + str(i+1)  + "," + "{:.2f}".format((segments[i] * 100) / total) + "\n"
    Output.write(newline)
Output.close()

Percentages = pd.read_csv(outfilePath)
Percentages_table= Percentages.groupby(['Data-Type','Segment'],sort=False)['RoundedPercentages'].sum().unstack('Segment')
Plot = sns.heatmap(Percentages_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngOutPath)
plt.clf()

df = pd.read_csv(outfilePath)
df = df.sort_values(['RoundedPercentages'],ascending=False)
df_table= df.groupby(['Data-Type','Segment'],sort=False)['RoundedPercentages'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=True,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(pngOutPath[:-4]+"_sorted.png")
plt.clf()