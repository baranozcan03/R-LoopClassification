import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

transformation_data = sys.argv[1]
emission_data = sys.argv[2]
transformation_csv = sys.argv[3]
emssion_csv = sys.argv[4]
trnasformation_png = sys.argv[5]
emission_png = sys.argv[6]

transformations = []

with open( transformation_data,"r") as inputs:
    next(inputs)
    for line in inputs:
        line = line.strip().split("\t")
        transformations.append(line[1:]) 


for i in range(0,10):
    for k in range(0,10):
        transformations[i][k] = float(transformations[i][k]) 


Output = open(transformation_csv , "w")

Output.write("From,To,Probability"+"\n")
for i in range(0,10):
    for k in range(0,10):
        newline= "E"+str(i+1)+","+"E"+str(k+1)+","+str(transformations[i][k])+"\n"
        Output.write(newline)

Output.close()


df = pd.read_csv(transformation_csv)
df_table= df.groupby(['From','To'],sort=False)['Probability'].sum().unstack('To')
Plot = sns.heatmap(df_table,square=True,annot=False,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(trnasformation_png)
plt.clf()

df = pd.read_csv(transformation_csv)
df = df.sort_values(['Probability'],ascending=False)
df_table= df.groupby(['From','To'],sort=False)['Probability'].sum().unstack('To')
Plot = sns.heatmap(df_table,square=True,annot=False,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(trnasformation_png[:-4]+"_sorted.png")
plt.clf()


emissions = []


with open( emission_data,"r") as inputs:
    tags = inputs.readline()
    tags = tags.strip().split("\t")
    tags = tags[1:]
    for l in tags:
        emissions.append(l)
    lines_list=[]
    for line in inputs:
            line = line.strip().split("\t")
            lines_list.append(line)

emissions_new=[]
length=len(emissions)
for k in range(0,length):   
    for l in range(0,10):
        emissions_new.append(str(emissions[k])+","+"E"+str(l+1)+","+str(lines_list[l][k+1]))        


print(emissions_new)




Output = open(emssion_csv , "w")
Output.write("CellLine_Regulator,Segment,Probability"+"\n")
for i in emissions_new:
    Output.write(i+"\n") 
Output.close()

df = pd.read_csv(emssion_csv)
df_table= df.groupby(['CellLine_Regulator','Segment'],sort=False)['Probability'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=False,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(emission_png)
plt.clf()

df = pd.read_csv(emssion_csv)
df = df.sort_values(['Probability'],ascending=False)
df_table= df.groupby(['CellLine_Regulator','Segment'],sort=False)['Probability'].sum().unstack('Segment')
Plot = sns.heatmap(df_table,square=True,annot=False,cmap="Blues",cbar_kws={'shrink': 0.1,'pad':0.04,'aspect': 4})
Figure= Plot.get_figure()
Figure.savefig(emission_png[:-4]+"_sorted.png")
plt.clf()
