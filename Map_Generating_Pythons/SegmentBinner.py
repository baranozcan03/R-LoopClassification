import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


raw = sys.argv[1]
out= sys.argv[2]



Output = open(out , "w")

with open( raw , "r" ) as inputs:
     for line in inputs:
        line = line.strip().split("\t")
        times = int((int(line[2])-int(line[1])) /5000)
        for m in range(0,times):
            start = (int(line[1])) + (5000*m)
            end= start + 5000
            Output.write(line[0] + "\t" + str(start) + "\t" + str(end) + "\t" + line[3] + "\n" )

Output.close()