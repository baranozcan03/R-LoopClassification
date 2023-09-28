#!/bin/bash
#
# -= Resources =-
#SBATCH --job-name=ft  
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
##SBATCH --mem 240G 
#SBATCH --account=mdbf
#SBATCH --qos=short_mdbf
#SBATCH --partition=short_mdbf
#SBATCH --time=0-2:00:00
#SBATCH --output=model1.out 
#SBATCH --mail-type=ALL
##SBATCH --mail-user=baran.ozcan@sabanciuniv.edu 


cd
java -mx4000M -jar ChromHMMapp/ChromHMM/ChromHMM.jar BinarizeBed \
        -b 5000 \
        -f 2 \
        ChromHMMapp/ChromHMM/CHROMSIZES/hg38.txt \
        ChromHMM/Data_finished_simple \
        ChromHMM/cellmarkfiletable.txt \
        ChromHMM/Binary
       #-p possionthreshold  

java -mx4000M -jar ChromHMMapp/ChromHMM/ChromHMM.jar LearnModel \
    -b 5000 \
    -s 316 \
    -r 3000 \
    -holdcolumnorder \
    -init random \
    -printstatebyline \
    -printposterior \
    ChromHMM/Binary \
    ChromHMM/Conclusion \
    10\
    hg38
    #5-10 bin 
    #10 state hatta 5 bile düşünülebilir. (en son 10 a karar verildi)


