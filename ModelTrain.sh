#add relevant SBATCH commands here


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


