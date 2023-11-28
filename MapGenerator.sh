#add relevant SBATCH commands here


##Settings##
Create_Bins="no" #create a bed file that has all the bins ID'ed for use. (yes/no)
Create_Intersections="no" #Create Intersections for other things to work. (yes/no)
Create_Percentage="no" #Do you want percentage heat map with current segments file. (yes/no)
Create_ATAC="no" #Do you want ATAC-seq Rpkm heat map with current segments file. (yes/no)
Create_RNA="no" #Do you want RNA-seq Rpkm heat map with current segments file. (yes/no)
Create_RLOOP="yes" #Do you want R-loop Rpkm heat map with current segments file. (yes/no)
Create_XrDS="no"  #Do you want repair/damage heat maps with current segments file. (yes/no)
Create_Regulator="no" # Do you want to create Regulator RPKM map with current segments file. (yes/no)
Create_TransitionEmission="no" #Do you want to create ChromHMM maps? (yes/no)
Create_Histograms="no" #do you want to create RNA-seq histograms per bin? (yes/no)
Create_Segment_Regulator_RloopIntersectionPercentage="no" #do you want to create For all segments, per bin regulator-rLoop intercestion percentage chart.  



##FileNames
Percentage_Filename="Percentages"
ATACseq_Filename="ATAC-seq"
RNAseq_Filename="RNA-seq"
RLOOP_Filename="R-Loops"
XR_Filename="XR"
DS_Filename="DS"
XR_Relative="XR_Relative"
DS_ratio="DS_ratio"
ATR_Filename="ATR"
BRCA1_Filename="BRCA1"
DDX21_Filename="DDX21"
FIP1L1_Filename="FIP1L1"
NFAT5_Filename="NFAT-5"
PRMT5_Filename="PRMT5"
RAD51_Filename="RAD51"
RPA1_Filename="RPA1"
SETX_Filename="SETX"
SIRT7_Filename="SIRT7"
SMC3_Filename="SMC3"
SRPK2_Filename="SRPK2"
SRSF1_Filename="SRSF1"
XRN2_Filename="XRN2"
Transformation_Filename="Transformation_matrix"
Emission_Filename="Emission_matrix"
Bin_Histogram="RNAseq_Bin_Histogram"


#BedFileLocations
#Please Check all bed files in the explorer and Regulator files down below
##Non_regulator
Intersections_file_path="/cta/users/baran.ozcan/HeatMaps/Intersections"
Segment_file_path_raw="/cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/genome_10_316_segments.bed"
Segment_file_path="/cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/genome_binned_segments.bed" 
ATACseq_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/ATACseq_HeLa_hg38.bed"
RNAseq_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/RNA_seq.bed"
RLOOP_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/RL_regions_strassigned_uncommon_rmdup_hybrids_chr_mod.bed"
XR_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/HXACA4_TGACCA_hg38_XR.bed"
XR_sim_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/HXACA4_TGACCA_hg38_XR_sim.bed"
DS_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/HDACA6_GCCAAT_hg38_DS.bed"
DS_sim_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Non-Regulator/HDACA6_GCCAAT_hg38_DS_sim.bed"
transformation_file_path="/cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/transitions_10_316.txt"
emission_file_path="/cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/emissions_10_316.txt"
##Regulator UPDATE DOWN BELOW IF YOU ADD 
ATR_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/ATR.bed"
BRCA1_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/BRCA1.bed"
DDX21_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/DDX21.bed"
FIP1L1_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/FIP1L1.bed"
NFAT5_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/NFAT-5.bed"
PRMT5_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/PRMT5.bed"
RAD51_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/RAD51.bed"
RPA1_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/RPA1.bed"
SETX_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/SETX.bed"
SIRT7_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/SIRT7.bed"
SMC3_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/SMC3.bed"
SRPK2_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/SRPK2.bed"
SRSF1_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/SRSF1.bed"
XRN2_file_path="/cta/users/baran.ozcan/HeatMaps/Beds/Regulator/XRN2.bed"



###Code
if [ ${Create_Bins} == "yes" ] 
then
    python /cta/users/baran.ozcan/HeatMaps/Pythons/SegmentBinner.py ${Segment_file_path_raw} ${Segment_file_path}
fi
if [ ${Create_Intersections}  == "yes" ]
then
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${ATACseq_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${ATACseq_Filename}_intersection.bed -c
    bedtools intersect  -a ${Segment_file_path} -b ${RNAseq_file_path} \
    > /cta/users/baran.ozcan/HeatMaps/Intersections/${RNAseq_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${RLOOP_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${RLOOP_Filename}_intersection.bed -c
    #bedtools intersect -F 0.50 -a ${Segment_file_path} -b ${XR_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${XR_Filename}_intersection.bed -c
    #bedtools intersect -F 0.50 -a ${Segment_file_path} -b ${XR_sim_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${XR_Filename}_sim_intersection.bed -c
    #bedtools intersect -F 0.50 -a ${Segment_file_path} -b ${DS_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${DS_Filename}_intersection.bed -c
    #bedtools intersect -F 0.50 -a ${Segment_file_path} -b ${DS_sim_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${DS_Filename}_sim_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${ATR_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${ATR_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${BRCA1_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${BRCA1_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${DDX21_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${DDX21_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${FIP1L1_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${FIP1L1_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${NFAT5_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${NFAT5_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${PRMT5_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${PRMT5_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${RAD51_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${RAD51_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${RPA1_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${RPA1_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${SETX_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${SETX_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${SIRT7_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${SIRT7_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${SMC3_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${SMC3_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${SRPK2_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${SRPK2_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${SRSF1_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${SRSF1_Filename}_intersection.bed -c
    #bedtools intersect -F 0.30 -a ${Segment_file_path} -b ${XRN2_file_path} \
    #> /cta/users/baran.ozcan/HeatMaps/Intersections/${XRN2_Filename}_intersection.bed -c
fi
if [ ${Create_Percentage} == "yes" ]
then
    python /cta/users/baran.ozcan/HeatMaps/Pythons/Percentage.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${Percentage_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${Percentage_Filename}.png 
fi
if [ ${Create_ATAC} == "yes" ]
then
    ATAC_Line_Num="$(grep -c "^" ${ATACseq_file_path})"
    python /cta/users/baran.ozcan/HeatMaps/Pythons/ATACseqRPKM.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${ATACseq_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${ATACseq_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${ATACseq_Filename}.png \
    $ATAC_Line_Num
fi
if [ ${Create_RNA} == "yes" ]
then
    RNA_Line_Num="$(grep -c "^" ${RNAseq_file_path})"
    python /cta/users/baran.ozcan/HeatMaps/Pythons/RNAseqRPKM.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RNAseq_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${RNAseq_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${RNAseq_Filename}.png \
    $RNA_Line_Num
fi
if [ ${Create_RLOOP} == "yes" ]
then
    RLOOP_Line_Num="$(grep -c "^" ${RLOOP_file_path})"
    python /cta/users/baran.ozcan/HeatMaps/Pythons/RloopRPKM.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RLOOP_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${RLOOP_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${RLOOP_Filename}.png \
    $RLOOP_Line_Num
fi
if [ ${Create_XrDS} == "yes" ]
then
    repair_line1="$(grep -c "^" ${XR_file_path})"
    repair_line2="$(grep -c "^" ${XR_sim_file_path})"
    repair_line3="$(grep -c "^" ${DS_file_path})"
    repair_line4="$(grep -c "^" ${DS_sim_file_path})"
    python /cta/users/baran.ozcan/HeatMaps/Pythons/XRDS.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${XR_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${XR_Filename}_sim_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${DS_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${DS_Filename}_sim_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${XR_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${XR_Filename}_sim.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${DS_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${DS_Filename}_sim.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${XR_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${XR_Filename}_sim.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${DS_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${DS_Filename}_sim.png \
    ${repair_line1} ${repair_line2} ${repair_line3} ${repair_line4} \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${XR_Relative}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${DS_ratio}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${XR_Relative}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${DS_ratio}.png
fi
if [ ${Create_Regulator} == "yes" ]
then
    regulator_line1="$(grep -c "^" ${ATR_file_path})" 
    regulator_line2="$(grep -c "^" ${BRCA1_file_path})"
    regulator_line3="$(grep -c "^" ${DDX21_file_path})" 
    regulator_line4="$(grep -c "^" ${FIP1L1_file_path})"
    regulator_line5="$(grep -c "^" ${NFAT5_file_path})"
    regulator_line6="$(grep -c "^" ${PRMT5_file_path})"
    regulator_line7="$(grep -c "^" ${RAD51_file_path})"
    regulator_line8="$(grep -c "^" ${RPA1_file_path})"
    regulator_line9="$(grep -c "^" ${SETX_file_path})"
    regulator_line10="$(grep -c "^" ${SIRT7_file_path})"
    regulator_line11="$(grep -c "^" ${SMC3_file_path})"
    regulator_line12="$(grep -c "^" ${SRPK2_file_path})"
    regulator_line13="$(grep -c "^" ${SRSF1_file_path})"
    regulator_line14="$(grep -c "^" ${XRN2_file_path})"
    regulator_line15="$(grep -c "^" ${ATACseq_file_path})"
    regulator_line16="$(grep -c "^" ${RNAseq_file_path})"
    python /cta/users/baran.ozcan/HeatMaps/Pythons/RegulatorRPKM.py ${Segment_file_path} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${ATR_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${BRCA1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${DDX21_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${FIP1L1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${NFAT5_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${PRMT5_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RAD51_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RPA1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SETX_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SIRT7_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SMC3_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SRPK2_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SRSF1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${XRN2_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${ATR_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${BRCA1_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${DDX21_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${FIP1L1_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${NFAT5_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${PRMT5_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${RAD51_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${RPA1_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${SETX_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${SIRT7_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${SMC3_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${SRPK2_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${SRSF1_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${XRN2_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/${ATR_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${BRCA1_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${DDX21_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${FIP1L1_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${NFAT5_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${PRMT5_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${RAD51_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${RPA1_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${SETX_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${SIRT7_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${SMC3_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${SRPK2_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${SRSF1_Filename}.png \
    /cta/users/baran.ozcan/HeatMaps/Maps/${XRN2_Filename}.png \
    ${regulator_line1} ${regulator_line2} ${regulator_line3} ${regulator_line4} ${regulator_line5} \
    ${regulator_line6} ${regulator_line7} ${regulator_line8} ${regulator_line9} ${regulator_line10} \
    ${regulator_line11} ${regulator_line12} ${regulator_line13} ${regulator_line14} \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/total.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/total.png \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${ATACseq_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${ATACseq_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/total_ATAC_normalized.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/total_ATAC_normalized.png \
    ${regulator_line15} \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RNAseq_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/${RNAseq_Filename}.csv \
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/total_RNA_normalized.csv \
    /cta/users/baran.ozcan/HeatMaps/Maps/total_RNA_normalized.png \
    ${regulator_line16}
fi
if [ ${Create_TransitionEmission} == "yes" ]
then
    python /cta/users/baran.ozcan/HeatMaps/Pythons/ChromHMMplotter.py \
        /cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/transitions_10_316.txt \
        /cta/users/baran.ozcan/HeatMaps/ChromHMMstuff/emissions_10_316.txt \
        /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/transformation.csv \
        /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/emission.csv \
        /cta/users/baran.ozcan/HeatMaps/Maps/transformation.png \
        /cta/users/baran.ozcan/HeatMaps/Maps/emission.png 
fi
if [ ${Create_Histograms} == "yes" ]
then
    python /cta/users/baran.ozcan/HeatMaps/Pythons/Histograms.py ${Intersections_file_path}\
    /cta/users/baran.ozcan/HeatMaps/Maps/_Histograms\
    /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/_Histograms \
    ${Segment_file_path}
fi
if [ ${Create_Segment_Regulator_RloopIntersectionPercentage} == "yes" ]
then
    python /cta/users/baran.ozcan/HeatMaps/Pythons/StackPercentages_Rloop_Regulator.py \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RLOOP_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${ATR_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${BRCA1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${DDX21_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${FIP1L1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${NFAT5_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${PRMT5_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RAD51_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${RPA1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SETX_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SIRT7_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SMC3_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SRPK2_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${SRSF1_Filename}_intersection.bed \
    /cta/users/baran.ozcan/HeatMaps/Intersections/${XRN2_Filename}_intersection.bed \
    14 /cta/users/baran.ozcan/HeatMaps/CSV_Outcomes/_StackPercentages ${Segment_file_path}\
    /cta/users/baran.ozcan/HeatMaps/Maps/_StackPercentages \


fi

echo "ended"
#Rna-seq doğru omuş mu ayarla.
#CHROM HMM'den direkt mapleri tek regulatorde birleştir ve sırala
#median heapla rpkm
#peak call bin binarisation
