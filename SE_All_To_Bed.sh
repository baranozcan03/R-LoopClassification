#add relevant SBATCH commands here



###Modify with the corresponding number of SRA number count
Sample_FileName="XRN2_ChIP-seq_HeLa_SRR442116_Sample_noRep"
#name = (regName_ChIP-seq_cellType_SRAID1_sampleOrInput_repNumber)
Sample_Path="/cta/users/baran.ozcan/ChromHMM/Data_trimmed_all/${Sample_FileName}"
SRA_ID1="a"
SRA_ID2="a"
SRA_ID3="a"
SRA_ID4="a"

prefetch -O ${Sample_Path} ${SRA_ID1}
fastq-dump ${Sample_Path}/${SRA_ID1}/${SRA_ID1}.sra
prefetch -O ${Sample_Path} ${SRA_ID2}
fastq-dump ${Sample_Path}/${SRA_ID2}/${SRA_ID2}.sra
prefetch -O ${Sample_Path} ${SRA_ID3}
fastq-dump ${Sample_Path}/${SRA_ID3}/${SRA_ID3}.sra
prefetch -O ${Sample_Path} ${SRA_ID4}
fastq-dump ${Sample_Path}/${SRA_ID4}/${SRA_ID4}.sra

cd ${Sample_Path}

cat ${Sample_Path}/${SRA_ID1}/${SRA_ID1}.fastq \
    ${Sample_Path}/${SRA_ID2}/${SRA_ID2}.fastq\
    ${Sample_Path}/${SRA_ID3}/${SRA_ID3}.fastq \
    ${Sample_Path}/${SRA_ID4}/${SRA_ID4}.fastq \
    > ${Sample_Path}/${Sample_FileName}.fastq

fastp -i ${Sample_Path}/${Sample_FileName}.fastq -o ${Sample_Path}/${Sample_FileName}_trimmed.fq
fastqc ${Sample_Path}/${Sample_FileName}.fastq
fastqc ${Sample_Path}/${Sample_FileName}_trimmed.fq


cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.1.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.2.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.3.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.4.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.rev.1.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/bowtie2_index.rev.2.bt2 ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/Homo_sapiens.GRCh38.dna.toplevel.fa ${Sample_Path}
cp ./GeneAssemblies/Human/GRCH38/Homo_sapiens.GRCh38.dna.toplevel.fa.fai ${Sample_Path}
echo "GeneAssemblies Downloaded"

cd ${Sample_Path}
bowtie2 --very-fast-local -x bowtie2_index -U ${Sample_FileName}_trimmed.fq -S ${Sample_FileName}.sam
samtools view -q 20 -S -b ${Sample_FileName}.sam > ${Sample_FileName}_unsorted.bam
samtools sort ${Sample_FileName}_unsorted.bam > ${Sample_FileName}.bam
java -jar /cta/apps/opt/spack/linux-ubuntu16.04-skylake/gcc-9.2.0/picard-2.20.8-hehlgdlwak4s5uew4tuuwi353lx36mz3/bin/picard.jar MarkDuplicates  \
    REMOVE_DUPLICATES=true I=${Sample_FileName}.bam \
    O=${Sample_FileName}_picard.bam M=marked_dup_metrics.txt
bedtools bamtobed -i ${Sample_FileName}_picard.bam > ${Sample_FileName}_unsorted.bed
bedtools sort -i ${Sample_FileName}_unsorted.bed > ${Sample_FileName}.bed
awk '{print "chr"$0}' ${Sample_FileName}.bed > ${Sample_FileName}_labeled.bed 


