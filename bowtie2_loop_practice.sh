
#make an output directory for the aligned sequences
mkdir alignment_output
# assigning the variable "reference" for use in the following for loop
reference="${Tcongo_genome/Tcongo_reference}"

# loop through all FASTQ files in the 'reads/' directory
for file1 in prac/*_1.fq.gz; do
    file2=${file1/_1.fq.gz/_2.fq.gz} # Generate the corresponding reverse file name

    filename=$(basename "$file1")
    filename="${filename%_1.fq.gz}"

    # align the paired-end sequences using Bowtie2
    bowtie2 -x Tcongo_genome/Tcongo_reference -1 ${file1} -2 ${file2} -S alignment_output/${filename}.sam
    # convert .sam to .bam using samtools
    #samtools view -bS "alignment_output/${filename}.sam" > "alignment_output/${filename}.bam"
done
