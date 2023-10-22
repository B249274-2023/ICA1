#!/bin/bash

# making project directory and changing it to our PWD
mkdir ICA1_outputs
cd ICA1_outputs

# OPTIONAL – removing previous data from this folder if the script has previously been run – use with caution
#rm -r ICA1/*

#1 run quality check on raw sequence data
# copy the fastq files into PWD
cp -r /localdisk/data/BPSM/ICA1/fastq . 

# make directory for fastqc output files
mkdir fastqc_output
# run qc on all the files in the fastq file and output them in fastqc_output
fastqc -o fastqc_output fastq/*


#2 assess the number and quality of the raw sequenced data
# unzip all the zipped output files from qc
cd fastqc_output
unzip \*.zip
cd ..

# OPTIONAL – removing the .zip and .html files
#rm fastqc_output/*.zip
#rm fastqc_output/*.html

# define the output directories and files
output_dir="fastqc_output"
pass_file="pass_output.txt"
warn_file="warn_output.txt"
fail_file="fail_output.txt"
summary_file="summary.txt"

# create or clear the output files
> "${pass_file}"
> "${warn_file}"
> "${fail_file}"

# define the threshold for number of fails
fail_threshold=4

# loop through FastQC output folders
for folder in "${output_dir}"/*_fastqc/; do
    # read the summary file and add information to the appropriate output file for each line
    while IFS=$'\t' read -r line; do
        outcome_status=$(echo "${line}" | awk '{print $1}')
        if [ "${outcome_status}" == "PASS" ]; then
            echo "${line}" >> "${output_dir}/${pass_file}"
        elif [ "${outcome_status}" == "WARN" ]; then
            echo "${line}" >> "${output_dir}/${warn_file}"
        elif [ "$outcome_status" == "FAIL" ]; then
            echo "${line}" >> "${output_dir}/${fail_file}"
        fi
    done < "${folder}/summary.txt"

    # check if the summary file exists and contains 3 or more fails
    if [ -e "${folder}/${summary_file}" ]; then
        fails=$(grep -o "FAIL" "$folder/$summary_file" | wc -l)
        # if number of fails is equal to or  higher than 3, remove the file
        if [ "${fails}" -ge "${fail_threshold}" ]; then
            # extract the sequence name from the folder name
            sequence_name= $(basename "${folder}")
            sequence_name="${sequence_name%_fastqc}"

            # remove the sequence from Tco2.fqfiles
            sed -i "/${sequence_name}/d" Tco2.fqfiles
            # tells user that the sequence has been removed
            echo '"${sequence_name}" was removed' >> removed_sequences.txt
            # removes the sequence files
            rm -r "fastq/$sequence_name"*
        fi
    fi
done



#3 align read pairs to reference genome

# copy the reference genome into the PWD
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .
# make the reference genome from the .fasta file copied from the above folder
bowtie2-build Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz Tcongo_reference
# move files into the Tcongo_reference_genome folder (not neccessary but makes it neater)
mv Tcongo_reference.* Tcongo_genome

#make an output directory for the aligned sequences
mkdir alignment_output
# assigning the variable "reference" for use in the following for loop
reference="${Tcongo_genome/Tcongo_reference}"

# loop through all FASTQ files in the 'reads/' directory
for file1 in fastq/*_1.fq.gz; do
    file2=${file1/_1.fq.gz/_2.fq.gz} # generate the corresponding reverse file name

    filename=$(basename "$file1") # extracts the base filename from file1 by removing the path 
    filename="${filename%_1.fq.gz}" # removes the "_1.fq.gz" from filename leaving the sequence number

    # align the paired-end sequences using Bowtie2
    bowtie2 --local -x Tcongo_genome/Tcongo_reference -1 ${file1} -2 ${file2} -S alignment_output/${filename}.sam
    # convert .sam to .bam using samtools
    samtools view -bS "alignment_output/${filename}.sam" > "alignment_output/${filename}.bam"
    # sort alignments by leftmost coordinates
    samtools sort "alignment_output/${filename}.bam" -o "alignment_output/${filename}.sorted.bam"
done


#4 generate counts data
# copy the bed file into PWD
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed .

#make a new directory to output the information from bedtools
mkdir counts_data


for file in alignment_output/*.sorted.bam; do
	filename=$(basename "${file}")
        filename="${filename%.sorted.bam}"
        output_file="${filename}_counts.txt"

        bedtools coverage -counts -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file} > "counts_data/${output_file}"
done



#5 generate output files of the statistical mean of the counts per gene for each group

# separate the groups into different conditions
input_dir="${PWD}/fastq"
input_file="Tco2.fqfiles"
mkdir combinations


# read the header to get column names
read -r header < "${input_dir}/${input_file}"

declare -A output_files

# read through each line of the file
while IFS='' read -r line; do
    SampleName=$(echo "$line" | cut -f1)
    SampleType=$(echo "$line" | cut -f2)
    Time=$(echo "$line" | cut -f4)
    Treatment=$(echo "$line" | cut -f5)

    # skip the header line
    if [[ "$SampleName" == "SampleName" ]]; then
        continue
    fi

    # specifiy the conditions
    condition="${SampleType}_${Time}_${Treatment}"

    # add the line to the corresponding output file
    output_files["$condition"]+="$line\n"
done < <(tail -n +2 "${input_dir}/${input_file}")

# create the output files
for condition in "${!output_files[@]}"; do
    filename="${condition}.txt"
    echo -e "${header}\n${output_files["$condition"]}" > "combinations/${filename}"
done





#6 generate fold change data
