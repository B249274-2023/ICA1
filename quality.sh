#!/usr/bin/bash

#Unzip the fastqc output files
#unzip /*.zip

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
fail_threshold=3

# loop through FastQC output folders
for folder in "${output_dir}"/*_fastqc/; do
    # read the summary file and add information to the appropriate output file for each line
    while IFS=$'\t' read -r line; do
        outcome_status=$(echo "${line}" | awk '{print $1}')
        if [ "${outcome_status}" == "PASS" ]; then
            echo "${line}" >> "${output_dir}"/"${pass_file}"
        elif [ "${outcome_status}" == "WARN" ]; then
            echo "${line}" >> "${output_dir}"/"${warn_file}"
        elif [ "$outcome_status" == "FAIL" ]; then
            echo "${line}" >> "${output_dir}"/"${fail_file}"
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

