#!/bin/bash
# separate the groups into different conditions
input_dir="${PWD}/fastq"
input_file="Tco2.fqfiles"
mkdir combinations
combos="${PWD}/combinations"

# creating files that have the sample names for every different combination
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
    echo -e "${header}\n${output_files["$condition"]}" > "${combos}/${filename}"
done


for combination_file in "${combos}"/*; do
    # Create a temporary file to store the modified content
    tmp_file=$(mktemp)
    # Use awk to keep only the first column
    awk 'BEGIN{FS="\t"} NR>1 {print $1}' "${combination_file}" > "$tmp_file"
    
    # Replace the original file with the modified content
    mv "$tmp_file" "$combination_file"
done

for combo_file in "$combos"/*; do
    #Change subject names like Tco258 to Tco-258
   sed -i 's/Tco\([0-9]\+\)/Tco-\1/' "$combo_file"
done



# add the counts information for each 
# Create the output directory if it doesn't exist
counts_dir="${PWD}/counts_data"


# Iterate through files in the combinations directory
for combination_file in "${combos}"/*; do
    condition=$(basename "$combination_file")

    while IFS="\n" read -r sample_name; do
        echo "Processing sample: ${sample_name}"

    # get corresponding counts file   
    counts_data="$counts_dir/${sample_name}_counts.txt"
   
	if [ -e "$counts_data" ]; then
	   # create a temporary file for modified counts data
	   tmp_counts_data=$(mktemp)

    	   # extract field 6 and save it to the temporary file
	   awk -F'\t' -v OFS='\t' 'NR==FNR{a[NR]=$6; next} {print $0, a[FNR]}' "$counts_data" "$output_file" > "${output_file}.tmp"

	   # create a new output file for each condition
   	   output_file="${counts_dir}/${condition}_counts.txt"
   	   touch "${output_file}"

   	   # append counts data as a new column
   	   paste "${counts_data}" "${tmp_counts_data}" > "${output_file}"

   	   # remove temporary counts data file
   	   rm "${tmp_counts_data}"
	else
   	   echo "Counts data file not found for sample: ${sample_name}"
	fi
   done < "${combination_file}"
done



