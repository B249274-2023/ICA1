!#/bin/bash


# Define the folder containing the files
folder="your_folder_path_here"

# Define the output file for mean values
output_file="mean_values.txt"

# Define the file containing columns 3 and 4
additional_file="path_to_additional_file.txt"

# Initialize the output file
> "$output_file"

# Loop through each file in the folder
for file in "${combinations}"/*_counts.txt; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        
	# Define output file with a similar name as input file
        output_file="${file%.*}_mean.txt"

	# Calculate mean for each row in the file
        awk '{
            sum = 0;
            for (i = 1; i <= NF; i++) {
                sum += $i;
            }
            mean = sum / NF;
            print mean;
        }' "${file}" >> temp_file.txt
        
	# Combine with columns 3 and 4 from the additional file
        paste temp_file.txt <(awk '{print $3, $4}' "$additional_file") >> "$output_file"
        rm temp_file.txt
    fi
done



