#!/bin/bash


mkdir random_data

# Generate random data and save it to the file
for ((i = 0; i < 10; i++)); do
    echo $((RANDOM % 20 + 1)) $((RANDOM % 20 + 1)) $((RANDOM % 20 + 1))
done > random_data/random_file1.txt
# generate random file 2
for ((i = 0; i < 10; i++)); do
    echo $((RANDOM % 20 + 1)) $((RANDOM % 20 + 1)) $((RANDOM % 20 + 1))
done > random_data/random_file2.txt

#calculating the averages
bed_file="${TriTrypDB-46_TcongolenseIL3000_2019.bed}"
combos="${PWD}/combinations"

# loop through each file in the folder
for file in random_data/*.txt; do
        # define the output file with a similar name as input file
        output_file="${file%.*}_mean.txt"
        touch "${output_file}"
	# calculate mean for each row in the file
        awk '{
            sum = 0;
            for (i = 1; i <= NF; i++) {
                sum += $i;
            }
            mean = sum / NF;
            print mean;
        }' "${file}" >> "${output_file}"
done

