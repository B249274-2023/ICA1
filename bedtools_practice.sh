#!/bin/bash

mkdir counts_data

for file in alignment_output/*.sorted.bam; do
	filename=$(basename "${file}")
	filename="${filename%.soorted.bam}"
	output_file="${filename}_counts.txt"

	bedtools coverage -header -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file} > "counts_data/${output_file}" 
done






