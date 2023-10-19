

ICAdir=${"/home/s2506305/Exercises/ICA1"}


#copying the fastq files into ICA folder (PWD)
cp -r /localdisk/data/BPSM/ICA1/fastq . 

#1 running quality check
#make directory for fastqc output files
mkdir ${ICAdir}/fastqc_output

#running qc on all the files in the fastq file
fastqc -o fastqc_output fastq/*


#2 assessing the number and quality of the raw sequenced data
#unzipping all the zipped output files from qc
unzip \*.zip



#3 aligning read pairs to reference genome
#copying the reference genome into the folder
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/ .




#4 generating counts data
cp /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed .





#5 generating output files of the statistical mean of the counts per gene for each group







#6 generating fold change data
