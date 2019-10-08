
### A pipeline for the analysis of paired-end RNA-Seq  reads using fastqc, bowtie2, samtools and bedtools 

## Assumptions: 
# raw sequence files are in fq.gz format and stored in one folder
# a mapping file is located within this folder 
# a fasta format genome sequence is provided in the specified genome_dir  
# threading is used at an approproiate level for the active computer - please check the defaults before running. 

## Global defaults 

working_directory=BPSM_assignment1/Assignment1
threads=3

## OPTIONAL: check number of CPU cores and free memory and set to higher of these divided by 250 (fastqc threading)

### FastQC pipeline
## Pipeline defaults 


fastqc_data_dir="$working_directory/fastq"

#fastqc_cutoffs= 


# As fastqc doesn't create directories we want to create a folder to store the output in if it doesn't already exist:
# check the default folder name (FastQC_output) doesn't exist in working directory and if not, create it.

if [ ! -d "$working_directory/FastQC_output" ]; then
    mkdir "$working_directory/FastQC_output"
fi 


if [ ! -d "$working_directory/FastQC_summary" ]; then
    mkdir "$working_directory/FastQC_summary"
fi 
fastqc_out_dir="$working_directory/FastQC_output"

# for each file in the fastqc_data_dir
# if the file is a fasta file (excluding the map file which has no file extension)
for file in $(ls "$fastqc_data_dir/"*.fq.gz)
do 
# run fastqc on the data and output the results to the fastqc_out_dir
fastqc $file -o $fastqc_out_dir -t $threads
done

### fastqc outputs the data to FastQC_output/filename and creates many files in a zip file and a html file. 
## On extraction of the zip files, mutiple folders and files exist: 
##  fastqc_data.txt  fastqc.fo  fastqc_report.html  Icons  Images  summary.txt
## In summary.txt and pass/warn/fail report is generated with columns (tab seperated) of 'Pass  Test    File'

for outfile in $(find "$fastqc_out_dir/"*.zip -type f -printf "%f\n")
do
#unzip -p $fastqc_out_dir/216_L8_1_fastqc.zip *summary.txt >  "$working_directory/FastQC_summary/216_L8_1_fastqc.txt"
unzip -p "$fastqc_out_dir/$outfile" *summary.txt >  "$working_directory/FastQC_summary/$outfile.txt"
done

## When we have completed the analysis of all reads, we should go into this report for each tally and 
## create a composite output result for each test if the result was not pass 
## Pull out any column2 where column1 is 'WARN' or 'FAIL' and store FileName ($3) and Test ($2) and Result ($1) in full_summary.txt
## print test_summay to screen 

for file in $(ls "$working_directory/FastQC_summary"); 
do 
while IFS=$'\t' read result testtype filename; 
do  
if [ $result != "PASS" ]; then
echo "$filename $testtype $result" >> "$working_directory/FastQC_summary/full_summary.txt";
fi; 
done < "$working_directory/FastQC_summary/$file"; 
done


================================ GOT TO HERE WORKING ============================

# OPTIONAL: Display QC scores below cuttoffs to screen 
# OPTIONAL: Display statistics on QC scores
# OPTIONAL: Look into using MultiQC after FastQC but this may be cheating as not specified in assignment



for file

# QUESTION: How can you look at the html or png files? 
# QUESTION: Do we read in files paired, in groups or one by one 

### Bowtie2 pipeline
## Pipeline defaults

genome_dir = 



