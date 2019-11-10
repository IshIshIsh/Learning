git = https://github.com/IshIshIsh/Learning.git

#1) copy folder Assignment1 to home space (resets permissions):  
#2) If needing to retest, delete contents of ./BPSM_assignment1 with 
cp /localdisk/data/BPSM/Assignment1/* ./BPSM_assignment1/ -r
rm -f -r Assignment1


### A pipeline for the analysis of paired-end RNA-Seq  reads using fastqc, bowtie2, samtools and bedtools 

## Assumptions: 
# raw sequence files are in fq.gz format and stored in one folder
# a mapping file is located within this folder 
# a fasta format genome sequence is provided in the specified genome_dir  


## Global defaults 

working_directory=BPSM_assignment1

# threading is used at an approproiate level for the active computer - please check the defaults before running. 
# TO GET CPU INFO: cat /proc/cpuinfo
## Select a thread count based on the output (in a power of two to make use of alignmment)
threads=16 

bedfile="$working_directory/Tbbgenes.bed"
fastqc_data_dir="$working_directory/fastq"

### -------------- Fastqc Pipeline  --------------- ### 

# As fastqc doesn't create directories we want to create a folder to store the output in if it doesn't already exist:
# check the default folder name (fastqc_output) doesn't exist in working directory and if not, create it.

if [ ! -d "$working_directory/fastqc_output" ]; then
    mkdir "$working_directory/fastqc_output"
fi 
fastqc_out_dir="$working_directory/fastqc_output"

# for each file in the fastqc_data_dir
# if the file is a fasta file (excluding the map file which has no file extension)
filesforfastqc=''
for file in $(ls "$fastqc_data_dir/"*.fq.gz)
do 
filesforfastqc="$filesforfastqc $file"
# run fastqc on the data and output the results to the fastqc_out_dir
done
fastqc $filesforfastqc -o $fastqc_out_dir -t $threads

### fastqc outputs the data to fastqc_output/filename and creates many files in a zip file and a html file. 
## On extraction of the zip files, mutiple folders and files exist: 
##  fastqc_data.txt  fastqc.fo  fastqc_report.html  Icons  Images  summary.txt
## In summary.txt and pass/warn/fail report is generated with columns (tab seperated) of 'Pass  Test    File'

for outfile in $(find "$fastqc_out_dir/"*.zip -type f -printf "%f\n")
do
#unzip -p $fastqc_out_dir/216_L8_1_fastqc.zip *summary.txt >  "$working_directory/fastqc_output/216_L8_1_fastqc.txt"
unzip -p "$fastqc_out_dir/$outfile" *summary.txt >  "$working_directory/fastqc_output/$outfile.txt"
done

## When we have completed the analysis of all reads, we should go into this report for each tally and 
## create a composite output result for each test if the result was not pass 
## Pull out any column2 where column1 is 'WARN' or 'FAIL' and store FileName ($3) and Test ($2) and Result ($1) in full_summary.txt
## print test_summay to screen 

pass=0
warn=0
fail=0 

rm -f "$working_directory/fastqc_output/full_summary.txt"
rm -f "$working_directory/fastqc_output/testing_summary.tsv" 

for file in $(find "$working_directory/fastqc_output/"*.txt -type f -printf "%f\n"); 
do 
while IFS=$'\t' read result testtype filename; 
do  
if [ "$result" == "PASS" ]; then
pass=$((pass+1)) 
fi; 
if [ "$result" == "WARN" ]; then
warn=$((warn+1))
fi; 
if [ "$result" == "FAIL" ]; then
fail=$((fail+1))
fi; 
echo -e "$filename\t$testtype\t$result" >> "$working_directory/fastqc_output/full_summary.txt"
done < "$working_directory/fastqc_output/$file"
done
total=$((pass+warn+fail))
echo "Pass: $pass/$total"
if [ $warn != 0 ] || [ $fail != 0 ]; then 
echo "WARNING: failed tests: $fail/$total, warnings: $warn/$total"
# print the number of files which failed or passed each test and write the filenames to a file in the format 'test', 'error', 'filelist'
# this would allow (if a very important test failed, to pull out an array of those filenames to filter your data for further analysis)
awk '{FS="\t"; if($3 == "FAIL") failedtest[$2]++} END {for (k in failedtest) print k, "\t", failedtest[k] }' "$working_directory/fastqc_output/full_summary.txt" 
awk '{FS="\t"; if($3 == "WARN") warntest[$2]++} END {for (k in warntest) print k, "\t", warntest[k] }'  "$working_directory/fastqc_output/full_summary.txt" 
awk '{FS="\t"; if($3 == "FAIL") failedtest[$2]++} {filenames=$1 "," filenames} END {for (k in failedtest) print k, "\t", "FAIL", "\t", filenames }' "$working_directory/fastqc_output/full_summary.txt"  >> "$working_directory/fastqc_output/testing_summary.tsv" 
awk '{FS="\t"; if($3 == "WARN") warntest[$2]++} {filenames=$1 "," filenames} END {for (k in warntest) print k, "\t", "WARN" , "\t", filenames}'  "$working_directory/fastqc_output/full_summary.txt"  >> "$working_directory/fastqc_output/testing_summary.tsv"
fi; 

# OPTIONAL: Look into using MultiQC after FastQC but this may be cheating as not specified in assignment
# QUESTION: How can you look at the html or png files? 
# QUESTION: Do we read in files paired, in groups or one by one 

### -------------- Bowtie Pipeline --------------- ### 

genome="$working_directory/Tbb_genome"
zipfilename=Tb927_genome.fasta.gz
fasta=Tb927_genome.fasta
genomename=Tb927_genome


### Creating 'GenomeSequence.bt2'
if [ ! -d "$working_directory/bowtie" ]; then
    mkdir "$working_directory/bowtie"
fi 
genomeseq="$working_directory/bowtie"


## Unzip the files in the reference genome
# OPTIONAL: Fina a way to pipe into -build so that we don't take up extra space? 
gunzip "$genome/$zipfilename" 

## Create a Bowtie2 index genome from the fasta files from the reference genome
bowtie2-build "$genome/$fasta" "$working_directory/bowtie/$genomename"
================================ GOT TO HERE WORKING ============================

# -x is the index file, -p is the number of threads (preset globally), --m is memory mapping index 
pair1s=$(find "$fastqc_data_dir/"*1.fq.gz -type f -printf "$fastqc_data_dir/%f,")
pair2s=$(find "$fastqc_data_dir/"*2.fq.gz -type f -printf "$fastqc_data_dir/%f,")
bowtie2 -x "$working_directory/bowtie/$genomename" -p $threads --mm -1 $pair1s -2 $pair2s -S "$working_directory/bowtie/$genomename"

#### this works before folder realignment 
# bioinfmsc5:~/BPSM_assignment1/Assignment1/genomeseq$ bowtie2 -x $genomename -p $threads --mm -1 ../fastq/216_L8_1.fq.gz -2 ../fastq/216_L8_2.fq.gz -S OutputFile
####


## Pipeline defaults

