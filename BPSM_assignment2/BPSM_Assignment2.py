"""
Assignment 2: Identifying Protein Sequences. 

Thoughts & Notes: 
    -Allow obtaining of DNA or just Protein Levels? If we have time at the end, 
    using a fasta of DNA sequences & converting to protein using a PAM matrix (PAM250?) OR BLOSUM matrix
    for alignment could be a useful add on
        - To select PAM or BLOSUM, do some form of logic on how related the groups of sequences are? 


## Recomended test set: genefamily = glucose-6-phosphatase taxid = Aves
Aves[Organism:exp] AND (glucose-6-phosphatase[Gene Name])
ie: 
txid8782[Organism:exp] AND G6pc[Gene Name]
"""


 ### TO DO LIST ### 
"""
 - ADD: flexibility to find txid and gene name from inputs....
 - ADD: method of limiting sequence number BEFORE retrival (default 10,000)
 - ADD: method for checking species number for dataset with user prompt to continue
 - ADD: Default binning method (100) max (250) and user prompt select
 - ADD: Determine level of conservation between sequences
 - ADD: Plot level of conservation between sequences for selected bin(s)
 - ADD: Method to run as verbose (no print statements) so we can print to file
 - ADD: Method to run as silent (no print or user inputs) so we can run as script and leave 
 - ADD: Choice for conservation selection???
 - ADD: Scan prosite databse for known motifs in protein subset (bin). Note: if multiple bins selected, run in background and limit processing power
 - ADD: Add motif record containing motif, sequence, position (+ link??)
 - ADD: Kwargs dict option for subprocess, allowing var: varentry (key, value pair) for more flexibility in bash commands not programmed in the main python script (Note: Bash doesnt have dicts so this option would only work if running script from python). Alternative method would be to add a string onto end (but that's less fun)
 - CHECK: Other analysis (emboss) 
 - CHECK: timeout for subprocess


EMBOSS Packages which could be useful: 
- Pepstats
- pepwindow
- pepinfo


"""

### ---------- Import Statements ---------- ###
import os
from pathlib import Path 
import edirect 
import re
import subprocess
from datetime import datetime
import re

### ---------- TO RUN AS SCRIPT ---------- ### 
def __main__():
    run_protein_ident(genefamily, taxid, filtered_partial = True, filtered_predicted = True, filepath= None, foldername=None, foldertime = False, iteration=30,
    overwite_max_seqno =False)

### ---------- Master Function ---------- ###

def run_protein_ident(genefamily, taxid, filtered_partial = True, filtered_predicted = True, filepath= None, foldername=None, foldertime = False, iteration=30,
    overwite_max_seqno =False):
    """    
    Parent Function wrapper for processing Protein Alignments & Statistics 
    Requires: 


    1) Before starting any analysis we check the inputs are of the correct type. If they are not we break or return to default.
    2) Check query is not > hundreds of lines long & replace any spaces with + (will break ncbi request)
    3) Create a folder for the output
    4) Retrives Fasta Sequences from NCBI
    5) Generates a dictionary with protein stats for each fasta (and a pepstats file with extra info)
    6) Creates an alignment of the fasta sequences using Clustalo
    7) Creates a consensus sequence from the alignment file using cons package
    8) 

    
    """
    check_edirect_installation()
    check_type(genefamily, str, 'genefamily', default = None)
    check_type(taxid, str, 'taxid', default = None)
    genefamily = replace_non_alphanumeric_chars(genefamily)
    pathtaxid = replace_non_alphanumeric_chars(taxid)
    query_filename = genefamily+'_'+pathtaxid+'_'
    filepath = check_type(filepath, str, 'filepath', str(Path().absolute()))
    if foldername != None: 
        foldername =  check_type(foldername, str,'foldername', query_filename)
        analysispath = create_folder_path(foldername, filepath, foldertime)
    else: 
        foldername = check_type(query_filename, str, 'query_filename')
        analysispath = create_folder_path(foldername, filepath, foldertime)
    fastafiles = get_from_ncbi(genefamily, taxid, query_filename, analysispath)
    split_fasta = check_read_no_by_split(fastafiles, '>')
    fastano = len(split_fasta) + 1
    if fastano > 10000:
        print('WARNING: Sequence number is greater than 10000. For further analysis this may be cut off at 10000 reads unless specified by the user')
    protein_dict = check_protein_basic_stats(split_fasta)
    protein_stats = get_protein_stats(fastafiles)
    clustalo_files = align_with_clustalo(fastafiles, iteration)
    consensus = get_consensus_from_alignment(clustalo_files)
    blastdb = create_blastdb(fastafiles)
    query = query_choice(analysispath) 
    blastp_output = query_blastdb(blastdb, query, analysispath)

### ---------- Utility Functions ---------- ###


def check_type(inputquery, expected, varname, default = None, firm = True):
	"""
	Checks the input against an expected type and either throws an error or returns default value if passed in 
	"""
	if isinstance(inputquery, expected) == True:
		return inputquery
	else:
		realtype = type(inputquery)
		if default == None:
			raise ValueError('The input: '+str(inputquery)+' for var:'+varname+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
		else: 
			print('WARNING: The input: '+str(inputquery)+'  for var:'+varname+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
			print('WARNING: Attempting to return to default value if applicable:'+str(default))
			return default 


def replace_non_alphanumeric_chars(input_string, replacement_char = '_'): 
    """
    Utility function used to replace any non alphanumeric charecter 
    """
    output = re.sub('[^0-9a-zA-Z]+', replacement_char, input_string)
    return output 


def run_nonpython_process(query, timeout = 6000):
    """
    Function wich uses subprocess rather than os due to increased security.
    Note TO SELF: Timeout doesn't seem to actually cause a timeout error???
    """
    try: 
        print(query)
        response = subprocess.check_output(query, shell = True)
        #process = subprocess.Popen(query, stdout = subprocess.PIPE , stderr = subprocess.PIPE, shell = True)
        #response, error = process.communicate(timeout = timeout)
        response = response.decode('ascii')
        #error = error.decode('ascii')
        #print(error)
        return response 
    except Exception as e: 
        raise RuntimeError('Error running edirect: '+str(e))


def check_edirect_installation():
    """
    Check that edirect is installed or install
    User Prompt for update included. 
    """
    print ('Checking edirect is installed in home space...')
    if os.path.isdir(os.path.expanduser('~/edirect')):
        updateinput = input('Edirect is already installed, to update please input \'update\' or press any key to continue with analysis')
        if updateinput.strip().lower() == 'update': 
            print('Attempting to update Edirect now')
            run_nonpython_process('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"')
            print('Edirect has been sucessfully updated. Continuing with analysis')
    else: 
        print('Edirect is not currently installed, trying to install now')
        run_nonpython_process('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"')
        print('Edirect has been sucessfully installed. Continuing with analysis')


def threads_from_cpu(process):
    """
    Finds the cpu number and asks user prompt for thread number. 
    Defaults to a third of available threads if a valid integer is not entered.
    Valid integers are between 1 and cpu number 
    """
    cpu = os.cpu_count()
    threads = input('Your cpu count is: '+str(cpu)+' please enter a integer of the threads you wish to use for the process:'+process)
    default_threads = int(cpu/3) 
    threads = check_type(threads, int, 'threads', default_threads)
    if threads < 1 or threads > cpu:
        raise ValueError('Thread number cannot be less than 1 or greater than the max cpu count: '+str(cpu_count))
    else: 
        return threads 


def return_file_contents(file):
    with open(file) as f:
        data = f.read()
    return data 


def check_read_no_by_split(input_file, splitchar): 
    data = return_file_contents(input_file)
    datasplit = data.split(splitchar)
    if isinstance(datasplit, list) == False: 
        raise ValueError('input file could not be split by specified charecter:'+str(splitchar))
    datasplit = list(filter(None, datasplit))
    return(datasplit)


### ---------- Python File Analysis ---------- ###

def check_protein_basic_stats(fasta_list):
    if isinstance(fasta_list, list) == False:
        raise ValueError('Fasta_list input was not of type list. No dictionary created from fasta files.')
    short_name_regex = re.compile('^[\S]*')
    protein_dict = {}
    for f in fasta_list:
        name = f[0:f.find('\n')]
        label = short_name_regex.search(name)
        seq = f[f.find('\n'):]
        seqlength = len(seq.replace('\n', ''))
        protein_dict[label.group(0)] = {'full_name': name, 'sequence': seq, 'length': seqlength}
    return protein_dict 
    

#def add_pepstats_info_to_dict(protein_dict, pepstatsfile):
#    data = return_file_contents(pepstatsfile)
#    for k in protein_dict:
#        proteindict[k]['Mol_Weight'] = 
#        proteindict[k]['Residues'] =
#        proteindict[k]['Av_residue_weight'] =
#        proteindict[k]['Charge'] =
#        proteindict[k]['Isoeelectric_point'] =

            


def check_alignment_basic_stats():
    pass


### ---------- File outputing & File locations ---------- ### 


def create_folder_path(name, filepath, foldertime = False):
    """
    Checks if the inputted folder (name  + prefix + filepath) exists 
    if not it will create it. Both cases will result in a print statement 
    and return a string of the folder path to store results in.

    name: If no user defined name, it will default to family_taxid
    filepath: string, if no user defined input it will default to current working directory 

    Note: assignment1 feedback - 'what if we do the same analysis tomorrow?' - 
        assuming if we run analysis multiple times in one day, we'd like to overwrite so default time = False.

    """
    prefix = 'Protein_Analysis_'
    if foldertime == False: 
        datesuffix = datetime.today().strftime('%Y_%m_%d')
    else: 
        datesuffix = datetime.today().strftime('%Y_%m_%d_time_%H_%M')
    folder = filepath+'/'+prefix+name+datesuffix
    if not os.path.exists(folder):
        os.mkdir(folder)
        print(folder+' was created and will be used for storing the results of the analysis')
    else:
        print(folder+' already exists and will be used for storing the results of the analysis')
    return folder 



### ---------- NCBI Data Collection ---------- ### 
"""
Create BLASTdb (makeblastdb) or select existing - better to access remotely 
take query and run blastn, blastx or tblastn 

- We use subprocess module to use Unix edirect package inside the python terminal (safer than os?) with a timeout and error catching on communicate 
- We have defaults to filter both/either/none partial or predicted sequences.  
"""
def get_from_ncbi(family, taxid, query_filename, filepath, filtered_partial = True, filtered_predicted = True):
    print('Attempting to retrive data from NCBI for '+family+' and '+taxid+'...')
    if filtered_partial == True and filtered_predicted == True: 
        print('Note: (Defaults) Partial or Predicted sequences will be excluded')
    elif filtered_partial == False and filtered_predicted == True: 
        print('Note: (Default) Predicted sequences will be excluded')
    elif filtered_partial == True and filtered_predicted == False: 
        print('Note: (Default) Partial sequences will be excluded')
    else: 
        print('Note: Partial and Predicted sequences will be included in output')
    edirect_query_task = taxid+"[Organism:exp] AND "+family+"[Gene Name]"
    if filtered_partial == True:
        query_search_term = edirect_query_task+' NOT partial[All Fields]'
    if filtered_predicted == True:
        query_search_term = query_search_term+' NOT predicted[All Fields]'
    print('Query Searched is:'+str(query_search_term))
    save_location = filepath+'/'+query_filename+'fasta.fasta'
    edirect_query = "~/edirect/esearch -db protein -query \""+query_search_term+"\" | ~/edirect/efetch -format fasta" 
    edirect_reponse = run_nonpython_process(edirect_query)
    print(edirect_reponse)
    with open(save_location, 'w') as result_file:
        result_file.write(edirect_reponse)
    print('Data Retrival sucessful, fasta files stored in '+str(save_location))
    return save_location

### --------- Clustalo Alignment ---------  ####

def align_with_clustalo(fastafile_location, iteration, overwite_max_seqno =False):
    """
    Uses the clustalo package through bash (through subprocess) to generate a multiple sequence alignment file
    fastafile_location: The full path of the fasta file location 
    iteration: the number of iterations for the clustalo tree
    overwrite_max_seqno: default value False, if true removes limit on the number of sequences WARNING: DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!!!
    """
    print('Attempting Clustalso Alignment...')
    clustal_output = fastafile_location.replace('fasta.fasta', 'clustal.msf')
    threads = threads_from_cpu('clustalo alignment')
    noreads = check_read_no_by_split(fastafile_location, '>')
    if overwite_max_seqno == False: 
        maxnumseq = 10000
        if len(noreads) > maxnumseq:
            print('WARNING: Number of reads inputted is greater than the maxnumseq allowed. If absouletly required use overwrite_max_seqno = True and re-analyse your data')
            print('WARNING: Changing overwrite_max_seqno to True is a dangerous option. DO NOT DO THIS UNLESS ABSOLUTELY REQUIRED')
        edirect_reponse = run_nonpython_process("clustalo -i \""+fastafile_location+"\" -t Protein --outfmt=msf --maxnumseq "+str(maxnumseq)+" --full --threads "+str(threads)) 
    else: 
        edirect_reponse = run_nonpython_process("clustalo -i \""+fastafile_location+"\" -t Protein --outfmt=msf --full --threads "+str(threads)) 
    print(edirect_reponse)
    with open(clustal_output, 'w') as result_file:
         result_file.write(edirect_reponse)
    print('Alignment sucessful, clustal files stored in '+str(clustal_output))
    return clustal_output




def get_consensus_from_alignment(clustal_output):
    """
    Uses EMBOSS Cons to get a consensus sequence from the Clustalo alignment files. 
    """
    output = clustal_output.replace('clustal.msf', 'consensus.fa' )
    edirect_reponse = run_nonpython_process("cons -sequence \""+clustal_output+"\" -outseq \""+output+"\" -sprotein1 true -sid1 \"consensus\" -osformat2 \"fasta\"")
    print('Consensus creation sucessful. Consensus file stored in: '+ str(output))
    return output 

### ---------- BLAST ---------- #### 


def create_blastdb(fastafile_location):
    """
    Uses makeblastdb to make a blast-searchable database from fasta files and save them to an output for use in blastp processing. 
    """
    print('Attempting to create a blastdb from fasta files')
    output = fastafile_location.replace('fasta.fasta', 'blastdb')
    edirect_reponse = run_nonpython_process("makeblastdb -in \""+fastafile_location+"\" -dbtype prot -out \""+output+"\"")
    print('blast database creation sucessful. Blastdb file stored in: '+ str(output))
    return output     


def query_choice(filepath):
    fastas = [fastafile for fastafile in os.listdir(os.path.expanduser(filepath)) if fastafile.endswith(".fasta") == True or fastafile.endswith(".fa") == True ]
    print('You have the following files with valid file extension for querying against blasdb:')
    count = 0
    for f in fastas:
        print(str(count) +':'+ f)
        count = count + 1
    selectedfile = input('Please select the file you wish to use by entering the name or index. If your file is not present, enter your filename instead \\n Please note: Files must be present in the directory:'+str(filepath))
    try:
        selectedfileindex = int(selectedfile)
        if selectedfileindex < count: 
            selectedfile = fastas[selectedfileindex]
        else: 
            raise ValueError('Selected Index:'+selectedfile+'is not valid it must be in the range: 0-'+str(len(fastas-1)))
    except ValueError: 
        pass 
    if selectedfile not in fastas:
        if os.path.exists(selectedfile):
            pass 
        else: 
            raise ValueError('Selected file:'+str(selectedfile)+' either does not exist or is not a valid option')
    print('You have selected file: '+str(selectedfile)+'to query against blastdb')
    reads = check_read_no_by_split(filepath+'/'+selectedfile, '>')
    if len(reads) > 1:
        print(len(reads))
        print('WARNING: The file you have selected contains multiple sequences... creating a dictionary of sequence choices')
        fasta_dict = check_protein_basic_stats(reads)
        selected_record = input('Please enter the accession number of the fasta sequence you wish to continue with...\\n'+str(fasta_dict.keys()))
        if selected_record == '':
            raise ValueError('Entry was None, invalid response. Terminating program')
        if selected_record in fasta_dict.keys():
            selectedfile = selected_record+'_query_record.fasta'
            with open(filepath+'/'+selectedfile, 'w') as f:
                f.write(fasta_dict[selected_record]['full_name'] + '\n' + fasta_dict[selected_record]['sequence'])
        else: 
            raise ValueError('Selected accession number is not in the list of accession numbers from file: '+selectedfile)
    print('File selection sucessful:'+str(selectedfile))
    return filepath+'/'+selectedfile


def query_blastdb(blastdb, query_file, filepath):
    """
    Uses blastp query a blast-searchable database with selected query_file. 
    """
    print('Attempting to run blastp for selected file:'+ str(query_file)+' against blastdb:'+str(blastdb))
    output = query_file.replace('.fasta','_blastp_output.txt')
    edirect_reponse = run_nonpython_process("blastp -db \""+blastdb+"\" -query \""+query_file+"\" -outfmt 7 > \""+output+"\"")
    print(edirect_reponse)
    print('blastp sucessful. Output file stored in: '+ str(output))
    return output 


### ---------- Plotting Routines ---------- ### 


### ---------- Protein Stats ---------- ### 
def get_protein_stats(fastafile_location):
    """
    Uses EMBOSS pepstats to get protein statistics from the fasta files and save them to an output. 
    """
    print('Attempting to get protein statistics from pepstats')
    output = fastafile_location.replace('fasta.fasta', 'stats.pepstats')
    edirect_reponse = run_nonpython_process("pepstats -sequence \""+fastafile_location+"\" -outfile \""+output+"\"")
    print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
    return output 


### ---------- Prosite Data Collection ---------- ### 




run_protein_ident('G6pc', 'txid8782')


"/localdisk/home/ifarquha/Protein_Analysis_G6pc_txid8782_2019_11_10/G6pc_txid8782_consensus.fa"
