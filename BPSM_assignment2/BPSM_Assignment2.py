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
- pepinfo
- pepscan 


"""

### ---------- Import Statements ---------- ###

import os
from pathlib import Path 
import re 
import subprocess
from datetime import datetime
import math
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

### ---------- TO RUN AS SCRIPT ---------- ### 

def __main__():
    run_protein_ident(genefamily, taxid, filtered_partial = True, filtered_predicted = True, filepath= None, foldername=None, foldertime = False, iteration=30,
    overwite_max_seqno =False)


### ---------- Master Function ---------- ###

def run_protein_ident(genefamily, taxid, filtered_partial = True, filtered_predicted = True, filepath= None, foldername=None, foldertime = False, iteration=30,
    overwite_max_seqno =False, bin_no = 100):
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
    print('Number of fasta sequences is: '+str(fastano))
    if fastano > 10000:
        print('WARNING: Sequence number is greater than 10000. For further analysis this may be cut off at 10000 reads unless specified by the user')
    fasta_dict = check_protein_basic_stats(split_fasta)
    protein_stats = get_protein_stats_from_pepstats(fastafiles)
    fasta_dict = add_pepstats_info_to_dict(fasta_dict, protein_stats)
    clustalo_files = align_with_clustalo(fastafiles, iteration)
    create_hydropathy_plot(clustalo_files)
    create_plotcon(clustalo_files)
    print('Printing clustalo_files variable')
    print(clustalo_files)
    print("\n")
    consensus = get_consensus_from_alignment(clustalo_files)
    print('Consensus Sequence from multiple alignment:')
    print(consensus)
    blastdb = create_blastdb(fastafiles)
    query = query_choice(analysispath) 
    blastp_output = query_blastdb(blastdb, query, analysispath)
    bins = bin_results(blastp_output, bin_no)
    bin_labels = bins.keys()
    selected_bin_labels = bin_selection(bin_labels, 'Prosite Motif search')
    os.mkdir(analysispath+'/motifs')
    if isinstance(selected_bin_labels, int):
        for selected_bin_label in selected_bin_labels:
            os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
            print('Creating individual fasta records for:'+str(selected_bin_label))
            fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
            for fasta in fasta_list_from_bin:
                create_new_fasta_record(fasta, fasta_dict[fasta]['sequence'], analysispath+'/motifs/bin_'+str(selected_bin_label))
            search_for_motifs(fasta_list_from_bin, analysispath+'/motifs/bin_'+str(selected_bin_label))
            print('Prosite scan completed for bin:'+selected_bin_label)
    else: 
        selected_bin_label = selected_bin_labels
        os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
        print('Creating individual fasta records for:'+str(selected_bin_label))
        fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
        for fasta in fasta_list_from_bin:
            create_new_fasta_record(fasta, fasta_dict[fasta]['sequence'], analysispath+'/motifs/bin_'+str(selected_bin_label))
        search_for_motifs(fasta_list_from_bin, analysispath+'/motifs/bin_'+str(selected_bin_label))
        print('Prosite scan completed for bin:'+selected_bin_label)
    for folder in os.listdir(os.path.expanduser(analysispath+'/motifs/')):
        report_motif_results(folder)
    print('For more information on motifs per sequence, please see relevant files in folder:'+str(analysispath+'/motifs/'))


### ---------- Utility Python Functions ---------- ###


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
    try: 
        default_threads = int(threads)
    except: 
        'WARNING: Threads was not a valid entry: '+threads+'returning to default value of:'+str(default_threads)
    if default_threads < 1 or default_threads > cpu:
        raise ValueError('Thread number cannot be less than 1 or greater than the max cpu count: '+str(cpu))
    else: 
        return default_threads 


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


def bin_results(input_file, bin_no, ignore_lines = '#'):
    result_bin = {}
    data = return_file_contents(input_file)
    lines = [line for line in data.split('\n')]
    validlines = [line for line in lines if line.startswith(ignore_lines) == False]
    validlines = list(filter(None, validlines))
    ## As lists have defined order we can parse this using the first as the top hits
    totalno = len(validlines)
    no_bins = math.ceil(totalno/bin_no)
    rank_counter = 0 
    for i in range(0, no_bins-1):
        result_bin[i] = []
        for j in range(0,bin_no-1):
            if j < totalno:
                result_bin[i].append(convert_blastline_to_dict(validlines[j], rank_counter))
                rank_counter = rank_counter + 1
            else: 
                pass
    bins = result_bin.keys()
    bins_created = len(bins)
    print(bins_created+ ' were created:' +str(bins))
    return result_bin


def convert_blastline_to_dict(line, rank_counter):
    entry = line.split('\t')
    if len(entry) != 12:
        raise ValueError('Blast output expected to be length 12 but is actually length:'+str(len(entry)))
    else:
        entry_dict = {
            "rank": rank_counter,
            "query_acc": entry[0],
            "subject_acc": entry[1], 
            "percent_identity": entry[2], 
            "alignment_length": entry[3],
            "mismatches": entry[4],
            "gap_opens": entry[5],
            "query_start": entry[6], 
            "query_end": entry[7],
            "subject_start": entry[8],
            "subject_end": entry[9],
            "evalue": entry[10],
            "bit_score": entry[11]
            }
    return entry_dict 


def bin_selection(valid_bin_list, method):
    default_bin = False
    print('Bins for analysis: '+str(method)+' are '+str(valid_bin_list))
    input_bins = input('Please select which bin number to continue with, or enter a list seperated by \',\'. Default value is the top scoring group (0) ')
    if len(input_bins.split(',')) > 1: 
        input_bins = list(filter(None, input_bins.split(',')))
        selected_bin = []
        for input_bin in input_bins:
            if input_bin.isdigit():
                if int(input_bin) not in valid_bin_list:       
                    print ('WARNING: Selected bin: '+str(input_bin)+' from list of bins:'+str(input_bins)+' does not exist in list of valid bins:'+str(valid_bin_list)+'. Returning to default value: Group 0')
                    default_bin = True
                else: 
                    selected_bin.append(int(input_bin))
            else:
                print('WARNING: Non integer bin:'+str(input_bin)+' found in list of input bins:'+str(input_bins)+' Returning to default value of: Bin 0')
                default_bin = True    
    else:
        if input_bins.isdigit():
            if int(input_bins) not in valid_bin_list:       
                print ('WARNING: Selected bin: '+str(input_bins)+' does not exist in list of valid bins:'+str(valid_bin_list)+'. Returning to default value: Group 0')
                default_bin = True
            else:
                selected_bin = input_bins
        else: 
            print('WARNING: Non integer bin:'+str(input_bins)+' selected. Returning to default value of: Bin 0')
            default_bin = True    
    if default_bin == True:
        selected_bin = 0        
    return selected_bin


def create_pd_from_blastp(inputfile, ignored = '#'):
    headers=['query_acc.', 'subject_acc.', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
    blastdb = pd.read_csv(inputfile, comment = ignored, sep ='\t', names = headers)
    return blastdb


short_name_regex = re.compile('^[\S]*')
mol_weight = re.compile('Molecular weight = (\S*)')
isoelectric = re.compile('Isoelectric Point = (\S*)')
residues = re.compile('Residues = (\S*)')
charge = re.compile('Charge   = (\S*)')
Av_res_weight = re.compile('Average Residue Weight  = (\S*)')


def check_protein_basic_stats(fasta_list):
    if isinstance(fasta_list, list) == False:
        raise ValueError('Fasta_list input was not of type list. No dictionary created from fasta files.')
    fasta_dict = {}
    for f in fasta_list:
        name = f[0:f.find('\n')]
        label = short_name_regex.search(name)
        seq = f[f.find('\n'):]
        seqlength = len(seq.replace('\n', ''))
        fasta_dict[label.group(0)] = {'full_name': name, 'sequence': seq, 'length': seqlength}
    return fasta_dict 
    

def add_pepstats_info_to_dict(fasta_dict, pepstatsfile):
    data = return_file_contents(pepstatsfile)
    entries = list(filter(None, data.split('PEPSTATS of ')))
    for entry in entries: 
        short_name = short_name_regex.search(entry).group(1)
        if short_name in fasta_dict.keys():
            fasta_dict[short_name]['Mol_Weight'] = mol_weight.search(entry).group(1)
            fasta_dict[short_name]['Residues'] =residues.search(entry).group(1)
            fasta_dict[short_name]['Av_residue_weight'] =short_name_regex.search(entry).group(1)
            fasta_dict[short_name]['Charge'] =charge.search(entry).group(1)
            fasta_dict[short_name]['Isoeelectric_point'] =isoelectric.search(entry).group(1)
        else: 
            raise ValueError('Short Label:'+str(short_name)+' was found in pepstats output but not in fasta_dict')
    return fasta_dict


def check_alignment_basic_stats(alignment_file, protein_dict):
    data = check_read_no_by_split(alignment_file, '>')
    for d in data: 
        name = d[0:d.find('\n')]
        label = short_name_regex.search(name)
        seq = d[d.find('\n'):]
        if label.group() in protein_dict.keys():
            protein_dict[label.group()]['alignment_seq'] = seq
        else: 
            raise ValueError('Protein dict does not contain key:'+str(label.group())+' something has gone wrong processing fastas or alignment file')
    return protein_dict


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
    clustal_output = fastafile_location.replace('fasta.fasta', 'clustal.fasta')
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
    run_nonpython_process("cons -sequence \""+clustal_output+"\" -outseq \""+output+"\" -sprotein1 true -sid1 \"consensus\" -osformat2 \"fasta\"")
    print('Consensus creation sucessful. Consensus file stored in: '+ str(output))
    return output 


### ---------- File Conversion ---------- #### 

def convert_filetypes_seqret(msf_infile, fasta_outfile):
    """
    Uses makeblastdb to make a blast-searchable database from fasta files and save them to an output for use in blastp processing. 
    """
    print('Attempting to convert msf alignment to fasta')
    run_nonpython_process("seqret -sequence \""+msf_infile+"\" -outseq  \""+fasta_outfile+"\"")
    print('Consensus msf file converted to fasta sucessfully. New file stored in: '+ str(fasta_outfile))


### ---------- BLAST ---------- #### 

def create_blastdb(fastafile_location):
    """
    Uses makeblastdb to make a blast-searchable database from fasta files and save them to an output for use in blastp processing. 
    """
    print('Attempting to create a blastdb from fasta files')
    output = fastafile_location.replace('fasta.fasta', 'blastdb')
    run_nonpython_process("makeblastdb -in \""+fastafile_location+"\" -dbtype prot -out \""+output+"\"")
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
    print('You have selected file: '+str(selectedfile)+' to query against blastdb')
    reads = check_read_no_by_split(filepath+'/'+selectedfile, '>')
    if len(reads) > 1:
        print(len(reads))
        print('WARNING: The file you have selected contains multiple sequences... creating a dictionary of sequence choices')
        fasta_dict = check_protein_basic_stats(reads)
        selected_record = input('Please enter the accession number of the fasta sequence you wish to continue with...\\n'+str(fasta_dict.keys())+"\\n")
        if selected_record == '':
            raise ValueError('Entry was None, invalid response. Terminating program')
        if selected_record in fasta_dict.keys():
            selectedfile = selected_record+'_query_record.fasta'
            with open(filepath+'/'+selectedfile, 'w') as f:
                f.write('>'+fasta_dict[selected_record]['full_name'] + fasta_dict[selected_record]['sequence'])
        else: 
            raise ValueError('Selected accession number is not in the list of accession numbers from file: '+selectedfile)
    print('File selection sucessful:'+str(selectedfile))
    return filepath+'/'+selectedfile


def query_blastdb(blastdb, query_file, filepath):
    """
    Uses blastp query a blast-searchable database with selected query_file. 

    Note: FAILURE - MOST HAVE EVALULE 0.0. This is means there is a high probability that these sequences are related (as a measure of how likely they'd align by chance)
        however as we're asked to plot the similarity level, we'd like to have a measure we can plot. The Bit score is plotable, but this is not the most reliable method 
        of relatedness (and is used to calculate evalue. Continuing with assignment regardless & will come back if time at the end).
    """
    print('Attempting to run blastp for selected file:'+ str(query_file)+' against blastdb:'+str(blastdb)+'with sequence:')
    print(return_file_contents(query_file))
    outputfile = query_file.replace('.fasta', '.fa').replace('.txt', '.fa')
    output = outputfile.replace('.fa','_blastp_output.txt')
    threads = str(threads_from_cpu('blastp query'))
    edirect_reponse = run_nonpython_process("blastp -db \""+blastdb+"\" -query \""+query_file+"\" -num_threads \""+threads+"\" -outfmt 7 > \""+output+"\"")
    print(edirect_reponse)
    print('blastp sucessful. Output file stored in: '+ str(output))
    return output 


### ---------- Protein Stats ---------- ### 

def get_protein_stats_from_pepstats(fastafile_location):
    """
    Uses EMBOSS pepstats to get protein statistics from the fasta files and save them to an output. 
    """
    print('Attempting to get protein statistics from pepstats')
    output = fastafile_location.replace('fasta.fasta', 'stats.pepstats')
    run_nonpython_process("pepstats -sequence \""+fastafile_location+"\" -outfile \""+output+"\"")
    print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
    return output 


### ---------- Prosite Data Collection ---------- ### 

def get_prosite_db_scratch(ftp_address = "http://ftp.ebi.ac.uk/pub/databases/prosite/" ): 
    """
    Downloads the prosite.dat and prosite.doc files using wget into home directory
    Used for prosextract db compile (and then for prosite scan)
    """
    directory = os.path.expanduser('~')
    print('Checking if prosite db is present in home directory...')
    if os.path.exists(directory+"/prosite/prosite.doc") and  os.path.exists(directory+"/prosite/prosite.dat"):
        print('Prosite directory is already populated. Continuing with analysis')
    else: 
        print('Prosite db was not found: Attempting to download prosite database (.doc & .dat files) to home directory...')
        if not os.path.exists(directory+"/prosite"):
            os.mkdir(directory+"/prosite")
        try: 
            run_nonpython_process("wget "+ftp_address+"/prosite.doc -O "+directory+"/prosite/prosite.doc")
            print('prosite.doc has been downloaded')
            run_nonpython_process("wget "+ftp_address+"/prosite.dat -O "+directory+"/prosite/prosite.dat")
            print('prosite.dat has been downloaded')
        except: 
            print('ftp request failed. Retrying now')
            run_nonpython_process("wget "+ftp_address+"/prosite.doc -O "+directory+"/prosite/prosite.doc")
            print('prosite.doc has been downloaded')
            run_nonpython_process("wget "+ftp_address+"/prosite.dat -O "+directory+"/prosite/prosite.dat")
            print('prosite.dat has been downloaded')
    return directory+"/prosite/"


def process_prositedb_with_prosextract(prosite_dir):
    """
    Uses prosextract to process prosite db for use in prosite search  
    """
    print('Attempting to process prosite db using prosextract')
    run_nonpython_process("prosextract -prositedir "+prosite_dir)
    print('Prosextract method sucessful. Continuing with analysis')
    return prosite_dir    


def get_prosite_db(prosite_db_directory = '/localdisk/software/EMBOSS-6.6.0/share/EMBOSS/data/PROSITE/prosite.lines'):
    if os.path.exists(prosite_db_directory):
        print('Prosite file \"prosite.lines\" was found. Continuing with analysis')
    else: 
        print('Prosite file \"prosite.lines\" was not found. Attempting to create by downloading .dat & .doc files from EMB=EBI')
        prosite_db_directory = get_prosite_db_scratch()
        prosite_db_directory = process_prositedb_with_prosextract(prosite_db_directory)
    return prosite_db_directory


def search_for_motifs(fasta_file_name, directory): 
    """
    Uses prosite to search motifs in a given fasta file 
    """
    f = directory+'/'+fasta_file_name
    run_nonpython_process("patmatmotifs -sequence "+f+" -outfile "+f+".motif")
    print('Prosite search sucessful for sequence:'+str(fasta_file_name)+' output stored in directory:')
           

def get_selected_seq(blast_bin_no, blast_bin, fasta_dict, directory):
    ## Get the acc from bin dict
    ## Look up sequence in fasta_dict 
    directory = os.mkdir(directory+'/blast_bin'+str(blast_bin_no))
    fasta_in_bin = []
    for blast_result in blast_bin:
        name = blast_bin[blast_result]['subject_acc']
        if name not in fasta_dict.keys():
            raise ValueError('Fasta dict does not contain key:'+str(name)+' something has gone wrong processing fastas or binned blast results')
        else: 
            create_new_fasta_record(name, fasta_dict[name]['sequence'], directory)
            fasta_in_bin.append(name)
    return fasta_in_bin


def search_motifs_from_list(fasta_list, directory):
    """
    Wraps the motif finder using prosite for each fasta file in a list. 
    Useful for searching bins but kept as a wrapper in case user wishes to search specific single sequence
    """
    print('Attempting to search for motifs using prosite in sequences:'+str(fasta_list))
    for fasta_file_name in fasta_list:
        search_for_motifs(fasta_file_name, directory)


def create_new_fasta_record(name, sequence, directory):
    """
    Creates a new fasta file for the specified identifier and sequence and writes to given directory
    """
    with open(directory+'/'+str(name)+'.fasta') as f:
        f.write('>'+str(name)+sequence)
    print(name+' has been made into an individual fasta file')

find_motif_no = re.compile('# HitCount:(.*)\n')
find_motif_id = re.compile('Motif = (.*)\n')

def report_motif_results(directory):
    motif_files = [mfile for mfile in os.path.expanduser(directory) if mfile.endswith(".motif")]
    motif_dict = {}

    for mfile in motif_files:
        with open(os.path.expanduser(directory)+'/'+mfile) as f:
            data = f.read()
            motif_no = find_motif_no.search(data).group(1)
            motifs = find_motif_id.findall(data)
            motif_dict[mfile] = {'no_motifs': motif_no, 'motifs:': motifs}
    print('Files Analysed:'+str(len(motif_files))+' motifs found in: '+str(len(motif_dict.keys())))
    for k in motif_dict.keys():
        print('File:'+str(k)+': '+str(motif_dict[k]['no_motifs'])+' motifs found with sequences: '+str(motif_dict[k]['motifs']))      
    print('For more information, please see relevant .motif file in directory:'+str(directory))      


### ---------- Plotting Functions ---------- ### 

def create_barplot_hor(x, y, data, figname, col = sns.colour_palette('Blues')):
    sns.barplot(x = x, y=y, data=data, colour = col)
    plt.title(figname)
    plt.tight_layout()
    plt.savefig(figname+'.png')
    plt.show()


def create_unique_res_line(x, y, data, figname, col = sns.colour_palette('GnBu')):
    pass


def create_stacked_bar_percent(x, y, data, figname, colourlist = sns.colour_palette('GnBu'), ignore_x = False):
    # Using multiple alignment file, create total numbers each AA at each residue (of the consensus sequence)
    # get the % of each AA at this position and plot. 
    pass


def create_plotcon(alignment_file):
    """
    Uses emboss plotcon to make a similarity plot from the clustalo alignment file. 
    """
    print('Attempting to create a similarity plot from clustalo alignment file')
    output = alignment_file.replace('_clustal.msf', '_similarityplot')
    run_nonpython_process("plotcon -sequences \""+alignment_file+"\" -graph svg -goutfile \""+output+"\"")
    print('Similarity Plot created sucessfully. Plotcon output file stored as: '+ str(output)+'.svg')


def plot_motif_counts():
    pass


def create_hydropathy_plot(alignment_file):
    print('Attempting to create Kyte-Doolittle Hydropathy plot for alignment')
    output = alignment_file.replace('clustal.msf', 'hydropathy_alignment')
    run_nonpython_process("pepwindowall \""+alignment_file+"\" -graph svg -gxtitle \'Residue\' -gtitle2 \"Kyte-Doolittle Hydropathy for Alignment\" -goutfile \""+output+"\"")
    print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))

#def create_pepinfo_plots(fasta_file, fasta_file_label):
#    print('Attempting to create pep-info plot for selected fasta file'+str(fasta_file_label))
#    output = fasta_file.replace('.fasta', 'pepinfo')
#    ## Makes graph in file called pepinfo.ps 
#    ## pepinfo "/localdisk/home/ifarquha/testing_patterns/XP_030822190.1" -graph cps -gtitle XP_030822190.1 #
#
  

"""

run_protein_ident('G6pc', 'txid8782', foldertime = True, bin_no = 10)

/localdisk/home/ifarquha/prosite 

"/localdisk/home/ifarquha/Protein_Analysis_G6pc_txid8782_2019_11_10/G6pc_txid8782_consensus.fa"

"/localdisk/home/ifarquha/Protein_Analysis_G6pc_txid8782_2019_11_11_time_13_27/G6pc_txid8782_clustal.fasta"

"/localdisk/home/ifarquha/Protein_Analysis_G6pc_txid8782_2019_11_11/G6pc_txid8782_clustal.msf"



example if motif found: 
########################################
# Program: patmatmotifs
# Rundate: Mon 11 Nov 2019 15:27:53
# Commandline: patmatmotifs
#    -sequence testing_patterns/OPJ74548.1
#    -outfile testing_patterns/OPJ74548.1.motif
# Report_format: dbmotif
# Report_file: testing_patterns/OPJ74548.1.motif
########################################

#=======================================
#
# Sequence: OPJ74548.1     from: 1   to: 358
# HitCount: 1
#
# Full: No
# Prune: Yes
# Data_file: /localdisk/software/EMBOSS-6.6.0/share/EMBOSS/data/PROSITE/prosite.lines
#
#=======================================

Length = 4
Start = position 139 of sequence
End = position 142 of sequence

Motif = AMIDATION

ILSAAAGKKQSRTL
     |  |
   139  142


#---------------------------------------
#---------------------------------------

example if motif not found:





########################################
# Program: patmatmotifs
# Rundate: Mon 11 Nov 2019 15:27:51
# Commandline: patmatmotifs
#    -sequence testing_patterns/XP_030822190.1
#    -outfile testing_patterns/XP_030822190.1.motif
# Report_format: dbmotif
# Report_file: testing_patterns/XP_030822190.1.motif
########################################

#=======================================
#
# Sequence: XP_030822190.1     from: 1   to: 359
# HitCount: 0
#
# Full: No
# Prune: Yes
# Data_file: /localdisk/software/EMBOSS-6.6.0/share/EMBOSS/data/PROSITE/prosite.lines
#
#=======================================


#---------------------------------------
#---------------------------------------


"""
