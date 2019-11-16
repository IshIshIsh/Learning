"""
BPSM Assignment2: Protein Analysis 

A pipeline to analyse groups of related proteins:

1) Basic checks for inputs 
2) Retrival of sequences from NCBI
3) Basic Report (No.Sequences, No.Genes, No.Species)
4) Clustalo alignment & Consensus creation
5) Blastdb creation
6) Blastp analysis
7) Binning results
8) Motif Reporting 
9) Protein Analysis
10) Graphing  

The pipeline can be used as a script, using a master function or using each 
task's master function individually (as long as formatting of inputs is correct)

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
import sys  
import json 


### ---------- Parameter Dict ---------- ### 

parameter_dict = {
	'name': [str, ''],
	'working_directory': [str, str(Path().absolute())],
	'filtered_partial': [bool, True],
	'filtered_partial': [bool, True], 
	'save_summary': [bool, True],
	'silent': [bool, True], 
	'bin_no': [int, 250], 
	}

### ---------- Master Function [whole script] ---------- ### 

def master_analysis_requests(gene_family, taxid, defaults_requests, parameter_dict=parameter_dict):
	"""
	Used as master function when parameters are entered via user prompt
	"""
	check_edirect_installation()
	gene_family = replace_non_alphanumeric_chars(gene_family)
	pathtaxid = replace_non_alphanumeric_chars(taxid)
	name = gene_family+'_'+pathtaxid+'_'
	parameter_dict['name'] = [str, name]
	if defaults_requests == True: 
		for k,v in parameter_dict.items():
			newvalue = input('Please enter a parameter value for '+k+' of type: '+str(v[0])+' or press enter to use default value:'+str(v[1]))
			
	pass

def master_analysis(gene_family, taxid, parameter_dict):
	"""
	Used as master function when parameters are passed to function
	"""
	pass

### ---------- Master Functions [each task seperately] ---------- ### 

def master_seq_retrieve(working_directory, gene_family, taxid, ncbi_file_name, filtered_partial, filtered_predicted, silent, save_summary):
	ncbi_data_path = get_from_ncbi(gene_family, taxid, ncbi_file_name, working_directory, filtered_partial, filtered_predicted, silent)
	get_species_no_from_taxid(gene_family)
	ncbi_data_summary(working_directory, ncbi_data_path, gene_family, taxid, filtered_partial, filtered_predicted, save_summary)
	return ncbi_data_path


def master_clustalo(fasta_path, iteration, thread_process_dict, silent = True):
	alignment = align_with_clustalo(fasta_path, iteration, thread_process_dict, silent)   
	consensus = get_consensus_from_alignment(alignment)
	return alignment, consensus 


def master_protein_analysis(fasta_path):
	split_fasta = check_read_no_by_split(fasta_path, '>')
	fasta_dict = check_protein_basic_stats(split_fasta)
	protein_stats = get_protein_stats_from_pepstats(split_fasta)
	fasta_dict = add_pepstats_info_to_dict(fasta_dict, protein_stats)
	return fasta_dict    


def master_blast(working_directory, outputname, fastafile_path, fasta_dict, bin_no, thread_process_dict, silent):
	blastdb_filepath = create_blastdb(working_directory, outputname, fastafile_path)
	query_filepath = query_choice(fastafile_path, fasta_dict)
	blast_output = query_blastdb(blastdb_filepath, query_filepath, thread_process_dict,	silent)
	bins = bin_results(blast_output, bin_no)
	return blast_output, bins


def master_motifs(analysispath, fasta_dict, bins):
	"""
	Parent Function for searching fasta sequences from selected bins using Prosite and creating a summary of results: 

	Vars:
		analysispath: (str) path to directory containing the analysis output files
		fasta_dict: (dict) dictionary containing various useful information on the fasta inputs
		selected_bin_labels: (int, list[int]) the labels of the bins to continue analysis on
		bins: (dict) containing the bin and fasta_id
	Returns:
		fasta_dict: (dict) of key(fasta_id) with values(summary dict of analysis results) for further processing  
	"""
	bin_labels = bins.keys()
	selected_bin_labels = bin_selection(bin_labels, 'Prosite Motif search')
	get_prosite_db()
	if not os.path.isdir(os.path.expanduser(analysispath+'/motifs')):
		os.mkdir(analysispath+'/motifs')
		if isinstance(selected_bin_labels, list):
			for selected_bin_label in selected_bin_labels:
				if not os.path.isdir(os.path.expanduser(analysispath+'/motifs/bin_'+str(selected_bin_label))):
					os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
				print('Creating individual fasta records for: bin_'+str(selected_bin_label))
				fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
				search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label))
				print('Prosite scan completed for bin:'+str(selected_bin_label))
				print(' ')
		else: 
			selected_bin_label = selected_bin_labels
			if not os.path.isdir(os.path.expanduser(analysispath+'/motifs/bin_'+str(selected_bin_label))):
				os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
			print('Creating individual fasta records for: bin_'+str(selected_bin_label))
			fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
			for fasta in fasta_list_from_bin:
				create_new_fasta_record(fasta, fasta_dict[fasta]['sequence'], analysispath+'/motifs/bin_'+str(selected_bin_label))
			search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label))
			print('Prosite scan completed for bin:'+str(selected_bin_label))
			print(' ')
		for folder in os.listdir(os.path.expanduser(analysispath+'/motifs/')):
			report_motif_results(analysispath+'/motifs/'+folder)
		print('For more information on motifs per sequence, please see relevant files in folder:'+analysispath+'/motifs/')
	return fasta_dict 


def master_graphs(clustalo_files):
	create_hydropathy_plot(clustalo_files)
	create_plotcon(clustalo_files)


### ---------- Utility Functions ---------- ###

def user_input_formating(input_string):
	"""
	Internal function to parse inputs into true, false, default or quit if given 
	in format described in conversion dict. Otherwise returns the input string for
	further parsing.
	"""
	conversion_dict = { True : ['y', 'yes', 'true'],
						False : ['n', 'no', 'false'],
						'default' : ['.', ''],
						'terminate' : ['q', 'quit']}
	input_string = input_string.lower().strip()
	for k, v in conversion_dict.items():
		if input_string in v and k != 'terminate':
			return k 
		elif input_string in v and k == 'terminate': 
			print('INTERUPT: User terminated process with input: '+input_string)
			sys.exit(0)
		else:
			return input_string


def check_type(input_string, expected_type, parameter_name, parameter_default = None, user_input = True):
	"""
	Internal function: Checks the input against an expected type and either throws an error or returns default value 

	Vars: 
		input_string: (ANY) value of the parameter to check the type of
		expected_type: (ANY) the type the parameter should be
		parameter_name: (str) the name of the parameter which is being checked 
		parameter_default: (ANY/None) the value to return if the type is wrong.  
	Returns: 
		input_string IF value of parameter was of correct type
		default IF value of parameter was not of the correct type but a valid default value was passed in 
		ValueError IF value of parameter was not correct type and no acceptable default value for parameter was passed in 
	"""
	error = False
	error_message = ('WARNING: The input: '+str(input_string)+'  for parameter:'+parameter_name+' was not of the expected type: '+
						str(expected_type))'.
	if input_string == 'default':
		return parameter_default
	if expected_type == int:
		if input_string.isnumeric():
			return int(input_string)
		else: 
			error = True
	elif expected_type == dict or expected_type == list: 
		try: 
			input_object = json.loads(input_string)
			if expected_type == dict and isinstance(input_object, dict):
				return input_object
			elif expected_type = list and isinstance(input_object, list):
				return input_object
			else: 
				error = True
		except:
			error = True
	elif expected_type == string and isinstance(input_string, str):
		return input_string
	if error == True and parameter_default == None: 
		raise ValueError(error_message)
	elif error == True and parameter_default != None: 
		print (error_messgae)
		print('Attempting to return to default value:'+str(parameter_default))
		return parameter_default
	else: 
		raise Exception('Check_type function should not reach this point. An internal error has occured')


def parameter_input(parameter_name, parameter_default, task_name, expected_type):
	"""
	Internal function to enable auto-format of parameter user input request
	"""
	print(task_name + " can have a parameter for: "+parameter_name+" or will use default value: "+str(parameter_default))
	input_string = input('Please enter a value or press enter to use default value')
	input_string = user_input_formating(input_string)
	input_string = check_type(input_string, expected_type, task_name, parameter_default)
	return 


def replace_non_alphanumeric_chars(input_string, replacement_char = '_', ignore_char = 'a'): 
	"""
	Utility function used to replace any non alphanumeric charecter with an option to ignore specified string

	Vars: 
		input_string: (str) the string to check/replace non-alphanumeric charecters in
		replacement_char: (str) the string to replace any non-alphanumeric charecter with
		ignore_char: (str) [default: 'a'] a string to add to the regex to allow an non-alphanumeric charecter to 
			be ingnored. If default value ('a') no non-alphanumeric charecters are allowed.
	Returns:
		output:(str) the input_string with charecters replaced with replacement_char if they were not allowed in 
			function call
	"""
	if ignore_char == 'a':
		output = re.sub('[^0-9a-zA-Z]+', replacement_char, input_string)
	else: 
		output = re.sub('[^0-9a-zA-Z'+ignore_char+']+', replacement_char, input_string)
	return output 


def check_thread_format(thread_input, cpu, default_threads):
	"""
	Internal function for processing thread input and printing errors or warnings
	Vars
	"""
	try: 
		input_threads = int(thread_input)
	except: 
		print('WARNING: Threads was not a valid integer: '+thread_input+' returning to default value of:'+str(default_threads))
	if input_threads < 1 or input_threads > cpu:
		print('WARNING: Thread number cannot be less than 1 or greater than the max cpu count: '+str(cpu)+' but was:'+str(input_threads)+'  returning to default value of:'+str(default_threads))
	else: 
		default_threads = input_threads
	return default_threads 


def auto_thread():
	"""
	Internal function which calculates and returns the available cpu count and default thread number
	"""    
	cpu = os.cpu_count()
	default_threads = int(cpu/3)
	return cpu, default_threads


def thread_entry(thread_processes):
	"""
	Internal function. Calculates the max cpu threads available and default values (an integer of 1/3rd available cpu)
	Lists processes which use threading and asks user to either set the thread number for all or individually
	Vars: 
		thread_processes: (list) containing all processes which can use a threading input
	Returns:
		thread_process: (dict) mapping thread_number to process
	"""
	thread_process_dict = {}
	print('The following processes can use threading: '+str(thread_processes))
	input_string = input('Would you like to use the same threading for all?')
	input_string = user_input_formating(input_string)
	input_string = check_type(input_string, bool, 'thread_number', True)
	cpu, default_threads = auto_thread()
	print('Your cpu count is: '+str(cpu)+' valid thread number range is: 1:'+str(cpu)+'. Default value is: '+str(default_threads))
	if input_string == False: 
		for process in thread_processes:
			input_string = ('Please enter a valid integer to use as thread number to use for process: '+str(process)+' or press enter to use default')
			input_string = user_input_formating(input_string)
			threads = check_thread_format(input_string, cpu, default_threads)
			thread_process_dict[process] = threads
	else: 
		for process in thread_processes: 
			thread_process_dict[process] = threads
	return thread_process_dict 


def run_nonpython_process(query, timeout = 6000, silent = True):
	"""
	Function wich uses subprocess rather than os due to increased security.

	Vars:
		query: (str) a string containing a valid bash command for processing through subprocess module
		timeout: (int) the number of seconds to wait for a response ( Note TO SELF: Timeout doesn't seem to actually cause a timeout error???) 
		silent: (boolean) [default: True] if false prints the query string 
	returns
		response from subprocess module call of query command 
	"""
	try: 
		if silent != True:
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
	Checks that edirect is accessible/installed in the users home directory
	If accesible installation found allows user prompt to update(reinstall) or continue
	if no installation found. Runs a bash command through subprocess module to install it in the users home directory 
	"""
	print ('Checking edirect is installed in home space...')
	if os.path.isdir(os.path.expanduser('~/edirect')):
		updateinput = input('Edirect is already installed, to update please input \'update\' or press any key to continue with analysis')
		if updateinput.strip().lower() == 'update': 
			print('Attempting to update Edirect now')
			run_nonpython_process('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"')
			print('Edirect has been sucessfully updated. Continuing with analysis...')
			print(' ')
	else: 
		print('Edirect is not currently installed, trying to install now')
		run_nonpython_process('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"')
		print('Edirect has been sucessfully installed. Continuing with analysis')
		print(' ')


### ---------- File and Directory Processing ---------- ###  

def return_file_contents(file):
	"""
	Internal function: Basic file opener to return file contents
	"""
	with open(file) as f:
		data = f.read()
	return data 


def check_read_no_by_split(input_file, splitchar): 
	"""
	Internal Function: Processes an input file into lines based on splitchar
	Vars:
		input_file: (str) a valid file with file extension (if not in current directory)
		splitchar: (int/str) a charecter to split the data by
	Returns: 
		datasplit: (list) contains the data contained in the file split by specified charecter
	"""
	data = return_file_contents(input_file)
	datasplit = data.split(splitchar)
	if isinstance(datasplit, list) == False: 
		raise ValueError('input file could not be split by specified charecter:'+str(splitchar))
	datasplit = list(filter(None, datasplit))
	return(datasplit)


def create_folder_path(name, filepath, foldertime = False):
	"""
	Checks if the inputted folder (name + prefix + filepath) exists 
	if not it will create it. Both cases will result in a print statement 
	and return a string of the folder path to store results in.

	Vars: 
		name: (str/None) If no user defined name, it will default to family_taxid
		filepath: (str) if no user defined input it will default to current working directory 
	Returns:
		folder: (str) the path to directory used for analysis output
	"""
	if name == None:
		prefix = 'Protein_Analysis_'
	else: 
		prefix = name+'_'
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


def convert_filetypes_seqret(msf_infile, fasta_outfile):
	"""
	TODO: THIS ISNT USED 
	Converts a msf clustal alignment to fasta record version for processing differently
	"""
	print('Attempting to convert msf alignment to fasta')
	run_nonpython_process("seqret -sequence \""+msf_infile+"\" -outseq  \""+fasta_outfile+"\"")
	print('msf file converted to fasta sucessfully. New file stored in: '+ str(fasta_outfile))
	print(' ')


### ---------- NCBI Data Collection ---------- ### 

hit_count = re.compile('<Count>([0-9]*)<')
species_fasta = re.compile('>.*\[(.*)\]')


def rescursive_ncbi_string(query_string, entry_list, join, entry_type):
	"""
	Internal Function: Recursivly adds each entry in the entry list to a query string,
	adding the join string (if recursion has occured previously) and the entry type and
	returns updated string. 
	"""
	counter = 0
	for entry in entry_list:
		if counter > 0:
			query_string = query_string + join
		query_string = query_string+str(entry)+entry_type
	return query_string


def ncbi_string_processing(gene_family, taxid, filtered_partial = True, filtered_predicted = True):
	"""
	Internal function: Processes inputs to create an NCBI valid query string for data retrival 
	"""
	query_string = ''
	if isinstance(taxid, list): 
		query_string = rescursive_ncbi_string(query_string, taxid, ' | ', '[Organism:exp]')
	elif isinstance(taxid, str): 
			query_string = query_string+taxid+'[Organism:exp]'
	query_string = query_string + ' AND '
	if isinstance(gene_family, list):
		query_string = rescursive_ncbi_string(query_string, gene_family, ' | ', '[Gene Name]')
	if filtered_partial == True and filtered_predicted == True:
		print('WARNING: filtered_partial is True and filtered_predicted is True so Partial or Predicted sequences will be excluded')
		query_string = query_string+' NOT partial[All Fields] NOT predicted[All Fields]'
	elif filtered_partial == False and filtered_predicted == True: 
		print('WARNING: filtered_predicted is True so Predicted sequences will be excluded')
		query_string = query_string+' NOT predicted[All Fields]'
	elif filtered_partial == True and filtered_predicted == False: 
		print('WARNING: filtered_partial is True so Partial sequences will be excluded')
		query_string = query_string+' NOT partial[All Fields] NOT predicted[All Fields]'
	else: 
		print('Note: Partial and Predicted sequences will be included in output')
	return query_string


def get_from_ncbi(gene_family, taxid, file_name, working_directory, filtered_partial = True, filtered_predicted = True, silent = True):
	"""
	Requests data in fasta format from the ncbi protein database and if silent is not True pipes this to screen
	Saves the output as a single output file in the working_directory as the file_name.fasta
	Prints a sucess message to screen containing file path to data if no errors encountered. 
	"""
	query_string = ncbi_string_processing(gene_family, taxid, filtered_partial, filtered_predicted)
	print('Attempting to retrive data from ncbi for query')
	print('Query searched is '+str(query_string))
	save_location = working_directory+'/'+file_name+'.fasta'
	edirect_query = "~/edirect/esearch -db protein -query \""+query_string+"\" | ~/edirect/efetch -format fasta" 
	edirect_response = run_nonpython_process(edirect_query)
	if silent != True:
		print(edirect_response)
	with open(save_location, 'w') as result_file:
		result_file.write(edirect_response)
	print('Data Retrival sucessful, fasta files stored in '+str(save_location))
	print(' ')
	return save_location   


def ncbi_data_summary(working_directory, ncbi_data_path, gene_family, taxid, filtered_partial, filtered_predicted, save_summary = False):
	"""
	Prints a basic report to screen with the option to save to the working directory 
	Includes: 
		Number of sequences
		Number of genes
		Number of species 
		Filter options used 
	"""
	fasta_data = return_file_contents(ncbi_data_path)
	no_sequences = fasta_data.count()
	species = species_fasta.findall(fasta_data)  
	no_species = len(species)
	if isinstance(gene_family, list):
		no_genes = len(gene_family)
	else:
		no_genes = 1
	print('Analysis includes:'+str(no_sequences)+' protein sequences from '+str(no_species)+' including '+str(no_genes)+' genes')
	input_string = input('Would you like to continue with the analysis?')
	if user_input_formating(input_string) != True:
		print('Processing Terminated: Interupt by user')
		sys.exit()
	if no_genes > 10000:
		input_string = input("WARNING: Number of sequences is greater than 10,000. Would you like to continue? This will make things slow")
		if input_string != True:
			print('Processing Terminated: Interupt by user')
			sys.exit()
	if save_summary == True:
		with open(working_directory+'/NCBI_retrival_summary.txt', 'w') as f:
			print('Protein Sequences:'+str(no_sequences), file = f)
			print('Species Number:'+str(no_species), file = f)
			print('Gene Number:'+str(no_genes), file = f) 
			print('Search options: Filtering Partial is '+str(filtered_partial)+', Filtering predicted is '+str(filtered_predicted), file = f)
			print('', file = f)
			print('Species included:'+str(species), file = f)
			print('Taxonomic Id(s) searched:'+str(gene_family), file = f)
			print('Gene(s) searched:'+str(gene_family), file = f)
		print('NCBI_retrival_summary.txt was saved in working directory')


def get_species_no_from_taxid(gene_family):
	edirect_query = "~/edirect/esearch -db protein -query \""+gene_family+"[Subtree]"+"\""
	edirect_response = run_nonpython_process(edirect_query)
	counts = hit_count.search(edirect_response)
	print('INFO: Total number of subtree taxonomy groups in taxid:'+str(gene_family)+' is '+str(counts))
	edirect_query = "~/edirect/esearch -db protein -query \""+gene_family+"[Subtree] AND species[rank]"+"\""
	edirect_response = run_nonpython_process(edirect_query)
	counts = hit_count.search(edirect_response)
	print('INFO: Total number of subtree species in taxid:'+str(gene_family)+' is '+str(counts))


### --------- Clustalo Alignment ---------  ####

def align_with_clustalo(fasta_path, iteration, thread_process_dict, silent = True):
	"""
	Uses the clustalo package through bash (through subprocess) to generate a multiple sequence alignment file
	fastafile_location: The full path of the fasta file location 
	iteration: the number of iterations for the clustalo tree
	overwrite_max_seqno: default value False, if true removes limit on the number of sequences WARNING: DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING!!!
	"""
	print('Attempting Clustalso Alignment...')
	clustal_output = fasta_path.replace('.fasta', '_clustal.msf')
	threads = thread_process_dict['clustalo alignment']
	fasta_data = return_file_contents(fasta_path)
	if fasta_data.count() > 10000: 
		input_string = input('WARNING: Number of reads inputted is greater than maxseqno allowed (10,000). Are you certain you want to continue? It will be very slow.')
		if user_input_formating(input_string) != True:
			print('Processing Terminated: Interupt by user')
			sys.exit()
	edirect_reponse = run_nonpython_process("clustalo -i \""+fasta_path+"\" -t Protein --outfmt=msf --full --threads "+str(threads)) 
	if silent != True:
		print(edirect_reponse)
	with open(clustal_output, 'w') as result_file:
		 result_file.write(edirect_reponse)
	print('Alignment sucessful, clustal files stored in '+str(clustal_output))
	print(' ')
	return clustal_output


def get_consensus_from_alignment(clustal_output):
	"""
	Uses EMBOSS Cons to get a consensus sequence from the Clustalo alignment files. 
	"""
	print('Attempting to create a consensus sequence from Clustalo output')
	output = clustal_output.replace('clustal.msf', 'consensus.fa' )
	run_nonpython_process("cons -sequence \""+clustal_output+"\" -outseq \""+output+"\" -sprotein1 true -sid1 \"consensus\" -osformat2 \"fasta\"")
	print('Consensus creation sucessful. Consensus file stored in: '+ str(output))
	print('Consensus Sequence from multiple alignment:')
	print(output)
	print(' ')
	return output 


### ---------- BLAST ---------- #### 

def create_blastdb(working_directory, outputname, fastafile_path):
	"""
	Uses makeblastdb to make a blast-searchable database from a file containing multiple fastas and save them to an output for use in blastp processing. 
	"""
	print('Attempting to create a blastdb from fasta files')
	output = working_directory+'/'+outputname+'blastdb'
	run_nonpython_process("makeblastdb -in \""+fastafile_path+"\" -dbtype prot -out \""+output+"\"")
	print('blast database creation sucessful. Blastdb file stored in: '+ str(output))
	print(' ')
	return output     


def query_choice(filepath, fasta_dict):
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
	print(' ')
	return filepath+'/'+selectedfile


def query_blastdb(blastdb_filepath, query_filepath, thread_process_dict,  silent = True):
	"""
	Uses blastp query a blast-searchable database with selected query_file. 
	"""
	print('Attempting to run blastp for selected file:'+ str(query_filepath)+' against blastdb:'+str(blastdb_filepath)+'with sequence:')
	print(return_file_contents(query_filepath))
	outputfile = query_filepath.replace('.fasta', '.fa').replace('.txt', '.fa')
	output = outputfile.replace('.fa','_blastp_output.txt')
	threads = str(thread_process_dict['blastp query'])
	edirect_reponse = run_nonpython_process("blastp -db \""+blastdb_filepath+"\" -query \""+query_filepath+"\" -num_threads \""+threads+"\" -outfmt 7 > \""+output+"\"")
	if silent != True:
		print(edirect_reponse)
	print('blastp sucessful. Output file stored in: '+ str(output))
	print(' ')
	return output 


def bin_results(input_file, bin_no, ignore_lines = '#', print_bin_contents = False):
	print('Binning results from Blastp')
	result_bin = {}
	data = return_file_contents(input_file)
	lines = [line for line in data.split('\n')]
	validlines = [line for line in lines if line.startswith(ignore_lines) == False]
	validlines = list(filter(None, validlines))
	## As lists have defined order we can parse this using the first as the top hits
	totalno = len(validlines)
	no_bins = math.ceil(totalno/bin_no)
	item_no = 0
	rank_counter = 0 
	for i in range(0, no_bins-1):
		result_bin[i] = []
		for j in range(0,bin_no-1):
			if item_no < totalno:
				result_bin[i].append(convert_blastline_to_dict(validlines[item_no], rank_counter))
				rank_counter = rank_counter + 1
				item_no =  item_no + 1
			else: 
				pass
	bins = result_bin.keys()
	bins_created = len(bins)
	print(str(bins_created)+ ' bins were created')
	if print_bin_contents != False:
		for bin_no in result_bin: 
			print(bin_no)
			bins_in_no = []
			for i in result_bin[bin_no]:
				bins_in_no.append(i['subject_acc'])
			print('contains: '+str(bins_in_no))
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
	print('Bins for analysis: '+str(method)+' are '+str(list(valid_bin_list)))
	input_bins = input('Please select which bin number to continue with, or enter a list seperated by \',\'. Default value is the top scoring group (0) ')
	input_bins = replace_non_alphanumeric_chars(input_bins, '_', ignore_char= ',').replace('_','')
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
	"""
	TODO NOT USED.
	"""
	headers=['query_acc.', 'subject_acc.', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
	blastpd = pd.read_csv(inputfile, comment = ignored, sep ='\t', names = headers)
	return blastpd

### ---------- Protein Stats ---------- ### 
short_name_regex = re.compile('^[\S]*')
short_name_regex_pepstats = re.compile('^([\S]*) from')
mol_weight = re.compile('Molecular weight = (\S*)')
isoelectric = re.compile('Isoelectric Point = (\S*)')
residues = re.compile('Residues = (\S*)')
charge = re.compile('Charge   = (\S*)')
Av_res_weight = re.compile('Average Residue Weight  = (\S*)')


def get_protein_stats_from_pepstats(fastafile_location):
	"""
	Uses EMBOSS pepstats to get protein statistics from the fasta files and save them to an output. 
	"""
	print('Attempting to get protein statistics from pepstats')
	output = fastafile_location.replace('.fasta', 'stats.pepstats')
	run_nonpython_process("pepstats -sequence \""+fastafile_location+"\" -outfile \""+output+"\"")
	print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
	print(' ')
	return output 


def check_protein_basic_stats(fasta_list):
	if isinstance(fasta_list, list) == False:
		raise ValueError('Fasta_list input was not of type list. No dictionary created from fasta files.')
	fasta_dict = {}
	for f in fasta_list:
		name = f[0:f.find('\n')]
		label = short_name_regex.search(name)
		seq = f[f.find('\n'):]
		seqlength = len(seq.replace('\n', ''))
		fasta_dict[label.group()] = {'full_name': name, 'sequence': seq, 'length': seqlength}
	return fasta_dict 


def add_pepstats_info_to_dict(fasta_dict, pepstatsfile):
	data = return_file_contents(pepstatsfile)
	entries = list(filter(None, data.split('PEPSTATS of ')))
	for entry in entries:
		short_name = short_name_regex_pepstats.search(entry).group(1)
		if short_name in fasta_dict.keys():
			fasta_dict[short_name]['Mol_Weight'] = mol_weight.search(entry).group(1)
			fasta_dict[short_name]['Residues'] =residues.search(entry).group(1)
			fasta_dict[short_name]['Av_residue_weight'] =Av_res_weight.search(entry).group(1)
			fasta_dict[short_name]['Charge'] =charge.search(entry).group(1)
			fasta_dict[short_name]['Isoeelectric_point'] =isoelectric.search(entry).group(1)
		else: 
			raise ValueError('Short Label:'+str(short_name)+' was found in pepstats output but not in fasta_dict')
	return fasta_dict


def check_alignment_basic_stats(alignment_file, protein_dict):
	"""
	TODO THIS ISN'T USED
	"""
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


def search_for_motifs(fasta_file_name, directory, silent = True): 
	"""
	Uses prosite to search motifs in a given fasta file 
	"""
	f = directory+'/'+fasta_file_name+'.fasta'
	fout = f.replace('.fasta', '.motif')
	run_nonpython_process("patmatmotifs -sequence "+f+" -auto Y -outfile "+fout)
	if silent != True:
		print('Prosite search sucessful for sequence:'+str(fasta_file_name)+' output stored in directory:')
		print(' ')
		   

def get_selected_seq(blast_bin_no, blast_bin, fasta_dict, directory):
	## Get the acc from bin dict
	## Look up sequence in fasta_dict 
	fasta_in_bin = []
	for blast_result in blast_bin:
		name = blast_result['subject_acc']
		if name not in fasta_dict.keys():
			raise ValueError('Fasta dict does not contain key:'+str(name)+' something has gone wrong processing fastas or binned blast results')
		else: 
			create_new_fasta_record(name, fasta_dict[name]['sequence'], directory)
			fasta_in_bin.append(name)
	print('New fasta files were created for: bin'+str(blast_bin_no))
	return fasta_in_bin


def search_motifs_from_list(fasta_list, group_name, directory, keep_fastas = False):
	"""
	Wraps the motif finder using prosite for each fasta file in a list. 
	Useful for searching bins but kept as a wrapper in case user wishes to search specific single sequence
	"""
	print('Attempting to search for motifs using prosite in sequences:'+str(group_name))
	for fasta_file_name in fasta_list:
		search_for_motifs(fasta_file_name, directory)
	if keep_fastas == False: 
		for fasta_file_name in fasta_list: 
			f = directory+'/'+fasta_file_name+'.fasta'
			os.remove(f)


def create_new_fasta_record(name, sequence, directory):
	"""
	Creates a new fasta file for the specified identifier and sequence and writes to given directory
	"""
	with open(directory+'/'+str(name)+'.fasta', 'w+')  as f:
		f.write('>'+str(name)+sequence)


find_motif_no = re.compile('# HitCount:(.*)\n')
find_motif_id = re.compile('Motif = (.*)\n')


def report_motif_results(directory):
	"""
	Prints to screen a summary report on the motifs found from PROSITE data in the given directory
	Does not include reports where no motifs were found
	eg [protein_id] has [INT] motifs found with prosite name [PROSITE ENTRY NAME] 
	"""

	motif_files = [mfile for mfile in os.listdir(directory) if mfile.endswith(".motif")]
	motif_dict = {}
	for mfile in motif_files:
		with open(os.path.expanduser(directory)+'/'+mfile) as f:
			data = f.read()
			motif_no = find_motif_no.search(data).group(1)
			motifs = find_motif_id.findall(data)
			if motif_no > 0:
				motif_dict[mfile] = {'no_motifs': motif_no, 'motifs': motifs}
	print('Files Analysed:'+str(len(motif_files))+' motifs found in: '+str(len(motif_files)))
	for k in motif_dict.keys():
		print('File:'+str(k)+': '+str(motif_dict[k]['no_motifs'])+' motifs found with Prosite motif entry name(s): '+str(motif_dict[k]['motifs']))      
	print('For more information, please see relevant .motif file in directory:'+str(directory))      
	print(' ')


### ---------- Plotting Functions ---------- ### 

def create_hydropathy_plot(alignment_file):
	print('Attempting to create Kyte-Doolittle Hydropathy plot for alignment')
	output = alignment_file.replace('clustal.msf', 'hydropathy_alignment')
	run_nonpython_process("pepwindowall \""+alignment_file+"\" -graph svg -gxtitle \'Residue\' -gtitle2 \"Kyte-Doolittle Hydropathy for Alignment\" -goutfile \""+output+"\"")
	print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
	print(' ')


def create_plotcon(alignment_file):
	"""
	Uses emboss plotcon to make a similarity plot from the clustalo alignment file. 
	"""
	print('Attempting to create a similarity plot from clustalo alignment file')
	output = alignment_file.replace('_clustal.msf', '_similarityplot')
	run_nonpython_process("plotcon -sequences \""+alignment_file+"\" -graph svg -goutfile \""+output+"\"")
	print('Similarity Plot created sucessfully. Plotcon output file stored as: '+ str(output)+'.svg')
	print(' ')


### ---------- TO RUN AS SCRIPT ---------- ### 

"""
To run as a script in command line with the example genefamily = G6pc and taxid = txid8782 with prompts for options: 

python3 Protein_Analysis.py G6pc txid8782

To run as a script in command line with the same example but using all default values: 

python3 Protein_Analysis.py G6pc txid8782 defaults

To run as a script in command line with multiple entries for genefamily or taxid use a comma seperated string (NO SPACES)

python3 Protein_Analysis.py G6pc,COX1 txid8782

Note: This limits the input arguments to a single genefamily and taxid and results in user 
prompts for all other variables. For more complicated usage please import the module and call 
the master function or each tasks master function.  
"""

if __name__ == "__main__":
	gene_family = sys.argv[1]
	taxid = sys.argv[2]
	if len(sys.argv) > 3:
		if sys.argv[3] == 'defaults':
			defaults_requests = True
		else: 
			raise Exception('Third argument passed to script was not recognised. Terminating process')
	else: 
		defaults_requests = False
	master_analysis(gene_family, taxid, defaults_requests)


"""
Things to include if time allows: 


def get_child_taxid():
esearch -db taxonomy -query "vertebrata[orgn]" | efetch -db taxonomy -format docsum | xtract -pattern DocumentSummary -if Rank -equals family -element Id,Division,ScientificName,CommonName | more

def get_gene_alias():
esearch -db gene -query "Liver cancer AND Homo sapiens" | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element Name OtherAliases OtherDesignations
glucose-6-phosphatase AND Aves 
"""


