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
"""
Method of storing parameter_name, expected_type and default values to be parsed for user inputs.
Additional vars can be added in and treated the same way easily
Global default. New one can be passed in. 
	- TODO checking for k, v pairs in passed in dict 
"""
global_parameter_dict = {
	'gene_family': ['StrList', None],
	'taxid': ['StrList', None], 
	'name': [str, ''],
	'working_directory': [str, str(Path().absolute())],
	'filtered_partial': [bool, True],
	'filtered_predicted': [bool, True], 
	'save_summary': [bool, True],
	'silent': [bool, True], 
	'bin_no': [int, 250], 
	'folder_time': [bool, False],
	'clustalo_iteration': [int, 30],
	'keep_fastas': [bool, True]
	}

### ---------- Master Function [whole script] ---------- ### 

def master_analysis(gene_family, taxid, defaults_requests = True, parameter_dict=global_parameter_dict, **kwargs):
	"""
	Used as master function for automated analysis. 
	Three use modes: 
		- defautlt_requests = True: User prompt for each parameter value, with defaults and type checking
		- default_requests = False: automatically uses default values without prompting
		- **kwargs = dict to pass in values to overwrite defaults
	Vars: 
		gene_family (str) a string of either a string or list (with appropriate escaped " or pairs of "') of gene names to search 
		taxid (str) a string of either a string or list (with appropriate escaped " or pairs of "') of taxonomic IDs to search in format txid0000
		defaults_requests [bool] [true] whether to prompt user for parameter values or use defaults
		parameter_dict (dict) [global_parameter_dict] dict containinng parameter_name: [parameter expected value, parameter default]
		**kwargs (dict) containing k,v arguments for overwriting global_parameter_dict
	Produces: 

	Returns: 
		fasta_dict [dict], a dict containing basic protein info, sequence_id and sequence
		fastadb: [pandas dataframe] a dataframe of fasta_dict for plotting vars
		blastdb: [pandas dataframe] a dataframe of blast results for plotting vars
	"""
	if defaults_requests == True and len(kwargs.items()) > 0:
		print('WARNING: Default requests is true but key word args were entered, ignoring key word arguments and prompting for var entry')
	check_edirect_installation()
	pathgene_family = replace_non_alphanumeric_chars(gene_family)
	pathtaxid = replace_non_alphanumeric_chars(taxid)
	name = pathgene_family+'_'+pathtaxid+'_'
	parameter_dict['name'] = [str, name]
	if defaults_requests == True: 
		for k,v in parameter_dict.items():
			if k == 'gene_family':
				gene_family = check_type(gene_family, v[0],'gene_family', v[1])
				print(type(gene_family))
				parameter_dict['gene_family'] = [v[0], gene_family]
			elif k == 'taxid': 
				taxid = check_type(taxid, v[0],'taxid', v[1])
				parameter_dict['taxid'] = [v[0], taxid]
			else:
				newvalue = user_input_formating(input('Please enter a parameter value for '+k+' of type: '+str(v[0])+' or press enter to use default value:'+str(v[1])))
				newvalue = check_type(newvalue, v[0], k, v[1])
				parameter_dict[k] = [v[0], newvalue]
	else:
		entry_type, default_val = parameter_dict['gene_family']
		gene_family = check_type(gene_family, entry_type, 'gene_family', default_val)
		parameter_dict['gene_family'] = [entry_type, gene_family]
		entry_type, default_val = parameter_dict['taxid']
		parameter_dict['taxid'] = [entry_type, check_type(taxid, entry_type, 'taxid', default_val)]
		for k, v in kwargs.items():
			if k in parameter_dict.items():
				newvalue = check_type(v, parameter_dict[k][0], parameter_dict[k][1])
				parameter_dict[k] = [v, newvalue]
			else: 
				raise ValueError("Kwards argument for k:"+str(k)+" was not regonised as a valid keyword argument. Proceess terminated.")
	print('')		 
	thread_process_dict = thread_entry(['clustalo alignment', 'blastp query', 'similarity matrix'])
	working_dir0 = create_folder_path(parameter_dict['name'][1], parameter_dict['working_directory'][1], parameter_dict['folder_time'][1])
	ncbi_data_path = master_seq_retrieve(working_dir0, gene_family, taxid, parameter_dict['name'][1], parameter_dict['filtered_partial'][1], parameter_dict['filtered_predicted'][1], parameter_dict['silent'][1], parameter_dict['save_summary'][1])
	alignment_file_path, consensus = master_clustalo(ncbi_data_path, parameter_dict['clustalo_iteration'][1], thread_process_dict, parameter_dict['silent'][1])
	fasta_dict = master_protein_analysis(ncbi_data_path)
	blast_data_path, bins = master_blast(working_dir0, parameter_dict['name'][1], ncbi_data_path, fasta_dict, parameter_dict['bin_no'][1], thread_process_dict, parameter_dict['silent'][1])
	fasta_dict = master_motifs(working_dir0, fasta_dict, bins, parameter_dict['keep_fastas'][1],parameter_dict['save_summary'][1])
	fastadb, blastdb = master_graphs(alignment_file_path,ncbi_data_path, thread_process_dict,  blast_data_path, fasta_dict)
	print('Anaylsis Complete')
	return fasta_dict, fastadb, blastdb
	

### ---------- Master Functions [each task seperately] ---------- ### 

def master_seq_retrieve(working_directory, gene_family, taxid, ncbi_file_name, filtered_partial= True, filtered_predicted= True, silent= True, save_summary = True):
	"""
	Parent Function for sequence retrival from NCBI
	Vars:
		working_directory: (str, path) path to directory to store results in
		gene_family (str) a string of either a string or list (with appropriate escaped " or pairs of "') of gene names to search 
		taxid (str) a string of either a string or list (with appropriate escaped " or pairs of "') of taxonomic IDs to search in format txid0000
		ncbi_file_name: (str) prefix name for seqs.fasta output
		filtered_partial (bool) [True] whether to add or exclude partial sequences from the fasta results
		filtered_predicted (bool) [True] whether to add or exclude predicted sequences from the fasta results
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
		save_summary (bool) [True] whether to save print summaries to text files.
	Returns: 
		ncbi_data_path: (str, path) path to the multi-sequence fasta_file
	"""
	ncbi_data_path = get_from_ncbi(gene_family, taxid, ncbi_file_name, working_directory, filtered_partial, filtered_predicted, silent)
	get_species_no_from_taxid(taxid)
	ncbi_data_summary(working_directory, ncbi_data_path, gene_family, taxid, filtered_partial, filtered_predicted, save_summary)
	return ncbi_data_path


def master_clustalo(fasta_path, iteration, thread_process_dict, silent = True):
	"""
	Parent function for clustalo alignment 
	Vars: 
		fasta_path: (str, path) path to the multi-sequence fasta file
		iteration: (int) [30] number of trees-iteratations for alignment method
		thread_process_dict (dict) contains the threads used for each method
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
	Returns: 
		alignment (str, path) path to the clustal alignment output file
		consensus (str, path) path to the consensus fasta file
	"""
	alignment = align_with_clustalo(fasta_path, iteration, thread_process_dict, silent)   
	consensus = get_consensus_from_alignment(alignment)
	return alignment, consensus 


def master_protein_analysis(fasta_path):
	"""
	Parent function to create a fasta_dict and populate it with sequence_ids, sequences and basic protein statistics
	Also creates a file containing in depth statitics about each sequence using pepstats. 
	Vars: 
		fasta_path: (str, path) path to the multi-sequence fasta file
	Returns: 
		fasta_dict: (dict) containg key (sequence_id), value (long name, sequence, basic stats) pairs on fasta interpretation results
	"""
	split_fasta = check_read_no_by_split(fasta_path, '>')
	fasta_dict = check_protein_basic_stats(split_fasta)
	protein_stats = get_protein_stats_from_pepstats(fasta_path)
	fasta_dict = add_pepstats_info_to_dict(fasta_dict, protein_stats)
	return fasta_dict    


def master_blast(working_directory, outputname, fastafile_path, fasta_dict, bin_no, thread_process_dict, silent):
	"""
	Parent function for blastp analysis. Option to blast any sequence (includig from multi-sequence fasta) or use consensus
	to give relatedness ranking. Creates bins a dictionary of bin, rank & sequence_ids
	Vars: 
		working_directory: (str, path) path to directory to store results in
		outputname: (str) prefix name for blastdb database output
		fastafile_path: (str, path) path to the multi-sequence fasta file
		fasta_dict: (dict) containg key (sequence_id), value (long name, sequence, basic stats) pairs on fasta interpretation results
		bin_no: (int) [250] the number of sequences to put in each bin 
		thread_process_dict (dict) contains the threads used for each method
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
	Returns:
		blast_output:(str, path) path to the blastp output file
		bins: (dict) of bin_number: sequence_id, rank
	"""
	blastdb_filepath = create_blastdb(working_directory, outputname, fastafile_path)
	query_filepath = query_choice(working_directory, fasta_dict)
	blast_output = query_blastdb(blastdb_filepath, query_filepath, thread_process_dict,	silent)
	bins = bin_results(blast_output, bin_no)
	return blast_output, bins


def master_motifs(analysispath, fasta_dict, bins, keep_fastas, save_summary):
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
			search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label), keep_fastas)
			print('Prosite scan completed for bin:'+str(selected_bin_label))
			print(' ')
	else: 
		selected_bin_label = selected_bin_labels
		if not os.path.isdir(os.path.expanduser(analysispath+'/motifs/bin_'+str(selected_bin_label))):
			os.mkdir(analysispath+'/motifs/bin_'+str(selected_bin_label))
		print('Creating individual fasta records for: bin_'+str(selected_bin_label))
		print('printing bins:')
		for k, v in bins.items():
			print(k)
			print(v)
			print()
		print(bins[selected_bin_label])
		print('printing fastadict')
		for k, v in fasta_dict.items():
			print(k)
			print(v)
			print()
		fasta_list_from_bin = get_selected_seq(selected_bin_label, bins[selected_bin_label], fasta_dict, analysispath+'/motifs/bin_'+str(selected_bin_label))
		for fasta in fasta_list_from_bin:
			create_new_fasta_record(fasta, fasta_dict[fasta]['sequence'], analysispath+'/motifs/bin_'+str(selected_bin_label))
		search_motifs_from_list(fasta_list_from_bin, selected_bin_label, analysispath+'/motifs/bin_'+str(selected_bin_label), keep_fastas)
		print('Prosite scan completed for bin:'+str(selected_bin_label))
		print(' ')
	for folder in os.listdir(os.path.expanduser(analysispath+'/motifs')):
		report_motif_results(analysispath+'/motifs/'+folder, save_summary)
	print('For more information on motifs per sequence, please see relevant files in folder:'+analysispath+'/motifs/')
	return fasta_dict 


def master_graphs(clustalo_files, fasta_path, thread_process_dict, blast_data_path, fasta_dict):
	"""
	Parent function to create some basic graphs. Return pandas database format output for blast & fasta 
	informtion to pass into further plotting routinues (to be added)
	Vars: 
		clustalo_files (str, path) path to the clustal alignment output file
		fasta_path: (str, path) path to the multi-sequence fasta file	
		thread_process_dict (dict) contains the threads used for each method
		blast_data_path:(str, path) path to the blastp output file
		fasta_dict: (dict) containg key (sequence_id), value (long name, sequence, basic stats) pairs on fasta interpretation results
	Output: 
		fastadb: [pandas dataframe] a dataframe of fasta_dict for plotting vars
		blastdb: [pandas dataframe] a dataframe of blast results for plotting vars
	"""
	create_hydropathy_plot(clustalo_files)
	create_plotcon(clustalo_files)
	matrix = similarity_matrix_clustalo(fasta_path, thread_process_dict)
	if matrix != False:
		plot_similarity_matrix_heatmap(matrix)
	blastdb = create_pd_from_blastp(blast_data_path)
	fastadb = create_pd_from_fasta_dict(fasta_dict)
	return blastdb, fastadb

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
	return input_string


def check_type(input_string, expected_type, parameter_name, parameter_default = None):
	"""
	Internal function: Checks the input against an expected type and either throws an error or returns default value 
	More complicated due to parsing input() information as this is always a string

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
						str(expected_type))
	if input_string == 'default':
		return parameter_default
	if expected_type == 'StrList':
		try:
			input_string=input_string.replace("'",'"')
			input_object = json.loads(input_string)
			if isinstance(input_object, list):
				return input_object
			elif isinstance(input_object, str):
				return input_object
			elif input_object.isnumeric == True:
				return str(input_object)
			else: 
				error = True
		except: 
			try:
				json.loads('"'+input_string+'"')
				return input_string
			except:
				error = True 
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
			elif expected_type == list and isinstance(input_object, list):
				return input_object
			else: 
				error = True
		except:
			error = True
	elif expected_type == bool: 
		if input_string == True or input_string == False:
			return input_string
		else:
			error = True
	elif expected_type == str and isinstance(input_string, str):
		return input_string
	if error == True and parameter_default == None: 
		raise ValueError(error_message)
	elif error == True and parameter_default != None: 
		print (error_message)
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
	while '__' in output:
		output = output.replace('__', '_')
	output = output.strip('_')
	return output 


def check_thread_format(thread_input, cpu, default_threads):
	"""
	Internal function for processing thread input and printing errors or warnings
	Vars:
		thread_input: [str None ] the (processed) input string for parsing into int or default
		cpu: [int] the available cpu threads
		default_threads [int] the default value to use for threading
	Returns: 
		default_threads: [int] either correct input integer or default value
	"""
	if thread_input == 'default' or thread_input == '' or thread_input == None:
		return default_threads
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
	as ints 
	"""    
	cpu = os.cpu_count()
	default_threads = int(cpu/3)
	return cpu, default_threads


def thread_entry(thread_processes):
	"""
	Internal function. Calculates the max cpu threads available and default values (an integer of 1/3rd available cpu)
	Lists processes which use threading and asks user to either set the thread number for all or individually and stores 
	output in a dictionary of thread_processes
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
			input_string = input('Please enter a valid integer to use as thread number to use for process: '+str(process)+' or press enter to use default')
			input_string = user_input_formating(input_string)
			threads = check_thread_format(input_string, cpu, default_threads)
			thread_process_dict[process] = threads
	else: 
		input_string = input('Please enter a valid integer to use as thread number to use for all processess or press enter to use default')
		input_string = user_input_formating(input_string)
		threads = check_thread_format(input_string, cpu, default_threads)
		for process in thread_processes: 
			thread_process_dict[process] = threads
	print('')
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
	print('')
	return folder 


def convert_filetypes_seqret(msf_infile, fasta_outfile):
	"""
	TODO: THIS ISNT USED 
	Converts a msf clustal alignment to fasta record version for processing differently
	Can be helpful for further analysis. 
	"""
	print('Attempting to convert msf alignment to fasta')
	run_nonpython_process("seqret -sequence \""+msf_infile+"\" -outseq  \""+fasta_outfile+"\"")
	print('msf file converted to fasta sucessfully. New file stored in: '+ str(fasta_outfile))
	print(' ')


### ---------- NCBI Data Collection ---------- ### 

hit_count = re.compile('<Count>([0-9]*)<')
species_fasta = re.compile('>.*\[(.*)\]')



def ncbi_string_processing(gene_family, taxid, filtered_partial = True, filtered_predicted = True):
	"""
	Internal function: Processes inputs to create an NCBI valid query string for data retrival 
	Vars:
		gene_family (str) a string of either a string or list (with appropriate escaped " or pairs of "') of gene names to search 
		taxid (str) a string of either a string or list (with appropriate escaped " or pairs of "') of taxonomic IDs to search in format txid0000
		filtered_partial (bool) [True] whether to add or exclude partial sequences from the fasta results
		filtered_predicted (bool) [True] whether to add or exclude predicted sequences from the fasta results
	Returns: 
		query_string (str) complete string query for searching ncbi
	"""
	query_string = ''
	if isinstance(taxid, list):
		join_string = '[Organism:exp] OR ' 
		entry_string = join_string.join(taxid)
		query_string = query_string + '('+entry_string+'[Organism:exp])'
	elif isinstance(taxid, str): 
		query_string = query_string+taxid+'[Organism:exp]'
		query_string = query_string + ' AND '
	if isinstance(gene_family, list):
		join_string = '[Gene Name] OR ' 
		entry_string = join_string.join(gene_family)
		query_string = query_string + '('+entry_string+'[Gene Name])'
	elif isinstance(gene_family, str):
		query_string = query_string+gene_family+'[Gene Name]'
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
	Vars:
		gene_family (str) a string of either a string or list (with appropriate escaped " or pairs of "') of gene names to search 
		taxid (str) a string of either a string or list (with appropriate escaped " or pairs of "') of taxonomic IDs to search in format txid0000
		file_name: (str) prefix name for seqs.fasta output
		working_directory: (str, path) path to directory to store results in
		filtered_partial (bool) [True] whether to add or exclude partial sequences from the fasta results
		filtered_predicted (bool) [True] whether to add or exclude predicted sequences from the fasta results
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
	Returns: 
		save_location: (str, path) path to the multi-sequence fasta_file
	"""
	query_string = ncbi_string_processing(gene_family, taxid, filtered_partial, filtered_predicted)
	print('Attempting to retrive data from ncbi for query')
	print('Query searched is '+str(query_string))
	save_location = working_directory+'/'+file_name+'seqs.fasta'
	edirect_query = "~/edirect/esearch -db protein -query \""+query_string+"\" | ~/edirect/efetch -format fasta" 
	edirect_response = run_nonpython_process(edirect_query)
	if silent != True:
		print(edirect_response)
	with open(save_location, 'w') as result_file:
		result_file.write(edirect_response)
	print('Data Retrival sucessful, fasta files stored in '+str(save_location))
	return save_location   


def ncbi_data_summary(working_directory, ncbi_data_path, gene_family, taxid, filtered_partial, filtered_predicted, save_summary = False):
	"""
	Prints a basic report to screen with the option to save to the working directory 
	Includes: 
		Number of sequences
		Number of genes
		Number of species 
		Filter options used 
	Vars:
		working_directory: (str, path) path to directory to store results in
		ncbi_data_path: (str, path) path to the multi-sequence fasta_file
		gene_family (str) a string of either a string or list (with appropriate escaped " or pairs of "') of gene names to search 
		taxid (str) a string of either a string or list (with appropriate escaped " or pairs of "') of taxonomic IDs to search in format txid0000
		filtered_partial (bool) [True] whether to add or exclude partial sequences from the fasta results
		filtered_predicted (bool) [True] whether to add or exclude predicted sequences from the fasta results
		save_summary (bool) [True] whether to save print summaries to text files.
	"""
	fasta_data = return_file_contents(ncbi_data_path)
	no_sequences = fasta_data.count('>')
	species = species_fasta.findall(fasta_data)  
	no_species = len(set(species))
	if isinstance(gene_family, list):
		no_genes = len(gene_family)
	else:
		no_genes = 1
	if isinstance(taxid, list):
		no_taxid = len(taxid)
	else:
		no_taxid = 1
	print('Analysis includes:'+str(no_sequences)+' protein sequences from '+str(no_species)+' species. Results for '+str(no_genes)+' gene(s) and '+str(no_taxid)+' taxid(s)')
	input_string = user_input_formating(input('Would you like to continue with the analysis?'))
	if input_string != True and input_string !='default':
		print('Processing Terminated: Interupt by user with input: '+str(input_string))
		sys.exit()
	if no_genes > 10000:
		input_string = input("WARNING: Number of sequences is greater than 10,000. Would you like to continue? This will make things slow")
		if input_string != True:
			print('Processing Terminated: Interupt by user')
			sys.exit()
	if save_summary == True:
		with open(working_directory+'/NCBI_retrival_summary.txt', 'w') as f:
			print('Number of Protein Sequences:'+str(no_sequences), file = f)
			print('Number of Unique Sdentified Species :'+str(no_species), file = f)
			print('Gene Number:'+str(no_genes), file = f) 
			print('Taxonomic Id(s) searched:'+str(gene_family), file = f)
			print('Gene(s) searched:'+str(taxid), file = f)
			print('Search options: Filtering Partial is '+str(filtered_partial)+', Filtering predicted is '+str(filtered_predicted), file = f)
			print('', file = f)
			print('Species included:'+str(list(set(species))), file = f)
		print('NCBI_retrival_summary.txt was saved in working directory')
	print(' ')


def get_species_no_from_taxid(taxid):
	"""
	Internal function to provide more information on number of sequences searched in specified taxid
	"""
	if isinstance(taxid, list):
		for gf in taxid:
			edirect_query = "~/edirect/esearch -db protein -query \""+gf+"[Subtree]"+"\""
			edirect_response = run_nonpython_process(edirect_query)
			counts = hit_count.search(edirect_response)
			counts = counts.group(1)
			print('INFO: Total number of subtree taxonomy groups searched in taxid:'+str(gf)+' is '+str(counts))
			edirect_query = "~/edirect/esearch -db protein -query \""+gf+"[Subtree] AND species[rank]"+"\""
			edirect_response = run_nonpython_process(edirect_query)
			counts = hit_count.search(edirect_response)
			counts = counts.group(1)
			print('INFO: Total number of subtree species searched in taxid:'+str(gf)+' is '+str(counts))
	elif isinstance(taxid, str):
		edirect_query = "~/edirect/esearch -db protein -query \""+taxid+"[Subtree]"+"\""
		edirect_response = run_nonpython_process(edirect_query)
		counts = hit_count.search(edirect_response)
		counts = counts.group(1)
		print('INFO: Total number of subtree taxonomy groups searched in taxid:'+str(taxid)+' is '+str(counts))
		edirect_query = "~/edirect/esearch -db protein -query \""+taxid+"[Subtree] AND species[rank]"+"\""
		edirect_response = run_nonpython_process(edirect_query)
		counts = hit_count.search(edirect_response)
		counts = counts.group(1)
		print('INFO: Total number of subtree species searched in taxid:'+str(taxid)+' is '+str(counts))	
	else:
		raise ValueError('Internal Error: Type formatting was incorrectly passed to get_species_no_from_taxid')	


### --------- Clustalo Alignment ---------  ####

def align_with_clustalo(fasta_path, iteration, thread_process_dict, silent = True):
	"""
	Uses the clustalo package through bash (through subprocess) to generate a multiple sequence alignment file
	Vars:
		fasta_path: (str, path) path to the multi-sequence fasta file
		iteration: (int) [30] number of trees-iteratations for alignment method
		thread_process_dict (dict) contains the threads used for each method
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
	Returns: 
		clustal_output (str, path) path to the clustal alignment output file
	"""
	print('Attempting Clustalso Alignment...')
	clustal_output = fasta_path.replace('seqs.fasta', 'clustal.msf')
	if not clustal_output.endswith('clustal.msf'):
		clustal_output = clustal_output + 'clustal.msf'
	threads = thread_process_dict['clustalo alignment']
	fasta_data = return_file_contents(fasta_path)
	if fasta_data.count('>') > 10000: 
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


def similarity_matrix_clustalo(fasta_path, thread_process_dict):
	"""
	Creates a pairwise distance matrix using clustalo aligment for use in heatmap plotting
	Requires user input to continue as large number of sequences will take a long time
	Vars: 
		fasta_path: (str, path) path to the multi-sequence fasta file
		thread_process_dict (dict) contains the threads used for each method
	Returns: 
		matrix_output: (str, path) path to matrix output file. 
	"""
	threads = thread_process_dict['similarity matrix']
	input_string = user_input_formating(input('WARNING: Pairwise distance matrix are very slow with large sequence numbers. Are you sure you want to continue?'))
	if input_string != True: 
		return False
	print('Attempting to create similarity matrix')
	matrix_output = fasta_path.replace('seqs.fasta', 'matrix.txt')
	align_output = fasta_path.replace('seqs.fasta', 'align.txt')
	if not matrix_output.endswith('matrix.txt'):
		matrix_output = matrix_output + 'matrix.txt'
	query = "clustalo -i \""+fasta_path+"\" --threads "+str(threads)+" --full -o \""+align_output+"\" -v --distmat-out=\""+matrix_output+"\""
	print(query)
	run_nonpython_process(query)
	print('Similarity matrix produced')	
	return matrix_output


def get_consensus_from_alignment(clustal_output):
	"""
	Uses EMBOSS Cons to get a consensus sequence from the Clustalo alignment files.
	Prints the consensus to screen.
	Vars: 
		clustal_output (str, path) path to the clustal alignment output file
	Returns: 
		output (str, path) path to the consensus fasta file
	"""
	print('Attempting to create a consensus sequence from Clustalo output')
	output = clustal_output.replace('clustal.msf', 'consensus.fa' )
	if not output.endswith('consensus.fa'):
		output = output + 'consensus.fa'
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
	Vars: 
		working_directory: (str, path) path to directory to store results in
		outputname: (str) prefix name for blastdb database output
		fastafile_path: (str, path) path to the multi-sequence fasta file
	Returns: 
		output (str, path) path to blastdb files 
	"""
	print('Attempting to create a blastdb from fasta files')
	output = working_directory+'/'+outputname+'blastdb'
	run_nonpython_process("makeblastdb -in \""+fastafile_path+"\" -dbtype prot -out \""+output+"\"")
	print('blast database creation sucessful. Blastdb file stored in: '+ str(output))
	print(' ')
	return output     


def query_choice(filepath, fasta_dict):
	"""
	Internal function to allow selection of valid sequence files to use as a query in the blastp analysis
	Gives list of valid file (based on extension) and user prompt to select. 
	If user selects multi-fasta file, allows user to select a single fasta from list of ids contained in multi-sequence fasta
	then will create a single-sequence fasta to use in blastp. 

	For measure of general similarity, using the consensus file created from cons in the clustalo sub-task is recomended.
	"""
	fastas = [fastafile for fastafile in os.listdir(os.path.expanduser(filepath)) if fastafile.endswith(".fasta") == True or fastafile.endswith(".fa") == True ]
	print('You have the following files with valid file extension for querying against blasdb:')
	count = 0
	for f in fastas:
		print(str(count) +':'+ f)
		count = count + 1
	selectedfile = user_input_formating(input('Please select the file you wish to use by entering the name or index. If your file is not present, enter your filename instead \\n Please note: Files must be present in the directory:'+str(filepath)))
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
		print('Selected file contains more than one sequence:' +str(len(reads)))
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
	Vars: 
		blastdb_filepath (str, path) path to custom blastdb 
		query_filepath (str, path) path to fasta file(s) to blast against db 
		thread_process_dict (dict) contains the threads used for each method
		silent: (bool) [True] whether to print outputs to screen (note: summary reports should still be printed regardless)
	Returns:
		output (str, path) pth to the blastp output file	
	"""
	print('Attempting to run blastp for selected file:'+ str(query_filepath)+' against blastdb:'+str(blastdb_filepath)+'with sequence:')
	print(return_file_contents(query_filepath))
	outputfile = query_filepath.replace('.fasta', '.fa').replace('.txt', '.fa')
	output = outputfile.replace('.fa','_blastp_output.txt')
	if not output.endswith('_blastp_output.txt'):
		output = output + '_blastp_output.txt'
	threads = str(thread_process_dict['blastp query'])
	edirect_reponse = run_nonpython_process("blastp -db \""+blastdb_filepath+"\" -query \""+query_filepath+"\" -num_threads \""+threads+"\" -max_target_seqs 10000 -outfmt 7 > \""+output+"\"")
	if silent != True:
		print(edirect_reponse)
	print('blastp sucessful. Output file stored in: '+ str(output))
	print(' ')
	return output 


def bin_results(input_file, bin_no, ignore_lines = '#', print_bin_contents = False):
	"""
	Function to split results into specified bins: 
	Vars: 
		input_file:(str, path) path to the blastp output file
		bin_no: (int) [250] the number of sequences to put in each bin 
		ignore_lines: (str) ['#'] marker for comment lines to ignore
		print_bin_contents: (bool), [False] useful for troubleshooting. Prints sequence_ids 
		contained in bin. 
	Returns: 
		result_bin (dict) containing bin_number as key with blast_output, rank as value
	"""
	print('Binning results from Blastp')
	result_bin = {}
	data = return_file_contents(input_file)
	lines = [line for line in data.split('\n')]
	validlines = [line for line in lines if line.startswith(ignore_lines) == False]
	print(len(validlines))
	validlines = list(filter(None, validlines))
	## As lists have defined order we can parse this using the first as the top hits
	totalno = len(validlines)
	print('Total results to bin:'+str(totalno))
	no_bins = math.ceil(totalno/bin_no)
	print('Bins expected:'+str(no_bins))
	item_no = 0
	rank_counter = 0 
	for i in range(0, no_bins):
		result_bin[i] = []
		for j in range(0,bin_no):
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
	"""
	Internal Function Creates a dictionary of blast results including rank and blast results
	"""
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
	"""
	Internal function to allow one or more bins to be selected for further analysis (motif search & plotting)
	will default to only using the top bin (0) if non valid user prompt entry.
	"""
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
				selected_bin = int(input_bins)
		else: 
			print('WARNING: Non integer bin:'+str(input_bins)+' selected. Returning to default value of: Bin 0')
			default_bin = True    
	if default_bin == True:
		selected_bin = 0        
	return selected_bin



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
	Vars: 
		fastafile_location: (str, path) path to the multi-sequence fasta file
	Returns: 
		output: (str, path) path to the pepstats output file
	"""
	print('Attempting to get protein statistics from pepstats')
	print(fastafile_location)
	output = fastafile_location.replace('seqs.fasta', 'stats.pepstats')
	if not output.endswith('stats.pepstats'):
		output = output + 'stats.pepstats'
	run_nonpython_process("pepstats -sequence \""+fastafile_location+"\" -outfile \""+output+"\"")
	print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
	print(' ')
	return output 


def check_protein_basic_stats(fasta_list):
	"""
	Internal function to create the fasta_dict and populate it with basic info 
	from the fasta file. 
	"""
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
	"""
	Internal function to populate fasta_dict with the basic results of the pepstats search
		-TODO can be extended to process more complex results from pepstats to plot protein 
			of intrest in more detail
	"""
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
	"""
	Checks if the prosite database exists in an accesible directory. If not downloads the required file
	and compiles them into a patmatmotifs searchable format. 
	"""
	if os.path.exists(prosite_db_directory):
		print('Prosite file \"prosite.lines\" was found. Continuing with analysis')
	else: 
		print('Prosite file \"prosite.lines\" was not found. Attempting to create by downloading .dat & .doc files from EMB=EBI')
		prosite_db_directory = get_prosite_db_scratch()
		prosite_db_directory = process_prositedb_with_prosextract(prosite_db_directory)


def search_for_motifs(fasta_file_name, directory, silent = True): 
	"""
	Uses prosite to search motifs in a given fasta file 
		fasta_file_name (str) the name of the fasta file without .fasta extension
		directory: (str, path) path to the individual-sequence fasta files storage location
		silent: (bool) [True] whether to print outputs to screen 
	"""
	f = directory+'/'+fasta_file_name+'.fasta'
	fout = f.replace('.fasta', '.motif')
	run_nonpython_process("patmatmotifs -sequence "+f+" -auto Y -outfile "+fout)
	if silent != True:
		print('Prosite search sucessful for sequence:'+str(fasta_file_name)+' output stored in directory:')
		print(' ')
		   

def get_selected_seq(blast_bin_no, blast_bin, fasta_dict, directory):
	"""
	Internal function. Checks the blast_bin results for sequence_id and compares this
	with the fasta_dict. If the sequence_id is found, it create an individual fasta record in
	that directory from the sequence in fasta_dict[sequence_id] as individual fasta required for motif searching
	"""
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
	Vars:
		fasta_list: (list) list of fasta_ids to use in motif searching
		group_name: (str int) name of bin group to use in motif searching 
		directory:  (str path) name of directory contaning fasta sequences
		keep_fastas:  (bool) [False] whether to store the single sequence fastas for each fasta_id in fasta_list
	"""
	print('Attempting to search for motifs using prosite in sequences: bin_'+str(group_name))
	for fasta_file_name in fasta_list:
		search_for_motifs(fasta_file_name, directory)
	if keep_fastas == False: 
		for fasta_file_name in fasta_list: 
			f = directory+'/'+fasta_file_name+'.fasta'
			os.remove(f)


def create_new_fasta_record(name, sequence, directory):
	"""
	Creates a new fasta file for the specified identifier and sequence and writes to given directory
	Vars:
		Name: (str) unique string id
		Sequence: (str) protein sequence 
		directory: Directory to store fasta results in (motifs/bin_x)
	"""
	with open(directory+'/'+str(name)+'.fasta', 'w+')  as f:
		f.write('>'+str(name)+sequence)


find_motif_no = re.compile('# HitCount:(.*)\n')
find_motif_id = re.compile('Motif = (.*)\n')


def report_motif_results(directory, save_summary):
	"""
	Prints to screen a summary report on the motifs found from PROSITE data in the given directory
	Does not include reports where no motifs were found
	eg [protein_id] has [INT] motifs found with prosite name [PROSITE ENTRY NAME] 
	Vars:
		directory: (str) of where the motif files are stored
		save_summary: (bool) wheather to save the summary report to motif_summary.txt file
	"""

	motif_files = [mfile for mfile in os.listdir(directory) if mfile.endswith(".motif")]
	motif_dict = {}
	for mfile in motif_files:
		with open(os.path.expanduser(directory)+'/'+mfile) as f:
			data = f.read()
			motif_no = find_motif_no.search(data).group(1)
			motifs = find_motif_id.findall(data)
			if int(motif_no) > 0:
				motif_dict[mfile] = {'no_motifs': motif_no, 'motifs': motifs}
	print('Files Analysed:'+str(len(motif_files))+' motifs found in: '+str(len(motif_files)))
	for k in motif_dict.keys():
		print('File:'+str(k)+': '+str(motif_dict[k]['no_motifs'])+' motifs found with Prosite motif entry name(s): '+str(motif_dict[k]['motifs']))      
	if save_summary == True: 
		with open(directory+'/motif_summary.txt', 'w') as f:
			for k in motif_dict.keys():
				print('File:'+str(k)+': '+str(motif_dict[k]['no_motifs'])+' motifs found with Prosite motif entry name(s): '+str(motif_dict[k]['motifs']), file = f)
		print('motif_summary.txt was saved in motif bin directory')
	print('For more information, please see relevant .motif file in directory:'+str(directory))      
	print(' ')


### ---------- Plotting Functions ---------- ### 

def create_hydropathy_plot(alignment_file):
	"""
	Uses the EMBOSS package pepwindowall to create a hydropathy_alignment 
	from the input alignment file path
	Vars: 
		alignment_file (str, path) path to the clustal alignment output file
	"""
	print('Attempting to create Kyte-Doolittle Hydropathy plot for alignment')
	output = alignment_file.replace('clustal.msf', 'hydropathy_alignment')
	if not output.endswith('hydropathy_alignment'):
		output = output + 'hydropathy_alignment'
	run_nonpython_process("pepwindowall \""+alignment_file+"\" -graph svg -gxtitle \'Residue\' -gtitle2 \"Kyte-Doolittle Hydropathy for Alignment\" -goutfile \""+output+"\"")
	print('Protein statistics creation sucessful. Protein Statistics file stored in: '+ str(output))
	print(' ')


def create_plotcon(alignment_file):
	"""
	Uses emboss plotcon to make a similarity plot from the clustalo alignment file. 
	Vars: 
		alignment_file (str, path) path to the clustal alignment output file
	"""
	print('Attempting to create a similarity plot from clustalo alignment file')
	output = alignment_file.replace('_clustal.msf', '_similarityplot')
	run_nonpython_process("plotcon -sequences \""+alignment_file+"\" -graph svg -goutfile \""+output+"\"")
	print('Similarity Plot created sucessfully. Plotcon output file stored as: '+ str(output)+'.svg')
	print(' ')


def plot_similarity_matrix_heatmap(matrix_filepath):
	"""
	Creates a heatmap from the similarity matrix calculated by Clustalo
	Vars: 
		matrix_filepath (str, path) the path to clustalo output similarity matrix file 
	"""
	data = return_file_contents(matrix_filepath)
	labels = []
	matrix = []
	for record_index, record_data in enumerate(filter(None, data.split('\n'))):
		if record_index != 0 and record_data:
			matrix.append([])
			for column_index, column_data in enumerate(filter(None, record_data.split(' '))):
				if column_index == 0:
					labels.append(column_data)
				else: 
					matrix[record_index-1].append(float(column_data))
	df = pd.DataFrame(matrix, columns = labels, index = labels)
	ax = sns.heatmap(df)
	plt.title('Heatmap of Similarity Matrix for sequence_id')
	plt.tight_layout()
	plt.savefig('similarity_matrix_heatmap.png')
	plt.show()


def create_pd_from_blastp(inputfile, ignored = '#'):
	"""
	Internal function, creates a blast dataframe using expected headers, split on \t from a blastp_outputfile
	Ignores comment lines
	Vars: 
		inputfile:(str, path) path to the blastp output file
	Returns: 
		blastdb: [pandas dataframe] a dataframe of blast results for plotting vars
	"""
	headers=['query_acc.', 'subject_acc.', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
	blastpd = pd.read_csv(inputfile, comment = ignored, sep ='\t', names = headers)
	return blastpd


def create_pd_from_fasta_dict(fasta_dict):
	"""
	Simple dataframe creator from fasta
	Vars: 
		fasta_dict: (dict) containg key (sequence_id), value (long name, sequence, basic stats) pairs on fasta interpretation results
	Output: 
		fastadb: [pandas dataframe] a dataframe of fasta_dict for plotting vars
	"""
	df = pd.DataFrame.from_dict(fasta_dict, orient = 'index')
	return df

def create_pepinfo_plots(fasta_file_path, fasta_file_label):
	"""
	Creates pepinfo plots (hydropathy and intrest units)
	TODO: NOT Currently used. 
	TODO: Not Currently working
	Vars: 
		fasta_file_path (str, path) a path to a single-sequence fasta file
		fasta_file_label (str) the label id for fasta file
	"""
	print('Attempting to create pep-info plot for selected fasta file'+str(fasta_file_label))
	## Makes graph in file called pepinfo.ps 
	run_nonpython_process("pepinfo \""+fasta_file_path+"\" -graph cps -gtitle \""+fasta_file_label+"\"")


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



def create_catagorical_line_plot(x, y, title, fasta_df, motif_bin_filepath):
	bin_paths = [f.path for f in os.scandir(folder) if f.is_dir()]
	bins_ids = {}
	count = 0
	for folder in bin_paths: 
		bin_ids[count] = []
		for file_name in os.listdir(folder):
			if file_name.replace('.fasta', '').replace('.motif') not in bin_ids[count]:
				bin_ids[count].append([])
	if len(bin_paths) == 1:
		plt.plot(subdf.index, subdf['Mol_Weight'], 'r')
		plt.plot(subdf.index, subdf['Charge'], 'b')
		plt.plot(subdf.index, subdf['length'], 'p')		
	else: 
		fig, ax = plt.subplots(1,len(bin_paths))
		ax
	
def subdf_for_bin(id_list, df):
	
	df = df[df.index.isin(id_list)]
	return df
['CAP10123.1', 'CAP10110.1', 'BAE72078.1', 'BAE46847.1']


"""

"""
Test cases: 
glucose-6-phosphatase proteins from Aves (G6pc txid8782)/(G6pc,COX1 txid8782)


ABC transporters in mammals (txid40674)

kinases in rodents

adenyl cyclases in vertebrates

"""

