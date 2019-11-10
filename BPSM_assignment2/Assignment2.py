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

 - need to add flexibility to find txid and gene name from inputs....
get_from_ncbi(G6pc, txid8782, 'testing')
"""

### ---------- Import Statements ---------- ###
import os
from pathlib import Path 
import edirect 
import re
import subprocess



### ---------- Master Function ---------- ###

def run_protein_ident(genefamily, taxid, filtered_partial = True, filtered_predicted = True, filepath= None, foldername=None)
    """    
    1) Before starting any analysis we check the inputs are of the correct type. If they are not we break or return to default.
    2) Check query is not > hundreds of lines long & replace any spaces with + 
    3) Create a folder for the output to be saved to 
    
    """
    check_type(genefamily, str, default = None)
    check_type(taxid, str, default = None)
    genefamily = replace_non_alphanumeric_chars(genefamily)
    taxid = replace_non_alphanumeric_chars(taxid)
    query_filename = genefamily+'_'+taxid
    filepath = check_type(filepath, str, Path().absolute())
    filepath = check_type(foldername, str, query_filename)
    
    




### ---------- Utility Functions ---------- ###

def check_type(input, expected, default, firm = True):
	"""
	Checks the input against an expected type and either throws an error
	if firm != True will allow continuation and return 'default'
	"""
	if isinstance(input, expected) == True:
		return 'pass'
	else:
		realtype = type(input)
		if default == None:
			raise ValueError('The input: '+str(input)+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
		else: 
			print('WARNING: The input: '+str(input)+' was not of the expected type: '+
						str(expected)+ ' instead was of type: '+str(realtype))
			print('WARNING: Attempting to return to default value if applicable')
			return default 

def replace_non_alphanumeric_chars(input_string, replacement_char = '_'): 
    output = re.sub('[^0-9a-zA-Z]+', replacement_char, input_string)
    return output 

def run_nonpython_process(query, timeout = 6000):
    try: 
        response = subprocess.check_output(query, shell = True)
        #process = subprocess.Popen(query, stdout = subprocess.PIPE , stderr = subprocess.PIPE, shell = True)
        #response, error = process.communicate(timeout = timeout)
        response = response.decode('ascii')
        #error = error.decode('ascii')
        #print(error)
        return response 
    except Exception as e: 
        print('Error running edirect: '+str(e))
    
def check_edirect_installation():
    """
    check that edirect is installed or install it with: 
    """
    print ('Checking edirect is installed in home space...')
    if exists ~/edirect:
        print('Edirect is already installed, continuing with analysis')
    else: 
        print('Edirect is not currently installed, trying to install now')
        run_nonpython_process('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"')
        print('Edirect has been sucessfully installed, continuing with analysis')



### ---------- File outputing & File locations ---------- ### 

def create_folder_path(name, filepath):
    """
    Checks if the inputted folder (name  + prefix + filepath) exists 
    if not it will create it. Both cases will result in a print statement 
    and return a string of the folder path to store results in.

    name: If no user defined name, it will default to family_taxid
    filepath: string, if no user defined input it will default to current working directory 
    """
    prefix = 'Protein_Analysis_'
    folder = filepath+'/'+prefix+name
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
"""
def get_from_ncbi(family, taxid, query_filename, filepath, filtered_partial = True, filtered_predicted = True):
    edirect_query_task = taxid+"[Organism:exp] AND "+family+"[Gene Name]"
    if filtered_partial == True:
        query_search_term = edirect_query_task+' NOT partial[All Fields]'
    if filtered_predicted == True:
        query_search_term = query_search_term+' NOT predicted[All Fields]'
    save_location = filepath+'/'+query_filename
    edirect_query = "~/edirect/esearch -db protein -query query_search_term | ~/edirect/efetch -format fasta" 
    edirect_reponse = run_nonpython_process(edirect_query)
    print(edirect_reponse)
    with open(save_location, 'w') as result_file:
        result_file.write(edirect_reponse)






### ---------- Plotting Routines ---------- ### 

### ---------- Prosite Data Collection ---------- ### 

