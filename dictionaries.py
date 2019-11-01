ctionaries Lecture Notes 


# Interrogating the users

def check_inputs(question, inputvalue):
	check_dict={'name': str,
				'age': int,
				'colour': str,
				'python_preference': bool,
				'sanity_check': bool}
				
	if check_dict[question] == int:
		try:
			int(inputvalue)
		except ValueError:
			print ('Entered value type for '+str(question)+' does not match expected type:'+str(check_dict[question]))
	elif check_dict[question] == bool:
		if inputvalue not in ['True', 'False', 'Yes', 'No']:
			print (inputvalue)
			print (type(inputvalue))
			raise ValueError('Entered value type for '+str(question)+' does not match expected type:'+str(check_dict[question]))
	else:
		if isinstance(inputvalue, check_dict[question]) == False:
			raise ValueError('Entered value type for '+str(question)+' does not match expected type:'+str(check_dict[question]))
		

def interrogation():
	user_info_dict = {}
	user_info_dict['name'] = input('What is your name?')
	check_inputs('name',user_info_dict['name'])
	user_info_dict['age'] = input('How old are you?')
	check_inputs('age',user_info_dict['age'])
	user_info_dict['age'] = int(user_info_dict['age'])
	if user_info_dict['age'] >= 100:
		print('Really? I don\'t think you are...')
	user_info_dict['colour'] = input('What is your favorite colour?')
	check_inputs('colour',user_info_dict['colour'])
	user_info_dict['python_preference'] = input('Do you like Python?')
	print(user_info_dict['python_preference'].lower())
	if user_info_dict['python_preference'].lower() == 'no':
		print('Your crazy...')
	user_info_dict['sanity_check'] = input('The world is flat: True or False')
	check_inputs('sanity_check',user_info_dict['sanity_check'])
	if user_info_dict['sanity_check'] == True:
		print('Oh no a flat-earther was discovered!')
		

## DNA translation

def check_translation_dict(translation_dict):
	lengths = []
	for c in translation_dict.keys():
		lengths.append(len(translation_dict[c]))
	if len(set(lengths)) > 1:
		raise ValueError('Translation dictionary has differing lengths for different codons. This is not compatible with method')
	else: 
		return lengths[0]

def translation_dictionary(coding_dict, codon):
	gencode = {	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
					'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
					'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
					'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
					'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
					'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
					'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
					'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
					'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
					'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
					'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
					'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
					'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
					'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
					'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
					'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
	if coding_dict != 'default':
		if isinstance(coding_dict, dict) == False:
			raise ValueError('Coding_dict variable is not a valid dictionary')
	else: 
		coding_dict = gencode
	
	
	trans_length = check_translation_dict(coding_dict)
	if codon not in coding_dict.keys():
		codonvalue = entry_codons(codon, trans_length)
		coding_dict[codon] = codonvalue
	return coding_dict[codon]

def input_quit_check(query):
	if query.lower() == 'q':
		raise Exception('User interupt. No values added & translation aborted')
	else:
		return query
		
	
def entry_codons(codon, trans_length):
	query = input_quit_check(input('Codon: '+codon+' is not in translation dictionary. Enter a translation or press \'Q\' to exit'))
	if len(query) != trans_length:
		query = input_quit_check(input('WARNING: Unequal translation code between user entry and existing dictionary. Please enter a translation with length: '+str(trans_length)+' or press \'Q\' to exit'))
	if len(query) != trans_length:
		raise ValueError('Incorrect length translation code entered for the second time. Aborting Translation.')
	return query 
		
		
		
		
def translator(DNAseq, coding_dict = 'default', startframe = 'all', direction = 'both', n = '?'):
	if isinstance(DNAseq, str) == False:
		raise ValueError('DNAseq value must be a string. Aborting translation')
	product_dict = {}
	if startframe == 'all':
		start = [1, 2, 3]
	else:
		if isinstance(startframe, int) == True:
			if startframe in [1,2,3]: 
				start = [startframe]
			else: 
				raise ValueError('Non default value for startframe is either not an integer or not in the range 1-3')
		else: 
			raise ValueError('Non default value for startframe is either not an integer or not in the range 1-3')
	
	for s in start:
		name = s
		product_dict[str(name) + '_F'] = {'codons':[], 'translation': []}
		end = 3 + s-1
		for i in range(int(len(DNAseq)/3)):		
			codon = DNAseq[s-1:end]
			print(codon)
			aa = translation_dictionary(coding_dict,n)	
			s = s +3
			end = end + 3
			product_dict[str(name)+'_F']['codons'].append(codon)
			product_dict[str(name)+'_F']['translation'] = product_dict[str(name) + '_F']['translation']+aa
	for k in product_dict.keys():
		print(k)
		print(product_dict[k])
	
		

	
