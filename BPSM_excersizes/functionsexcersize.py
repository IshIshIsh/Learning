
ns Lecture: 

# Amino Acid %

def amino_acid_percent(protein_seq, amino_acid_code): 
	amino_acid_code = amino_acid_code.upper()
	protein_seq = protein_seq.upper()
	amino_acid_count = protein_seq.count(amino_acid_code)
	protein_length = len(protein_seq)
	amino_acid_percent = amino_acid_count/protein_length * 100
	print (amino_acid_percent)
	return amino_acid_percent
	
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", "M")) == round(5)
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", "r")) == round(10)
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", "L")) == round(50)
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", "Y")) == round(0)

# Amino Acid % part 2 

def amino_acid_percent(protein_seq, amino_acid_codes = []): 
	amino_acid_dict = {}
	protein_seq = protein_seq.upper()
	total_percent = 0
	if isinstance(amino_acid_codes, list) == True: 
		if amino_acid_codes == []:
			amino_acid_codes = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
		else:
			pass
		for aa in amino_acid_codes: 
			aa = aa.upper()
			amino_acid_count = protein_seq.count(aa)
			protein_length = len(protein_seq)
			amino_acid_percent = amino_acid_count/protein_length * 100	
			total_percent = total_percent + amino_acid_percent
			amino_acid_dict[aa] = amino_acid_percent
		print ('Total percentage of protein sequence for amino acids given: '+str(total_percent))
		for aa in amino_acid_dict.keys(): 
			print (aa + ' count is '+ str(amino_acid_dict[aa]))
	else: 
		raise ValueError('amino_acid_codes variable was not entered as a list or left blank for default')
	return total_percent

assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", ["M"])) == 5
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L'])) == 70
assert round(amino_acid_percent("MSRSLLLRFLLFLLLLPPLP")) == 65

# Amino Acid % part 2 Contracted: 
def amino_acid_percent(protein_seq, amino_acid_codes = []): 
	amino_acid_dict = {}
	total_percent = 0
	if isinstance(amino_acid_codes, list) == True: 
		if amino_acid_codes == []:
			amino_acid_codes = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
		for aa in amino_acid_codes: 
			amino_acid_count = protein_seq.upper().count(aa.upper())
			amino_acid_dict[aa] = amino_acid_count/len(protein_seq) * 100	
			total_percent = total_percent + amino_acid_dict[aa]
		print ('Total percentage of protein sequence for amino acids given: '+str(total_percent))
		for aa in amino_acid_dict.keys(): 
			print (aa + ' count is '+ str(amino_acid_dict[aa]))
	else: 
		raise ValueError('amino_acid_codes variable was not entered as a list or left blank for default')
	return total_percent
	
# Base counter 

Write a function that will take a DNA sequence along with an optional threshold 
and return True or False to indicate whether the DNA sequence contains a high 
proportion of undetermined bases (i.e not A, T, G or C).

def base_counter(DNAseq, threshold, ignoredbases = ['A', 'T', 'G', 'C']): 
	base_count_dict = {}
	total_count = 0
	for charecter in DNAseq.upper(): 
		if charecter not in ignoredbases and charecter not in base_count_dict.keys():
			count = DNAseq.upper().count(charecter.upper())
			percentage = count/len(DNAseq) * 100
			base_count_dict[charecter] = {'percentage':percentage, 'boolean':threshold <= percentage}
			total_count = total_count + count
	total_percent = total_count/len(DNAseq) * 100
	if total_percent <= threshold:
		statement = '. Sequence contains a low proportion based on threshold'
	else: 
		statement = '. Sequence contains a high proportion based on threshold'
	print('Total of undetermined bases is '+str(round(total_percent))+'%'+statement)
	for charecter in base_count_dict:
		if base_count_dict[charecter]['boolean'] == False:
			print('Sequemce contains a low proportion '+str(charecter)+'('+str(round(base_count_dict[charecter]['percentage']))+'%) based on threshold')
		else: 
			print('Sequemce contains a high proportion '+str(charecter)+'('+str(round(base_count_dict[charecter]['percentage']))+'%) based on threshold')
	return threshold <= total_percent
	
assert (base_counter('ABBGCT', 50)) == False
assert (base_counter('ABBBBGCT', 50)) == True
assert (base_counter('ABBfFGC', 50)) == True
assert (base_counter('ABFhGCt', 50)) == False
assert (base_counter('ABFhGCt1', 50)) == True

## Kmer counting ## 
Write a function that, given any DNA sequence, will print all the k-mers (e.g. 4-mers) that occur more than n times.

E.g. with dna="ATGCATCATG", kmersize=2 and minfrequency=2 will output
AT
because the kmers (kmersize) are 2 bases long, and there are three (minfrequency was 2) instances of "AT"

def kmer_counting(dnaseq, kmersize=2, minfrequency=2):
	kmer_dict = {}
	kmer_list = []
	if kmersize > dnaseq or isinstance(kmersize, int) == False:
		print('WARNING: kmersize inappropriate resetting to default of 2')
		kmersize = 2
	if isinstance(minfrequency, int) == False: 
		print('WARNING: minfrequency inappropriate resetting to default of 2')
		kmersize = 2		
	for i in range(len(dnaseq)):
		kmer =(dnaseq[i:i+kmersize])
		if len(kmer) < kmersize:
			pass 
		else:
			if dnaseq.count(kmer) > 2:
				kmer_dict[kmer] = dnaseq.count(kmer) 
	for kmer in kmer_dict:
		kmer_list.append(str(kmer)+':'+str(kmer_dict[kmer]))
		print('kmer '+str(kmer)+' occurance: '+str(kmer_dict[kmer]))
	if kmer_list != []
		print('List of kmers found in sequence: '+(str(kmer_list)))
	else:
		print('No kmers found')

		
kmer_counting('ATGCATCATG')
kmer_counting('AATTGGCCAATTGCGCAATGGCGGC', 3, 2)
	



