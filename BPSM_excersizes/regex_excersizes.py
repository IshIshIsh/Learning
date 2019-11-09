ccession Numbers ###
accs = ['xkn59438', 'yhdck2', 'eihd39d9', 'chdsye847', 'hedle3455', 'xjhd53e', '45da', 'de37dp']

def sort_accs(accs):
	fives = []
	d_or_e = []
	de = []
	d_e = []
	d_and_e = []
	start_x_or_y = []
	start_x_or_y_end_e = []
	contains_3numbers = []
	ends_in_d_arp = []
	for a in accs: 
		a=a.lower()
		if re.search(r'5', a): 
			fives.append(a) 
		if re.search(r'[de]', a):
			d_or_e.append(a)
		if re.search(r'de', a): 
			de.append(a)
		if re.search(r'd.{1}e', a): 
			d_e.append(a)
		if re.search(r'd', a) and re.search(r'e', a): 
			d_and_e.append(a)
		if re.search(r'^x|^y', a):
			start_x_or_y.append(a)
		if re.search(r'^x|^y', a): and re.search(r'e$', a)
			start_x_or_y_end_e.append(a)
		if re.search(r'[1-9]{3,3}',a):
			contains_3numbers.append(a)
		if re.search(r'da$|dr$|dp$', a):
			ends_in_d_arp.append(a)
		
	print('The following accession numbers contains the number 5:')
	print(fives)
	print('The following accession numbers contain the letter d or e:')
	print(d_or_e)
	print('The following accession numbers contain the letters d and e in that order:')
	print(de)
	print('The following accession numbers contain the letters d and e in that order with a single letter between them:')
	print(d_e)
	print('The following accession numbers contain both the letters d and e in any order:')
	print(d_and_e)
	print('The following accession numbers start with x or y:')
	print(start_x_or_y)
	print('The following accession numbers start with x or y and end with e:')
	print(start_x_or_y_end_e)
	print('The following accession numbers contain three or more numbers in a row:')
	print(contains_3numbers)
	print('The following accession numbers end with d followed by either a, r or p:')
	print(ends_in_d_arp)
	
	
	
## Restriction Digests ### 

file = '//localdisk//home//ifarquha//BPSMLectureExcercises//Python//Regex//long_dna.txt'
site1 = 'ANT/AAT'
site2 = 'GCRW/TG'

code_dict ={'N': '[GCTA]',
			'R': '[AG]',
			'Y': '[CT]',
			'S': '[GC]',
			'W': '[AT]',
			'K': '[GT]',
			'M': '[AC]',
			'B': '[CGT]',
			'D': '[AGT]',
			'H': '[ACG]'}

def convert_sites(site, code):
	for k in code.keys():
		if k in site: 
			site = site.replace(k, code[k])
		else: 
			pass
	return site 


	
def digest(dnafile, site, code=False):
	site_dig_lens = []
	positions = {}
	site_dict = {}
	if code!= False: 
		site = convert_sites(site, code)		
	with open(dnafile, 'r') as file:
		dna = file.read()
		dna = dna.upper()
		print('Total DNA length: '+str(len(dna)))
		partialsites = site.split('/')
		for i in partialsites:
			site_dig_lens.append(len(i)
		site = site.replace('/', '')
		print(site)
		found_sites = re.findall(site, dna)
		if found_sites:
			print('Sites found: '+ str(found_sites))
		else: 
			print('')
		for fs in found_sites:
			position = re.search(site, dna)
			print('site: '+str(fs))
			print('starts at '+str(position.start())+' and ends at '+str(position.end()))
			positions.append(position)
			
