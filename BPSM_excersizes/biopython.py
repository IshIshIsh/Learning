eq import Seq
from Bio import SeqIO
from Bio import Entrez


Entrez.email = "s123456@ed.ac.uk"
Entrez.api_key ="f8230825b624062e53bb5ab709b586364c07"

## How many complete COX1 protein records are there for mammals

result_cox1 = Entrez.read(Entrez.esearch(db='protein', term = 'COX1[gene] AND 40674[taxid] AND NOT partial[All Fields]' ))
number = result_cox['Count']

## What is their average length (the proteins that is, not the mammals!)?


## Write a function that will answer the question for any gene name and any taxonomic group.

def find_gene_info(genename, searchtaxid, partial=False): 
	if partial == False: 
		result = Entrez.read(Entrez.esearch(db='protein', term = 'genename[gene] AND searchtaxid, [taxid] AND NOT partial[All Fields]' ))
	else: 
		result = Entrez.read(Entrez.esearch(db='protein', term = 'genename[gene] AND searchtaxid, [taxid] AND NOT partial[All Fields]' ))
	number = result_cox['Count']
	print('Number of search results for:'+str(genename)+','+str(searchtaxid)+' is: '+str(number))

## What else could you sensibly include in the function (in terms of retrieving possibly useful data)?
