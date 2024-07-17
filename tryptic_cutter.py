import search

import re

# komentar #4
# tryptic_cutter treba imati aminokiselinsku sekvencu kao ulaz i niz sekvenci
# kao izlaz

def tryptic_cutter (accession_number):
	tryptic_list=[]
	tryptic_list_final=[]
	protein=""
	sequence = search.canonical_sequence(accession_number)
	protein = re.sub(r'(?<=[RK])(?=[^P])','\n', sequence, re.DOTALL)
	tryptic_list = protein.split()
	return (tryptic_list)
	


