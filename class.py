import os

import xml.etree.ElementTree as et

import re

base_path = os.path.dirname(os.path.realpath("search.py"))
xml_file = os.path.join(base_path, "gli.xml")
tree = et.parse(xml_file)
root = tree.getroot()

ns={"uniprot":"http://uniprot.org/uniprot"}

class Protein_searcher:				

	def __init__(self, accession_number):
		
		self.accession_number = accession_number
	
	def protein_name (self):
		
		name=""
		
		for entry in root.findall("uniprot:entry",ns):
			if self.accession_number == entry.find("uniprot:accession",ns).text :
				name = entry.find("uniprot:name",ns).text
				return (name)
	
	
	def canonical_sequence (self):
		
		sequence=""
		
		for entry in root.findall("uniprot:entry",ns):
			if self.accession_number == entry.find("uniprot:accession",ns).text :
				sequence = entry.find("uniprot:sequence",ns).text
				sequence = ''.join(sequence.split())
				return (sequence)
				
	def glycosilation_sites (self):
		
		glycosilation_list = []
		location = ""
		glycosilation_site = {}
		
		for entry in root.findall("uniprot:entry",ns):
			if self.accession_number == entry.find("uniprot:accession",ns).text :
				for feature in entry.findall("uniprot:feature[@type='glycosylation site']" ,ns ):
					location = feature.find("uniprot:location",ns)
					glycosilation_site = location.find("uniprot:position",ns).attrib
					glycosilation_site = int(glycosilation_site["position"])
					glycosilation_list.append(glycosilation_site)
		
		return glycosilation_list
		

	def alternative_seq_ids_dictionary (self):
		
		isoform_number=""
		ref_id={}
		iso_dict={}
		
		for entry in root.findall("uniprot:entry",ns):
			if self.accession_number == entry.find("uniprot:accession",ns).text :
				for comment in entry.findall("uniprot:comment[@type='alternative products']",ns):
					for isoform in comment.findall("uniprot:isoform",ns):
						isoform_number=isoform.find("uniprot:name",ns).text
						ref_id = isoform.find("uniprot:sequence",ns).attrib
						if "ref" in ref_id:
							iso_dict[isoform_number]= ref_id["ref"].split()    
		
		return(iso_dict)

		
	def isoform_data_dictionary (self):
		
		ids_dict={}
		ids = ""
		original = ""
		variation = ""
		index_st = 0
		index_end = 0
		
		for entry in root.findall("uniprot:entry",ns):
			if self.accession_number == entry.find("uniprot:accession",ns).text :	
				for feature in entry.findall("uniprot:feature[@type='splice variant']",ns):
						feat = feature.attrib
						ids = feat["id"]
						original=feature.find("uniprot:original",ns)
						
						if original is None:
							original = ""
						else:
							original = feature.find("uniprot:original",ns).text
						
						variation = feature.find("uniprot:variation",ns)
						
						if variation is None:
							variation = ""
						else:
							variation = feature.find("uniprot:variation",ns).text
						
						location = feature.find("uniprot:location",ns)
						
						index_st = location.find("uniprot:begin",ns)
						
						if index_st is None:
							index_st = location.find("uniprot:position",ns).attrib
							index_st = int(index_st["position"])
						else:
							index_st = location.find("uniprot:begin",ns).attrib
							index_st = int(index_st["position"])

						index_end=location.find("uniprot:end",ns)
						
						if index_end is None:
							index_end = location.find("uniprot:position",ns).attrib
							index_end = int(index_end["position"])
						else:
							index_end = location.find("uniprot:end",ns).attrib
							index_end = int(index_end["position"])

						ids_dict[ids] = (original, variation, index_st, index_end)
		
		return(ids_dict)
		
	
def isoform_creator (accession_number):

	seqt = Protein_searcher.canonical_sequence (self)
	ids_dict = Protein_searcher.isoform_data_dictionary (self)
	iso_dict = Protein_searcher.alternative_seq_ids_dictionary (self)
	
	isoforms_dict = dict()
	original=""
	variation=""
	index_st=0
	index_end=0
	ref_id=[]
	seq_temp=""
	tmp_index=0
	idovi=[]
	
	idovi=ids_dict.keys()
	ref_id=iso_dict.keys()
	print(ref_id)
	
	for ref_id,idovi  in iso_dict.items():
		seq_temp = seqt
		tmp_index = 0
		for ids in idovi:
			original, variation, index_st, index_end = ids_dict[ids]
			seq_temp = '%s%s%s'%(seqt[:index_st-1 + tmp_index],variation,seqt[index_end+tmp_index:])
			if len(original) == index_end - index_st:
				tmp_index += len(variation)-len(original)
			else:
				tmp_index += index_end - index_st

		isoforms_dict[ref_id] = seq_temp	
		
	return(isoforms_dict)
	
def tryptic_cutter (accession_number):

	tryptic_list=[]	
	tryptic_list_final=[]
	protein=""
	sequence = Protein_searcher.canonical_sequence (self)
	
	protein = re.sub(r'(?<=[RK])(?=[^P])','\n', sequence, re.DOTALL)
	tryptic_list = protein.split()
	
	return (tryptic_list)	
	
	def mass_calculator_canonical (self):
		
		sequence = Protein_searcher.canonical_sequence (self)
		
		a = sequence.count("A")
		r = sequence.count("R")
		n = sequence.count("N")	
		d = sequence.count("D")	
		c = sequence.count("C")	
		q = sequence.count("Q")	
		e = sequence.count("E")	
		g = sequence.count("G")	
		h = sequence.count("H")	
		i = sequence.count("I")	
		l = sequence.count("L")	
		k = sequence.count("K")	
		m = sequence.count("M")	
		f = sequence.count("F")	
		p = sequence.count("P")	
		s = sequence.count("S")
		t = sequence.count("T")
		w = sequence.count("W")
		y = sequence.count("Y")
		v = sequence.count("V")
			
		carbon=a*3 + r*6 + n*4 + d*4 + c*3 + 5*q + 5*e + 2*g + 6*h + 6*i + 6*l + 6*k + 5*m + f*9 + p*5 + s*3 + t*4 + w*11 + y*9 + v*5 
		hydrogen= 2 + 5*a + 12*r + 6*n + 5*d + 5*c + 8*q + 7*e + 3*g + 7*h + 11*i + 11*l + 12*k + 9*m + 9*f + 7*p + 5*s + 7*t + 10*w + 9*y + 9*v
		nitrogen= a + 4*r + 2*n + d + c + 2*q + e + g + 3*h + i + l + 2*k + m + p + f + s + t + 2*w + y + v
		sulfur= m + c
		oxygen= 1 + a + c + 2*n + 3*d + c + 2*q + 3*e + g + h + i + l + k + m + f + p + 2*s + 2*t + w + 2*y + v
		total_mass= carbon * 12.0116 + oxygen * 15.99977 + nitrogen * 14.00643 + sulfur * 32.059 + hydrogen * 1.00784
		
		return total_mass  


def mass_calculator_tryptic (accession_number):
	tryptic_dictionary = {}
	
	sequences = Protein_searcher.tryptic_cutter(self)
	
	for sequence in sequences:
	
		a = sequence.count("A")
		r = sequence.count("R")
		n = sequence.count("N")	
		d = sequence.count("D")	
		c = sequence.count("C")	
		q = sequence.count("Q")	
		e = sequence.count("E")	
		g = sequence.count("G")	
		h = sequence.count("H")	
		i = sequence.count("I")	
		l = sequence.count("L")	
		k = sequence.count("K")	
		m = sequence.count("M")	
		f = sequence.count("F")	
		p = sequence.count("P")	
		s = sequence.count("S")
		t = sequence.count("T")
		w = sequence.count("W")
		y = sequence.count("Y")
		v = sequence.count("V")
			
		carbon = a*3 + r*6 + n*4 + d*4 + c*3 + 5*q + 5*e + 2*g + 6*h + 6*i + 6*l + 6*k + 5*m + f*9 + p*5 + s*3 + t*4 + w*11 + y*9 + v*5 
		hydrogen = 2 + 5*a + 12*r + 6*n + 5*d + 5*c + 8*q + 7*e + 3*g + 7*h + 11*i + 11*l + 12*k + 9*m + 9*f + 7*p + 5*s + 7*t + 10*w + 9*y + 9*v
		nitrogen = a + 4*r + 2*n + d + c + 2*q + e + g + 3*h + i + l + 2*k + m + p + f + s + t + 2*w + y + v
		sulfur = m + c
		oxygen = 1 + a + c + 2*n + 3*d + c + 2*q + 3*e + g + h + i + l + k + m + f + p + 2*s + 2*t + w + 2*y + v
		total_mass= carbon * 12.0116 + oxygen * 15.99977 + nitrogen * 14.00643 + sulfur * 32.059 + hydrogen * 1.00784
	
		tryptic_dictionary[sequence] = total_mass	
		
		return tryptic_dictionary	

	
def mass_calculator_isoforms (accession_number):
	
	isoform_numbers = []
	isoforms_dictionary = {}
	i =	1
	sequences_dictionary = Protein_searcher.isoform_creator (self)
	isoform_numbers = sequences_dictionary.keys()
	
	for isoform_number in isoform_numbers:
		
		seq = sequences_dictionary[isoform_number]
		
		a = seq.count("A")
		r = seq.count("R")
		n = seq.count("N")	
		d = seq.count("D")	
		c = seq.count("C")	
		q = seq.count("Q")	
		e = seq.count("E")	
		g = seq.count("G")	
		h = seq.count("H")	
		i = seq.count("I")	
		l = seq.count("L")	
		k = seq.count("K")	
		m = seq.count("M")	
		f = seq.count("F")	
		p = seq.count("P")	
		s = seq.count("S")
		t = seq.count("T")
		w = seq.count("W")
		y = seq.count("Y")
		v = seq.count("V")
			
		carbon = a*3 + r*6 + n*4 + d*4 + c*3 + 5*q + 5*e + 2*g + 6*h + 6*i + 6*l + 6*k + 5*m + f*9 + p*5 + s*3 + t*4 + w*11 + y*9 + v*5 
		hydrogen = 2 + 5*a + 12*r + 6*n + 5*d + 5*c + 8*q + 7*e + 3*g + 7*h + 11*i + 11*l + 12*k + 9*m + 9*f + 7*p + 5*s + 7*t + 10*w + 9*y + 9*v
		nitrogen = a + 4*r + 2*n + d + c + 2*q + e + g + 3*h + i + l + 2*k + m + p + f + s + t + 2*w + y + v
		sulfur = m + c
		oxygen = 1 + a + c + 2*n + 3*d + c + 2*q + 3*e + g + h + i + l + k + m + f + p + 2*s + 2*t + w + 2*y + v
		total_mass= carbon * 12.0116 + oxygen * 15.99977 + nitrogen * 14.00643 + sulfur * 32.059 + hydrogen * 1.00784
		
		isoforms_dictionary[isoform_number] = total_mass
		
		return isoforms_dictionary