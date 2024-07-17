#IMPORT 

import os

import re

import xml.etree.cElementTree as et

base_path = os.path.dirname(os.path.realpath("glycoproteins.py"))
xml_file = os.path.join(base_path, "gli.xml")
tree = et.parse(xml_file)
root = tree.getroot()

ns={"uniprot":"http://uniprot.org/uniprot"}



#FUNTION WHICH CREATES LISTS FOR SEQUENCES AND ACCESSION NUMBERS AND CREATES DICTIONARIES FOR ISOFORMS AND THEIR SPECIFIC IDENTIFIERS

iso_dict= {}
ids_dict= {}
accession_list=[]
sequence_list=[]

for entry in root.findall("uniprot:entry",ns):
	sequence = entry.find("uniprot:sequence",ns).text
	sequence = ''.join(sequence.split())
	sequence_list.append(sequence)
	accession=entry.find("uniprot:accession",ns).text
	accession_list.append(accession)
	for comment in entry.findall("uniprot:comment[@type='alternative products']",ns):
		for isoform in comment.findall("uniprot:isoform",ns):
			iso = accession + '-' + isoform.find("uniprot:name",ns).text
			ref_id=isoform.find("uniprot:sequence",ns).attrib
			
#CREATION OF A DICTIONARY ACCESSION NUMBER-ISFORM NUMBER : VSP LIST	
				
			if "ref" in ref_id:
				iso_dict[iso]= ref_id["ref"].split()      
	for feature in entry.findall("uniprot:feature[@type='splice variant']",ns):
		feat = feature.attrib
		ids = feat["id"]
		
# GETTING THE ORIGINAL SEQUENCE, "" IS ADDED FOR DELETIONS	
	
		original=feature.find("uniprot:original",ns)
		if original is None:
			original= ""	
		else:
			original=feature.find("uniprot:original",ns).text
			
# GETTING THE VARIATION SEQUENCE, "" IS ADDED FOR DELETIONS	
			
		variation=feature.find("uniprot:variation",ns)
		if variation is None:
			variation= ""	
		else:
			variation=feature.find("uniprot:variation",ns).text	
		location = feature.find("uniprot:location",ns)
		index_st=location.find("uniprot:begin",ns)
		
# RETURNING THE VARIANT START POSTION, TO MARK THE STARTING POSITION FOR THE CUT IN CANONICAL SEQUENCE
		
		if index_st is None:	
			index_st=location.find("uniprot:position",ns).attrib
			index_st=int(index_st["position"])	
		else:
			index_st=location.find("uniprot:begin",ns).attrib	
			index_st=int(index_st["position"])
		index_end=location.find("uniprot:end",ns)
		
# RETURNING THE SEQUENCE END INDEX, TO MARK THE STARTING POSITION FOR THE CUT IN CANONICAL SEQUENCE	
			
		if index_end is None:		
			index_end=location.find("uniprot:position",ns).attrib			
			index_end=int(index_end["position"])			
		else:		
			index_end=location.find("uniprot:end",ns).attrib			
			index_end=int(index_end["position"])
			
#FORMATION OF THE VSP DICTIONARY, CONTAINING FORMER VALUES

		ids_dict[ids] = (original, variation, index_st, index_end)
        
		
#FUNCTION WHICH TAKES SEQUENCES FROM THEIR RESPECTIVE LISTS AND THEIR VSPS (ref_id variable) AND USES THEM TO CREATE AN UNIQUE LIST OF SEQUENCE NAMES AND THEIR ISOFORMS		


original, variation, index_st, index_end = ids_dict[ids]
ref_id=iso_dict[iso]
isoforms_dictionary = {}	
#SEARCHING VSP KEYS	
for ref_id,idovi  in iso_dict.items():	
#SETTING A TEMPORARY SEQUENCE	
	for seqt in sequence_list:		
		temporary_sequence = seqt			
		tmp_index = 0			
			for ids in idovi:
			
	#REPLACING ORIGINAL WITH THE VARIANT AND CORRECTING THE LENGTH		
	
				original, variation, index_st, index_end = ids_dict[ids]
				temporary_sequence = '%s%s%s'%(seqt[:index_st-1 + tmp_index],variation,seqt[index_end+tmp_index:])
				if len(original) == index_end - index_st:				
					tmp_index += len(variation)-len(original)					
				else:				
					tmp_index += index_end - index_st					
			isoforms_dictionary[ref_id] = temporary_sequence			
	return(isoforms_dictionary)
	

#SEQUENCE MASS CALCULATIONS

mass_list =	 {}
for seq in sequence list
	for acs in accession_list:
		
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
			
			carbon=a*3 + r*6 + n*4 + d*4 + c*3 + 5*q + 5*e + 2*g + 6*h + 6*i + 6*l + 6*k + 5*m + f*9 + p*5 + s*3 + t*4 + w*11 + y*9 + v*5 
			hydrogen= 2 + 5*a + 12*r + 6*n + 5*d + 5*c + 8*q + 7*e + 3*g + 7*h + 11*i + 11*l + 12*k + 9*m + 9*f + 7*p + 5*s + 7*t + 10*w + 9*y + 9*v
			nitrogen= a + 4*r + 2*n + d + c + 2*q + e + g + 3*h + i + l + 2*k + m + p + f + s + t + 2*w + y + v
			sulfur= m + c
			oxygen= 1 + a + c + 2*n + 3*d + c + 2*q + 3*e + g + h + i + l + k + m + f + p + 2*s + 2*t + w + 2*y + v
			total_mass= carbon * 12.0107 + oxygen * 15.9994 + nitrogen * 14.0067 + sulfur * 32.065 + hydrogen * 1.0079
			
			mass_list[acs]=total_mass
			

	#TRYPTIC DIGESTION OF PROTEINS

tryptic_list =	 {}

for seq in sequence_list:
	for acs in accession_list:
		protein = re.sub(r'(?<=[RK])(?=[^P])','\n', seq, re.DOTALL)
		tryptic_list  =	protein=protein.split()
			


#GLYCOSILATION POSITION


glycosilation_position = {}
for entry in root.findall("uniprot:entry",ns):
	i=0
	for feature in entry.findall("uniprot:feature[@type='glycosylation site']" ,ns ):
			
		i=i+1
		acs=entry.find("uniprot:accession",ns).text
		location = feature.find("uniprot:location",ns)
		glycosilation_site=location.find("uniprot:position",ns).attrib
		glycosilation_site=int(glycosilation_site["position"])
		acs = acs + ':' + str(i)
		glycosilation_position[acs]=glycosilation_site




	