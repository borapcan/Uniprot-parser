import os

import xml.etree.ElementTree as et

base_path = os.path.dirname(os.path.realpath("search.py"))
xml_file = os.path.join(base_path, "gli.xml")
tree = et.parse(xml_file)
root = tree.getroot()

ns={"uniprot":"http://uniprot.org/uniprot"}

def canonical_sequence (accession_number):
	sequence=""
	for entry in root.findall("uniprot:entry",ns):
		if accession_number == entry.find("uniprot:accession",ns).text :
			sequence = entry.find("uniprot:sequence",ns).text
			sequence = ''.join(sequence.split())
			return (sequence)
			
def glycosilation_sites (accession_number):
	glycosilation_list = []
	location = ""
	glycosilation_site = {}
	for entry in root.findall("uniprot:entry",ns):
		if accession_number == entry.find("uniprot:accession",ns).text :
			for feature in entry.findall("uniprot:feature[@type='glycosylation site']" ,ns ):
				location = feature.find("uniprot:location",ns)
				glycosilation_site = location.find("uniprot:position",ns).attrib
				glycosilation_site = int(glycosilation_site["position"])
				glycosilation_list.append(glycosilation_site)
	return glycosilation_list
	

def alternative_seq_ids_dictionary (accession_number):
	iso_dict={}
	isoform_number=""
	ref_id={}
	for entry in root.findall("uniprot:entry",ns):
		if accession_number == entry.find("uniprot:accession",ns).text :
			for comment in entry.findall("uniprot:comment[@type='alternative products']",ns):
				for isoform in comment.findall("uniprot:isoform",ns):
					isoform_number=isoform.find("uniprot:name",ns).text
					ref_id = isoform.find("uniprot:sequence",ns).attrib
					if "ref" in ref_id:
						iso_dict[isoform_number]= ref_id["ref"].split()    
	return(iso_dict)

	
def isoform_data_dictionary(accession_number):
	ids_dict={}
	ids = ""
	original = ""
	variation = ""
	index_st = 0
	index_end = 0
	for entry in root.findall("uniprot:entry",ns):
		if accession_number == entry.find("uniprot:accession",ns).text :	
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
	
