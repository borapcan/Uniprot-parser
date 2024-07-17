import search


# komentar #5
# Teško mi je pratiti isoform_creator. Molim te preimenuj varijable opisnije.

# komentar #6
# Samo kao podsjetnik:
# Za svaku izoformu modificiraš kanonsku sekvencu od kraja prema naprijed kako
# ne bi imao problema s pozicijama opisanih varijanti.
# Dakle skupiš sve splice varijante za izoformu, te ih sortiraš prema startnoj
# poziciji unatrag i zatim ih uklopiš u kanonsku sekvencu jednu po jednu.

isoforms_dict = dict()
def isoform_creator (accession_number):
	
	seqt = search.canonical_sequence(accession_number)
	ids_dict = search.isoform_data_dictionary(accession_number)
	iso_dict = search.alternative_seq_ids_dictionary (accession_number)
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




		

