import search

import tryptic_cutter

import isoforms

# komentar #1:
# Sve funkcije imaju istu logiku, dio koda je potpuno identičan.
# Ako trebaš isti kod na više mjesta, izdvoji ga u zasebnu funkciju.

# komentar #3:
# Računanje mase peptida ovisi samo o sekvenci peptida, pa bi to trebao biti i
# jedini argument.

def mass_calculator_canonical (accession_number):
	
	sequence = search.canonical_sequence(accession_number)
	

        # komentar #2:
        # Ovako brojanje elemenata nije čitljivo, pa time niti jednostavno za
        # održavanje.  Ono što bi trebao napraviti jest for petlju po sequence
        # ili for petlju po aminokiselinama.

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
	sequences = tryptic_cutter.tryptic_cutter(accession_number)
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
	sequences_dictionary = isoforms.isoform_creator (accession_number)
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
			
	
