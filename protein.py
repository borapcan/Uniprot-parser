import re
import os
from abc import ABC, abstractmethod
import xml.etree.ElementTree as ET


class AbcPeptide(ABC):

	@abstractmethod
	def sequence(self):
		pass

	@property
	def mass(self):
		"""This function calculates protein mass.

			Parameters
			----------
			sequence: str

			Returns
			-------
			total_mass: int
				Returns protein mass.

			"""

		amino_atoms = {
			'A': (3, 5, 1, 1, 0, 0), 'C': (3, 5, 1, 1, 1, 0), 'D': (4, 5, 3, 1, 0, 0), 'E': (5, 7, 3, 1, 0, 0),
			'F': (9, 9, 1, 1, 0, 0), 'G': (2, 3, 1, 1, 0, 0), 'H': (6, 7, 1, 3, 0, 0), 'I': (6, 11, 1, 1, 0, 0),
			'K': (6, 12, 1, 2, 0, 0), 'L': (6, 11, 1, 1, 0, 0), 'M': (5, 9, 1, 1, 1, 0), 'N': (4, 6, 2, 2, 0, 0),
			'P': (5, 7, 1, 1, 0, 0), 'Q': (5, 8, 2, 2, 0, 0), 'R': (6, 12, 1, 4, 0, 0), 'S': (3, 5, 2, 1, 0, 0),
			'T': (4, 7, 2, 1, 0, 0), 'V': (5, 9, 1, 1, 0, 0), 'W': (11, 10, 1, 2, 0, 0), 'Y': (9, 9, 2, 1, 0, 0),
			'U':  (3, 5, 1, 1, 0, 1)
		}

		carbons = 0
		hydrogens = 2
		oxygens = 1
		nitrogens = 0
		sulfurs = 0
		seleniums = 0

		sequence = self.sequence.upper()

		for letter in sequence:
			carbon, hydrogen, oxygen, nitrogen, sulfur, selenium = amino_atoms[letter]
			carbons += carbon
			hydrogens += hydrogen
			oxygens += oxygen
			nitrogens += nitrogen
			sulfurs += sulfur
			seleniums += selenium

		total_mass = carbons * 12.0116 + oxygens * 15.99977 + nitrogens * 14.00643 + sulfurs * 32.059 + hydrogens * 1.00784 + seleniums*78.971

		return total_mass

	@property
	def tryptic_peptides(self):
		"""This function makes trypsine cuts in a protein for all relevant cuts ( longer than 2 amino acids).

		Parameters
		----------
		sequence: str

		Returns
		-------
		total_peptide: dict
			Returns dictionary of a order number of a trypsine cut with the cut, start and end position of a cut and the
			mass of the cut.

		"""
		tryptic_peptides = []
		pattern = re.compile('([KR]?[^P].*?[KR](?!P))')
		tryptic_list = pattern.split(self.sequence)
		k = 1
		for cut in tryptic_list:
			cut_start = k
			cut_end = k + len(cut) - 1
			if len(cut) > 0:
				peptide_cut = Peptide(cut)
				peptide_cut.start, peptide_cut.end = cut_start, cut_end
				tryptic_peptides.append(peptide_cut)
				k += len(cut)


		return tryptic_peptides

	@property
	def avg_hydrophobicity (self):
		# Amino acid scale: Hydration potential (kcal/mole) at 25°C
		# Author(s): Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F.
		# authors haven't provided hydrophobicity for Sec/U so it was designated with a 0.000
		amino_atoms = {
			'A': (1.940), 'C': (-1.240), 'D': (-10.950), 'E': (-10.200),
			'F': (-0.760), 'G': (2.390), 'H': (-10.270), 'I': (2.150),
			'K': (-9.520), 'L': (2.280), 'M': (-1.480), 'N': (-9.680),
			'P': (0.000), 'Q': (-9.380), 'R': (-19.920), 'S': (-5.060),
			'T': (-4.880), 'V': (1.990), 'W': (-5.880), 'Y': (-6.110),
			'U': (0.000)
		}

		hydrophobicity = 0.000
		sequence = self.sequence.upper()

		for letter in sequence:
			hydrophobicity += amino_atoms[letter]

		return hydrophobicity


class Peptide(AbcPeptide):

	def __init__(self, input_sequence):
		self.input_sequence = input_sequence

	@property
	def sequence(self):
		return self.input_sequence


class AnnotatedPeptide(Peptide):
	"""Class that uses XML document to return data relevant to a single UNIPROT protein entry using an accession number.
		Accession number must be valid in order for entry to be found (https://www.uniprot.org/help/accession_numbers).
		Class returns proteins name, canonical sequence, glycosylation sites, isoforms.

		Parameters
		----------
		accession_number : str
				as described above
		path : str
		Attributes
		----------
				Location of XML file database of proteins
		ns : str
				XML namespace used for parsing
		entry :obj:`str`
				object of a single protein entry used for data mining

		Examples
		--------
		>>> protein = AnnotatedPeptide("gli.xml", "P52798")
		>>> protein.name
		'Ephrin-A4'

		>>> protein.sequence
		'MRLLPLLRTVLWAAFLGSPLRGGSSLRHVVYWNSSNPRLLRGDAVVELGLNDYLDIVCPHYEGPGPPEGPETFALYMVDWPGYESCQAEGPRAYKRWVCSLPFGHVQFSEKIQRFTPFSLGFEFLPGETYYYISVPTPESSGQCLRLQVSVCCKERKSESAHPVGSPGESGTSGWRGGDTPSPLCLLLLLLLLILRLLRIL'

		>>> protein.glycosylation
		{33: 'N-linked (GlcNAc...) asparagine'}

		>>> protein.isoforms[2].sequence
		'MRLLPLLRTVLWAAFLGSPLRGGSSLRHVVYWNSSNPRLLRGDAVVELGLNDYLDIVCPHYEGPGPPEGPETFALYMVDWPGYESCQAEGPRAYKRWVCSLPFGHVQFSEKIQRFTPFSLGFEFLPGETYYYISVPTPESSGQCLRLQVSVCCKERNLPSHPKEPESSQDPLEEEGSLLPALGVPIQTDKMEH'

	"""

	def __init__(self, path, accession_number):
		self.path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), path))
		self.ns = {"uniprot": "http://uniprot.org/uniprot"}
		self.accession_number = accession_number
		self.entry = self._entry_finder()
		self.protein = None

	def _entry_finder(self):
		parse = ET.iterparse(self.path, events=("start", "end"))
		parse = iter(parse)
		event, root = next(parse)

		for event, elem in parse:
			if event == "end" and elem.tag.endswith("}entry"):
				accession_elem = elem.findall("{" + self.ns['uniprot'] + "}" + "accession")
				if self.accession_number in [i.text for i in accession_elem]:
					return elem
					root.remove(elem)

	@property
	def name(self):
		"""Returns name of the protein.

		"""
		protein_names = self.entry.find("uniprot:protein", self.ns)
		recommended_name = protein_names.find("uniprot:recommendedName", self.ns)
		full_name = recommended_name.find("uniprot:fullName", self.ns).text
		return full_name

	@property
	def sequence(self):
		"""Returns canonical sequence of the protein.

		"""

		sequence = self.entry.find("uniprot:sequence", self.ns).text
		sequence = ''.join(sequence.split())
		return sequence

	@property
	def glycosylation(self):
		"""Returns dictionary of glycosilation sites and their types.

		"""
		try:
			glycosylation_dictionary = {}
			for feature in self.entry.findall("uniprot:feature[@type='glycosylation site']", self.ns):
				glycosylation_feature = feature.attrib
				glycosylation_type = glycosylation_feature["description"]
				location = feature.find("uniprot:location", self.ns)
				glycosylation_site = location.find("uniprot:position", self.ns).attrib
				glycosylation_site = int(glycosylation_site["position"])
				glycosylation_dictionary[glycosylation_site] = glycosylation_type
		except AttributeError:
			glycosylation_dictionary = {}
		return glycosylation_dictionary

	def _isoform_numbers(self):

		isoform_numbers = []
		comment = self.entry.find("uniprot:comment[@type='alternative products']", self.ns)
		try:
			for isoform in comment.findall("uniprot:isoform", self.ns):
				isoform_number = isoform.find("uniprot:name", self.ns).text
				if int(isoform_number) == 1:
					pass
				else:
					isoform_numbers.append(int(isoform_number))
		except AttributeError:
			isoform_numbers = []
		return isoform_numbers

	@property
	def isoforms(self):
		isoforms = {}
		if self._isoform_numbers() == []:
			return AttributeError("No isoforms found")
		else:
			for isoform_number in self._isoform_numbers():
				isoform = Isoform.__new__(Isoform)
				isoform.entry, isoform.isoform_number, isoform.ns = self.entry, isoform_number, self.ns
				isoforms[isoform_number] = isoform
			return isoforms


class Isoform(AnnotatedPeptide):

	def __init__(self, path, accession_number, isoform_number):
		super().__init__(path, accession_number)
		self.isoform_number = isoform_number

	@property
	def name(self):
		return super().name + "-" + str(self.isoform_number)

	def _alternative_products(self):

		isoform_ids = {}
		comment = self.entry.find("uniprot:comment[@type='alternative products']", self.ns)
		try:
			for isoform in comment.findall("uniprot:isoform", self.ns):
				isoform_number = int(isoform.find("uniprot:name", self.ns).text)
				ref_id = isoform.find("uniprot:sequence", self.ns).attrib
				if "ref" in ref_id:
					isoform_ids[isoform_number]=(ref_id["ref"].split())
		except AttributeError:
			isoform_ids = {}
		return isoform_ids

	def _vsps(self):

		ids_information = {}
		for feature in self.entry.findall("uniprot:feature[@type='splice variant']", self.ns):
			feat = feature.attrib

			id_ = feat["id"]

			variation = feature.find("uniprot:variation", self.ns)

			if variation is None:
				variation = ""
			else:
				variation = variation.text

			location = feature.find("uniprot:location", self.ns)

			try:
				begin = int(location[0].get("position"))
				end = int(location[1].get("position"))
			except IndexError:
				position = int(location[0].get("position"))
				begin = position
				end = position

			ids_information[id_] = (variation, begin, end)
		return ids_information

	@property
	def sequence(self):
		id_list = self._alternative_products()[self.isoform_number]
		isoform_information = [self._vsps()[key] for key in id_list]
		isoform_information.sort(key=lambda x: x[1], reverse=True)
		isoform_sequence = super().sequence

		for single_variation in isoform_information:
			variation, index_start, \
			index_end = single_variation
			isoform_sequence = isoform_sequence[:index_start - 1] + variation + isoform_sequence[index_end:]

		return isoform_sequence

	def _isoform_glycosylation(self,isoform_number):
		"""Returns dictionary of a secific isoform glycosilation sites and their types.

		This method filters annotated glycosylation sites for a single isoform by checking 3 specific cases.
		"""
		isoform_glycosylation = {}
		id_list = self._alternative_products()[isoform_number]
		variation_information_collection = [self._vsps()[key] for key in id_list]
		appended_sites = []
		index_end_last = 0
		i = 0
		evariation, eindex_start, \
		eindex_end = variation_information_collection[-1]
		for single_variation in variation_information_collection:
			variation, index_start,\
			index_end = single_variation
			i2 = i
			i = index_end - index_start + 1 - len(variation) + i
			for site in super().glycosylation.keys():
				if site not in appended_sites:
					if index_end_last == 0:
						if site < index_start:
							isoform_glycosylation[site] = super().glycosylation[site]
							appended_sites.append(site)
					elif site > index_end_last:
						if site < index_start:
							appended_sites.append(site)
							sitem = site - i2
							isoform_glycosylation[sitem] = super().glycosylation[site]
					if index_end == eindex_end:
						if site > index_end:
							appended_sites.append(site)
							sitem = site - i
							isoform_glycosylation[sitem] = super().glycosylation[site]
			index_end_last = index_end
		return isoform_glycosylation

	@property
	def glycosylation(self):
		return self._isoform_glycosylation(self.isoform_number)

	@property
	def isoforms(self):
		return NotImplementedError()

# factory to create final product
def protein(sequence=None, isoform_number=None,  acc=None, xmlfile=None):
	if not xmlfile:
		xmlfile ='gli.xml'
	if not isoform_number:
		if acc:
				if re.match('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', acc):
					return AnnotatedPeptide(xmlfile, acc)
				else:
					raise ValueError('You have entered a non-valid accession number, please review your sequence before '
																					'trying to process it again.')
		elif sequence:
			sequence = sequence.upper()
			if re.match("[ACDEFGHIKLMNPQRSTUVWY]", sequence):
				return Peptide(sequence)
			else:
				raise ValueError('You have entered a non-valid sequence, please review your sequence before '
																							'trying to process it again.')
		else:
			raise ValueError('Sequence or accession number are required!')
	else:
		if acc:
			return Isoform(xmlfile, acc, isoform_number)
		else:
			return ValueError('Both isoform number and accession number are required!')


