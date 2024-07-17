from typing import NamedTuple, Dict, Optional
from pprint import PrettyPrinter
from textwrap import TextWrapper


pp = PrettyPrinter(indent=4)
wrapper = TextWrapper(width=81, initial_indent="\t", subsequent_indent="\t")


class ProteinOutput(NamedTuple):
	"""	Data structure holding information about a protein

	- used in Interface to facilitate outputing to stdout or a file
	- attributes can be accessed using 'dot notation'

	Example
	--------
	>>> p = Protein(
		name="FINC_HUMAN",
		sequence= "MLRGPGPGLLLLAVQCLGTAVPSTGASKSKRQAQQMVQPQS...QADREDSRE",
		mass=262630.05697,
		glycosylationSites: {279: 'O-linked (GalNAc...) threonine',...},
		trypticDigestion: ["MLR", "GPGPGLLLLAVQCLGTAVPSTGASK", "SK",...],
		isoforms: {"Isoform 2": "MLRGPGPGLLLLAVQCLGTAVPSTGASKSKRQAQQMVQPQSPVAVSQSK...",...},
		massUnit: "Da"  # defaults to "Da" no need to write
		)

	>>> print(p.name)
	FINC_HUMAN
	"""
	sequence: str
	mass: float
	trypticDigestion: list
	name: Optional[str] = "Custom input sequence"
	glycosylationSites: Optional[Dict] = {}
	isoforms: Optional[Dict] = {}
	massUnit: str = "Da"  # if kDa -> replace: massUnit="kDa"

	def _wrap_dictionary_helper(self, dictionary, wrapper, in_line=False):

		if in_line:
			for key, value in dictionary.items():
				print(f"\t{key}: {value}")
		else:
			for key, value in dictionary.items():
				print(f"{key}:\n{wrapper.fill(value)}")

	def _wrap_list_helper(self, list):
		for item in list:
			print(f"Cut Sequence: {item.sequence}, Cut start: {item.start}, Cut end: {item.end}")


	def formatted_output(self):
		""" formatted_output prepares and prints Protein in human readable form """

		print(
			f"Name: {self.name}\nSequence:\n{wrapper.fill(self.sequence)}\nMass: {self.mass} {self.massUnit}\n"
			f"Tryptic Digests: ")
		self._wrap_list_helper(list=self.trypticDigestion)
		print("Glycosylation Sites:")
		self._wrap_dictionary_helper(dictionary=self.glycosylationSites, wrapper=wrapper, in_line=True)
		print("Isoforms:")
		self._wrap_dictionary_helper(dictionary=self.isoforms, wrapper=wrapper, in_line=True)
		print("\n"+">"*40+"FINISH"+"<"*47+"\n")

	def json_output(self):
		""" json_output prepares Protein in json format """
		return_object = {
			"name": self.name,
			"sequence": self.sequence,
			"mass": self.mass,
			"tryptic product": self.trypticDigestion,
			"isoforms": self.isoforms,
			"glycosylation sites": self.glycosylationSites
		}
		return return_object

	def __str__(self):
		return f"Name: {self.name}\nMass: {self.mass} {self.massUnit}\nSequence:\n{wrapper.fill(self.sequence)}\n" \
				f"Number of tryptic products: {len(self.trypticDigestion)}\n" \
				f"Number of glycosylation sites: {len(self.glycosylationSites.keys())}\n" \
				f"Number of isoforms: {len(self.isoforms.keys())}"

	def __repr__(self):
		return pp.pformat(self._asdict())








