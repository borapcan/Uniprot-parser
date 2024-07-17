from utilities.types import ProteinOutput


def sequence_helper(protein, sequence_counter: int):
	""""sequence_helper receives sequence and returns:

	- protein name
	- protein sequence
	- protein mass
	- glycosylation sites (empty)
	- tryptic digest sites
	- isoforms (empty)

	Parameters
	----------
	sequence : str
			sequence of a custom protein
	sequence_counter : int
			custom sequence position in interface
	Returns
	-------
	output_variable : named tuple
			named tuple of protein data
	"""
	name = "Input " + "(" + str(sequence_counter) + ")"
	sequence = protein.sequence
	mass = protein.mass
	glycosylation_sites = {}
	tryptic_digestion = protein.tryptic_peptides
	p_isoforms = {}

	output_variable = ProteinOutput(
		name=name,
		sequence=sequence,
		mass=mass,
		glycosylationSites=glycosylation_sites,
		trypticDigestion=tryptic_digestion,
		isoforms=p_isoforms
	)
	return output_variable
