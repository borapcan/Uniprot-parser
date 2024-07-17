import re, os
import argparse
from contextlib import redirect_stdout
import logging
import json
import datetime

from utilities.helpers import sequence_helper
from protein import protein
from utilities.types import ProteinOutput


def Main():
	# reading config
	with open("config.json", "r") as read:
		cfg = json.load(read)

	# setting up argparser
	parser = argparse.ArgumentParser()
	parser.add_argument("-p", "--path", help="Path to XML file containing protein data you want to mine. Accession "
																					"numbers must be separated.", type=str)
	parser.add_argument("sequences", help="As input write the string of your sequence or your accession number",
																					type=str, nargs='+')
	parser.add_argument("--debug", help="Show debug messages to stdout", action="store_true")
	parser.add_argument("-q", "--quiet", help = "Show compact output", action="store_true")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-o", "--output", help="Output your search or result as a human readable file", action="store_true")
	group.add_argument("-j", "--json", help="Output your query results to json file(computer readable file)",
																				action="store_true")

	args = parser.parse_args()

	# setting up logger
	logger = logging.getLogger('MAIN')
	if args.debug:
		logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s][%(levelname)s]  %(message)s")
	else:
		logging.basicConfig(level=logging.INFO, format="[%(asctime)s][%(levelname)s]  %(message)s")

	sequence_counter = 1

	output = {}
	isoforms = {}

	if args.path:
		if args.path.endswith(".xml"):
			path = os.path.normpath(args.path)
		else:
			raise EnvironmentError(f"Initializing parsing failed with file: {args.path}")
	else:
		path = None

	for seq in args.sequences:
		t1 = datetime.datetime.timestamp(datetime.datetime.now())
		arg_sequence = seq.upper()  # arg_sequences upper cased for proper matching
		logger.debug(f"Working on sequence: {arg_sequence}")
		if re.match('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', arg_sequence):
			pprotein = protein(xmlfile=path, acc=arg_sequence)
			logger.debug(f"Looking for your entry in input file: {args.path}")
			logger.debug(f"MATCH FOUND FOR SEQUENCE: {arg_sequence} -- PROCESSING --")
			for protein_number in pprotein.isoforms.keys():
				isoforms[protein_number] = pprotein.isoforms[protein_number].sequence
			parsed_protein = ProteinOutput(
							name=pprotein.name,
							sequence=pprotein.sequence,
							mass=pprotein.mass,
							glycosylationSites=pprotein.glycosylation,
							trypticDigestion=pprotein.tryptic_peptides,
							isoforms=isoforms
							)
			if args.json:
				output[parsed_protein.name] = parsed_protein.json_output()
			else:
				output[parsed_protein.name] = parsed_protein

		elif re.match("[ACDEFGHIKLMNPQRSTUVWY]+$", arg_sequence):
			logger.debug("WORKING ON CUSTOM IMPORT SEQUENCE")
			pprotein = protein(xmlfile=path, sequence=arg_sequence)
			parsed_protein = sequence_helper(pprotein, sequence_counter)
			if args.json:
				output[parsed_protein.name] = parsed_protein.json_output()
			else:
				output[parsed_protein.name] = parsed_protein
		else:
			raise ValueError(
				'You have entered a non-existent sequence, please review your sequence before trying to process it again.')
		t2 = datetime.datetime.timestamp(datetime.datetime.now())
		logger.debug(f"Processing took: {t2-t1}s")
		sequence_counter += 1

	if args.output:
		logger.debug("Outputting data as text")
		with open(cfg["text_output"], 'w') as out:
			with redirect_stdout(out):
				for pprotein in output.values():
					pprotein.formatted_output()
	elif args.json:
		logger.debug("Outputting data as json")
		with open(cfg["json_output"], 'w') as outfile:
			json.dump(output, outfile)
	else:
		if args.quiet:
			for pprotein in output.values():
				print(pprotein)
		else:
			for pprotein in output.values():
				pprotein.formatted_output()


if __name__ == "__main__":
	Main()
