import os
import json
from tkinter import *
from tkinter import messagebox
from contextlib import redirect_stdout
from pprint import PrettyPrinter
from textwrap import TextWrapper


from utilities.helpers import sequence_helper
from protein import protein
from utilities.types import ProteinOutput

pp = PrettyPrinter(indent=4)
wrapper = TextWrapper(width=81, initial_indent="\t", subsequent_indent="\t")


class ProteinGui():
	def __init__(self):

		self.root = Tk()
		self.root.configure(bg="orange")
		self.root.title("Diplomski-rad: Protein queries")

		self.statusbar_text = StringVar()
		self.statusbar_text.set("Ready to begin your query")

		self.title_frame = Frame(self.root, bg="#dfdfdf")

		self.title = Label(self.title_frame, text="Diplomski-rad: Protein queries", bg="#dfdfdf")
		self.title.pack()


		self.canvas_area = Canvas(self.root, background="#cce0ff")
		self.canvas_area.grid(row=1, column=1)

		self.path_label = Label(self.canvas_area, text="Input your path", bg="#cce0ff")
		self.sequences_label = Label(self.canvas_area, text="Input your sequences*", bg="#cce0ff")
		self.path_entry = Entry(self.canvas_area, width=250)
		self.sequence_entry = Entry(self.canvas_area, width=250)

		self.path_label.grid(row=0, sticky=W+E)
		self.sequences_label.grid(row=30, sticky=W+E)

		self.path_entry.grid(row=0, column=12, sticky=W+E)
		self.sequence_entry.grid(row=30, column=12, sticky=W+E)

		self.text = Checkbutton(self.canvas_area, text="Output as text document", bg="#cce0ff")
		self.json = Checkbutton(self.canvas_area, text="Output as json document", bg="#cce0ff")

		self.text.grid(row=60, columnspan=20)
		self.json.grid(row=90, columnspan=20)

		self.run_button = Button(self.canvas_area, text="Run", command=lambda : self.Run_GUI(), bg="green")
		self.run_button.grid(row=150, column=0, columnspan=20, sticky=W)

		self.quit_button = Button(self.canvas_area, text="Quit", command=self.root.quit, bg="red")
		self.quit_button.grid(row=150, column=200, columnspan=20, sticky=W)

		self.status_frame = Frame(self.root)
		self.status = Label(self.status_frame, textvariable=self.statusbar_text, relief=SUNKEN, bd=2, anchor="w")
		self.status.pack(fill="both", expand=True)

		self.title_frame.grid(row=0, column=1, sticky="ew")
		self.canvas_area.grid(row=1, column=1, sticky="nsew")
		self.status_frame.grid(row=2, column=0, columnspan=2, sticky="ew")

		self.root.grid_rowconfigure(1, weight=1)
		self.root.grid_columnconfigure(1, weight=1)

		self.root.mainloop()

	def _wrap_dictionary_helper(self, dictionary, wrapper, in_line=False):
		output = ""
		if in_line:
			for key, value in dictionary.items():
				output = output + f"\t{key}: {value}\n"
		else:
			for key, value in dictionary.items():
				output = output + f"{key}:\n{wrapper.fill(value)}\n"
		return output

	def _wrap_list_helper(self, list):
		output = ""
		for item in list:
			output = output + f"Cut Sequence: {item.sequence}, Cut start: {item.start}, Cut end: {item.end}\n"
		return output

	def _output_helper(self, parsed_protein):

		formated_output = ""

		formated_output = formated_output + f"Name: {parsed_protein.name}\nSequence:\n{wrapper.fill(parsed_protein.sequence)}\nMass: {parsed_protein.mass} {parsed_protein.massUnit}\n" \
											f"Tryptic Digests: \n"
		tryptic_list = self._wrap_list_helper(list=parsed_protein.trypticDigestion)
		formated_output = formated_output + tryptic_list
		formated_output = formated_output + "Glycosylation sites:"
		glycosylation_sites = self._wrap_dictionary_helper(dictionary=parsed_protein.glycosylationSites,
														   wrapper=wrapper, in_line=True)
		formated_output = formated_output + "\n" + glycosylation_sites
		formated_output = formated_output + "Isoforms:"
		Isoforms = self._wrap_dictionary_helper(dictionary=parsed_protein.isoforms, wrapper=wrapper,
												in_line=True)
		formated_output = formated_output + Isoforms
		formated_output = formated_output + "\n" + ">"*40 + "FINISH"+"<"*47 + "\n"

		return formated_output

	def Run_GUI(self):
		import re
		# reading config
		with open("config.json", "r") as read:
			cfg = json.load(read)

		self.statusbar_text.set("Working on your query")

		isoforms = {}
		sequence_counter = 1
		output = {}
		total_output = ""

		path = str(self.path_entry.get())
		if path == "":
			path = "gli.xml"
		elif path.endswith(".xml"):
			path = os.path.normpath(path)
		else:
			messagebox.showerror("Invalid file type entered", "Oops your path is not valid, you need .xml file to make a query.")
			raise EnvironmentError(f"Initializing parsing failed with file: {path}")

		entered_sequence = str(self.sequence_entry.get())
		entered_sequence = entered_sequence.split()
		for seq in entered_sequence:

			if re.match('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', seq.upper()):
				pprotein = protein(xmlfile=path, acc=seq.upper())
				for protein_number in pprotein.isoforms.keys():
					isoforms[protein_number] = pprotein.isoforms[protein_number].sequence
				output_protein = ProteinOutput(
					name=pprotein.name,
					sequence=pprotein.sequence,
					mass=pprotein.mass,
					glycosylationSites=pprotein.glycosylation,
					trypticDigestion=pprotein.tryptic_peptides,
					isoforms=isoforms
				)

				total_output += self._output_helper(output_protein)
				if self.json:
					output[output_protein.name] = output_protein.json_output()

			elif re.match('[ACDEFGHIKLMNPQRSTVWUY]+$', seq.upper()):

				pprotein = protein(xmlfile=path, sequence=seq.upper())
				output_protein = sequence_helper(pprotein, sequence_counter)
				sequence_counter += 1


				total_output += self._output_helper(output_protein)
				if self.json:
					output[output_protein.name] = output_protein.json_output()

			else:
				messagebox.showerror("Wrong sequence inserted", f"Sequence {seq} please review it before trying again")




		if self.text:
			with open(cfg["text_output"], 'w') as out:
				with redirect_stdout(out):
					print(total_output)

		elif self.json:
			with open(cfg["json_output"], 'w') as outfile:
				json.dump(output, outfile)

		output_box = Text(self.canvas_area, width=200, height=55, wrap=WORD, background="white")
		output_box.grid(row=180, column=12, columnspan=2, sticky=W)
		output_box.delete(0.0, END)
		output_box.insert(END, total_output)

		self.statusbar_text.set("Ready to begin new query")

ProteinGui()