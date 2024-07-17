from protein import protein
import matplotlib.pyplot as plt
import seaborn as sns
import math
#from contextlib import redirect_stdout
import json

with open("config.json", "r") as read:
	cfg = json.load(read)


sns.set()
plasma_proteins = []
input = open("plasmagp.list", "r")

input_content = input.read()
input_list = input_content.split()

for plasma_protein in input_list:
	protein_instance = protein(acc=plasma_protein,
							   xmlfile="C:\programing\dip\diplomski-rad\human_plasma_analysis\plasmagp.xml")
	plasma_proteins.append(protein_instance)


def mean(values):
	return sum(values) * 1.0 / len(values)


def stanDev(values):
	length = len(values)
	m = mean(values)
	total_sum = 0

	for i in range(length):
		total_sum += (values[i] - m) ** 2

	root = total_sum * 1.0 / length

	return math.sqrt(root)


def length_set(instance_list):
	sequences = []

	for plasma_protein in instance_list:
		sequences.append(len(plasma_protein.sequence))

	return sequences


def avg_hydrophobicity_set(instance_list):
	sequences = []

	for plasma_protein in instance_list:
		sequences.append(plasma_protein.avg_hydrophobicity)

	return sequences


def mass_set(instance_list):
	sequences = []

	for plasma_protein in instance_list:
		sequences.append(plasma_protein.mass)

	return sequences

def sequence_set(instance_list):
	sequences = []

	for plasma_protein in instance_list:
		sequences.append(plasma_protein.sequence)

	print (sequences)

def output(accs, length, hydrophobicity, mass):
	i = 0
	output = ""

	for leng in length:
		output = output + f"{accs[i]}: {length[i]} {hydrophobicity[i]} {mass[i]}\n"
		i += 1

	print(output)


#with open(cfg["text_output"], 'w') as out:
#	with redirect_stdout(out):
		#output(input_list, length_set(plasma_proteins), avg_hydrophobicity_set(plasma_proteins), mass_set(plasma_proteins))
#		sequence_set(plasma_proteins)

print(str(mean(length_set(plasma_proteins))) + '+/-' + str(stanDev(length_set(plasma_proteins))))

print(str(mean(avg_hydrophobicity_set(plasma_proteins))) + '+/-' + str(stanDev(avg_hydrophobicity_set(plasma_proteins))))

print(str(mean(mass_set(plasma_proteins))) + '+/-' + str(stanDev(mass_set(plasma_proteins))))
l_bins = [0, 250, 500, 750,  1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
h_bins = [-30000, -27500, -25000, -22500, -20000, -17500,-15000, -12500, -10000, -7500, -5000, -2500, -1000, 0]
m_bins = [20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000, 180000, 200000, 220000, 240000, 260000, 280000,
		  300000, 320000, 340000, 360000, 380000, 400000, 420000, 440000, 460000, 480000, 500000, 520000, 540000, 560000, 580000,
		  600000, 620000, 640000, 660000, 680000, 700000, 720000, 740000,  760000, 780000, 1500000, 1600000]


plt.hist(length_set(plasma_proteins), l_bins, rwidth=0.9)
plt.ylabel("Number of glycoproteins")
plt.xlabel("Glycoprotein length")
plt.title("Glycoprotein length distribution")
plt.show()

plt.hist(mass_set(plasma_proteins), m_bins, rwidth=0.9)
plt.ylabel("Number of glycoproteins")
plt.xlabel("Glycoprotein mass")
plt.title("Glycoprotein mass distribution")
plt.show()

plt.hist(avg_hydrophobicity_set(plasma_proteins),h_bins, rwidth=0.9)
plt.ylabel("Number of glycoproteins")
plt.xlabel("Average glycoprotein hydrophobicity")
plt.title("Glycoprotein hydrophobicity distribution")
plt.show()



