from protein import protein
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math
import statistics

sns.set()
plt.style.use('seaborn-bright')
plasma_proteins = []
input = open("plasma.list", "r")

input_content = input.read()
input_list = input_content.split()

for plasma_protein in input_list:
	protein_instance = protein(acc=plasma_protein,
							   xmlfile="human_plasma_analysis/plasma.xml")
	plasma_proteins.append(protein_instance)

fig, axes = plt.subplots()

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
	hydrophobicities = []

	for plasma_protein in instance_list:
		hydrophobicities.append(plasma_protein.avg_hydrophobicity)

	return hydrophobicities


def mass_set(instance_list):
	masses = []

	for plasma_protein in instance_list:
		masses.append(plasma_protein.mass)

	return masses

def sequence_set(instance_list):
	sequences = []

	for plasma_protein in instance_list:
		sequences.append(plasma_protein.sequence)

	print(sequences)

def glycosylation_analysis(list):
	N = 0
	O = 0
	C = 0
	S = 0
	G = 0
	T = 0

	for protein in list:
		for glycosylation in protein.glycosylation.values():
			T += 1
			if 'N-linked' in glycosylation:
				N += 1
			elif 'O-linked' in glycosylation:
				O += 1
			elif 'C-linked' in glycosylation:
				C += 1
			elif 'S-linked' in glycosylation:
				S += 1
			elif 'Glycation' in glycosylation:
				G += 1

	glycosylations = [N, O, C, S, G, T]
	return glycosylations



def peptide_set_len(instance_list):
	peptides = []
	N_glyco_peptides = []
	O_glyco_peptides = []
	C_glyco_peptides = []

	for protein in instance_list:
		for glycosylation in protein.glycosylation.keys():
			peptides.append(glycosylation)
			for peptide in protein.tryptic_peptides:
				if int(glycosylation) > peptide.start-1 and glycosylation < peptide.end:
					if 'N-linked' in protein.glycosylation[glycosylation]:
						N_glyco_peptides.append(len(peptide.sequence))
					elif 'O-linked' in protein.glycosylation[glycosylation]:
						O_glyco_peptides.append(len(peptide.sequence))
					elif 'C-linked'in protein.glycosylation[glycosylation]:
						C_glyco_peptides.append(len(peptide.sequence))
				else:
					peptides.append(len(peptide.sequence))

	list_peptides = [N_glyco_peptides, O_glyco_peptides, C_glyco_peptides, peptides]
	return list_peptides

def peptide_set_hydro(instance_list):
	peptides = []
	N_glyco_peptides = []
	O_glyco_peptides = []
	C_glyco_peptides = []

	for protein in instance_list:
		for glycosylation in protein.glycosylation.keys():
			peptides.append(glycosylation)
			for peptide in protein.tryptic_peptides:
				if int(glycosylation) > peptide.start-1 and glycosylation < peptide.end:
					if 'N-linked' in protein.glycosylation[glycosylation]:
						N_glyco_peptides.append(peptide.avg_hydrophobicity)
					elif 'O-linked' in protein.glycosylation[glycosylation]:
						O_glyco_peptides.append(peptide.avg_hydrophobicity)
					elif 'C-linked'in protein.glycosylation[glycosylation]:
						C_glyco_peptides.append(peptide.avg_hydrophobicity)
				else:
					peptides.append(peptide.avg_hydrophobicity)

	list_peptides = [N_glyco_peptides, O_glyco_peptides, C_glyco_peptides, peptides]
	return list_peptides

def peptide_set_mass(instance_list):
	peptides = []
	N_glyco_peptides = []
	O_glyco_peptides = []
	C_glyco_peptides = []

	for protein in instance_list:
		for glycosylation in protein.glycosylation.keys():
			peptides.append(glycosylation)
			for peptide in protein.tryptic_peptides:
				if int(glycosylation) > peptide.start-1 and glycosylation < peptide.end:
					if 'N-linked' in protein.glycosylation[glycosylation]:
						N_glyco_peptides.append(peptide.mass)
					elif 'O-linked' in protein.glycosylation[glycosylation]:
						O_glyco_peptides.append(peptide.mass)
					elif 'C-linked'in protein.glycosylation[glycosylation]:
						C_glyco_peptides.append(peptide.mass)
				else:
					peptides.append(peptide.mass)

	list_peptides = [N_glyco_peptides, O_glyco_peptides, C_glyco_peptides, peptides]
	return list_peptides


#data_dict_len = {"N-glikozilirani":peptide_set_len(plasma_proteins)[0], "O-glikozilirani":peptide_set_len(plasma_proteins)[1],
#					"C-glikozilirani":peptide_set_len(plasma_proteins)[2], "Neglikozilirani":peptide_set_len(plasma_proteins)[3]}

#print(len(peptide_set_len(plasma_proteins)[0]), len(peptide_set_len(plasma_proteins)[1]), len(peptide_set_len(plasma_proteins)[2]), len(peptide_set_len(plasma_proteins)[3]))

#dfl = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in data_dict_len.items()]))
#len_plot = sns.violinplot(data=dfl, showfliers=False)
#plt.show()

ngl = peptide_set_len(plasma_proteins)[0]
ngl.append(20)

ogl = peptide_set_len(plasma_proteins)[1]*6
medstacko=[14]*225
ogl +=medstacko

cgl = peptide_set_len(plasma_proteins)[2]*45
medstackc=[25]*36
cgl += medstackc

i=0
negl=[]
neg = peptide_set_len(plasma_proteins)[3]
list.sort(neg)
while i < 222001:
	negl.append(neg[i])
	i += 124

dfl_1791 = [ngl, ogl, cgl, negl]

i=0
ngls=[]
ngl = peptide_set_len(plasma_proteins)[0]
list.sort(ngl)
while i < 1790:
	negl.append(neg[i])
	i += 850

ogl = peptide_set_len(plasma_proteins)[1]
ogl.append(14)

cgl = peptide_set_len(plasma_proteins)[2]*6
medstackc=[25]*28
cgl += medstackc

i=0
negl=[]
neg = peptide_set_len(plasma_proteins)[3]
list.sort(neg)
while i < 222001:
	negl.append(neg[i])
	i += 850

dfl_262 = [ngl, ogl, cgl, negl]

#medstackn=[20]*41
#n_len=peptide_set_len(plasma_proteins)[0]*124
#n_len+=medstackn

#medstacko=[14]*151
#o_len=peptide_set_len(plasma_proteins)[1]*850
#o_len+=medstacko

#medstackc=[25]*13
#c_len=peptide_set_len(plasma_proteins)[2]*5692
#c_len+=medstackc

#dfl = [n_len, o_len, c_len, peptide_set_len(plasma_proteins)[3]]

plt.violinplot(dfl_1791, showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.ylabel("Duljina (broj aminokiselina)")
#plt.xlabel("Mase glikoproteina")
#plt.title("Violinski graf  distribucije duljina svih analiziranih tipova glikozilacije peptida")
plt.show()

plt.violinplot(dfl_262, showmeans=True, showmedians=False, showextrema=False, widths=0.6)
plt.show()

#nlen_plot = sns.violinplot(data=peptide_set_len(plasma_proteins)[0])
#plt.show()
#axes.violinplot(peptide_set_len(plasma_proteins)[0], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()
#olen_plot = sns.violinplot(data=peptide_set_len(plasma_proteins)[1])
#plt.show()
#slen_plot = sns.violinplot(data=peptide_set_len(plasma_proteins)[2])
#plt.show()x
#plen_plot = sns.violinplot(data=peptide_set_len(plasma_proteins)[3])
#plt.show()
#plt.violinplot(peptide_set_len(plasma_proteins)[3], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()

#data_dict_mass = {"N-glikozilirani":peptide_set_mass(plasma_proteins)[0], "O-glikozilirani":peptide_set_mass(plasma_proteins)[1],
#					"C-glikozilirani":peptide_set_mass(plasma_proteins)[2], "Neglikozilirani":peptide_set_mass(plasma_proteins)[3]}

#dfm = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in data_dict_mass.items()]))
#mass_plot = sns.violinplot(data=dfm, showfliers=False)
#plt.show()

#mmedstackn=[2653.9392591564247]*41
#n_mass=peptide_set_mass(plasma_proteins)[0]*124
#n_mass+=mmedstackn

#mmedstacko=[2243.0259499233716]*151
#o_mass=peptide_set_mass(plasma_proteins)[1]*850
#o_mass+=mmedstacko

#mmedstackc=[3212.992697179487]*13
#c_mass=peptide_set_mass(plasma_proteins)[2]*5692
#c_mass+=mmedstackc

#dfm = [n_mass, o_mass, c_mass, peptide_set_mass(plasma_proteins)[3]]

#plt.violinplot(dfm, showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.ylabel("Masa(Da)")
#plt.xlabel("Mase glikoproteina")
#plt.title("Violinski graf  distribucije masa svih analiziranih tipova glikozilacije peptida")
#plt.show()


#nmass_plot = sns.violinplot(data=peptide_set_mass(plasma_proteins)[0])
#plt.show()
#axes.violinplot(peptide_set_mass(plasma_proteins)[0], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()
#omass_plot = sns.violinplot(data=peptide_set_mass(plasma_proteins)[1])
#plt.show()
#smass_plot = sns.violinplot(data=peptide_set_mass(plasma_proteins)[2])
#plt.show()
#pmass_plot = sns.violinplot(data=peptide_set_mass(plasma_proteins)[3])
#plt.show()
#axes.violinplot(peptide_set_mass(plasma_proteins)[3], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()

#data_dict_hydro = {"N-glikozilirani":peptide_set_hydro(plasma_proteins)[0], "O-glikozilirani":peptide_set_hydro(plasma_proteins)[1],
#						"C-glikozilirani":peptide_set_hydro(plasma_proteins)[2], "Neglikozilirani":peptide_set_hydro(plasma_proteins)[3]}

#dfh = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in data_dict_hydro.items()]))


#hmedstackn=[-92.0769217877095]*41
#n_hydro=peptide_set_hydro(plasma_proteins)[0]*124
#n_hydro += hmedstackn


#hmedstacko=[-72.512030651341]*151
#o_hydro=peptide_set_hydro(plasma_proteins)[1]*850
#o_hydro+=hmedstacko


#hmedstackc=[-106.60897435897436]*13
#c_hydro=peptide_set_hydro(plasma_proteins)[2]*5692
#c_hydro+=hmedstackc

#dfh = [n_hydro, o_hydro, c_hydro, peptide_set_hydro(plasma_proteins)[3]]

#hydro_plot = sns.violinplot(data=dfh, showfliers=False)
#plt.show()

#plt.violinplot(dfh, showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.ylabel("Hidrofobnost(Kcal pri 25 °C)")
#plt.xlabel("Mase glikoproteina")
#plt.title("Violinski graf  distribucije hidrofobnosti svih analiziranih tipova glikozilacije peptida")
#plt.show()

#nhydro_plot = sns.violinplot(data=peptide_set_hydro(plasma_proteins)[0])
#plt.show()
#axes.violinplot(peptide_set_hydro(plasma_proteins)[0], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()
#ohydro_plot = sns.violinplot(data=peptide_set_hydro(plasma_proteins)[1])
#plt.show()
#shydro_plot = sns.violinplot(data=peptide_set_hydro(plasma_proteins)[2])
#plt.show()
#phydro_plot = sns.violinplot(data=peptide_set_hydro(plasma_proteins)[3])
#plt.show()
#axes.violinplot(peptide_set_hydro(plasma_proteins)[3], showmeans=True, showmedians=False, showextrema=False, widths=0.6)
#plt.show()



#print(len(plasma_proteins))

#print(glycosylation_analysis(plasma_proteins))

#print(str(statistics.mean(length_set(plasma_proteins))) + '+/-' + str(statistics.stdev(length_set(plasma_proteins))))

#print(str(statistics.mean(avg_hydrophobicity_set(plasma_proteins))) + '+/-' + str(statistics.stdev(avg_hydrophobicity_set(plasma_proteins))))

#print(str(statistics.mean(mass_set(plasma_proteins))) + '+/-' + str(statistics.stdev(mass_set(plasma_proteins))))



l_bins = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500,
		  4750, 5000, 5250, 5500]

m_bins = [0, 20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000, 180000, 200000, 220000, 240000, 260000, 280000,
		  300000, 320000, 340000, 360000, 380000, 400000, 420000, 440000, 460000, 480000, 500000, 520000, 540000, 560000, 580000]

h_bins = [-21000, -20000, -19000, -18000, -17000, -16000, -15000, -14000, -13000, -12000, -11000, -10000, -9000, -8000,
		  -7000, -6000, -5000, -4000, -3000, -2000, -1000, 0]


#plt.hist(length_set(plasma_proteins), l_bins, rwidth=0.9)
#plt.ylabel("Broj glikoproteina")
#plt.xlabel("Duljine glikoproteina")
#plt.title("Histogram distribucije duljina glikoproteina humane plazme")
#plt.show()

#plt.hist(mass_set(plasma_proteins), m_bins, rwidth=0.9)
#plt.ylabel("Broj glikoproteina")
#plt.xlabel("Mase glikoproteina")
#plt.title("Histogram distribucije masa glikoproteina humane plazme")
#plt.show()

#plt.hist(avg_hydrophobicity_set(plasma_proteins), h_bins, rwidth=0.9)
#plt.ylabel("Broj glikoproteina")
#plt.xlabel("Prosječna hidrofobnost glikoproteinskog peptida")
#plt.title("Histogram distribucije prosječnih hidrofobnosti peptida glikoproteina humane plazme")
#plt.show()
