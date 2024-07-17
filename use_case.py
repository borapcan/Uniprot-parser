"""
User cases for diplomski-rad.

Importing necessary packages:

>>> from protein import protein


Protein initialization.
Form protein contains 3 optional arguments xmlfile (file to search for your protein of choice, if none is given default
one will be used), acc (protein accession number, which is used to find a protein of interest), sequence (for basic
custom sequence operations, mass calculation/trypsine digestion). If there are no acc and no sequence arguments user
will get an error message.
>>> pprotein = protein(xmlfile="gli.xml", acc="P52798")

Protein information stored:

Protein name:
>>> pprotein.name
'Ephrin-A4'

Protein cannonical sequence:
>>> pprotein.sequence
'MRLLPLLRTVLWAAFLGSPLRGGSSLRHVVYWNSSNPRLLRGDAVVELGLNDYLDIVCPHYEGPGPPEGPETFALYMVDWPGYESCQAEGPRAYKRWVCSLPFGHVQFSEKIQRFTPFSLGFEFLPGETYYYISVPTPESSGQCLRLQVSVCCKERKSESAHPVGSPGESGTSGWRGGDTPSPLCLLLLLLLLILRLLRIL'

Protein glycosylation sites :
>>> pprotein.glycosylation
{33: 'N-linked (GlcNAc...) asparagine'}

Protein isoforms (stored as a dictionary and can be accessed as such):

>>> pprotein.isoforms[2].sequence
'MRLLPLLRTVLWAAFLGSPLRGGSSLRHVVYWNSSNPRLLRGDAVVELGLNDYLDIVCPHYEGPGPPEGPETFALYMVDWPGYESCQAEGPRAYKRWVCSLPFGHVQFSEKIQRFTPFSLGFEFLPGETYYYISVPTPESSGQCLRLQVSVCCKERNLPSHPKEPESSQDPLEEEGSLLPALGVPIQTDKMEH'

>>> pprotein.mass
22386.357689999997

>>> pprotein.tryptic_peptides[1].sequence
'LLPLLR'

>>> pprotein.tryptic_peptides[1].start
3

>>> pprotein.tryptic_peptides[1].end
8

Isoform instance creation:
Creating and isoform instance ( requires isoform number argument, if the user doesn't know which are available they can
 always use .isoform_numbers to return a list of numbers)
>>> Isoform_2 = pprotein.isoforms[2]

Isoform class contains the same options as protein class with an exception of .isoforms
>>> Isoform_2.name
'Ephrin-A4-2'

"""
