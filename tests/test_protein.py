import pytest

from protein import protein

@pytest.fixture(scope="module")
def entry():
	entry = protein(xmlfile="gli.xml", acc="Q9UKJ0")
	return entry


@pytest.fixture(scope="module")
def Q9UKJ0_1():
	Q9UKJ0_1_sequence = """
			MGRPLLLPLLLLLQPPAFLQPGGSTGSGPSYLYGVTQPKHLSASMGGSVEIPFSFYYPWE
			LAIVPNVRISWRRGHFHGQSFYSTRPPSIHKDYVNRLFLNWTEGQESGFLRISNLRKEDQ
			SVYFCRVELDTRRSGRQQLQSIKGTKLTITQAVTTTTTWRPSSTTTIAGLRVTESKGHSE
			SWHLSLDTAIRVALAVAVLKTVILGLLCLLLLWWRRRKGSRAPSSDF
			"""
	Q9UKJ0_1_sequence = ''.join(Q9UKJ0_1_sequence.split())
	return Q9UKJ0_1_sequence


@pytest.fixture(scope="module")
def Q9UKJ0_2():
	Q9UKJ0_2_sequence = """
			MGRPLLLPLLLLLQPPAFLQPGLCEPALSELDRGSGERLPQDLKPAEGGPVCVFLPSRAG
			HPEIREAAVAVHQGDQTHHHPGCHNHHHLEAQQHNHHSRPQGHRKQRALRIMAPKSGHCH
			QGCIGCRCAQNCHFGTAVPPPPVVEEKER
			"""
	Q9UKJ0_2_sequence = ''.join(Q9UKJ0_2_sequence.strip().split())
	return Q9UKJ0_2_sequence


@pytest.fixture(scope="module")
def Q9UKJ0_3():
	Q9UKJ0_3_sequence = """
			MGRASYGGSPVGHQSSSDHPRAKTCRSPVRDGPRTLCGVPTVALSLHFLREASSGSRTCG
			RRTSLCTSAKSSWTYRSGRLSWQSIKGTHLTITQALRQPLHRAPLLPGQLCWSPRPLEKN
			KAMGRPLLLPLLLLLQPPAFLQPGLCEPALSELDRGSGERLPQDLKPAEGGPVCVFLPSR
			AGHPEIREAAVAVHQGDQTHHHPGCHNHHHLEAQQHNHHSRPQGHRKQRALRIMAPKSGH
			CHQGCIGCRCAQNCHFGTAVPPPPVVEEKER
			"""
	Q9UKJ0_3_sequence = ''.join(Q9UKJ0_3_sequence.strip().split())
	return Q9UKJ0_3_sequence


@pytest.fixture(scope="module")
def trypsine_cut_1():
	trypsine_cut_1 = ['MGRPLLLPLLLLLQPPAFLQPGGSTGSGPSYLYGVTQPK', 1, 39]
	return trypsine_cut_1

def test_protein_name(entry):
	name = entry.name
	assert name == "Paired immunoglobulin-like type 2 receptor beta"


def test_canonical_sequence(Q9UKJ0_1,entry):
	sequence = entry.sequence
	assert sequence == Q9UKJ0_1


def test_glycosylation_sites(entry):
	glycosylation = entry.glycosylation
	assert glycosylation == {100: 'N-linked (GlcNAc...) asparagine'}



def test_alternative_sequence_isoform_changes_ids(entry):
	ids_data = entry.isoforms[3]._alternative_products()
	assert ids_data == {2: ['VSP_017504'], 3: ['VSP_017503', 'VSP_017504']}


def test_isoform_data_changes_information(entry):
	isoform_data = entry.isoforms[3]._vsps()
	assert isoform_data == {
		'VSP_017503': ('MGRASYGGSPVGHQSSSDHPRAKTCRSPVRDGPRTLCGVPTVALSLHFLREASSGSRTCGRRTSLCTSAKSSWTYRSGRLSWQSIKGTHLTITQALRQPLHRAPLLPGQLCWSPRPLEKNKAM', 1, 1),
		'VSP_017504': ('LCEPALSELDRGSGERLPQDLKPAEGGPVCVFLPSRAGHPEIREAAVAVHQGDQTHHHPGCHNHHHLEAQQHNHHSRPQGHRKQRALRIMAPKSGHCHQGCIGCRCAQNCHFGTAVPPPPVVEEKER', 23, 227)
	}

def test_single_vsp(Q9UKJ0_2,entry):
	assert entry.isoforms[2].sequence == Q9UKJ0_2

def test_multiple_vsp(Q9UKJ0_3,entry):
	assert entry.isoforms[3].sequence == Q9UKJ0_3

def test_isoform_glycosylation(entry):
	assert entry.isoforms[3].glycosylation == {}

def test_mass(entry):
	assert entry.mass == 25542.989100000003


def test_trypsine_cutter(trypsine_cut_1,entry):
	trypsine_cut_list = [entry.tryptic_peptides[0].sequence, entry.tryptic_peptides[0].start, entry.tryptic_peptides[0].end]
	assert trypsine_cut_list == trypsine_cut_1

def test_average_hydrophobacity(entry):
	assert entry.avg_hydrophobicity == -920.0899999999993


