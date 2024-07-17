import mass_calculator

def test_mass_calculator_canonical():
	sequence = mass_calculator.mass_calculator_canonical("P02751")
	assert sequence == 261628.8713
	
	
def test_mass_calculator_tryptic():
	dict = mass_calculator.mass_calculator_tryptic("P02751")
	assert dict == {'MLR': 402.55389999999994, 'GPGPGLLLLAVQCLGTAVPSTGASK': 2323.7006, 'QAQQMVQPQSPVAVSQSKPGCYDNGK': 2792.0584, 'HYQINQQWER': 1385.4823999999999, 'TYLGNALVCTCYGGSR': 1693.8941000000002, 'GFNCESKPEAEETCFDK': 1966.0589000000002, 'YTGNTYR': 857.9076, 'VGDTYERPK': 1048.147, 'DSMIWDCTCIGAGR': 1543.74, 'ISCTIANR': 877.0183, 'CHEGGQSYK': 1024.0630999999998, 'IGDTWR': 730.8097, 'RPHETGGYMLECVCLGNGK': 2080.3628, 'GEWTCKPIAEKCFDHAAGTSYVVGETWEKPYQGWMMVDCTCLGEGSGRITCTSRNRCNDQDTRTSYRIGDTWSKKDNRGNLLQCICTGNGRGEWKCERHTSVQTTSSGSGPFTDVRAAVYQPQPHPQPPPYGHCVTDSGVVYSVGMQWLKTQGNKQMLCTCLGNGVSCQETAVTQTYGGNSNGEPCVLPFTYNGRTFYSCTTEGRQDGHLWCSTTSNYEQDQKYSFCTDHTVLVQTRGGNSNGALCHFPFLYNNHNYTDCTSEGRRDNMKWCGTTQNYDADQKFGFCPMAAHEEICTTNEGVMYRIGDQWDKQHDMGHMMRCTCVGNGRGEWTCIAYSQLRDQCIVDDITYNVNDTFHKRHEEGHMLNCTCFGQGRGRWKCDPVDQCQDSETGTFYQIGDSWEKYVHGVRYQCYCYGRGIGEWHCQPLQTYPSSSGPVEVFITETPSQPNSHPIQWNAPQPSHISKYILRWRPKNSVGRWKEATIPGHLNSYTIKGLKPGVVYEGQLISIQQYGHQEVTRFDFTTTSTSTPVTSNTVTGETTPFSPLVATSESVTEITASSFVVSWVSASDTVSGFRVEYELSEEGDEPQYLDLPSTATSVNIPDLLPGRKYIVNVYQISEDGEQSLILSTSQTTAPDAPPDTTVDQVDDTSIVVRWSRPQAPITGYRIVYSPSVEGSSTELNLPETANSVTLSDLQPGVQYNITIYAVEENQESTPVVIQQETTGTPRSDTVPSPRDLQFVEVTDVKVTIMWTPPESAVTGYRVDVIPVNLPGEHGQRLPISRNTFAEVTGLSPGVTYYFKVFAVSHGRESKPLTAQQTTKLDAPTNLQFVNETDSTVLVRWTPPRAQITGYRLTVGLTRRGQPRQYNVGPSVSKYPLRNLQPASEYTVSLVAIKGNQESPKATGVFTTLQPGSSIPPYNTEVTETTIVITWTPAPRIGFKLGVRPSQGGEAPREVTSDSGSIVVSGLTPGVEYVYTIQVLRDGQERDAPIVNKVVTPLSPPTNLHLEANPDTGVLTVSWERSTTPDITGYRITTTPTNGQQGNSLEEVVHADQSSCTFDNLSPGLEYNVSVYTVKDDKESVPISDTIIPAVPPPTDLRFTNIGPDTMRVTWAPPPSIDLTNFLVRYSPVKNEEDVAELSISPSDNAVVLTNLLPGTEYVVSVSSVYEQHESTPLRGRQKTGLDSPTGIDFSDITANSFTVHWIAPRATITGYRIRHHPEHFSGRPREDRVPHSRNSITLTNLTPGTEYVVSIVALNGREESPLLIGQQSTVSDVPRDLEVVAATPTSLLISWDAPAVTVRYYRITYGETGGNSPVQEFTVPGSKSTATISGLKPGVDYTITVYAVTGRGDSPASSKPISINYRTEIDKPSQMQVTDVQDNSISVKWLPSSSPVTGYRVTTTPKNGPGPTKTKTAGPDQTEMTIEGLQPTVEYVVSVYAQNPSGESQPLVQTAVTNIDRPKGLAFTDVDVDSIKIAWESPQGQVSRYRVTYSSPEDGIHELFPAPDGEEDTAELQGLRPGSEYTVSVVALHDDMESQPLIGTQSTAIPAPTDLKFTQVTPTSLSAQWTPPNVQLTGYRVRVTPKEKTGPMKEINLAPDSSSVVVSGLMVATKYEVSVYALKDTLTSRPAQGVVTTLENVSPPRRARVTDATETTITISWRTKTETITGFQVDAVPANGQTPIQRTIKPDVRSYTITGLQPGTDYKIYLYTLNDNARSSPVVIDASTAIDAPSNLRFLATTPNSLLVSWQPPRARITGYIIKYEKPGSPPREVVPRPRPGVTEATITGLEPGTEYTIYVIALKNNQKSEPLIGRKKTDELPQLVTLPHPNLHGPEILDVPSTVQKTPFVTHPGYDTGNGIQLPGTSGQQPSVGQQMIFEEHGFRRTTPPTTATPIRHRPRPYPPNVGEEIQIGHIPREDVDYHLYPHGPGLNPNASTGQEALSQTTISWAPFQDTSEYIISCHPVGTDEEPLQFRVPGTSTSATLTGLTRGATYNVIVEALKDQQRHKVREEVVTVGNSVNEGLNQPTDDSCFDPYTVSHYAVGDEWERMSESGFKLLCQCLGFGSGHFRCDSSRWCHDNGVNYKIGEKWDRQGENGQMMSCTCLGNGKGEFKCDPHEATCYDDGKTYHVGEQWQKEYLGAICSCTCFGGQRGWRCDNCRRPGGEPSPEGTTGQSYNQYSQRYHQRTNTNVNCPIECFMPLDVQADREDSRE': 242584.60039999997}
	
		
def test_mass_calculator_isoform():
	dictionary = mass_calculator.mass_calculator_isoforms("P02751")
	assert dictionary == {'2': 77005.8539, '3': 258236.23400000003, '4': 248808.01240000004, '5': 255486.25780000002, '6': 239657.65759999998, '7': 258236.23400000003, '8': 251879.30940000003, '9': 258083.18930000003, '10': 248598.70759999997, '11': 261442.62370000003, '12': 261628.8713, '13': 248598.70759999997, '14': 258192.04970000003, '15': 271292.5156, '16': 74852.3404, '17': 258270.32439999998}