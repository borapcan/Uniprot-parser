import tryptic_cutter

def test_tryptic_cutter ():
	tryptic_list_final = tryptic_cutter.tryptic_cutter("P02751")
	assert tryptic_list_final == ['MLR', 'GPGPGLLLLAVQCLGTAVPSTGASK', 'QAQQMVQPQSPVAVSQSKPGCYDNGK', 'HYQINQQWER', 'TYLGNALVCTCYGGSR', 'GFNCESKPEAEETCFDK', 'YTGNTYR', 'VGDTYERPK', 'DSMIWDCTCIGAGR', 'ISCTIANR', 'CHEGGQSYK', 'IGDTWR', 'RPHETGGYMLECVCLGNGK', 'GEWTCKPIAEKCFDHAAGTSYVVGETWEKPYQGWMMVDCTCLGEGSGRITCTSRNRCNDQDTRTSYRIGDTWSKKDNRGNLLQCICTGNGRGEWKCERHTSVQTTSSGSGPFTDVRAAVYQPQPHPQPPPYGHCVTDSGVVYSVGMQWLKTQGNKQMLCTCLGNGVSCQETAVTQTYGGNSNGEPCVLPFTYNGRTFYSCTTEGRQDGHLWCSTTSNYEQDQKYSFCTDHTVLVQTRGGNSNGALCHFPFLYNNHNYTDCTSEGRRDNMKWCGTTQNYDADQKFGFCPMAAHEEICTTNEGVMYRIGDQWDKQHDMGHMMRCTCVGNGRGEWTCIAYSQLRDQCIVDDITYNVNDTFHKRHEEGHMLNCTCFGQGRGRWKCDPVDQCQDSETGTFYQIGDSWEKYVHGVRYQCYCYGRGIGEWHCQPLQTYPSSSGPVEVFITETPSQPNSHPIQWNAPQPSHISKYILRWRPKNSVGRWKEATIPGHLNSYTIKGLKPGVVYEGQLISIQQYGHQEVTRFDFTTTSTSTPVTSNTVTGETTPFSPLVATSESVTEITASSFVVSWVSASDTVSGFRVEYELSEEGDEPQYLDLPSTATSVNIPDLLPGRKYIVNVYQISEDGEQSLILSTSQTTAPDAPPDTTVDQVDDTSIVVRWSRPQAPITGYRIVYSPSVEGSSTELNLPETANSVTLSDLQPGVQYNITIYAVEENQESTPVVIQQETTGTPRSDTVPSPRDLQFVEVTDVKVTIMWTPPESAVTGYRVDVIPVNLPGEHGQRLPISRNTFAEVTGLSPGVTYYFKVFAVSHGRESKPLTAQQTTKLDAPTNLQFVNETDSTVLVRWTPPRAQITGYRLTVGLTRRGQPRQYNVGPSVSKYPLRNLQPASEYTVSLVAIKGNQESPKATGVFTTLQPGSSIPPYNTEVTETTIVITWTPAPRIGFKLGVRPSQGGEAPREVTSDSGSIVVSGLTPGVEYVYTIQVLRDGQERDAPIVNKVVTPLSPPTNLHLEANPDTGVLTVSWERSTTPDITGYRITTTPTNGQQGNSLEEVVHADQSSCTFDNLSPGLEYNVSVYTVKDDKESVPISDTIIPAVPPPTDLRFTNIGPDTMRVTWAPPPSIDLTNFLVRYSPVKNEEDVAELSISPSDNAVVLTNLLPGTEYVVSVSSVYEQHESTPLRGRQKTGLDSPTGIDFSDITANSFTVHWIAPRATITGYRIRHHPEHFSGRPREDRVPHSRNSITLTNLTPGTEYVVSIVALNGREESPLLIGQQSTVSDVPRDLEVVAATPTSLLISWDAPAVTVRYYRITYGETGGNSPVQEFTVPGSKSTATISGLKPGVDYTITVYAVTGRGDSPASSKPISINYRTEIDKPSQMQVTDVQDNSISVKWLPSSSPVTGYRVTTTPKNGPGPTKTKTAGPDQTEMTIEGLQPTVEYVVSVYAQNPSGESQPLVQTAVTNIDRPKGLAFTDVDVDSIKIAWESPQGQVSRYRVTYSSPEDGIHELFPAPDGEEDTAELQGLRPGSEYTVSVVALHDDMESQPLIGTQSTAIPAPTDLKFTQVTPTSLSAQWTPPNVQLTGYRVRVTPKEKTGPMKEINLAPDSSSVVVSGLMVATKYEVSVYALKDTLTSRPAQGVVTTLENVSPPRRARVTDATETTITISWRTKTETITGFQVDAVPANGQTPIQRTIKPDVRSYTITGLQPGTDYKIYLYTLNDNARSSPVVIDASTAIDAPSNLRFLATTPNSLLVSWQPPRARITGYIIKYEKPGSPPREVVPRPRPGVTEATITGLEPGTEYTIYVIALKNNQKSEPLIGRKKTDELPQLVTLPHPNLHGPEILDVPSTVQKTPFVTHPGYDTGNGIQLPGTSGQQPSVGQQMIFEEHGFRRTTPPTTATPIRHRPRPYPPNVGEEIQIGHIPREDVDYHLYPHGPGLNPNASTGQEALSQTTISWAPFQDTSEYIISCHPVGTDEEPLQFRVPGTSTSATLTGLTRGATYNVIVEALKDQQRHKVREEVVTVGNSVNEGLNQPTDDSCFDPYTVSHYAVGDEWERMSESGFKLLCQCLGFGSGHFRCDSSRWCHDNGVNYKIGEKWDRQGENGQMMSCTCLGNGKGEFKCDPHEATCYDDGKTYHVGEQWQKEYLGAICSCTCFGGQRGWRCDNCRRPGGEPSPEGTTGQSYNQYSQRYHQRTNTNVNCPIECFMPLDVQADREDSRE']