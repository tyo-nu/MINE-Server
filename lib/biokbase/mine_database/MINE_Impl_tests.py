__author__ = 'JGJeffryes'
from Impl import mineDatabaseServices

test_db = '1GenEcoCyc'
glucose = {u'Inchi_key': u'WQZGKKKJIJFFOK-UHFFFAOYSA-N', u'KEGG_code': u'C00031', u'Mass': 180.063388104,
            u'Names': [u'D-Glucose', u'Grape sugar', u'Dextrose', u'Glucose'], u'Formula': u'C6H12O6',
            u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13'}

class Options():
    def __init__(self):
        self.adduct_file = "lib/All Adducts.txt"
        self.kbase_db = 'KBase'
        self.test_db = '1GenEcoCyc'

        self.positive_adduct_file = "lib/Positive Adducts full.txt"
        self.negative_adduct_file = "lib/Negative Adducts full.txt"
        self.adduct_list = ['M+H', 'M-H']
        self.tolerance = 0.002
        self.ppm = False
        self.ms_file = "lib/Unknowns.txt"
        self.halogens = True

config = Options()
services = mineDatabaseServices(None)

def test_quick_search():
    assert services.quick_search(config.test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']
    assert services.quick_search(config.test_db, 'C00031') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']
    assert services.quick_search(config.test_db, 'Glucose') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']


def test_database_query():
    assert services.database_query('admin', '', '', False) == ['Illegal query']
    assert services.database_query(test_db, 'KEGG_code', 'C00031', False) == [glucose]
    assert services.database_query(test_db, 'Names', 'Grape', True) == [glucose]


def test_get_comps():
    print services.get_comps(config.test_db, ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13'])


def test_get_models():
    meh = services.get_models()[0]
    assert meh[2] == (u'kb|fm.3375', u'Escherichia coli 97.0264')
    assert len(meh) == 2356


def test_get_adducts():
    meh = services.get_adducts()[0]
    assert len(meh[0]) == 33
    assert len(meh[1]) == 30
    assert meh[0][2] == 'M+Na '


def test_adduct_db_search():
    meh = services.adduct_db_search(test_db, 164.0937301, 0.002, ['M+H'], [], False, True, False)
    assert len(meh) == 4
    assert len(meh[1]) == 3


def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert len(meh) == 1
    assert meh[0] == ['C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rab0dd7bc1c91b88c6f6ba90362413cb31fe00a42',
                         u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R4cfa8ce3f06297e2282b42ad69356815ee18d94f',
                         u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True)) == 9

def test_similarity_search():
    meh = services.similarity_search('EcoCycexp', 'O=C1CC(OC1COP(=O)(OP(=O)(O)O)O)n1cc(C)c(nc1=O)O', 0.8)
    print meh