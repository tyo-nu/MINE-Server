__author__ = 'JGJeffryes'
from Impl import mineDatabaseServices


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
services = mineDatabaseServices(config)


class DB_query_params():
    def __init__(self, db, collection, field, value, regex=False):
        self.db = db
        self.collection = collection
        self.field = field
        self.value = value
        self.regex = regex


class Adduct_search_params():
    def __init__(self, db, mz, tolerance, adducts, charge):
        self.db = db
        self.mz = mz
        self.tolerance = tolerance
        self.adduct_list = adducts
        self.models = []
        self.ppm = False
        self.charge = charge
        self.halogens = False

class Pathway_query_params():
    def __init__(self, db, start, end, length, all_path):
        self.db = db
        self.start_comp = start
        self.end_comp = end
        self.len_limit = length
        self.all_paths = all_path
        self.np_min = -3
        self.gibbs_cap = 100

def test_quick_search():
    assert services.quick_search(config.test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']
    assert services.quick_search(config.test_db, 'C00031') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']
    assert services.quick_search(config.test_db, 'Glucose') == ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']


def test_database_query():
    assert services.database_query(DB_query_params('admin', 'users', '', '')) == [['Illegal query']]
    assert services.database_query(DB_query_params(config.test_db, 'compounds', 'KEGG_code', 'C00031')) == \
           [['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']]
    assert services.database_query(DB_query_params(config.test_db, 'compounds', 'Names', 'Grape', True)) == \
           [['Cb5b3273ab083d77ed29fbef8f7e464929af29c13']]


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
    meh = services.adduct_db_search(Adduct_search_params(config.test_db, 164.0937301, 0.002, ['M+H'], True))[0]
    assert len(meh) == 4
    assert len(meh[1]) == 3
    print len(meh[1][2])


def test_pathway_search():
    meh = services.pathway_search(Pathway_query_params(config.test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False))[0]
    assert len(meh) == 1
    assert meh[0] == ['C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rab0dd7bc1c91b88c6f6ba90362413cb31fe00a42',
                         u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R4cfa8ce3f06297e2282b42ad69356815ee18d94f',
                         u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(Pathway_query_params(config.test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                         'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True))[0]) == 9


