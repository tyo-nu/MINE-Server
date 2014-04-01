__author__ = 'JGJeffryes'
from Client import mineDatabaseServices

services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
test_db = '1GenEcoCyc'
glucose = [{u'Inchi_key': u'WQZGKKKJIJFFOK-UHFFFAOYSA-N', u'KEGG_code': u'C00031', u'Mass': 180.063388104,
            u'Names': [u'D-Glucose', u'Grape sugar', u'Dextrose', u'Glucose'], u'Formula': u'C6H12O6',
            u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13'}]


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
    assert services.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N') == glucose
    assert services.quick_search(test_db, 'C00031') == glucose
    assert services.quick_search(test_db, 'Glucose') == glucose


def test_database_query():
    assert services.database_query('admin', '', '', False) == ['Illegal query']
    assert services.database_query(test_db, 'KEGG_code', 'C00031', False) == glucose
    assert services.database_query(test_db, 'Names', 'Grape', True) == glucose


def test_get_comps():
    print services.get_comps(test_db, ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13'])

def test_get_models():
    meh = services.get_models()
    assert len(meh[0]) == 2
    assert len(meh) == 2356


def test_get_adducts():
    meh = services.get_adducts()
    assert len(meh[0]) == 33
    assert len(meh[1]) == 15
    assert meh[0][2] == 'M+Na '


def test_adduct_db_search():
    meh = services.adduct_db_search(test_db, 164.0937301, 0.002, ['M+H'], [], False, True, False)
    assert len(meh) == 4
    assert len(meh[1]) == 3
    print len(meh[1][2])


def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert len(meh) == 1
    assert meh[0] == ['C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rab0dd7bc1c91b88c6f6ba90362413cb31fe00a42',
                         u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R4cfa8ce3f06297e2282b42ad69356815ee18d94f',
                         u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True)) == 9

