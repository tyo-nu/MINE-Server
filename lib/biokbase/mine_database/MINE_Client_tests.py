__author__ = 'JGJeffryes'
from Client import mineDatabaseServices

services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
test_db = 'EcoCycexp'
glucose = {u'Formula': u'C6H12O6', u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13',
           u'Names': [u'D-Glucose', u'Grape sugar', u'Dextrose', u'Glucose'], u'Model_SEED': 17398}


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
    assert services.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N') == [glucose]
    #assert services.quick_search(test_db, 'C00031') == [glucose]
    assert glucose in services.quick_search(test_db, 'Glucose')


def test_database_query():
    assert services.database_query('admin', '', '', False) == ['Illegal query']
    assert services.database_query(test_db, 'DB Links.PubChem', '3333', False) == [glucose]
    assert services.database_query(test_db, 'Names', 'Grape', True) == [glucose]


def test_get_comps():
    meh = services.get_comps(test_db, ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13'])
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()


def test_get_rxns():
    meh = services.get_rxns(test_db, ['R2e28f382545b6b00b88bcd5b1bb927bc480c2711'])
    assert len(meh) == 1
    assert 'Operators' in meh[0].keys()


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
    assert len(meh) == 3
    assert isinstance(meh[1]['isomers'], list)


def test_batch_ms_adduct_search():
    result = services.batch_ms_adduct_search(test_db, "164.0937301", "form", 0.002, ['M+H'], [], False, True, False)[0]
    meh = result['adducts']
    assert len(meh) == 3
    assert isinstance(meh[1]['isomers'], list)


def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert len(meh) == 1
    assert meh[0] == [u'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rbbc40c762b05d59890c196c949522e3ee6ca08c6',
                      u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R0e87f4c178bcb78b9190938b3413b7889b3fbad4',
                      u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True)) == 4


def test_similarity_search():
    assert len(services.similarity_search('EcoCycexp', 'OCC1OC(O)C(C(C1O)O)O', 0.9, "FP2")) == 28
    assert len(services.similarity_search('EcoCycexp', 'OCC1OC(O)C(C(C1O)O)O', 0.9, "FP4")) == 46
