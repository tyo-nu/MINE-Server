__author__ = 'JGJeffryes'
from lib.biokbase.mine_database.Client import mineDatabaseServices

services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
test_db = 'EcoCycexp2'
test_molfile = open("./scripts/xanthine.mol", "r").read()
glucose = {u'Formula': u'C6H12O6', u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13',
           u'Names': [u'Glucose', u'beta-D-Galactose', u'D-gulo-Hexose', u'D-Mannose', u'D-Idose', u'D-talo-Hexose',
                      u'D-Gulose', u'1,4-beta-D-Mannooligosaccharide', u'beta-D-Mannose', u'D-Aldose2', u'D-Allose',
                      u'D-Aldose1', u'D-Hexose', u'alpha-D-Glucose', u'D-Glucose', u'beta-D-Glucose', u'Grape sugar',
                      u'Mannose', u'Dextrose', u'D-Altrose', u'Seminose', u'D-ido-Hexose', u'alpha-D-Galactose',
                      u'D-altro-Hexose', u'Carubinose', u'L-Galactose', u'L-Gulose', u'D-galactopyranose', u'D-Talose',
                      u'alpha-D-Mannose', u'D-glucopyranose', u'Hexose', u'D-Aldose', u'D-Galactose'], u'Model_SEED': 17398}

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
    assert services.quick_search(test_db, 'C00031') == [glucose]
    assert glucose in services.quick_search(test_db, 'Glucose')


def test_database_query():
    assert services.database_query('admin', '') == ['Illegal query']
    assert services.database_query(test_db, "{'DB_links.PubChem': '3333'}") == [glucose]
    assert services.database_query(test_db, "{'Names': 'Glucose'}") == [glucose]


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
    assert len(meh[0]) == 32
    assert len(meh[1]) == 15
    assert meh[0][2] == 'M+Na'


def test_batch_ms_adduct_search():
    result = services.batch_ms_adduct_search(test_db, "181.071188116\n0.0", "form", 2.0, ['M+H'], ['kb|fm.2944'], False,
                                             True, False)
    assert len(result) == 2
    meh = result[0]['adducts']
    assert len(meh) == 2
    assert isinstance(meh[1]['isomers'], list)
    assert result[0]['native_hit'] is False
    assert result[0]['min_steps'] == 0


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
    assert len(services.similarity_search(test_db, 'OCC1OC(O)C(C(C1O)O)O', 0.9, "FP2", 100)) == 28
    assert len(services.similarity_search(test_db, test_molfile, 0.8, 'FP4', 100)) == 9


def test_substructure_search():
    assert len(services.substructure_search('KEGGexp', 'cccccc', 100)) == 100
    assert isinstance(services.substructure_search('KEGGexp', 'Nc1ncnc2[nH]cnc12', 100)[0], dict)


def test_structure_search():
    assert services.structure_search(test_db, "mol", test_molfile)[0][u'_id'] == u'C84d297bb12c40a0996e449dfc54afd69ccc3dd54'
    assert services.structure_search(test_db, "smi", 'OCC1OC(O)C(C(C1O)O)O') == [glucose]

