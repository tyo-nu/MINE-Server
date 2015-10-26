__author__ = 'JGJeffryes'
from lib.biokbase.mine_database.Impl import mineDatabaseServices
import time

test_db = 'EcoCycexp2'
glucose = {u'SMILES': u'OCC1OC(O)C(C(C1O)O)O', u'Inchikey': u'WQZGKKKJIJFFOK-UHFFFAOYSA-N', u'Generation': 0.0,
           u'MINE_id': 19160, u'Mass': 180.063388104, u'Names': [u'Hexose', u'D-Idose', u'Glucose', u'Mannose',
            u'D-Gulose', u'D-Allose', u'D-Hexose', u'Dextrose', u'Seminose', u'L-Gulose', u'D-Talose', u'D-Aldose',
            u'D-Mannose', u'D-Aldose2', u'D-Aldose1', u'D-Glucose', u'D-Altrose', u'Carubinose', u'Grape sugar',
            u'L-Galactose', u'D-Galactose', u'D-ido-Hexose', u'D-gulo-Hexose', u'D-talo-Hexose', u'beta-D-Mannose',
            u'beta-D-Glucose', u'D-altro-Hexose', u'alpha-D-Glucose', u'alpha-D-Mannose', u'D-glucopyranose',
            u'beta-D-Galactose', u'alpha-D-Galactose', u'D-galactopyranose', u'1,4-beta-D-Mannooligosaccharide'],
           u'NP_likeness': 0, u'Formula': u'C6H12O6', u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13'}
test_molfile = open("./scripts/xanthine.mol", "r").read()


class MZParams():
    def __init__(self):
        self.db = 'KEGGexp2'
        self.tolerance = 2.0
        self.adducts = ['[M+H]+']
        self.models = ['Bacteria']
        self.ppm = False
        self.charge = True
        self.halogens = False


class Options():
    def __init__(self):
        self.adduct_file = "lib/All Adducts.txt"
        self.kbase_db = 'KBase'
        self.test_db = '1GenEcoCyc'
        self.positive_adduct_file = "lib/Positive Adducts full.txt"
        self.negative_adduct_file = "lib/Negative Adducts full.txt"
        self.adduct_list = ['M+H', 'M-H']
        self.tolerance = 2
        self.ppm = False
        self.ms_file = "lib/Unknowns.txt"
        self.halogens = True

config = Options()
services = mineDatabaseServices(None)

def test_quick_search():
    assert services.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N')[0] == [glucose]
    assert services.quick_search(test_db, 'C00031')[0] == [glucose]
    assert glucose in services.quick_search(test_db, 'Glucose')[0]


def test_database_query():
    assert services.database_query('admin', '', "", "")[0] == ['Illegal query']
    assert services.database_query(test_db, "{'MINE_id': 19160}", "", "")[0] == [glucose]
    assert services.database_query(test_db, "{'Names': 'Glucose'}", "", "")[0] == [glucose]


def test_get_comps():
    meh = services.get_comps(test_db, ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13'])[0]
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()
    meh = services.get_comps(test_db, [19160])[0]
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()


def test_get_rxns():
    meh = services.get_rxns(test_db, ['R2e28f382545b6b00b88bcd5b1bb927bc480c2711'])[0]
    assert len(meh) == 1
    assert 'Operators' in meh[0].keys()


def test_get_ops():
    meh = services.get_ops(test_db, ['2.7.1.a'])[0][0]
    assert meh['Reactions_predicted'] > 250


def test_get_adducts():
    meh = services.get_adducts()[0]
    assert len(meh[0]) == 33
    assert len(meh[1]) == 15
    assert meh[0][2] == '[M+NH4]+'


def test_ms_adduct_search():
    params = {'db': test_db, 'tolerance': 2.0, 'adducts': ['[M+H]+'], 'models': ['Bacteria'], 'ppm': False,
              'charge': True, 'halogens': False}
    result = services.ms_adduct_search("181.071188116\n0.0", "form", params)[0]
    assert len(result) == 31
    assert isinstance(result[0], dict)
    keys = [u'SMILES', u'NP_likeness', u'logP', u'adduct', u'maxKovatsRI', u'MINE_id', u'Inchikey', u'Generation',
             u'Formula', u'minKovatsRI', u'_id', u'peak_name']
    assert result[0].keys() == keys


def test_ms2_search():
    params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-'], 'models': ['Bacteria'], 'ppm': False,
              'charge': False, 'halogens': False, 'scoring_function': 'jacquard', 'energy_level': 1}
    result2 = services.ms2_search(open("./scripts/folate.mgf").read(), "mgf", params)[0]
    assert isinstance(result2[0], dict)
    assert result2
    print(result2[0])
    keys = [u'SMILES', u'NP_likeness', u'logP', u'adduct', u'maxKovatsRI', u'MINE_id', u'Inchikey', u'Generation',
            u'Spectral_score', u'Formula', u'minKovatsRI', u'_id', u'peak_name']
    assert result2[0].keys() == keys
    result2_2 = services.ms2_search(open("./scripts/folate_form.txt").read(), "form", params)[0]
    assert len(result2) == len(result2_2)


def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert meh[0] == ['Not yet implemented']
    """assert len(meh) == 1
    assert meh[0] == [u'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rbbc40c762b05d59890c196c949522e3ee6ca08c6',
                      u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R0e87f4c178bcb78b9190938b3413b7889b3fbad4',
                      u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True)) == 4"""


def test_similarity_search():
    print services.similarity_search(test_db, 'O=C1CC(OC1COP(=O)(OP(=O)(O)O)O)n1cc(C)c(nc1=O)O', 0.8, 'FP2', 100, "", "")
    print services.similarity_search(test_db, test_molfile, 0.8, 'FP4', 100, "", "")


def test_structure_search():
    print services.structure_search("EcoCycexp2", "smi", 'O=C1CC(OC1COP(=O)(OP(=O)(O)O)O)n1cc(C)c(nc1=O)O', "", "")
    print services.structure_search("EcoCycexp2", "mol", test_molfile, "", "")


def test_substructure_search():
    print services.substructure_search(test_db, 'O=C1CC(OC1COP(=O)(OP(=O)(O)O)O)n1cc(C)c(nc1=O)O', 20, "", "")
    print services.substructure_search('KEGGexp2', test_molfile, 20, "", "")


def test_model_search():
    assert services.model_search("human")[0] == [u'Eukaryotes', u'Animals', u'hsa', u'Mammals', u'Vertebrates',
                                                 u'Arthropods', u'Insects', u'phu']
