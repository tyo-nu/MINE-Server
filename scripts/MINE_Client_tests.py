__author__ = 'JGJeffryes'
from lib.biokbase.mine_database.Client import mineDatabaseServices

services = mineDatabaseServices(url='http://bio-data-1.mcs.anl.gov/services/mine-database')
test_db = 'EcoCycexp2'
test_molfile = open("scripts/xanthine.mol", "r").read()
glucose = {u'SMILES': u'OCC1OC(O)C(C(C1O)O)O', u'Inchikey': u'WQZGKKKJIJFFOK-UHFFFAOYSA-N', u'Generation': 0.0,
           u'MINE_id': 19160, u'Mass': 180.063388104, u'Names': [u'Hexose', u'D-Idose', u'Glucose', u'Mannose',
            u'D-Gulose', u'D-Allose', u'D-Hexose', u'Dextrose', u'Seminose', u'L-Gulose', u'D-Talose', u'D-Aldose',
            u'D-Mannose', u'D-Aldose2', u'D-Aldose1', u'D-Glucose', u'D-Altrose', u'Carubinose', u'Grape sugar',
            u'L-Galactose', u'D-Galactose', u'D-ido-Hexose', u'D-gulo-Hexose', u'D-talo-Hexose', u'beta-D-Mannose',
            u'beta-D-Glucose', u'D-altro-Hexose', u'alpha-D-Glucose', u'alpha-D-Mannose', u'D-glucopyranose',
            u'beta-D-Galactose', u'alpha-D-Galactose', u'D-galactopyranose', u'1,4-beta-D-Mannooligosaccharide'],
           u'NP_likeness': 0, u'Formula': u'C6H12O6', u'_id': u'Cb5b3273ab083d77ed29fbef8f7e464929af29c13'}


def test_quick_search():
    assert services.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N') == [glucose]
    assert services.quick_search(test_db, 'C00031') == [glucose]
    assert glucose in services.quick_search(test_db, 'Glucose')


def test_database_query():
    assert services.database_query('admin', '', "", "") == ['Illegal query']
    assert services.database_query(test_db, "{'MINE_id': 19160}", "", "") == [glucose]
    assert services.database_query(test_db, "{'Names': 'Glucose'}", "", "") == [glucose]


def test_get_ids():
    meh = services.get_ids(test_db, 'operators', None)
    assert len(meh) == 179
    assert u'1.1.-1.h' in meh
    assert len(services.get_ids(test_db, 'operators', "{'_id':'1.1.-1.h'}"))


def test_get_comps():
    meh = services.get_comps(test_db, ['Cb5b3273ab083d77ed29fbef8f7e464929af29c13'])
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()
    meh = services.get_comps(test_db, [19160])
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()


def test_get_rxns():
    meh = services.get_rxns(test_db, ['R2e28f382545b6b00b88bcd5b1bb927bc480c2711'])
    assert len(meh) == 1
    assert 'Operators' in meh[0].keys()


def test_get_ops():
    meh = services.get_ops(test_db, ['2.7.1.a'])[0]
    assert meh['Reactions_predicted'] > 250


def test_get_operator():
        meh = services.get_operator(test_db, '2.7.1.a')
        assert meh['Reactions_predicted'] > 250
        assert len(meh['Reaction_ids']) == meh['Reactions_predicted']


def test_get_adducts():
    meh = services.get_adducts()
    assert len(meh[0]) == 33
    assert len(meh[1]) == 15
    assert meh[0][3] == '[M+Na]+'


def test_ms_adduct_search():
    params = {'db': test_db, 'tolerance': 2.0, 'adducts': ['[M+H]+'], 'models': ['Bacteria'], 'ppm': False,
              'charge': True, 'halogens': False}
    result = services.ms_adduct_search("181.071188116\n0.0", "form", params)
    assert len(result) == 31
    assert isinstance(result[0], dict)


def test_ms2_search():
    params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-'], 'models': ['Bacteria'], 'ppm': False,
              'charge': False, 'halogens': False, 'scoring_function': 'dot_product', 'energy_level': 20}
    result2 = services.ms2_search(open("./scripts/folate.mgf").read(), "mgf", params)
    assert result2
    assert isinstance(result2[0], dict)
    print(result2[0])
    keys = [u'SMILES', u'NP_likeness', u'logP', u'adduct', u'maxKovatsRI', u'MINE_id', u'Inchikey', u'Generation',
            u'Formula', u'Spectral_score', u'minKovatsRI', u'_id', u'peak_name']
    assert u'Spectral_score' in result2[0].keys()
    result2_2 = services.ms2_search(open("./scripts/folate_form.txt").read(), "form", params)
    assert len(result2) == len(result2_2)
    result2_2 = services.ms2_search(open("./scripts/2870575.msp").read(), "msp", params)
    assert len(result2) == len(result2_2)



def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert meh == ['Not yet implemented']
    """assert len(meh) == 1
    assert meh[0] == [u'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1', u'Rbbc40c762b05d59890c196c949522e3ee6ca08c6',
                      u'C4d1c9d1a3841a799052b6e347f1a9553ed088092', u'R0e87f4c178bcb78b9190938b3413b7889b3fbad4',
                      u'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']
    assert len(services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, True)) == 4"""


def test_similarity_search():
    assert len(services.similarity_search(test_db, 'OCC1OC(O)C(C(C1O)O)O', 0.9, "FP2", 100, "", "")) == 28
    assert len(services.similarity_search(test_db, test_molfile, 0.8, 'FP4', 100, "", "")) == 7


def test_substructure_search():
    assert len(services.substructure_search(test_db, 'cccccc', 100, "", "")) == 100
    assert isinstance(services.substructure_search(test_db, 'Nc1ncnc2[nH]cnc12', 100, "", "")[0], dict)


def test_model_search():
    assert services.model_search("human") == [u'Animals', u'Eukaryotes', u'Mammals', u'Vertebrates',
                                              u'hsa', u'Arthropods', u'Insects', u'phu']


def test_structure_search():
    assert services.structure_search(test_db, "mol", test_molfile, "", "")[0][u'_id'] == u'C84d297bb12c40a0996e449dfc54afd69ccc3dd54'
    assert services.structure_search(test_db, "smi", 'OCC1OC(O)C(C(C1O)O)O', "", "") == [glucose]
