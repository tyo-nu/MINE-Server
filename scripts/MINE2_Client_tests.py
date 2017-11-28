from lib.biokbase.mine_database.Client import mineDatabaseServices, ServerError
from nose.tools import assert_raises

services = mineDatabaseServices(url='http://modelseed.org/services/mine-database')
test_db = 'EcoCycPickaxe'
test_molfile = open("scripts/xanthine.mol", "r").read()
glucose = {u'SMILES': u'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', u'Inchikey': u'WQZGKKKJIJFFOK-GASJEMHNSA-N',
           u'Generation': 1, u'MINE_id': 917030, u'NP_likeness': 2.62691337083175,
           u'Sources': [{u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C2ce917f1d3aaef7501595123894a68a7d786a9e7']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C2c3b118fb2dbd237cc8a4878b0643385860fb427',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Ca67dd0794223d284a9517566c3fd4107728e0808']},
                        {u'Operators': [u'3.1.1.b'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C71c9c2296f22b1aa1bd031d91e391fa415dc72fd']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C1f719e8cf42a266bc93fe185816e244f5498fc1e']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'Cd51c567a382bff0cba43ecc97d5807e0ffebb5f8',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'1.1.-1.h'], u'Compounds': [u'C20fa60742147ed7b7c328d7835b84d14036d7f36',
                                                                     u'C1b424c6768e3787b3453fa5d572b06106b8933d4',
                                                                     u'Cac741137c0d633e53486ed958334b206a85dcc03']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C8d3efae75c074cd20730fbc4752c0498d94ba40e',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'1.1.1.g'], u'Compounds': [u'C6d0b4d05908b6e7fd3f71b56ca0e81c7a9b9027d',
                                                                    u'C8e823fcc54dec371d81eebaa8e01dc4796e4add6']},
                        {u'Operators': [u'3.1.3.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cff41aa31cc2f8f2f13759ff3d4283726265072ee']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cd51de949ca3a35a5b1bdd3959225766ffd2fc027']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cfdd7350d4ca8405dd44990548ae18d5d56364003']},
                        {u'Operators': [u'1.14.13.d'], u'Compounds': [u'C20fa60742147ed7b7c328d7835b84d14036d7f36',
                                                                      u'C71c9c2296f22b1aa1bd031d91e391fa415dc72fd',
                                                                      u'C1aa818910461f5961959eb05dd06b8f7a4fbbb9f',
                                                                      u'Cac741137c0d633e53486ed958334b206a85dcc03']}],
           u'Mass': 180.063388104, u'Names': [u'Glucose', u'Dextrose', u'D-Glucose', u'Grape sugar', u'D-Glucopyranose'],
           u'Formula': u'C6H12O6', u'_id': u'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f'}


def test_quick_search():
    assert glucose in services.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N')
    assert services.quick_search(test_db, 'C00031') == [glucose]
    assert glucose in services.quick_search(test_db, 'Glucose')


def test_database_query():
    with assert_raises(ServerError):
        services.database_query('admin', '', "", "")
    assert services.database_query(test_db, "{'MINE_id': 917030}", "", "") == [glucose]
    assert services.database_query(test_db, "{'Inchikey': 'WQZGKKKJIJFFOK-GASJEMHNSA-N'}", "", "") == [glucose]


def test_get_ids():
    meh = services.get_ids(test_db, 'operators', None)
    assert len(meh) == 276
    assert 'dc8c6b00513bbc0a708527e99c138094798c8a74f3820618b9a3cc59ffea7881' in meh
    assert len(services.get_ids(test_db, 'reactions', '{"Products.c_id" : "Ca6efd3153bc763e9471785e4b38d0c227ab534fc"}'))


def test_get_comps():
    meh = services.get_comps(test_db, ['Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f'])
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()
    meh = services.get_comps(test_db, [19160])
    assert len(meh) == 1
    assert 'Reactant_in' in meh[0].keys()


def test_get_rxns():
    meh = services.get_rxns(test_db, ['999a273d9705dceef8e8ed7cd8a2f6a6a31a80a79ed3aa4707c9ee8a40dd5771'])
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
    assert len(result) == 101
    assert isinstance(result[0], dict)


def test_ms2_search():
    params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-'], 'models': ['Bacteria'], 'ppm': False,
              'charge': False, 'halogens': False, 'scoring_function': 'dot_product', 'energy_level': 20}
    result2 = services.ms2_search(open("./scripts/folate.mgf").read(), "mgf", params)
    assert result2
    assert isinstance(result2[0], dict)
    print(result2[0])
    assert u'Spectral_score' in result2[0].keys()
    result2_2 = services.ms2_search(open("./scripts/folate_form.txt").read(), "form", params)
    assert len(result2) == len(result2_2)
    result2_2 = services.ms2_search(open("./scripts/2870575.msp").read(), "msp", params)
    assert len(result2) == len(result2_2)


def test_pathway_search():
    meh = services.pathway_search(test_db, 'C1b443383bfb0f99f1afe6a37f3ff2dadc3dbaff1',
                                                       'C89b394fd02e5e5e60ae1e167780ea7ab3276288e', 3, False)
    assert meh == ['Not yet implemented']


def test_similarity_search():
    assert len(services.similarity_search(test_db, 'OCC1OC(O)C(C(C1O)O)O', 0.9, "MACCS", 100, "", "")) == 100
    assert len(services.similarity_search(test_db, test_molfile, 0.8, 'MACCS', 100, "", "")) == 9


def test_substructure_search():
    meh = services.substructure_search(test_db, 'CCC', 50, "", "")
    assert len(meh) == 50
    assert isinstance(meh[0], dict)


def test_model_search():
    assert services.model_search("human") == [u'Animals', u'Eukaryotes', u'Mammals', u'Vertebrates',
                                              u'hsa', u'Arthropods', u'Insects', u'phu']


def test_structure_search():
    assert services.structure_search(test_db, "smi", 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', "", "")[0][u'_id'] == u'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f'