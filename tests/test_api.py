"""Test all API functionality, ensuring that all responses are healthy and
contain data. The content of the responses is not tested, as the content is
generated by dependencies in minedatabase, which already have their own
tests."""

import json

from flask import url_for

from api.config import Config


def assert_response_fields(response, status_code=200):
    """Ensure that Flask response is healthy, in JSON format, and has data.

    Parameters
    ----------
    response : flask.Flask.Response
        Contains response data in JSON format.
    status_code : int, optional (default: 200)
        Expected response status code.
    """
    assert response
    assert response.status_code == status_code
    if status_code == 200:
        assert response.is_json
        assert response.json


def post_json(client, url, json_dict):
    """Send form data as a json to the specified url in a post request.

    Parameters
    ----------
    client : flask.Flask.test_client
        Client object that sends requests.
    url : str
        URL to send post request to.
    json_dict : dict
        Dictionary containing form data for post request.

    Returns
    -------
    response : flask.Flask.Response
        Contains response data in JSON format.
    """
    response = client.post(url, data=json.dumps(json_dict),
                           content_type='application/json')
    return response


def test_quick_search_api(client):
    """
    GIVEN a compound to query using quick search via the API
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.quick_search_api',
                  db_name='mongotest', query='cpd00348')
    response = client.get(url)
    assert_response_fields(response)


def test_similarity_search_api(client):
    """
    GIVEN a compound to query using similarity search via the API
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    smiles = r'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)' \
             r'[C@H]1O'
    url = url_for('mineserver_api.similarity_search_api', db_name='mongotest',
                  smiles=smiles)
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.similarity_search_api', db_name='mongotest',
                  smiles=smiles, min_tc=0.9)
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.similarity_search_api', db_name='mongotest',
                  smiles=smiles, limit=1)
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.similarity_search_api', db_name='mongotest',
                  smiles=smiles, min_tc=0.9, limit=1)
    response = client.get(url)
    assert_response_fields(response)


def test_structure_search_api(client):
    """
    GIVEN a structure in SMILES format
    WHEN that structure is queried against a MINE DB
    THEN make sure the response is healthy and contains response data
    """
    smiles = r'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)' \
             r'[C@H]1O'
    url = url_for('mineserver_api.structure_search_api', db_name='mongotest',
                  smiles=smiles)
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.structure_search_api', db_name='mongotest',
                  smiles=smiles, stereo=False)
    response = client.get(url)
    assert_response_fields(response)


def test_substructure_search_api(client):
    """
    GIVEN a substructure in SMILES format
    WHEN that substructure is queried against a MINE DB
    THEN make sure the response is healthy and contains response data
    """
    smiles = r'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)' \
             r'[C@H]1O'
    url = url_for('mineserver_api.substructure_search_api',
                  db_name='mongotest', smiles=smiles)
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.substructure_search_api',
                  db_name='mongotest', smiles=smiles, limit=2)
    response = client.get(url)
    assert_response_fields(response)


def test_database_query_api(client):
    """
    GIVEN a direct query (using Mongo syntax) to a MINE DB
    WHEN that query is processed via advanced search
    THEN make sure the response is healthy and contains response data
    """
    query = '{"ID": "cpd00348"}'
    url = url_for('mineserver_api.database_query_api', db_name='mongotest',
                  mongo_query=query)
    response = client.get(url)
    assert_response_fields(response)


def test_get_ids_api(client):
    """
    GIVEN a MINE DB and collection name within that DB
    WHEN all IDs are requested for that collection
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_ids_api', db_name='mongotest',
                  collection_name='compounds')
    response = client.get(url)
    assert_response_fields(response)

    query = '{"ID": "cpd00348"}'
    url = url_for('mineserver_api.get_ids_api', db_name='mongotest',
                  collection_name='compounds', query=query)
    response = client.get(url)
    assert_response_fields(response)


def test_get_comps_api(client):
    """
    GIVEN a MINE DB
    WHEN specified compounds are requested from that DB using a POST request
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_comps_api', db_name='mongotest')
    id_list = ["Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f",
               "C03e0b10e6490ce79a7b88cb0c4e17c2bf6204352"]
    response = post_json(client, url, {'id_list': id_list})
    assert_response_fields(response)


def test_get_rxns_api(client):
    """
    GIVEN a MINE DB
    WHEN specified reactions are requested from that DB
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_rxns_api', db_name='mongotest')
    id_list = ["4542c96f4bca04bfe2db15bc71e9eaee38bee5b87ad8a6752a5c4718ba19"
               "74c1"]
    response = post_json(client, url, {'id_list': id_list})
    assert_response_fields(response)


def test_get_ops_api(client):
    """
    GIVEN a MINE DB
    WHEN operators are requested from that DB
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_ops_api', db_name='mongotest')
    response = post_json(client, url, None)
    assert_response_fields(response)

    url = url_for('mineserver_api.get_ops_api', db_name='mongotest')
    op_id_list = ['2.7.1.a', 'test_op']
    response = post_json(client, url, {'id_list': op_id_list})
    assert_response_fields(response)


def test_get_op_w_rxns_api(client):
    """
    GIVEN a MINE DB and operator ID
    WHEN that operator and its associated reactions are requested
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_op_w_rxns_api', db_name='mongotest',
                  op_id='2.7.1.a')
    response = client.get(url)
    assert_response_fields(response)

    url = url_for('mineserver_api.get_op_w_rxns_api', db_name='mongotest',
                  op_id='invalid')
    response = client.get(url)
    assert_response_fields(response, status_code=400)


def test_get_adduct_names_api(client):
    """
    GIVEN a request to get the names of metabolomics adducts via the API
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.get_adduct_names_api',
                  adduct_type='negative')
    response = client.get(url)
    assert_response_fields(response)


def test_ms_adduct_search_api(client):
    """
    GIVEN a request with an MS1 adduct search query is made
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 1000,
        'charge_mode': 'Positive',
        'text': '161'
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms_adduct_search_api_error(client):
    """
    GIVEN a request with an invalid MS1 search query
    WHEN a response is received
    THEN make sure the response contains a bad request error code 400
    """
    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 10,
        'charge_mode': True,
        'text': '161',
    }
    for key in json_dict:
        new_dict = json_dict.copy()
        del new_dict[key]
        response = post_json(client, url, new_dict)
        assert_response_fields(response, status_code=400)


def test_ms_adduct_search_api_mgf(client):
    """
    GIVEN a request with an MS1 adduct search query is made in MGF format
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    mgf_text = "BEGIN IONS\nTITLE=Ion 1\nSCANS=7906\nRTINSECONDS=2999.0465\n" \
               "CHARGE=2+\nPEPMASS=188.3805\n80.1044 2.6\n81.3382 6.8\n" \
               "143.0018 19.6\nEND IONS\n\nBEGIN IONS\nTITLE=Ion 2\n" \
               "SCANS=7908\nRTINSECONDS=2999.8819\nCHARGE=2+\n" \
               "PEPMASS=95.95\n67.9357 3.7\nEND IONS\n\n"

    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 1000,
        'charge_mode': True,
        'text': mgf_text,
        'text_type': 'mgf'
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms_adduct_search_api_mzxml(client):
    """
    GIVEN a request with an MS1 adduct search query is made in mzXML format
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    config = Config()
    with open(config.TEST_DATA_DIR + '/mzxml_data.mzxml', 'r') as infile:
        mzxml_text = infile.read()

    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 1000,
        'charge_mode': True,
        'text': mzxml_text,
        'text_type': 'mzxml'
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms_adduct_search_api_msp(client):
    """
    GIVEN a request with an MS1 adduct search query is made in MSP format
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    msp_text = "NAME: test\nPRECURSORMZ: 153.0195\nPRECURSORTYPE: [M-H]-\n" \
               "RETENTIONTIME: 9.461333\nIONMODE: Negative\nNum Peaks: 5\n" \
               "82.1518	1000\n106.939	7500\n109.0055	11000000\n" \
               "109.8953	12000\n153.1177	12000\n\n"

    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 1000,
        'charge_mode': True,
        'text': msp_text,
        'text_type': 'msp'
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms_adduct_search_api_all_params(client):
    """
    GIVEN a request with an MS1 adduct search query with all optional params
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.ms_adduct_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 10,
        'charge_mode': False,
        'text': '259.02244262600003',
        'text_type': 'form',
        'adducts': "['[M-H]-']",
        'models': "['eco']",
        'ppm': False,
        # TODO: add kovats and logp data to mongotest MINE
        # 'kovats': "(0, 20000)",
        # 'logp': "(-35, 35)",
        'halogens': True,
        'verbose': True
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms2_search_api(client):
    """
    GIVEN a request with an MS2 search query
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.ms2_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 10,
        'charge_mode': True,
        'energy_level': 20,
        'scoring_function': 'dot product',
        'text': '261.037\n43.0189 1\n59.013 1\n96.970 10',
        'text_type': 'form'
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_ms2_search_api_error(client):
    """
    GIVEN a request with an invalid MS2 search query
    WHEN a response is received
    THEN make sure the response contains a bad request error code 400
    """
    url = url_for('mineserver_api.ms2_search_api', db_name='mongotest')
    json_dict = {
        'tolerance': 10,
        'charge_mode': True,
        'energy_level': 20,
        'scoring_function': 'dot product',
        'text': '261.037\n43.0189 1\n59.013 1\n96.970 10',
    }
    for key in json_dict:
        new_dict = json_dict.copy()
        del new_dict[key]
        response = post_json(client, url, new_dict)
        assert_response_fields(response, status_code=400)


def test_ms2_search_api_all_params(client):
    """
    GIVEN a request with an MS2 search query with all optional params
    WHEN a response is received
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.ms2_search_api', db_name='UniversalMINE')
    json_dict = {
        'tolerance': 10,
        'charge_mode': True,
        'energy_level': 20,
        'scoring_function': 'dot product',
        'text': '261.037\n43.0189 1\n59.013 1\n96.970 10',
        'text_type': 'form',
        'adducts': "['[M+H]+']",
        'models': "['hsa']",
        'ppm': True,
        # TODO: add kovats and logp data to mongotest MINE
        # 'kovats': (0, 1),
        # 'logp': (0, 1),
        'halogens': True,
        'verbose': True
    }
    response = post_json(client, url, json_dict)
    assert_response_fields(response)


def test_spectra_download_api(client):
    """
    GIVEN a request with a Mongo syntax query
    WHEN downloading spectra for compounds matching that query
    THEN make sure the response is healthy and contains response data
    """
    url = url_for('mineserver_api.spectra_download_api', db_name='mongotest')
    response = client.get(url)
    assert response
    assert response.status_code == 200
    assert response.data

    url = url_for('mineserver_api.spectra_download_api', db_name='mongotest',
                  mongo_query='{"ID": "cpd00348"}')
    response = client.get(url)
    assert response
    assert response.status_code == 200
    assert response.data

    url = url_for('mineserver_api.spectra_download_api', db_name='mongotest',
                  mongo_query='{"ID": "cpd00348"}', parent_filter='hsa',
                  putative=False)
    response = client.get(url)
    assert response
    assert response.status_code == 200
    assert response.data
