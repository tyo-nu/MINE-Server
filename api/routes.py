"""Here, routes are defined for all possible API requests. Note that nearly
all actual logic is imported from the minedatabase package."""

from ast import literal_eval

import requests
from flask import Blueprint
from flask import current_app as app
from flask import jsonify, request
from flask.helpers import send_from_directory
from minedatabase.metabolomics import (ms2_search, ms_adduct_search, read_adduct_names,
                                       spectra_download)
from minedatabase.utils import score_compounds

from api.config import Config
from api.database import mongo
from api.exceptions import InvalidUsage
from api.queries import (advanced_search, get_comps, get_ids, get_op_w_rxns, get_ops, get_rxns,
                         get_rxns_for_cpd, model_search, quick_search, similarity_search,
                         structure_search, substructure_search)

if Config.THERMO_ON:
    from api.database_thermo import mine_thermo

from api.utils import (extract_enzyme_pathway_from_kegg_data, get_extra_info,
                       get_smiles_from_mol_string)

# pylint: disable=invalid-name
mineserver_api = Blueprint('mineserver_api', __name__)
# pylint: enable=invalid-name


@mineserver_api.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    """Makes it so user can receive an informative error message rather than
    a default internal server error."""
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response


@mineserver_api.route('/quick-search/<db_name>/q=<query>')
def quick_search_api(db_name, query):
    """Perform a quick search and return results.

    .. :quickref: Quick search

    :param str db_name:
        Name of Mongo database to query against.
    :param str query:
        A MINE id, KEGG code, ModelSEED id, Inchikey, or Name.

    :return: JSON Documents matching query.
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]
    results = quick_search(db, core_db, query)
    json_results = jsonify(results)

    return json_results


# Routes for mol input
@mineserver_api.route('/similarity-search/<db_name>', methods=['POST'])
@mineserver_api.route('/similarity-search/<db_name>/<float:min_tc>',
                      methods=['POST'])
@mineserver_api.route('/similarity-search/<db_name>/<int:limit>',
                      methods=['POST'])
@mineserver_api.route('/similarity-search/<db_name>/<float:min_tc>'
                      '/<int:limit>', methods=['POST'])
# Routes for smiles input
@mineserver_api.route('/similarity-search/<db_name>/smiles=<smiles>')
@mineserver_api.route('/similarity-search/<db_name>/smiles=<smiles>'
                      '/<float:min_tc>')
@mineserver_api.route('/similarity-search/<db_name>/smiles=<smiles>'
                      '/<int:limit>')
@mineserver_api.route('/similarity-search/<db_name>/smiles=<smiles>'
                      '/<float:min_tc>/<int:limit>')
def similarity_search_api(db_name, smiles=None, min_tc=0.7, limit=100):
    """Perform a similarity search for a SMILES string and return results.

    .. :quickref: Compound; Structure similarity search

    Either a SMILES string or mol object string is required. SMILES is the
    recommended format. If GET request, uses smiles. If POST, uses mol.

    :param str db_name:
        Name of Mongo database to query against.
    :param str smiles:
        SMILES string describing molecular structure of query molecule. Either
        smiles or mol arg is required (smiles arg recommended over mol arg).
    :param str,optional mol:
        mol object in str format (should contain a lot of strange chars like
        "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from input
        structures, and cannot convert it to SMILES. Captured from form data.
        Defaults to None.
    :param float,optional min_tc:
        Minimum Tanimoto Coefficient required for similarity match. Defaults to
        0.7.
    :param int,optional limit:
        Maximum number of results (compounds) to return. By default, returns
        up to 100 results (limit=100).
    :param str,optional model:
        KEGG organism code (e.g. 'hsa'). Adds annotations to each compound
        based on whether it is in or could be derived from the KEGG compounds
        in this organism (provided in the 'Likelihood_score' field of each
        compound document). Defaults to None.

    :return: JSON Document of similar compounds.
    :rtype: flask.Response
    """
    json_data = request.get_json()

    if json_data and 'mol' in json_data:
        mol_str = str(json_data['mol'])
        smiles = get_smiles_from_mol_string(mol_str)

    if json_data and 'model' in json_data:
        model = str(json_data['model'])
    else:
        model = None

    model_db = mongo.cx[app.config['KEGG_DB_NAME']]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]

    results = similarity_search(db_name, core_db, smiles, min_tc=min_tc, limit=limit,
                                parent_filter=model, model_db=model_db)
    json_results = jsonify(results)

    return json_results


# Routes for mol input
@mineserver_api.route('/structure-search/<db_name>', methods=['POST'])
# Routes for smiles input
@mineserver_api.route('/structure-search/<db_name>/smiles=<smiles>')
def structure_search_api(db_name, smiles=None):
    """Perform an exact structure search and return results.

    .. :quickref: Compound; Exact structure search

    If GET request, uses smiles. If POST, uses mol.

    :param str db_name:
        Name of Mongo database to query against.
    :param str smiles:
        SMILES string describing molecular structure of query molecule.
    :param str,optional mol:
        mol object in str format (should contain a lot of strange chars like
        "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from input
        structures, and cannot convert it to SMILES. Captured from form data.
        Defaults to None.
    :param str,optional model:
        KEGG organism code (e.g. 'hsa'). Adds annotations to each compound
        based on whether it is in or could be derived from the KEGG compounds
        in this organism (provided in the 'Likelihood_score' field of each
        compound document). Defaults to None.

    :return: JSON Document of match (empty if no match).
    :rtype: flask.Response
    """
    json_data = request.get_json()

    if json_data and 'mol' in json_data:
        mol_str = str(json_data['mol'])
        smiles = get_smiles_from_mol_string(mol_str)

    if json_data and 'model' in json_data:
        model = str(json_data['model'])
    else:
        model = None

    model_db = mongo.cx[app.config['KEGG_DB_NAME']]

    db = mongo.cx[db_name]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]
    ref_db = mongo.cx[app.config['REF_DB_NAME']]
    results = structure_search(db, smiles, model_db=model_db, parent_filter=model)
    results = get_extra_info(db, core_db, ref_db, results)
    json_results = jsonify(results)

    return json_results


# Routes for mol input
@mineserver_api.route('/substructure-search/<db_name>', methods=['POST'])
@mineserver_api.route('/substructure-search/<db_name>/<int:limit>',
                      methods=['POST'])
# Routes for smiles input
@mineserver_api.route('/substructure-search/<db_name>/smiles=<smiles>')
@mineserver_api.route('/substructure-search/<db_name>/smiles=<smiles>'
                      '/<int:limit>')
def substructure_search_api(db_name, smiles=None, limit=100):
    """Perform a substructure search and return results.

    .. :quickref: Compound; Substructure search

    If GET request, uses smiles. If POST, uses mol.

    :param str db_name:
        Name of Mongo database to query against.
    :param str smiles:
        SMILES string describing molecular substructure to search for.
    :param str,optional mol:
        mol object in str format (may contain a lot of whitespace hex-chars
        like "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from
        input structures, and cannot convert it to SMILES. Captured from form
        data. Defaults to None.
    :param int limit:
        Maximum number of results (compounds) to return. By default, returns
        100 results.

    :return: JSON Documents of compounds containing given substructure.
    :rtype: flask.Response
    """
    json_data = request.get_json()

    if json_data and 'mol' in json_data:
        mol_str = str(json_data['mol'])
        smiles = get_smiles_from_mol_string(mol_str)

    if json_data and 'model' in json_data:
        model = str(json_data['model'])
    else:
        model = None

    model_db = mongo.cx[app.config['KEGG_DB_NAME']]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]

    db = mongo.cx[db_name]
    results = substructure_search(db, core_db, smiles, limit=limit, model_db=model_db,
                                  parent_filter=model)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/model-search/q=<query>')
def model_search_api(query):
    """Perform a model search and return results.

    .. :quickref: Model; KEGG model search

    :param str query:
        KEGG Org Code or Org Name of model(s) to search for (e.g. 'hsa' or
        'yeast'). Can provide multiple search terms by separating each term
        with a space.  TODO: change from space delimiter to something else

    :return: JSON Document with KEGG org codes matching query.
    :rtype: flask.Response
    """
    db = mongo.cx[app.config['KEGG_DB_NAME']]
    results = model_search(db, query)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/database-query/<db_name>/q=<mongo_query>')
def database_query_api(db_name, mongo_query):
    """Perform a direct query built with Mongo syntax.

    .. :quickref: MINE DB; Query MINE DB with Mongo syntax

    :param str db_name:
        Name of Mongo database to query against.
    :param str mongo_query:
        A valid Mongo query (e.g. .../q={"ID": "cpd00001"}).

    :return: JSON Documents matching provided Mongo query.
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    results = advanced_search(db, mongo_query)
    # TODO: add model to score_compounds (where None currently is)
    results = score_compounds(db, results, None)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-kegg-info/q=<kegg_id>')
def get_kegg_info(kegg_id):
    """Get EC numbers and Pathway names from KEGG API.

    .. :quickref: KEGG; Get KEGG info for compound

    :param str kegg_id:
        Compound KEGG ID (must start with C).
    """
    if not kegg_id.startswith('C'):
        return None
    uri = f'http://rest.kegg.jp/get/cpd:{kegg_id}'
    r = requests.get(uri)
    data = r.text
    enzymes, pathways = extract_enzyme_pathway_from_kegg_data(data)
    json_results = jsonify({'Enzymes': enzymes, 'Pathways': pathways})
    return json_results


@mineserver_api.route('/get-ids/<db_name>/<collection_name>')
@mineserver_api.route('/get-ids/<db_name>/<collection_name>/q=<query>')
def get_ids_api(db_name, collection_name, query=None):
    """Get Mongo IDs for a subset of a given database collection.

    .. :quickref: ID; Get Mongo IDs

    :param str db_name:
        Name of Mongo database to query against.
    :param str collection_name:
        Name of Mongo collection within database to query against.
    :param str,optional query:
        Specifies subset of collection to retrieve ids for. Formatted as a
        python dict as you would have in argument to db.collection.find().
        Defaults to None.

    :return: List of ids matching query in JSON format.
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    results = get_ids(db, collection_name, query)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-comps/<db_name>', methods=['GET', 'POST'])
def get_comps_api(db_name):
    """Get compounds for specified ids in database.

    .. :quickref: Compound; Get MINE compounds by ID

    :param str db_name:
        Name of Mongo database to query against.
    :param list id_list:
        List of compound ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs can be either MINE IDs or Mongo IDs (_id).
    :param list,optional return_extra_info:
        If provided get extra info on this compound from the core and reference
        databases based on the fields in the list. Possible fields are
        "spectra", "DB links", and "RDKit_fp".

    :return: List of compound JSON documents.
    :rtype: flask.Response
    """
    id_list = request.get_json()['id_list']

    if id_list == ['']:
        id_list = []

    return_extra_info = request.get_json()['return_extra_info']

    if not id_list:
        raise InvalidUsage('id_list must be specified in form data.')

    db = mongo.cx[db_name]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]
    ref_db = mongo.cx[app.config['REF_DB_NAME']]
    results = get_comps(db, id_list, core_db)

    if return_extra_info and results and all(results):
        results = get_extra_info(db, core_db, ref_db, results)

    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-rxns/<db_name>', methods=['POST'])
def get_rxns_api(db_name):
    """Get reactions for specified ids in database.

    .. :quickref: Reaction; Get MINE reactions by ID

    :param str db_name:
        Name of Mongo database to query against.
    :param list id_list:
        List of reaction ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs can be either MINE IDs or Mongo IDs (_id).

    :return: List of reaction JSON documents.
    :rtype: flask.Response
    """
    id_list = request.get_json()['id_list']

    db = mongo.cx[db_name]
    results = get_rxns(db, id_list)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-rxns-product-of/<db_name>/<cpd_id>')
def get_rxns_product_of_api(db_name, cpd_id):
    """Get reactions producing compound.

    .. :quickref: Reaction; Get MINE reactions producing a compound.

    :param str db_name:
        Name of Mongo database to query against.
    :param str cpd_id:
        Mongo ID of compound.

    :return: List of reaction JSON documents.
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    results = get_rxns_for_cpd(db, cpd_id, mode='product_of')
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-rxns-reactant-in/<db_name>/<cpd_id>')
def get_rxns_reactant_in_api(db_name, cpd_id):
    """Get reactions consuming compound.

    .. :quickref: Reaction; Get MINE reactions consuming a compound.

    :param str db_name:
        Name of Mongo database to query against.
    :param str cpd_id:
        Mongo ID of compound.

    :return: List of reaction JSON documents.
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    results = get_rxns_for_cpd(db, cpd_id, mode='reactant_in')
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-ops/<db_name>', methods=['POST'])
def get_ops_api(db_name):
    """Get operators for specified ids in database.

    .. :quickref: Operator; Get MINE operators by ID

    :param str db_name:
        Name of Mongo database to query against.
    :param list,optional id_list:
        List of operator ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs can be either operator ids (e.g. 1.1.-1.h) or Mongo IDs (_id). If
        not provided, all operators are returned (id_list=None).

    :return: List of operator JSON documents.
    :rtype: flask.Response
    """
    if request.get_json():
        id_list = request.get_json()['id_list']
    else:
        id_list = None

    db = mongo.cx[db_name]
    results = get_ops(db, id_list)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-op-w-rxns/<db_name>/<op_id>')
def get_op_w_rxns_api(db_name, op_id):
    """Get operator with all its associated reactions in selected database.

    .. :quickref: Operator; Get reactions for MINE operator

    :param str db_name:
        Name of Mongo database to query against.
    :param str op_id:
        Either operator id (e.g. 1.1.-1.h) or Mongo ID (_id) for operator.

    :return: Operator JSON document (including associated reactions).
    :rtype: flask.Response
    """
    db = mongo.cx[db_name]
    results = get_op_w_rxns(db, op_id)
    if results:
        json_results = jsonify(results)
        return json_results
    else:
        raise InvalidUsage('Operator with ID \"{}\" not found.'.format(op_id))


@mineserver_api.route('/get-op-image/<rule_name>')
def get_op_image_api(rule_name):
    """Get image of operator's reaction produced using ChemAxon software.

    .. :quickref: Operator; Get SMARTS image for operator

    :param str rule_name:
        Name of rule from MetaCyc generalized operator set (e.g. "rule0361").

    :return: Image of operator reaction (.jpg)
    :rtype: flask.Response
    """
    app.logger.info(f"Sending operator image: {rule_name}.jpg")
    return send_from_directory(app.config['OP_IMG_DIR'], f"{rule_name}.jpg")


@mineserver_api.route('/get-thermo-info/<c_id>')
def get_thermo_info(c_id):
    """Get dG formation for compound.

    .. :quickref: Compound; Get dG formation for compound

    :param str c_id:
        Mongo ID of compound.

    :return: JSON dict with dG formation.
    :rtype: flask.Response:
    """
    if app.config['THERMO_ON']:
        dG = mine_thermo.standard_dg_formation_from_cid(c_id)
        thermo_dict = {'dG': dG[0]}
    else:
        thermo_dict = {'dG': -9999}
    json_results = jsonify(thermo_dict)
    return json_results


@mineserver_api.route('/get-adduct-names')
@mineserver_api.route('/get-adduct-names/<adduct_type>')
def get_adduct_names_api(adduct_type='all'):
    """Get names of all adducts for the specified adduct type.

    .. :quickref: Adduct; Get names of metabolomics adducts

    :param str,optional adduct_type:
        Options are 'positive', 'negative', or 'all'. Defaults to 'all'.

    :return:
        JSON array of adduct names. If adduct_type == 'all', then this is an
        array of two arrays, with the first element being positive adducts and
        the second adduct being negative adducts.
    :rtype: flask.Response
    """

    if adduct_type.lower() == 'positive':
        results = read_adduct_names(app.config['POS_ADDUCT_PATH'])
    elif adduct_type.lower() == 'negative':
        results = read_adduct_names(app.config['NEG_ADDUCT_PATH'])
    elif adduct_type == 'all':
        pos_results = read_adduct_names(app.config['POS_ADDUCT_PATH'])
        neg_results = read_adduct_names(app.config['NEG_ADDUCT_PATH'])
        results = [pos_results, neg_results]
    else:
        raise InvalidUsage('URL param <adduct_type> must be "all", "pos", or '
                           '"neg".')

    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/ms-adduct-search/<db_name>', methods=['POST'])
def ms_adduct_search_api(db_name):
    """Search for commpound-adducts matching precursor mass(es).

    .. :quickref: Compound; Search MINE compounds with MS1 data

    Attach all arguments besides db_name as JSON data in POST request.

    :param str db_name:
        Name of Mongo database to query against.
    :param float tolerance:
        Specifies tolerance for m/z, in mDa by default. Can specify in ppm if
        ppm is set to True.
    :param bool charge:
        Positive or negative mode. (True for positive, False for negative).
    :param str text:
        Text as in metabolomics datafile for specific peak.
    :param str,optional text_type:
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines. Default is None.
    :param list,optional adducts:
        List of adducts to use. If not specified, uses all adducts
        (adducts=None).
    :param list,optional models:
        List of model _ids. If supplied, score compounds higher if present
        in metabolic model. Defaults to None.
    :param bool,optional ppm:
        Specifies whether tolerance is in ppm. Defaults to False.
    :param tuple,optional logp:
        Length 2 tuple specifying min and max logp to filter compounds (e.g.
        (-1, 2)). Defaults to None.
    :param bool,optional halogen:
        Specifies whether to filter out compounds containing F, Cl, or Br.
        Filtered out if set to True. Defaults to False.
    :param bool,optional verbose:
        If True, verbose output. Defaults to False.

    :return:
        JSON array of compounds that match m/z within defined tolerance and
        after passing other defined filters (such as logP).
    :rtype: flask.Response
    """
    json_data = request.get_json()

    if 'tolerance' in json_data:
        tolerance = float(json_data['tolerance'])
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge' in json_data:
        charge_pos = bool(json_data['charge'])
        if charge_pos:
            charge = "+"
        else:
            charge = "-"
    else:
        raise InvalidUsage('<charge> argument must be specified. "Positive" '
                           'for positive mode, "Negative" for negative mode.')

    if 'text' in json_data:
        text = json_data['text']
    else:
        raise InvalidUsage('<text> argument must be specified.')

    if 'text_type' in json_data:
        text_type = json_data['text_type']
    else:
        text_type = 'form'

    if 'adducts' in json_data:
        adducts = literal_eval(str(json_data['adducts']))
        assert isinstance(adducts, list)
    else:
        adducts = None

    if 'models' in json_data:
        models = literal_eval(str(json_data['models']))
        assert isinstance(models, list)
        if models == []:
            models = None
    else:
        models = None

    if 'ppm' in json_data:
        ppm = bool(json_data['ppm'])
    else:
        ppm = None

    if 'logP' in json_data:
        logp = tuple(literal_eval(str(json_data['logP'])))
    else:
        logp = None

    if 'halogens' in json_data:
        halogens = bool(json_data['halogens'])
    else:
        halogens = None

    if 'verbose' in json_data:
        verbose = bool(json_data['verbose'])
    else:
        verbose = False

    ms_params = {
        'tolerance': tolerance,
        'charge': charge,
        'adducts': adducts,
        'models': models,
        'ppm': ppm,
        'logp': logp,
        'halogens': halogens,
        'verbose': verbose
    }
    print(ms_params)

    db = mongo.cx[db_name]
    keggdb = mongo.cx[app.config['KEGG_DB_NAME']]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]

    results = ms_adduct_search(db, core_db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    if results:
        app.logger.info(f'MS Search successful ({len(results)} results found)')

    return json_results


@mineserver_api.route('/ms2-search/<db_name>', methods=['POST'])
def ms2_search_api(db_name):
    """Search for commpound-adducts matching precursor mass(es).

    .. :quickref: Compound; Search MINE compounds with MS2 data

    Attach all arguments besides db_name as form data in POST request.

    :param str db_name:
        Name of Mongo database to query against.
    :param float tolerance:
        Specifies tolerance for m/z, in mDa by default. Can specify in ppm if
        ppm is set to True.
    :param bool charge:
        Positive or negative mode. (True for positive, False for negative).
    :param int energy_level:
        Fragmentation energy level to use. May be 10, 20, or 40.
    :param str scoring_function:
        Scoring function to use. Can be either 'jaccard' or 'dot product'.
    :param str text:
        Text as in metabolomics datafile for specific peak.
    :param str,optional text_type:
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines. Default is None.
    :param list,optional adducts:
        List of adducts to use. If not specified, uses all adducts.
        (adducts=None)
    :param list,optional models:
        List of model _ids. If supplied, score compounds higher if present
        in metabolic model. Defaults to None.
    :param bool,optional ppm:
        Specifies whether tolerance is in ppm. Defaults to False.
    :param tuple,optional logp:
        Length 2 tuple specifying min and max logp to filter compounds (e.g.
        (-1, 2)). Defaults to None.
    :param bool,optional halogens:
        Specifies whether to filter out compounds containing F, Cl, or Br.
        Filtered out if set to True. Defaults to False.

    :return:
        JSON array of compounds that match m/z within defined tolerance and
        after passing other defined filters (such as logP).
    :rtype: flask.Response
    """
    json_data = request.get_json()

    if 'tolerance' in json_data:
        tolerance = float(json_data['tolerance'])
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge' in json_data:
        charge_pos = bool(json_data['charge'])
        if charge_pos:
            charge = "+"
        else:
            charge = "-"
    else:
        raise InvalidUsage('<charge> argument must be specified. "Positive" '
                           'for positive mode, "Negative" for negative mode.')

    if 'energy_level' in json_data:
        energy_level = int(json_data['energy_level'])
    else:
        raise InvalidUsage('<energy_level> argument must be specified. '
                           'Possible values are 10, 20, or 40.')

    if 'scoring_function' in json_data:
        scoring_function = json_data['scoring_function']
    else:
        raise InvalidUsage("<scoring_function> argument must be specified. "
                           "Possible values are 'jaccard' and 'dot product'.")

    if 'text' in json_data:
        text = json_data['text']
    else:
        raise InvalidUsage('<text> argument must be specified.')

    if 'text_type' in json_data:
        text_type = json_data['text_type']
    else:
        text_type = None

    if 'adducts' in json_data:
        adducts = literal_eval(str(json_data['adducts']))
        assert isinstance(adducts, list)
    else:
        adducts = None

    if 'models' in json_data:
        models = literal_eval(str(json_data['models']))
        assert isinstance(models, list)
        if models == []:
            models = None
    else:
        models = None

    if 'ppm' in json_data:
        ppm = bool(json_data['ppm'])
    else:
        ppm = None

    if 'logP' in json_data:
        logp = tuple(literal_eval(str(json_data['logP'])))
    else:
        logp = None

    if 'halogens' in json_data:
        halogens = bool(json_data['halogens'])
    else:
        halogens = None

    if 'verbose' in json_data:
        verbose = bool(json_data['verbose'])
    else:
        verbose = False

    ms_params = {
        'tolerance': tolerance,
        'charge': charge,
        'energy_level': energy_level,
        'scoring_function': scoring_function,
        'adducts': adducts,
        'models': models,
        'ppm': ppm,
        'logp': logp,
        'halogens': halogens,
        'verbose': verbose
    }

    db = mongo.cx[db_name]
    keggdb = mongo.cx[app.config['KEGG_DB_NAME']]
    core_db = mongo.cx[app.config['CORE_DB_NAME']]

    results = ms2_search(db, core_db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    return json_results

@mineserver_api.route('/spectra-download/<mongo_id>')
def spectra_download_api(mongo_id):
    """Download one or more spectra for compounds matching a given query.

    .. :quickref: Spectra; Get computationally predicted MS2 spectra

    :param str mongo_id:
        Mongo ID of compound to get spectra for.

    :return: Text of spectra for input compound.
    :rtype: flask.Response
    """
    core_db = mongo.cx[app.config['CORE_DB_NAME']]
    results = spectra_download(core_db, mongo_id)

    return app.response_class(results)
