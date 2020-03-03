"""Here, routes are defined for all possible API requests. Note that nearly
all actual logic is imported from the minedatabase package."""

from ast import literal_eval

from flask import Blueprint
from flask import current_app as app
from flask import jsonify, request

from api.database import mongo
from api.exceptions import InvalidUsage
from minedatabase.metabolomics import (ms2_search, ms_adduct_search,
                                       read_adduct_names, spectra_download)
from minedatabase.queries import (advanced_search, get_comps, get_ids,
                                  get_op_w_rxns, get_ops, get_rxns,
                                  model_search, quick_search,
                                  similarity_search, structure_search,
                                  substructure_search)
from minedatabase.utils import score_compounds, get_smiles_from_mol_string


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

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    query : str
        A MINE id, KEGG code, ModelSEED id, Inchikey, or Name.

    Returns
    -------
    json_results : flask.Response
        JSON Documents matching query.
    """
    db = mongo.cx[db_name]
    results = quick_search(db, query)
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
def similarity_search_api(db_name, smiles=None, min_tc=0.7, limit=-1):
    """Perform a similarity search for a SMILES string and return results.

    Either a SMILES string or mol object string is required. SMILES is the
    recommended format.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular structure of query molecule. Either
        smiles or mol arg is required (smiles arg recommended over mol arg).
    mol : str, optional (default: None)
        mol object in str format (should contain a lot of strange chars like
        "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from
        input structures, and cannot convert it to SMILES. Captured from form
        data.
    min_tc : float, optional (default: 0.7)
        Minimum Tanimoto Coefficient required for similarity match.
    limit : int, optional (default: -1)
        Maximum number of results (compounds) to return. By default, returns
        all results.
    model : str, optional (default: None)
        KEGG organism code (e.g. 'hsa'). Adds annotations to each compound
        based on whether it is in or could be derived from the KEGG compounds
        in this organism (provided in the 'Likelihood_score' field of each
        compound document).

    Returns
    -------
    json_results : flask.Response
        JSON Documents of similar compounds.
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
    results = similarity_search(db, smiles, min_tc=min_tc, limit=limit,
                                model_db=model_db, parent_filter=model)
    json_results = jsonify(results)

    return json_results


# Routes for mol input
@mineserver_api.route('/structure-search/<db_name>', methods=['POST'])
@mineserver_api.route('/structure-search/<db_name>/sterio=<sterio>',
                      methods=['POST'])
# Routes for smiles input
@mineserver_api.route('/structure-search/<db_name>/smiles=<smiles>')
@mineserver_api.route('/structure-search/<db_name>/smiles=<smiles>'
                      '/stereo=<stereo>')
def structure_search_api(db_name, smiles=None, stereo=True):
    """Perform an exact structure search and return results.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular structure of query molecule.
    mol : str, optional (default: None)
        mol object in str format (should contain a lot of strange chars like
        "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from
        input structures, and cannot convert it to SMILES. Captured from form
        data.
    stereo : bool, optional (default: True)
        If true, uses sterochemistry in finding exact match.
    model : str, optional (default: None)
        KEGG organism code (e.g. 'hsa'). Adds annotations to each compound
        based on whether it is in or could be derived from the KEGG compounds
        in this organism (provided in the 'Likelihood_score' field of each
        compound document).

    Returns
    -------
    json_results : flask.Response
        JSON Document of match (empty if no match).
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
    results = structure_search(db, smiles, stereo=stereo, model_db=model_db,
                               parent_filter=model)
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
def substructure_search_api(db_name, smiles=None, limit=-1):
    """Perform a substructure search and return results.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular substructure to search for.
    mol : str, optional (default: None)
        mol object in str format (may contain a lot of whitespace hex-chars
        like "%20" and end with "END"). Used only for the MINE website backend
        because MarvinJS on the front end only generates mol objects from
        input structures, and cannot convert it to SMILES. Captured from form
        data.
    limit : int (default: -1)
        Maximum number of results (compounds) to return. By default, returns
        all results.

    Returns
    -------
    json_results : flask.Response
        JSON Documents of compounds containing given substructure.
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
    results = substructure_search(db, smiles, limit=limit, model_db=model_db,
                                  parent_filter=model)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/model-search/q=<query>')
def model_search_api(query):
    """Perform a model search and return results.

    Parameters
    ----------
    query : str
        KEGG Org Code or Org Name of model(s) to search for (e.g. 'hsa' or
        'yeast'). Can provide multiple search terms by separating each term
        with a space.  TODO: change from space delimiter to something else
    """
    db = mongo.cx[app.config['KEGG_DB_NAME']]
    results = model_search(db, query)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/database-query/<db_name>/q=<mongo_query>')
def database_query_api(db_name, mongo_query):
    """Perform a direct query built with Mongo syntax.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    mongo_query : str
        A valid Mongo query (e.g. .../q={"ID": "cpd00001"}).

    Returns
    -------
    json_results : flask.Response
        JSON Documents matching provided Mongo query.
    """
    db = mongo.cx[db_name]
    results = advanced_search(db, mongo_query)
    # TODO: add model to score_compounds (where None currently is)
    results = score_compounds(db, results, None)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-ids/<db_name>/<collection_name>')
@mineserver_api.route('/get-ids/<db_name>/<collection_name>/q=<query>')
def get_ids_api(db_name, collection_name, query=None):
    """Get Mongo IDs for a subset of a given database collection.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    collection_name : str
        Name of Mongo collection within database to query against.
    query : str (default: None)
        Specifies subset of collection to retrieve ids for. Formatted as a
        python dict as you would have in argument to db.collection.find().

    Returns
    -------
    json_results : flask.Response
        List of ids matching query in JSON format.
    """
    db = mongo.cx[db_name]
    results = get_ids(db, collection_name, query)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-comps/<db_name>', methods=['POST'])
def get_comps_api(db_name):
    """Get compounds for specified ids in database.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    id_list : list
        List of compound ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs can be either MINE IDs or Mongo IDs (_id).

    Returns
    -------
    json_results : flask.Response
        List of compound JSON documents.
    """
    id_list = request.get_json()['id_list']

    db = mongo.cx[db_name]
    results = get_comps(db, id_list)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-rxns/<db_name>', methods=['POST'])
def get_rxns_api(db_name):
    """Get reactions for specified ids in database.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    id_list : list
        List of reaction ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs must be Mongo IDs (_id).

    Returns
    -------
    json_results : flask.Response
        List of reaction JSON documents.
    """
    id_list = request.get_json()['id_list']

    db = mongo.cx[db_name]
    results = get_rxns(db, id_list)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/get-ops/<db_name>', methods=['POST'])
def get_ops_api(db_name):
    """Get operators for specified ids in database.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    id_list : list, optional (default: None)
        List of operator ids. Attach as "dict" to POST request. For example,
        requests.post(<this_uri>, data="{'id_list': ['id1', 'id2', 'id3']}").
        IDs can be either operator ids (e.g. 1.1.-1.h) or Mongo IDs (_id). If
        not provided, all operators are returned.

    Returns
    -------
    json_results : flask.Response
        List of operator JSON documents.
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

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    op_id : str
        Either operator id (e.g. 1.1.-1.h) or Mongo ID (_id) for operator.

    Returns
    -------
    json_results : flask.Response
        Operator JSON document (including associated reactions).
    """
    db = mongo.cx[db_name]
    results = get_op_w_rxns(db, op_id)
    if results:
        json_results = jsonify(results)
        return json_results
    else:
        raise InvalidUsage('Operator with ID \"{}\" not found.'.format(op_id))


@mineserver_api.route('/get-adduct-names')
@mineserver_api.route('/get-adduct-names/<adduct_type>')
def get_adduct_names_api(adduct_type='all'):
    """Get names of all adducts for the specified adduct type.

    Parameters
    ----------
    adduct_type : str
        Options are 'positive', 'negative', or 'all'.

    Returns
    -------
    json_results : flask.Response
        JSON array of adduct names. If adduct_type == 'all', then this is an
        array of two arrays, with the first element being positive adducts and
        the second adduct being negative adducts.
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

    Attach all arguments besides db_name as JSON data in POST request.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    tolerance: float
        Specifies tolerance for m/z, in mDa by default. Can specify in ppm if
        ppm is set to True.
    charge: bool
        Positive or negative mode. (True for positive, False for negative).
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str, optional (default: None)
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines.
    adducts: list, optional (default: None)
        List of adducts to use. If not specified, uses all adducts.
    models: list, optional (default: None)
        List of model _ids. If supplied, score compounds higher if  present
        in metabolic model.
    ppm: bool, optional (default: False)
        Specifies whether tolerance is in ppm.
    logp: tuple, optional (default: None)
        Length 2 tuple specifying min and max logp to filter compounds (e.g.
        (-1, 2)).
    halogens: bool, optional (default: False)
        Specifies whether to filter out compounds containing F, Cl, or Br.
        Filtered out if set to True.
    verbose: bool, optional (default: False)
        If True, verbose output.

    Returns
    -------
    json_results : flask.Response
        JSON array of compounds that match m/z within defined tolerance and
        after passing other defined filters (such as logP).
    """
    json_data = request.get_json()

    if 'tolerance' in json_data:
        tolerance = float(json_data['tolerance'])
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge' in json_data:
        charge = bool(json_data['charge'])
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

    if 'logp' in json_data:
        logp = literal_eval(str(json_data['logp']))
        assert isinstance(logp, tuple)
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

    db = mongo.cx[db_name]
    keggdb = mongo.cx[app.config['KEGG_DB_NAME']]

    results = ms_adduct_search(db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    if results:
        app.logger.info(f'MS Search successful ({len(results)} results found)')

    return json_results


@mineserver_api.route('/ms2-search/<db_name>', methods=['POST'])
def ms2_search_api(db_name):
    """Search for commpound-adducts matching precursor mass(es).

    Attach all arguments besides db_name as form data in POST request.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    tolerance: float
        Specifies tolerance for m/z, in mDa by default. Can specify in ppm if
        ppm is set to True.
    charge: bool
        Positive or negative mode. (True for positive, False for negative).
    energy_level: int
        Fragmentation energy level to use. May be 10, 20, or 40.
    scoring_function: str
        Scoring function to use. Can be either 'jaccard' or 'dot product'.
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str, optional (default: None)
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines.
    adducts: list, optional (default: None)
        List of adducts to use. If not specified, uses all adducts.
    models: list, optional (default: None)
        List of model _ids. If supplied, score compounds higher if  present
        in metabolic model.
    ppm: bool, optional (default: False)
        Specifies whether tolerance is in ppm.
    logp: tuple, optional (default: None)
        Length 2 tuple specifying min and max logp to filter compounds (e.g.
        (-1, 2)).
    halogens: bool, optional (default: False)
        Specifies whether to filter out compounds containing F, Cl, or Br.
        Filtered out if set to True.

    Returns
    -------
    json_results : flask.Response
        JSON array of compounds that match m/z within defined tolerance and
        after passing other defined filters (such as logP).
    """
    json_data = request.get_json()

    if 'tolerance' in json_data:
        tolerance = float(json_data['tolerance'])
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge' in json_data:
        charge = bool(json_data['charge'])
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

    if 'logp' in json_data:
        logp = literal_eval(json_data['logp'])
        assert isinstance(logp, tuple)
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

    results = ms2_search(db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    return json_results


@mineserver_api.route('/spectra-download/<db_name>', methods=['GET', 'POST'])
@mineserver_api.route('/spectra-download/<db_name>/q=<mongo_query>',
                      methods=['GET', 'POST'])
def spectra_download_api(db_name, mongo_query=None):
    """Download one or more spectra for compounds matching a given query.

    Parameters
    ----------
    db : Mongo DB
        Contains compound documents to search.
    mongo_query : str, optional (default: None)
        A valid Mongo query as a literal string. If None, all compound spectra
        are returned.
    parent_filter : str, optional (default: None)
        If set to a metabolic model's Mongo _id, only get spectra for compounds
        in or derived from that metabolic model.
    putative : bool, optional (default: True)
        If False, only find known compounds (i.e. in Generation 0). Otherwise,
        finds both known and predicted compounds.

    Returns
    -------
    results : flask.Response
        Text of all matching spectra, including headers and peak lists.
    """
    parent_filter = None
    putative = None
    if request.method == 'POST':
        json_data = request.get_json()

        if 'parent_filter' in json_data:
            parent_filter = json_data['parent_filter']

        if 'putative' in json_data:
            putative = bool(json_data['putative'])

    db = mongo.cx[db_name]
    results = spectra_download(db, mongo_query=mongo_query,
                               parent_filter=parent_filter, putative=putative)

    return app.response_class(results)
