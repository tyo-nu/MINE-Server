"""Here, routes are defined for all possible API requests. Note that nearly
all actual logic is imported from the minedatabase package."""

from ast import literal_eval

from flask import jsonify, request

from app import app, mongo
from app.exceptions import InvalidUsage
from minedatabase.metabolomics import (ms_adduct_search, read_adduct_names,
                                       spectra_download)
from minedatabase.queries import (advanced_search, get_comps, get_ids,
                                  get_op_w_rxns, get_ops, get_rxns,
                                  quick_search, similarity_search,
                                  structure_search, substructure_search)
from minedatabase.utils import score_compounds


@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    """Makes it so user can receive an informative error message rather than
    a default internal server error."""
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response


@app.route('/mineserver/quick-search/<db_name>/q=<query>')
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


@app.route('/mineserver/similarity-search/<db_name>/smiles=<smiles>')
@app.route('/mineserver/similarity-search/<db_name>/smiles=<smiles>'
           '/<float:min_tc>')
@app.route('/mineserver/similarity-search/<db_name>/smiles=<smiles>'
           '/<int:limit>')
@app.route('/mineserver/similarity-search/<db_name>/smiles=<smiles>'
           '/<float:min_tc>/<int:limit>')
def similarity_search_api(db_name, smiles, min_tc=0.7, limit=-1):
    """Perform a similarity search and return results.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular structure of query molecule.
    min_tc : float (default: 0.7)
        Minimum Tanimoto Coefficient required for similarity match.
    limit : int (default: -1)
        Maximum number of results (compounds) to return. By default, returns
        all results.

    Returns
    -------
    json_results : flask.Response
        JSON Documents of similar compounds.
    """
    db = mongo.cx[db_name]
    results = similarity_search(db, smiles, min_tc=min_tc, limit=limit)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/structure-search/<db_name>/smiles=<smiles>')
@app.route('/mineserver/structure-search/<db_name>/smiles=<smiles>'
           '/stereo=<stereo>')
def structure_search_api(db_name, smiles, stereo=True):
    """Perform an exact structure search and return results.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular structure of query molecule.
    stereo : bool (default: True)
        If true, uses sterochemistry in finding exact match.

    Returns
    -------
    json_results : flask.Response
        JSON Document of match (empty if no match).
    """
    db = mongo.cx[db_name]
    results = structure_search(db, smiles, stereo=stereo)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/substructure-search/<db_name>/smiles=<smiles>')
@app.route('/mineserver/substructure-search/<db_name>/smiles=<smiles>'
           '/<int:limit>')
def substructure_search_api(db_name, smiles, limit=-1):
    """Perform a substructure search and return results.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    smiles : str
        SMILES string describing molecular substructure to search for.
    limit : int (default: -1)
        Maximum number of results (compounds) to return. By default, returns
        all results.

    Returns
    -------
    json_results : flask.Response
        JSON Documents of compounds containing given substructure.
    """
    db = mongo.cx[db_name]
    results = substructure_search(db, smiles, limit=limit)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/database-query/<db_name>/q=<mongo_query>')
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


@app.route('/mineserver/get-ids/<db_name>/<collection_name>')
@app.route('/mineserver/get-ids/<db_name>/<collection_name>/q=<query>')
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


@app.route('/mineserver/get-comps/<db_name>', methods=['POST'])
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
    id_list = literal_eval(request.form.get('id_list'))

    db = mongo.cx[db_name]
    results = get_comps(db, id_list)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/get-rxns/<db_name>', methods=['POST'])
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
    id_list = literal_eval(request.form.get('id_list'))

    db = mongo.cx[db_name]
    results = get_rxns(db, id_list)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/get-ops/<db_name>', methods=['POST'])
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
    if 'id_list' in request.form:
        id_list = literal_eval(request.form.get('id_list'))
    else:
        id_list = None

    db = mongo.cx[db_name]
    results = get_ops(db, id_list)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/get-op-w-rxns/<db_name>/<op_id>')
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


@app.route('/mineserver/get-adduct-names/<adduct_type>')
def get_adduct_names_api(adduct_type):
    """Get names of all adducts for the specified adduct type.

    Parameters
    ----------
    adduct_type : str
        Options are 'positive' or 'negative'.

    Returns
    -------
    json_results : flask.Response
        JSON array of adduct names.
    """

    if adduct_type.lower() == 'positive':
        results = read_adduct_names(app.config['POS_ADDUCT_PATH'])
    elif adduct_type.lower() == 'negative':
        results = read_adduct_names(app.config['NEG_ADDUCT_PATH'])
    else:
        raise InvalidUsage('URI argument <adduct_type> must be "positive" or '
                           '"negative".')

    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/ms-adduct-search/<db_name>', methods=['POST'])
def ms_adduct_search_api(db_name):
    """Search for commpound-adducts matching precursor mass(es).

    Attach all arguments besides db_name as form data in POST request.

    Parameters
    ----------
    db_name : str
        Name of Mongo database to query against.
    tolerance: float
        Specifies tolerance for m/z, in mDa by default. Can specify in ppm if
        ppm is set to True.
    charge_mode: bool
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
    kovats: tuple, optional (default: None)
        Length 2 tuple specifying min and max kovats retention index to filter
        compounds (e.g. (500, 1000)).
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
        after passing other defined filters (such as kovats or logP).
    """
    if 'tolerance' in request.form:
        tolerance = float(request.form.get('tolerance'))
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge_mode' in request.form:
        charge = int(request.form.get('charge_mode'))
    else:
        raise InvalidUsage('<charge_mode> argument must be specified. True '
                           'for positive mode, False for negative mode.')

    if 'text' in request.form:
        text = request.form.get('text')
    else:
        raise InvalidUsage('<text> argument must be specified.')

    if 'text_type' in request.form:
        text_type = request.form.get('text_type')
    else:
        text_type = None

    if 'adducts' in request.form:
        adducts = literal_eval(request.form.get('adducts'))
        assert isinstance(adducts, list)
    else:
        adducts = None

    if 'models' in request.form:
        models = literal_eval(request.form.get('models'))
        assert isinstance(models, list)
    else:
        models = None

    if 'ppm' in request.form:
        ppm = bool(request.form.get('ppm'))
    else:
        ppm = None

    if 'kovats' in request.form:
        kovats = literal_eval(request.form.get('kovats'))
        assert isinstance(kovats, tuple)
    else:
        kovats = None

    if 'logp' in request.form:
        logp = literal_eval(request.form.get('logp'))
        assert isinstance(logp, tuple)
    else:
        logp = None

    if 'halogens' in request.form:
        halogens = bool(request.form.get('halogens'))
    else:
        halogens = None

    ms_params = {
        'tolerance': tolerance,
        'charge': charge,
        'adducts': adducts,
        'models': models,
        'ppm': ppm,
        'kovats': kovats,
        'logp': logp,
        'halogens': halogens
    }

    db = mongo.cx[db_name]
    keggdb = mongo.cx[app.config['KEGG_DB_NAME']]

    results = ms_adduct_search(db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/ms2-search/<db_name>', methods=['POST'])
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
    charge_mode: bool
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
    kovats: tuple, optional (default: None)
        Length 2 tuple specifying min and max kovats retention index to filter
        compounds (e.g. (500, 1000)).
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
        after passing other defined filters (such as kovats or logP).
    """
    if 'tolerance' in request.form:
        tolerance = float(request.form.get('tolerance'))
    else:
        raise InvalidUsage('<tolerance> argument must be specified (in mDa).')

    if 'charge_mode' in request.form:
        charge = int(request.form.get('charge_mode'))
    else:
        raise InvalidUsage('<charge_mode> argument must be specified. True '
                           'for positive mode, False for negative mode.')

    if 'energy_level' in request.form:
        energy_level = int(request.form.get('energy_level'))
    else:
        raise InvalidUsage('<energy_level> argument must be specified. '
                           'Possible values are 10, 20, or 40.')

    if 'scoring_function' in request.form:
        scoring_function = request.form.get('scoring_function')
    else:
        raise InvalidUsage("<scoring_function> argument must be specified. "
                           "Possible values are 'jaccard' and 'dot product'.")

    if 'text' in request.form:
        text = request.form.get('text')
    else:
        raise InvalidUsage('<text> argument must be specified.')

    if 'text_type' in request.form:
        text_type = request.form.get('text_type')
    else:
        text_type = None

    if 'adducts' in request.form:
        adducts = literal_eval(request.form.get('adducts'))
        assert isinstance(adducts, list)
    else:
        adducts = None

    if 'models' in request.form:
        models = literal_eval(request.form.get('models'))
        assert isinstance(models, list)
    else:
        models = None

    if 'ppm' in request.form:
        ppm = bool(request.form.get('ppm'))
    else:
        ppm = None

    if 'kovats' in request.form:
        kovats = literal_eval(request.form.get('kovats'))
        assert isinstance(kovats, tuple)
    else:
        kovats = None

    if 'logp' in request.form:
        logp = literal_eval(request.form.get('logp'))
        assert isinstance(logp, tuple)
    else:
        logp = None

    if 'halogens' in request.form:
        halogens = bool(request.form.get('halogens'))
    else:
        halogens = None

    ms_params = {
        'tolerance': tolerance,
        'charge': charge,
        'energy_level': energy_level,
        'scoring_function': scoring_function,
        'adducts': adducts,
        'models': models,
        'ppm': ppm,
        'kovats': kovats,
        'logp': logp,
        'halogens': halogens
    }

    db = mongo.cx[db_name]
    keggdb = mongo.cx[app.config['KEGG_DB_NAME']]

    results = ms_adduct_search(db, keggdb, text, text_type, ms_params)
    json_results = jsonify(results)

    return json_results


@app.route('/mineserver/spectra-download/<db_name>')
@app.route('/mineserver/spectra-download/<db_name>/q=<mongo_query>')
def spectra_download_api(db_name, mongo_query):
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
    if 'parent_filter' in request.form:
        parent_filter = request.form.get('parent_filter')
    else:
        parent_filter = None

    if 'putative' in request.form:
        putative = bool(request.form.get('putative'))
    else:
        putative = None

    db = mongo.cx[db_name]
    results = spectra_download(db, mongo_query=mongo_query,
                               parent_filter=parent_filter, putative=putative)

    return app.response_class(results)
