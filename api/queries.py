"""Queries.py: Contains functions which power the API queries"""
import re
from ast import literal_eval
from typing import Dict, List

import pymongo
from rdkit.Chem import AllChem

from minedatabase.databases import MINE
from minedatabase.metabolomics import score_compounds

DEFAULT_PROJECTION = {
    "SMILES": 1,
    "InChI_key": 1,
    "Generation": 1,
    "MINE_id": 1,
    "KEGG_id": 1,
    "Formula": 1,
    "logP": 1,
    "Charge": 1,
    "Inchi": 1,
    "Mass": 1,
    "len_RDKit_fp": 1
}


def quick_search(
    db: MINE, core_db: MINE, query: str, search_projection: Dict[str, int] = DEFAULT_PROJECTION.copy()
) -> List:
    """This function takes user provided compound identifiers and attempts to
     find a related database ID.

    Parameters
    ----------
    db : MINE
        Database to search.
    core_db : MINE
        Core Database.
    query : str
        A MINE id, KEGG code, ModelSEED id, Inchikey or Name.
    search_projection : Dict[str, int]
        The fields which should be returned in the results.

    Returns
    -------
    results : List
        List of query results (documents in MINE database).
    """

    # Determine what kind of query was input (e.g. KEGG code, MINE id, etc.)
    # If it can't be determined to be an id or a key, assume it is the name
    # of a compound.
    if re.match(r"C\w{40}", query):
        query_field = "_id"
    elif re.match(r"C\d{5}", query):
        query_field = "KEGG_id"
    elif re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", query):
        query_field = "Inchikey"
        query = query.split("=", 1)[-1]
    elif query.isdigit():
        query_field = "MINE_id"
        query = int(query)
    else:
        query_field = None

    if query_field:
        results = [
            x
            for x in core_db.compounds.find({"$and": [{query_field: query}, {"MINES": db.name}]}, search_projection).limit(
                500
            )
            if x["_id"][0] == "C"
        ]
    else:
        results = []

    return results


def advanced_search(
    db: MINE,
    mongo_query: str,
    search_projection: Dict[str, int] = DEFAULT_PROJECTION.copy(),
) -> List:
    """Returns compounds in the indicated database which match the provided
    mongo query

    Parameters
    ----------
    db : MINE
        Database to search.
    mongo_query : str
        A query string with Mongo syntax.
    search_projection : Dict[str, int]
        The fields which should be returned in the results.

    Returns
    -------
    List
        List of query results (documents in MINE database).
    """
    # We don't want users poking around here
    if db.name == "admin" or not mongo_query:
        raise ValueError("Illegal query")
    # This transforms the string into a dictionary
    query_dict = literal_eval(mongo_query)
    return [x for x in db.compounds.find(query_dict, search_projection)]


def similarity_search(
    db: MINE,
    core_db: MINE,
    comp_structure: str,
    min_tc: float,
    limit: int,
    parent_filter: str = None,
    model_db: pymongo.database = None,
    search_projection: Dict[str, int] = DEFAULT_PROJECTION.copy(),
) -> List:
    """Returns compounds in the indicated database which have structural
     similarity to the provided compound.

    Parameters
    ----------
    db : MINE
        Database to search.
    comp_structure : str
        A molecule in molfile or SMILES format.
    min_tc : float
        Minimum Tanimoto score.
    limit : int
        The maximum number of compounds to return.
    parent_filter : str
        ID of the organism in KEGG to keep only compounds in this organism.
    model_db : pymongo.database
        MongoDB with KEGG organism codes and associated compounds.
    search_projection : Dict[str, int]
        The fields which should be returned in the results.

    Returns
    -------
    similary_search_results : List
        List of search results (documents in MINE database).
    """
    similarity_search_results = []
    fp_type = "RDKit_fp"
    # Create Mol object from Molfile (has newlines)
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    # Based on fingerprint type specified by user, get the finger print as an
    # explicit bit vector (series of 1s and 0s). Then, return a set of all
    # indices where a bit is 1 in the bit vector.
    query_fp = set(AllChem.RDKFingerprint(mol, fpSize=512).GetOnBits())

    len_fp = len(query_fp)
    # Return only id and fingerprint vector
    # search_projection[fp_type] = 1

    # Filter compounds that meet tanimoto coefficient size requirements
    search_projection['RDKit_fp'] = 1
    for x in core_db.compounds.find(
        {
            "$and": [
                {"len_" + fp_type: {"$gte": min_tc * len_fp}},
                {"len_" + fp_type: {"$lte": len_fp / min_tc}},
                {"MINES": db.name}
            ]
        },
        search_projection,
    ):
        # Put fingerprint in set for fast union (&) and intersection (|)
        # calculations
        test_fp = set(x[fp_type])
        # Calculate tanimoto coefficient
        tmc = len(query_fp & test_fp) / float(len(query_fp | test_fp))
        # If a sufficient tanimoto coefficient is calculated, append the
        # compound to the search results (until the limit is reached)
        #print(tmc)
        if tmc >= min_tc:
            similarity_search_results.append(x)
            if len(similarity_search_results) >= limit:
                break

    del search_projection['RDKit_fp']

    if parent_filter and model_db:
        similarity_search_results = score_compounds(
            similarity_search_results, parent_filter, core_db=core_db, mine_db=db,
            kegg_db=model_db, get_native=True
        )

    return similarity_search_results


def structure_search(
    db: MINE,
    core_db: MINE,
    comp_structure: str,
    parent_filter: str = None,
    model_db: pymongo.database = None,
    search_projection: Dict[str, int] = DEFAULT_PROJECTION.copy(),
) -> List:
    """Returns compounds in the indicated database which are exact matches to
    the provided structure.

    Parameters
    ----------
    db : MINE
        Database to search.
    comp_structure : str
        A molecule in molfile or SMILES format.
    parent_filter : str
        ID of the organism in KEGG to keep only compounds in this organism.
    model_db : pymongo.database
        MongoDB with KEGG organism codes and associated compounds.
    search_projection : Dict[str, int]
        The fields which should be returned in the results.

    Returns
    -------
    results : List
        List of search results (documents in MINE database).
    """
    # Create Mol object from Molfile (has newlines)
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    # Get InChI string from mol file (rdkit)
    inchi = AllChem.MolToInchi(mol)
    # Get InChI key from InChI string (rdkit)
    inchi_key = AllChem.InchiToInchiKey(inchi)
    results = [x for x in db.compounds.find(
               {"InChI_key": {"$regex": "^" + inchi_key.split("-")[0]}},
               search_projection)]

    if parent_filter and model_db:
        results = score_compounds(
            results, parent_filter, core_db=core_db, mine_db=db,
            kegg_db=model_db, get_native=True
        )

    return results


def substructure_search(
    db: MINE,
    core_db: MINE,
    sub_structure: str,
    limit: int,
    parent_filter: str = None,
    model_db: pymongo.database = None,
    search_projection: Dict[str, int] = DEFAULT_PROJECTION.copy(),
) -> List:
    """Returns compounds in the indicated database which contain the provided
    structure

    Parameters
    ----------
    db : MINE
        Database to search.
    sub_structure : str
        A compound's substructure in molfile or SMILES format.
    limit : int
        The maximum number of compounds to return.
    parent_filter : str
        ID of the organism in KEGG to keep only compounds in this organism.
    model_db : pymongo.database
        MongoDB with KEGG organism codes and associated compounds.
    search_projection : Dict[str, int]
        The fields which should be returned in the results.

    Returns
    -------
    substructure_search_results : List
        List of search results (documents in MINE database).
    """
    substructure_search_results = []
    # Create Mol object from Molfile (has newlines)
    if "\n" in sub_structure:
        mol = AllChem.MolFromMolBlock(str(sub_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(sub_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    # Based on fingerprint type specified by user, get the finger print as an
    # explicit bit vector (series of 1s and 0s). Then, return a set of all
    # indices where a bit is 1 in the bit vector.
    query_fp = list(AllChem.RDKFingerprint(mol, fpSize=512).GetOnBits())
    for x in core_db.compounds.find(
        {
            "$and": [
                {"RDKit_fp": {"$all": query_fp}},
                {"MINES": db.name}
            ]
        },
        search_projection):

        # Get Mol object from SMILES string (rdkit)
        comp = AllChem.MolFromSmiles(x["SMILES"])
        # Use HasSubstructMatch (rdkit) to determine if compound has a
        # specified substructure. If so, append it to the results (until
        # limit).
        if comp and comp.HasSubstructMatch(mol):
            substructure_search_results.append(x)
            if len(substructure_search_results) >= limit:
                break

    if parent_filter and model_db:
        substructure_search_results = score_compounds(
            substructure_search_results, parent_filter, core_db=core_db, mine_db=db,
            kegg_db=model_db, get_native=True
        )

    return substructure_search_results


def model_search(db: pymongo.database, query: str) -> List[str]:
    """Returns models that match a given KEGG Org Code query (e.g. 'hsa').

    Parameters
    ----------
    db : Mongo DB
        DB to search. Should contain a 'models' collection. '_id' and 'Name'
        fields should be indexed for text search.
    query : str
        KEGG Org Code or Org Name of model(s) to search for (e.g. 'hsa' or
        'yeast'). Can provide multiple search terms by separating each term
        with a space.  TODO: change from space delimiter to something else

    Returns
    -------
    model_ids : str
        List of KEGG Org codes (_id) found in the DB matching search query.
    """
    cursor = db.models.find(
        {"$text": {"$search": query}}, {"score": {"$meta": "textScore"}, "_id": 1, "Name": 1}
    )
    model_ids = [x["_id"] for x in cursor.sort([("score", {"$meta": "textScore"})])]
    return model_ids


def get_ids(db: pymongo.database, collection: str, query: str) -> List[str]:
    """Returns ids for documents in database collection matching query.

    Parameters
    ----------
    db : pymongo.database
        DB to search.
    collection : str
        Name of collection within db to query.
    query : str
        KEGG Org Code or Org Name of model(s) to search for (e.g. 'hsa' or
        'yeast'). Can provide multiple search terms by separating each term
        with a space.  TODO: change from space delimiter to something else

    Returns
    -------
    ids : List[str]
        List of KEGG Org codes (_id) found in the DB matching search query.
    """
    if query:
        query = literal_eval(query)
    else:
        query = {}
    ids = [x["_id"] for x in db[collection].find(query)]

    return ids


def get_comps(db: MINE, id_list: List[str], core_db: MINE) -> List:
    """Returns compounds with associated IDs from a Mongo database.

    Parameters
    ----------
    db : MINE
        DB to search.
    id_list : List[str]
        IDs to get compound documents for.

    Returns
    -------
    compounds : List
        List of compound documents with specified IDs.
    """
    compounds = []
    for cpd_id in id_list:
        if isinstance(cpd_id, int):
            cpd = core_db.compounds.find_one({"MINE_id": cpd_id})
            cpd_id = cpd['_id']
        cpd = db.compounds.find_one({"_id": cpd_id})
        # New MINEs won't have this precomputed
        if cpd and "Reactant_in" not in cpd and "Product_of" not in cpd:
            rxns_as_sub = db.reactions.find({"Reactants.c_id": cpd["_id"]})
            cpd["Reactant_in"] = [x["_id"] for x in rxns_as_sub]
            rxns_as_prod = db.reactions.find({"Products.c_id": cpd["_id"]})
            cpd["Product_of"] = [x["_id"] for x in rxns_as_prod]
        compounds.append(cpd)
    return compounds


def get_rxns(db: MINE, id_list: List[str]) -> List:
    """Returns reactions with associated IDs from a Mongo database.

    Parameters
    ----------
    db : MINE
        DB to search.
    id_list : List[str]
        IDs to get compound documents for.

    Returns
    -------
    reactions : List
        List of reaction documents with specified IDs.
    """
    reactions = []
    all_reactions = db.reactions.find({"_id": {"$in": id_list}})

    for rxn in all_reactions:
        for i, cpds in enumerate(rxn['SMILES_rxn'].split('=>')):
            if i == 0:  # reactants
                for reactant, reactant_array in zip(cpds.split('+'), rxn['Reactants']):
                    reactant = reactant.strip()
                    reactant_smiles = reactant.split(' ')[-1]  # removes stoich from beginning
                    reactant_array.append(reactant_smiles)
            elif i == 1:
                for product, product_array in zip(cpds.split('+'), rxn['Products']):
                    product = product.strip()
                    product_smiles = product.split(' ')[-1]  # removes stoich from beginning
                    product_array.append(product_smiles)
        reactions.append(rxn)

    return reactions


def get_rxns_for_cpd(db: MINE, cpd_id: str, mode: str = 'product_of') -> List:
    """Returns reactions either being producing or consuming a compound.

    Parameters
    ----------
    db : MINE
        DB to search.
    cpd_id : str
        Mongo ID of compound.
    mode : str
        If 'product_of', get reactions producing this compound. If
        'reactant_in', get reactions consuming this compound.

    Returns
    -------
    reactions : List
        List of reaction documents producing or consuming given compound.
    """
    reaction_ids = []
    reactions = []
    rxn_list_docs = db[mode].find({"c_id": cpd_id})
    field_name = mode.capitalize()
    for rxn_list_doc in rxn_list_docs:
        reaction_ids += rxn_list_doc[field_name]

    all_reactions = db.reactions.find({"_id": {"$in": reaction_ids}})

    for rxn in all_reactions:
        for i, cpds in enumerate(rxn['SMILES_rxn'].split('=>')):
            if i == 0:  # reactants
                for reactant, reactant_array in zip(cpds.split('+'), rxn['Reactants']):
                    reactant = reactant.strip()
                    reactant_smiles = reactant.split(' ')[-1]  # removes stoich from beginning
                    reactant_array.append(reactant_smiles)
            elif i == 1:
                for product, product_array in zip(cpds.split('+'), rxn['Products']):
                    product = product.strip()
                    product_smiles = product.split(' ')[-1]  # removes stoich from beginning
                    product_array.append(product_smiles)
        reactions.append(rxn)

    return reactions


def get_ops(db: MINE, operator_ids: List[str]) -> List:
    """Returns operators from a Mongo database.

    Parameters
    ----------
    db : MINE
        DB to search.
    operator_ids : List[str]
        IDs to get operator documents for (e.g. "1.1.-1.h").

    Returns
    -------
    operators : List
        List of operator documents with specified IDs.
    """
    if not operator_ids:
        operators = [op for op in db.operators.find()]
    else:
        operators = []
        for op_id in operator_ids:
            op = db.operators.find_one({"$or": [{"_id": op_id}, {"Name": op_id}]})
            operators.append(op)

    return operators


def get_op_w_rxns(db: MINE, operator_id: str, limit: int = 100) -> Dict:
    """Returns operator with all its associated reactions.

    Parameters
    ----------
    db : MINE
        DB to search.
    operator_id : str
        Mongo _id or operator name (e.g. rule0001).
    limit : int
        Max number of reaction _ids to return with operator.

    Returns
    -------
    operator : Dict
        Operator JSON document (with associated reactions).
    """
    operator = db.operators.find_one(
        {"$or": [{"_id": operator_id}, {"Name": operator_id}]}
    )
    if operator:
        op_rxns = db.reactions.find({"Operators": operator_id}).limit(limit)
        operator["Reaction_ids"] = list(set([op_rxn['_id'] for op_rxn in op_rxns]))
    else:
        return None

    return operator
