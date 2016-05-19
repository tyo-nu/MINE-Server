__author__ = 'JGJeffryes'

import platform
from pymongo import MongoClient

def establish_db_client():
    """This establishes a mongo database client for many environments"""
    try:
        # special case for working on a sapphire node
        if 'node' in platform.node():
            client = MongoClient(host='master')
        # special case for working on a SEED cluster
        elif ('bio' in platform.node()) or ('twig' in platform.node()):
            client = MongoClient(host='branch')
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        elif ('branch' in platform.node()):
            client = MongoClient()
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        # local database
        else:
            client = MongoClient()
    except:
        raise IOError("Failed to load database client. Please verify that mongod is running")
    return client


def quick_search(db, comp_data, search_projection={}):
    """This function takes user provided compound identifiers and attempts to find a related database ID"""
    # check if comp_data already is a _id
    if (len(comp_data) == 41) and (comp_data[0] == 'C'):
        query_field = '_id'
    elif (len(comp_data) == 6) and (comp_data[0] == 'C'):
        query_field = 'DB_links.KEGG'
    elif (len(comp_data) == 8) and (comp_data[0:2] == 'cpd'):
        query_field = 'DB_links.Model_SEED'
    elif len(comp_data.split('-')[0]) == 14 and comp_data.isupper():
        query_field = 'Inchikey'
        comp_data = comp_data.split('-')[0]
    elif comp_data.isdigit():
        query_field = "MINE_id"
        comp_data = int(comp_data)
    else:
        query_field = 'Names'

    if query_field == 'Inchikey':
        results = [x for x in db.compounds.find({query_field: {'$regex': '^'+comp_data}}, search_projection).limit(500)
                   if x['_id'][0] == "C"]
    elif query_field == 'Names':
        results = [x for x in db.compounds.find({"Names": {'$regex': '^'+comp_data+'$', '$options': 'i'}}, search_projection) if x['_id'][0] == "C"]
        if not results:
            cursor = db.compounds.find({"$text": {"$search": comp_data}}, {"score": {"$meta": "textScore"}, 'Formula': 1,
                                                                           'MINE_id': 1, 'Names': 1, 'Inchikey': 1,
                                                                           'SMILES': 1, 'Mass': 1})
            results.extend(x for x in cursor.sort([("score", {"$meta": "textScore"})]).limit(500) if x['_id'][0] == "C")
    else:
        results = [x for x in db.compounds.find({query_field: comp_data}, search_projection).limit(500)
                   if x['_id'][0] == "C"]
    if not results:
        raise ValueError("%s was not found in the database." % comp_data)

    return results


def score_compounds(db, compounds, model_id, parent_frac=0.5, reaction_frac=0.5):
    """This function validates compounds against a metabolic model, returning only the compounds which pass"""
    if not model_id:
        return compounds
    model = db.models.find_one({"_id": model_id})
    parents = set(model["Compound_ids"])
    operators = dict((x[0], x[1]) for x in model['Operators'])

    for comp in compounds:
        if comp['_id'] in parents:
            comp['Likelihood_score'] = parent_frac+reaction_frac
            continue
        elif comp['Generation'] == 0:
            comp['Likelihood_score'] = reaction_frac
            continue
        else:
            comp['Likelihood_score'] = 0.0

        for source in comp['Sources']:
            ls = reaction_frac
            for op in source['Operators']:
                if op in operators:
                    ls *= operators[op]
                else:
                    ls *= 0

            if source['Compound'] in parents:
                ls += parent_frac

            if ls > comp['Likelihood_score']:
                comp['Likelihood_score'] = ls

    return compounds


def print_sorted_dict(dict):
    list = [x for x in dict]
    for x in sorted(list, key=lambda x: dict[x], reverse=True):
        print("%s\t%s" % (x, dict[x]))


def uni_dict_2_string_dict(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary

def approximate_matches(list1, list2, epsilon=0.01):
    """
    Takes two list of tuples and searches for matches of tuples first value within the supplied epsilon. Emits tuples
    with the tuples second values where found. if a value in one dist does not match the other list, it is emitted alone.
    :param list1: first list of tuples
    :type list1: list
    :param list2: second list of tuples
    :type list2: list
    :param epsilon: maximum difference in
    :type epsilon: float
    :return: second values of tuples
    :rtype: generator
    """
    list1.sort()
    list2.sort()
    list1_index = 0
    list2_index = 0

    while list1_index < len(list1) or list2_index < len(list2):
        if list1_index == len(list1):
            yield (0, list2[list2_index][1])
            list2_index += 1
            continue
        if list2_index == len(list2):
            yield (list1[list1_index][1], 0)
            list1_index += 1
            continue

        list1_element = list1[list1_index][0]
        list2_element = list2[list2_index][0]

        difference = abs(list1_element - list2_element)

        if difference < epsilon:
            yield (list1[list1_index][1], list2[list2_index][1])
            list1_index += 1
            list2_index += 1
        elif list1_element < list2_element:
            yield (list1[list1_index][1], 0)
            list1_index += 1
        elif list2_element < list1_element:
            yield (0, list2[list2_index][1])
            list2_index += 1
        else:
            raise AssertionError('Unexpected else taken')