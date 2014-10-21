__author__ = 'JGJeffryes'

import platform
from pymongo import MongoClient

def establish_db_client():
    """This establishes a mongo database client for many environments"""
    try:
        #special case for working on a sapphire node
        if 'node' in platform.node():
            client = MongoClient(host='master')
        #special case for working on a SEED cluster
        elif ('bio' in platform.node()) or ('twig' in platform.node()):
            client = MongoClient(host='branch')
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        elif ('branch' in platform.node()):
            client = MongoClient()
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        #local database
        else:
            client = MongoClient()
    except:
        raise IOError("Failed to load database client. Please verify that mongod is running")
    return client


def quick_search(db, comp_data):
    """This function takes user provided compound identifiers and attempts to find a related database ID"""
    #check if comp_data already is a _id
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
        results = [x for x in db.compounds.find({query_field: {'$regex': '^'+comp_data}},
                                                {'Formula': 1, 'MINE_id': 1, 'Names': 1}) if x['_id'][0] == "C"]
    elif query_field == 'Names':
        cursor = db.compounds.find({"$text": {"$search": comp_data}}, {"score": {"$meta": "textScore"}, 'Formula': 1, 'MINE_id': 1, 'Names': 1})
        results = [x for x in cursor.sort([("score", {"$meta": "textScore"})]) if x['_id'][0] == "C"]
    else:
        results = [x for x in db.compounds.find({query_field: comp_data},
                                                {'Formula': 1, 'MINE_id': 1, 'Names': 1}) if x['_id'][0] == "C"]
    if not results:
        raise ValueError("%s was not found in the database." % comp_data)

    return results


def print_sorted_dict(dict):
    list = [x for x in dict]
    for x in sorted(list, key=lambda x: dict[x], reverse=True):
        print "%s\t%s" % (x, dict[x])


def uni_dict_2_string_dict(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary
