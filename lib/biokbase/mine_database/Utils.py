__author__ = 'JGJeffryes'

import platform
from pymongo import MongoClient
import re
import os


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
        query_field = 'Model_SEED'
    elif (len(comp_data.split('-')) == 3) or (len(comp_data) == 14):
        query_field = 'Inchikey'
        comp_data = comp_data.split('-')[0]
    else:
        query_field = 'Names'

    if query_field == 'Inchikey':
        results = [x for x in db.compounds.find({query_field: {'$regex': '^'+comp_data}},
                                                {'Formula': 1, 'Model_SEED': 1, 'Names': 1}) if x['_id'][0] == "C"]
    elif query_field == 'Names':
        results = [x['obj'] for x in db.command("text", "compounds", search=comp_data,
                                                   project={'Formula': 1, 'Model_SEED': 1, 'Names': 1})['results']
                   if x['obj']['_id'][0] == "C"]
    else:
        results = [x for x in db.compounds.find({query_field: comp_data},
                                                {'Formula': 1, 'Model_SEED': 1, 'Names': 1}) if x['_id'][0] == "C"]
    if not results:
        raise ValueError("%s was not found in the database." % comp_data)

    return results


def dict_push(dict, key, value):
    """This is a convenient way to build dictionaries of lists"""
    if not key in dict:
            dict[key] = [value]
    else:
        dict[key].append(value)


def print_sorted_dict(dict):
    list = [x for x in dict]
    for x in sorted(list, key=lambda x: dict[x], reverse=True):
        print "%s\t%s" % (x, dict[x])


def uni_dict_2_string_dict(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary
