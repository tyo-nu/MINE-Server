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
        elif 'bio' in platform.node():
            client = MongoClient(host='branch')
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        #local database
        else:
            client = MongoClient()
    except:
        raise IOError("Failed to load database client. Please verify that mongod is running")
    return client


def get_id(db, comp_data):
    """This function takes user provided compound identifiers and attempts to find a related database ID"""
    #check if comp_data already is a _id
    if (len(comp_data) == 41) and (comp_data[0] == 'C'):
        return comp_data
    elif (len(comp_data) == 6) and (comp_data[0] == 'C'):
        query_field = 'KEGG_code'
    elif (len(comp_data) == 8) and (comp_data[0:2] == 'cpd'):
        query_field = 'ModelSEED_id'
    elif (len(comp_data.split('-')) == 3) or (len(comp_data) == 14):
        query_field = 'Inchi_key'
        comp_data = comp_data.split('-')[0]
    else:
        query_field = 'Names'

    if query_field == 'Inchi_key':
        compound_id = db.compounds.find_one({query_field: {'$regex': comp_data}}, {'_id': 1})
    else:
        compound_id = db.compounds.find_one({query_field: comp_data}, {'_id': 1})
    if not compound_id:
        raise ValueError("%s was not found in the database." % comp_data)

    return compound_id['_id']


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
