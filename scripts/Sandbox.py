__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices
import numpy
import itertools

services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
pos_unks = open("Oce_Pos_MZ").read()
neg_unks = open("Oce_Neg_MZ").read()
test_db = "KEGGdb"
"""
test_compounds = open("/Users/JGJeffryes/Desktop/test comps.csv").read()
for line in test_compounds.split('\r')[:-1]:
    sl = line.split(',')
    result = services.batch_ms_adduct_search(test_db, sl[2], "form", 0.003, [sl[5]], [], False, sl[3] == "Positive", False)
    try:
        _ids = set([x['_id'] for x in services.quick_search(test_db, sl[6].strip())])
    except:
        print "%s not in db!" % sl[0]
        continue
    predicted_ids = set(itertools.chain.from_iterable(x['isomers'] for x in result[0]["adducts"]))
    if _ids & predicted_ids:
        print "Found correct id in list of %s possibilities" % result[0]['total_hits']
    else:
        print "missed %s!" % sl[0]

"""
pos_result = services.batch_ms_adduct_search(test_db, pos_unks, "form", 0.003, ['M+H', "M+Na", "M+NH4"], ['kb|fm.1697'], False, True, False)
neg_result = services.batch_ms_adduct_search(test_db, neg_unks, "form", 0.003, ['M-H', 'M+CH3COO'], ['kb|fm.1697'], False, False, False)
hits = []
for peak in pos_result+neg_result:
    if peak['total_hits']:
        hits.append(peak['total_hits'])
print len(hits), len(hits)/float(len(pos_result+neg_result))
print numpy.average(hits), numpy.median(hits)