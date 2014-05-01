__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices
import numpy


services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
pos_unks = open("PosUnk_MZ.txt").read()
test_db = "KEGGexp"

result = services.batch_ms_adduct_search(test_db, pos_unks, "form", 0.005, ['M+H'], [], False, True, False)
hits = []
for peak in result:
    if peak['total_hits']:
        hits.append(peak['total_hits'])
print len(hits), len(hits)/float(len(result))
print numpy.average(hits), numpy.median(hits)