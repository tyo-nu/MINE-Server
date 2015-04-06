__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices, ServerError
import numpy
import itertools

#test_compounds = open("/Users/JGJeffryes/Desktop/test comps.csv").read()
test_compounds = open('testMasses.csv').read()
services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
#pos_unks = open("Oce_Neg_MZ").read()
#neg_unks = open("Oce_Pos_MZ").read()


def known_compound_stats(test_db, test_compounds):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    found = []
    hits = []
    pos_params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M+H]+', '[M]+', '[M+Na]+'], 'models': [], 'ppm': False,
                  'charge': True, 'halogens': True}
    neg_params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-', '[M+CH3COO]-'], 'models': [], 'ppm': False,
                  'charge': False, 'halogens': True}
    for line in test_compounds.split('\n')[:-1]:
        sl = line.split(',')
        # try to find the compound by accurate mass
        if sl[3] == "Positive":
            result = services.ms_adduct_search(sl[2], "form", pos_params)
        else:
            result = services.ms_adduct_search(sl[2], "form", neg_params)
        # try to find the compound through Inchikey
        try:
            comp = services.quick_search(test_db, sl[6])[0]
        except ServerError:
            # found results w/ accurate mass but not Inchikey
            if result:
                result.sort(key=lambda x: x['steps_from_source'])
                if result[0]['steps_from_source']:
                    tp += 1
                fp += 1
                hits.append(len(result))
            # found nothing
            else:
                tn += 1
            continue
        #predicted_ids = set(y['_id'] for y in result)
        #if comp['_id'] in predicted_ids:
        result = [x['_id'] for x in sorted(result, key=lambda x: (x['steps_from_source'], x['NP_likeness']))]
        try:
            ind = result.index(comp['_id'])
            hits.append(len(result))
            tp += 1
            if not 'Name' in comp:
                found.append((sl[0], sl[6], ind))
        except ValueError:
            fn += 1
    print tp, fp, tn, fn
    print "P: %s" % (tp / float(tp + fp))
    print "Coverage: %s" % ((tp + fp) / float(tp + fp + tn + fn))
    print "Canidates\tMean: %s Median: %s" % (numpy.average(hits), numpy.median(hits))

    return found

for x in known_compound_stats("KEGGexp2", test_compounds):
    print(x)

"""

pos_result = services.batch_ms_adduct_search(test_db, pos_unks, "form", 0.003, ['M+H', "M+", "M+Na", "M+NH4"], ['kb|fm.1697'], False, True, False)
neg_result = services.batch_ms_adduct_search(test_db, neg_unks, "form", 0.003, ['M-H', 'M+CH3COO'], ['kb|fm.1697'], False, False, False)
hits = []
native = 0
total = float(len(pos_result+neg_result))
for peak in pos_result+neg_result:
    if peak['total_hits']:
        hits.append(peak['total_hits'])
    if peak['native_hit']:
        native += 1

print len(hits), len(hits)/total
print native, native/total
print numpy.average(hits), numpy.median(hits)
"""