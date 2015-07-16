__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices, ServerError
import numpy
import itertools

#test_compounds = open("/Users/JGJeffryes/Desktop/test comps.csv").read()
test_compounds = open('/Users/JGJeffryes/Documents/Research/Manuscripts/2015 MINE/MassBankTestSet.csv').read()
services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
#pos_unks = open("Oce_Neg_MZ").read()
#neg_unks = open("Oce_Pos_MZ").read()


def known_compound_stats(test_db, test_compounds):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    top_X_tot = 0
    top_X_new = 0
    cutoff = 5
    found = []
    hits = []
    pos_params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M+H]+', '[M]+', '[M+Na]+'], 'models': [], 'ppm': False,
                  'charge': True, 'halogens': True}
    neg_params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-', '[M+CH3COO]-'], 'models': [], 'ppm': False,
                  'charge': False, 'halogens': True}
    for i, line in enumerate(test_compounds.split('\n')[:-1]):
        print i + 1
        sl = line.split('\t')
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
                fp += 1
                hits.append(len(result))
            # found nothing
            else:
                tn += 1
            continue
        result = [x['_id'] for x in sorted(result, key=lambda x: (x['Generation'], x['NP_likeness']))]
        try:
            ind = result.index(comp['_id'])
            hits.append(len(result))
            tp += 1
            if ind < cutoff:
                top_X_tot += 1
            if not 'Names' in comp:
                found.append((sl[0], sl[6], ind))
                if ind < cutoff:
                    top_X_new += 1
        except ValueError:
            fn += 1
    print tp, fp, tn, fn
    print "P: %s" % (tp / float(tp + fp))
    print "Coverage: %s" % ((tp + fp) / float(tp + fp + tn + fn))
    print "Canidates\tMean: %s Median: %s" % (numpy.average(hits), numpy.median(hits))
    print "In top %s: %s total, %s new" % (cutoff, top_X_tot, top_X_new)

    return found

for x in known_compound_stats("KEGGexp2", test_compounds):
    pass

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