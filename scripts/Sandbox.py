__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices, ServerError
import numpy

services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')


def known_compound_stats(test_db, test_compounds):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    top_X_tot = 0
    top_X_new = 0
    cutoff = 3
    found = []
    hits = []
    params = {'db': test_db, 'tolerance': 10.0, 'adducts': ['[M+H]+'], 'models': ['Bacteria'], 'ppm': False,
              'charge': True, 'halogens': True, 'scoring_function': 'jacquard', 'energy_level': 1}
    neg_params = {'db': test_db, 'tolerance': 5.0, 'adducts': ['[M-H]-'], 'models': [], 'ppm': False,
                  'charge': False, 'halogens': True}
    for i, compound in enumerate(test_compounds.split('\n\n')[:-1]):
        print i + 1
        result = services.ms2_search(compound, "msp", params)
        ids = [x['Inchikey'] for x in result]
        scores = [x['Spectral_score'] for x in result]
        # try to find the compound through Inchikey
        try:
            code = compound.split('calculated InChI Key: ')[1].split('\n')[0]
            comp = services.quick_search(test_db, code)[0]
        except ServerError:
            # found results w/ accurate mass but not Inchikey
            if result:
                fp += 1
                hits.append(len(result))
            # found nothing
            else:
                tn += 1
            continue
        try:
            ind = ids.index(comp['Inchikey'])
            print scores[ind]
            hits.append(len(result))
            tp += 1
            if ind < cutoff:
                top_X_tot += 1
            if comp['Generation']:
                top_X_new += 1
        except ValueError:
            fn += 1
            print comp['Inchikey']
    print tp, fp, tn, fn
    print "P: %s" % (tp / float(tp + fp))
    print "Coverage: %s" % ((tp + fp) / float(tp + fp + tn + fn))
    print "Canidates\tMean: %s Median: %s" % (numpy.average(hits), numpy.median(hits))
    print "In top %s: %s total, %s new" % (cutoff, top_X_tot, top_X_new)

    return found

with open('/Users/JGJeffryes/Desktop/pos_MoNA_20eV.msp') as test_compounds:
    for x in known_compound_stats("KEGGexp2", test_compounds.read()):
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