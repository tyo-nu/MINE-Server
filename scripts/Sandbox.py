__author__ = 'JGJeffryes'

from lib.biokbase.mine_database.Client import mineDatabaseServices, ServerError
import numpy
import itertools

test_compounds = open("/Users/JGJeffryes/Documents/repository/MassBank/MassBankTestSet.csv").read()
services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
pos_unks = open("Oce_Neg_MZ").read()
neg_unks = open("Oce_Pos_MZ").read()


def known_compound_stats(test_db, test_compounds):
    fp = 0
    tn = 0
    fn = 0
    found = []
    hits = []
    for line in test_compounds.split('\n')[:-1]:
        sl = line.split('\t')
        if sl[3] == "Positive":
            result = services.batch_ms_adduct_search(test_db, sl[2], "form", 0.002, ['M+H', "M+", "M+Na"], [], False,
                                                     True, True)
        else:
            result = services.batch_ms_adduct_search(test_db, sl[2], "form", 0.002, ['M-H', 'M+CH3COO'], [], False,
                                                     False, True)
        try:
            _ids = [x['_id'] for x in services.quick_search(test_db, sl[6].strip())][0]
        except ServerError:
            if result[0]["total_hits"]:
                fp += 1
                hits.append(result[0]['total_hits'])
            else:
                tn += 1
            continue
        predicted_ids = list(itertools.chain.from_iterable(x['isomers'] for x in result[0]["adducts"]))
        if _ids in predicted_ids:
            hits.append(result[0]['total_hits'])
            found.append(predicted_ids.index(_ids)/float(len(predicted_ids)))
        else:
            fn += 1
    tp = len(found)
    print tp, fp, tn, fn
    print "P: %s" % (tp / float(tp + fp))
    print "Coverage: %s" % ((tp + fp) / float(tp + fp + tn + fn))
    print "Canidates\tMean: %s Median: %s" % (numpy.average(hits), numpy.median(hits))

    return found

print numpy.mean(known_compound_stats("KEGGexp", test_compounds))


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