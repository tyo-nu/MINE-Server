__author__ = 'JGJeffryes'

"""
This script reads a set of unknowns in a CSV, MGF, mzXML, or excel file and queries a mongo database for potential matches.
The matches are written back to a CSV, HTML doc, or printed to STDOUT
USAGE - BatchAdductQuery.py [options] run_name file_of_Unknowns
"""

import sys
from Utils import establish_db_client
from optparse import OptionParser
import xml.etree.ElementTree as ET
import numpy
import cPickle
import time
import string
import re
from collections import defaultdict

class Dataset():
    """A class containing all the information for a metabolomics data set"""
    def __init__(self, name, options):
        self.name = name
        self.options = options  # Object containing all options from the optparser
        dtype = numpy.dtype("S20, f8, f8")
        all_pos_adducts = numpy.loadtxt('./lib/biokbase/mine_database/Positive Adducts full.txt', dtype=dtype)
        all_neg_adducts = numpy.loadtxt('./lib/biokbase/mine_database/Negative Adducts full.txt', dtype=dtype)
        if hasattr(options, 'adducts'):
            self.pos_adducts = numpy.array(filter(lambda x: x[0] in options.adducts, all_pos_adducts), dtype=dtype)
            self.neg_adducts = numpy.array(filter(lambda x: x[0] in options.adducts, all_neg_adducts), dtype=dtype)
        else:
            if hasattr(options, 'positive_adduct_file'):
                self.pos_adducts = numpy.loadtxt(options.positive_adduct_file, dtype=dtype)
            else:
                self.pos_adducts = all_pos_adducts
            if hasattr(options, 'negative_adduct_file'):
                self.neg_adducts = numpy.loadtxt(options.negative_adduct_file, dtype=dtype)
            else:
                self.neg_adducts = all_neg_adducts
        if hasattr(options, 'kovats'):
            self.min_kovats = options.kovats[0]
            self.max_kovats = options.kovats[1]
        if hasattr(options, 'logP'):
            self.min_logP = options.logP[0]
            self.max_logP = options.logP[1]
        self.hit_projection = {'Formula': 1, 'MINE_id': 1, 'logP': 1, 'minKovatsRI': 1, 'maxKovatsRI': 1,
                               'NP_likeness': 1, 'Names': 1, 'SMILES': 1, 'Inchikey': 1, 'Generation': 1}
        self.known_peaks = []  # contains Peak objects for knowns
        self.unk_peaks = []  # contains Peak objects for unknowns
        self.clusters = []  # contains tuples of formula and list of matching peaks
        self.known_set = set()  # contains InChI key of all known compounds
        self.native_set = set()  # contains _ids of compounds in model set
        self.isomers = {}  # contains a list of all compounds that are isomers of the key formula
        self.total_formulas = 0
        self.total_hits = 0
        self.matched_peaks = 0

    def __str__(self):
        return self.name

    def find_db_hits(self, peak, db, adducts):
        """This function searches the database for matches of a peak given adducts and updates the peak and data set
           file with that information."""
        #find nominal mass for a given m/z for each adducts and the max and min values for db
        potential_masses = (peak.mz - adducts['f2']) / adducts['f1']
        if self.options.ppm:
            precision = (self.options.tolerance/100000.)*potential_masses
        else:
            precision = self.options.tolerance*0.001
        upper_bounds = potential_masses + precision
        lower_bounds = potential_masses - precision

        #search database for hits in the each adducts mass range that have no innate charge.
        hits = {}

        for i, adduct in enumerate(adducts):
            # build the query by adding the optional terms
            query_terms = [{"Mass": {"$gte": float(lower_bounds[i])}}, {"Mass": {"$lte": float(upper_bounds[i])}}, {'Charge': 0}]
            if hasattr(self, 'min_logP'):
                query_terms += [{"logP": {"$gte": self.min_logP}}, {"logP": {"$lte": self.max_logP}}]
            if hasattr(self, 'min_kovats'):
                query_terms += [{"maxKovatsRI": {"$gte": self.min_kovats}}, {"minKovatsRI": {"$lte": self.max_kovats}}]
            if adduct['f0'] == '[M]+':
                query_terms[2] = {'Charge': 1}
            for compound in db.compounds.find({"$and": query_terms}, self.hit_projection):
                #Filters out halogens if the flag is enabled by moving to the next compound before the current compound
                # is counted or stored.
                if not self.options.halogens:
                    if re.search('F[^e]|Cl|Br', compound['Formula']):
                        continue

                #update the total hits for the peak and make a note if the compound is in the native_set
                peak.total_hits += 1
                if compound['_id'] in self.native_set:
                    peak.native_hit = True
                    compound['native_hit'] = True
                if compound['Generation'] < peak.min_steps:
                    peak.min_steps = compound['Generation']

                peak.formulas.add((compound['Formula'], adduct['f0']))

                if compound['Formula'] not in self.isomers:
                    self.isomers[compound['Formula']] = []
                self.isomers[compound['Formula']].append(compound)

    def annotate_peaks(self, db):
        """ This function iterates the through the unknown peaks in the data set searches the database for compounds
            that match a peak m/z given the adducts permitted. Statistics on the annotated data set are printed"""
        for i, peak in enumerate(self.unk_peaks):
            if peak.charge == '+' or peak.charge == 'Positive' or peak.charge:
                self.find_db_hits(peak, db, self.pos_adducts)

            elif peak.charge == '-' or peak.charge == 'Negative' or not peak.charge:
                self.find_db_hits(peak, db, self.neg_adducts)

            else:
                raise ValueError("Invalid compound charge specification. Please use \"+\" or \"Positive\" for positive"
                                 " ions and \"-\" or \"Negative\"for negative ions")

            if peak.total_hits > 0:
                self.matched_peaks += 1
                self.total_hits += peak.total_hits
                self.total_formulas += len(peak.formulas)
            if self.options.verbose:
                print "%s percent of peaks processed" % int(float(i)/float(len(self.unk_peaks)) * 100)

        if self.options.verbose:
            print "Proposed matches for %s of %s peaks" % (self.matched_peaks, len(self.unk_peaks))
            try:
                print "Average hits per peak: %s" % (float(self.total_hits)/float(self.matched_peaks))
                print "Average formulas per peak: %s" % (float(self.total_formulas)/float(self.matched_peaks))
            except ZeroDivisionError:
                pass

    def print_ranked_isomers(self, formula):
        #print each isomer hit to stdout
        known_isomers = []
        native_isomers = []
        other_isomers = []
        for isomer in self.isomers[formula]:
            if isomer['Inchikey'].split('-')[0] in self.known_set:
                known_isomers.append(isomer)
            elif isomer['_id'] in self.native_set:
                native_isomers.append(isomer)
            else:
                other_isomers.append(isomer)
        if len(known_isomers) > 0:
            print "Matches with known/target compounds: %s" % len(known_isomers)
            print sorted(known_isomers, key=lambda x: x['NP Likeness'], reverse=True)
        if len(native_isomers) > 0:
            print "Matches with native compounds %s" % len(native_isomers)
            #print sorted(native_isomers, key=lambda x: x['NP Likeness'], reverse=True)
        if len(other_isomers) > 0:
            print "Matches with promiscuity products %s" % len(other_isomers)
            #print sorted(other_isomers, key=lambda x: x['NP Likeness'], reverse=True)


def get_modelSEED_comps(kb_db, models):
    seed_ids, _ids = set(), set()
    for x in models:
        for z in kb_db.models.find_one({'_id': x}, {'Compounds': 1})['Compounds']:
            seed_ids.add(z)
    for y in seed_ids:
        y = int(y.split('.')[1])
        for w in kb_db.compounds.find({'Model_SEED': y}, {'_id': 1}):
            _ids.add(w['_id'])
    return _ids


def get_KEGG_comps(db, kegg_db, models):
    kegg_ids, _ids = set(), set()
    for x in models:
        for z in kegg_db.models.find_one({'_id': x}, {'Compounds': 1})['Compounds']:
            kegg_ids.add(z)
    for y in kegg_ids:
        for w in db.compounds.find({'DB_links.KEGG': y}, {'_id': 1}):
            _ids.add(w['_id'])
    return _ids


def sort_NPLike(dic_list):
    return sorted(dic_list, key=lambda x: float(x['NP_likeness']), reverse=True)


def read_mgf(input_file):

    peaks = []
    with open(input_file, 'r') as infile:
        for line in infile:
            sl = line.strip(' \r\n').split('=')
            if sl[0] == "PEPMASS":
                mass = sl[1]
            if sl[0] == "TITLE":
                name = sl[1]
            if sl[0] == "RTINSECONDS":
                r_time = sl[1]
            if sl[0] == "END IONS":
                peaks.append(Peak(name, r_time, mass, "+", {}, "False"))
    return peaks


def read_mzXML(input_file):
    peaks = []
    tree = ET.parse(input_file)
    root = tree.getroot()
    prefix = root.tag.strip('mzXML')
    #we collect info about the instrument...'cus we can.
    msInst = {}
    for child in root.find('.//%smsInstrument' % prefix):
        try:
            msInst[child.attrib['category']] = child.attrib['value']
        except KeyError:
            msInst['msSoftware'] = child.attrib
    data.Instrument = msInst

    for scan in root.findall('.//%sscan' % prefix):
        #somewhat counter intuitively we will get the peak info from the second fragments precursor info.
        if scan.attrib['msLevel'] == '2':
            precursor = scan.find('./%sprecursorMz' % prefix)
            mz = precursor.text
            r_time = scan.attrib['retentionTime'][2:-1]
            name = "%s @ %s" % (mz, r_time)
            charge = scan.attrib['polarity']
            peaks.append(Peak(name, r_time, mz, charge, {}, "False"))
    return peaks


def read_csv(input_file):
    #import peaks from a tab or comma delimited file
    peaks = []
    with open(input_file, 'r') as infile:
        for line in infile:
            #This is some jury rigged code to do some automatic parsing of tab or comma delimited files.
            #It probably should be replaced with numpy automatic parsing
            sl = line.strip('\n').split('\t')
            #if tabs split a line up in to more than two chunks, there is a pretty good chance that's not by accident
            if len(sl) > 2:
                #If it has more than 4 fields that it must have data associated with being a known compound
                if len(sl) > 4:
                    peaks.append(Peak(sl[0], sl[1], sl[2], sl[3], {sl[5]: [sl[4]]}, sl[6]))
                else:
                    peaks.append(Peak(sl[0], sl[1], sl[2], sl[3], {}, "False"))
            #parse comma separated values
            else:
                sl = line.strip('\n').split(',')
                if len(sl) > 2:
                    if len(sl) > 4:
                        peaks.append(Peak(sl[0], sl[1], sl[2], sl[3], {sl[5]: [sl[4]]}, sl[6]))
                    else:
                        peaks.append(Peak(sl[0], sl[1], sl[2], sl[3], {}, "False"))
                #Congradualtions! You broke the jank-ass parser.
                else:
                    sys.exit("Failed to parse input file")

    return peaks


class Peak:
    """A class holding information about an unknown peak"""
    def __init__(self, name, r_time, mz, charge, inchi_key):
        self.name = name
        self.r_time = float(r_time)  # retention time
        self.mz = float(mz)  # mass to charge ratio
        self.charge = charge  # polarity of charge
        self.inchi_key = inchi_key  # the id of the peak if known, as an Inchikey
        self.formulas = set()
        self.total_hits = 0
        self.native_hit = False
        self.min_steps = 99

    def __str__(self):
        return self.name

    def print_formulas(self):
        #this sends the results to stdout in a "hopefully" human readable form"
        print "%s: %s formulas and %s structures" % (self.name, len(self.formulas), self.total_hits)
        print ''
        for form_tup in self.formulas:
            print string.center(form_tup[1], 25)
            print string.ljust("Formula", 18), string.center("Isomers", 7)
            print string.ljust(form_tup[0], 18), string.center(repr(len(data.isomers[form_tup[0]])), 7)
            print ''


######################################################################################
#                                   Main code                                        #
######################################################################################
if __name__ == '__main__':
    tstart = time.time()
    #This block handles user flags and arguments. For more information see the optparse API documentation

    usage = "usage: %prog [options] run_name"
    parser = OptionParser(usage)
    parser.add_option("-p", "--positive_adducts", dest="positive_adduct_file", default="Batch Adduct Query/Positive "
                      "Adducts.txt", help="The path to the desired positive adducts file")
    parser.add_option("-n", "--negative_adducts", dest="negative_adduct_file", default="Batch Adduct Query/Negative "
                      "Adducts.txt", help="The path to the desired negative adducts file")
    parser.add_option("-d", "--database", dest="database", default="1GenKEGG",
                      help="The name of the database to search")
    parser.add_option("-m", "--modelSEED", dest="modelSEED", default="null",
                      help="The model SEED id of the organism for the data set")
    parser.add_option("-k", "--known", dest="known_file", default="null",
                      help="The path to a file containing known or targeted peaks")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float", default=2,
                      help="The m/z tolerance(precision) for peak searching")
    parser.add_option("--ppm", dest="ppm", action="store_true", default=False,
                      help="Tolerance is in Parts Per Million")
    parser.add_option("-x", "--halogens", dest="halogens", action="store_true", default=False,
                      help="include compounds containing F, Cl, and Br")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="include compounds containing F, Cl, and Br")

    (options, args) = parser.parse_args()
    #Ensure that the right number of arguments are passed to the script.
    if len(args) < 2:
        sys.exit("Missing required arguments Usage: [options] run_name, path_to_unk_peaks")
    if len(args) > 2:
        sys.exit("Too many arguments submitted. Usage: [options] run_name, path_to_unk_peaks")

    data = Dataset(args[0], options)
    unknowns_file = args[1]

    #if the user specifies known peaks import them for use in clustering and ordering of isomers
    if options.known_file != 'null':
        print "Loading known peaks"
        data.known_peaks = read_csv(options.known_file)
        for peak in data.known_peaks:
            data.known_set.add(peak.known)

    #detect unknown compounds file type and parse accordingly
    if ('.txt' or '.csv') in unknowns_file:
        data.unk_peaks = read_csv(unknowns_file)
    elif '.mgf' in unknowns_file:
        data.unk_peaks = read_mgf(unknowns_file)
    elif '.mzXML' in unknowns_file:
        data.unk_peaks = read_mzXML(unknowns_file)
    else:
        sys.exit("Unknown file type not recognised. Please use .mgf, .xlsx, .txt or .csv file")

    client = establish_db_client()
    db = client[options.database]
    if not db.compounds.count():
        sys.exit('No compounds in supplied database')

    if options.modelSEED != 'null':
        print "Loading native compounds"
        kbase_db = client['KBase']
        data.native_set = get_modelSEED_comps(kbase_db, [options.modelSEED])

    data.annotate_peaks(db)

    for i, peak in enumerate(data.unk_peaks):
        peak.print_formulas()

    tend = time.time()
    print "BatchAdductQuery.py completed in %s seconds" %(tend-tstart)
