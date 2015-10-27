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
import re
import math
import Utils

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
                               'NP_likeness': 1, 'Names': 1, 'SMILES': 1, 'Inchikey': 1, 'Generation': 1,
                               'CFM_spectra': 1}
        self.known_peaks = []  # contains Peak objects for knowns
        self.unk_peaks = []  # contains Peak objects for unknowns
        self.clusters = []  # contains tuples of formula and list of matching peaks
        self.known_set = set()  # contains InChI key of all known compounds
        self.native_set = set()  # contains _ids of compounds in model set
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
                peak.formulas.add(compound['Formula'])
                compound['adduct'] = adduct['f0']
                compound['peak_name'] = peak.name
                peak.isomers.append(compound)

    def annotate_peaks(self, db):
        """ This function iterates the through the unknown peaks in the data set searches the database for compounds
            that match a peak m/z given the adducts permitted. Statistics on the annotated data set are printed"""
        for i, peak in enumerate(self.unk_peaks):
            if peak.charge == '+' or peak.charge == 'Positive' or (peak.charge and isinstance(peak.charge, bool)):
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


def read_mgf(input_string, charge):
    peaks = []
    ms2 = []
    for line in input_string.split('\n'):
        sl = line.strip(' \r\n').split('=')
        if sl[0] == "PEPMASS":
            mass = sl[1]
        elif sl[0] == "TITLE":
            name = sl[1]
        elif sl[0] == "RTINSECONDS":
            r_time = sl[1]
        elif sl[0] == "END IONS":
            peaks.append(Peak(name, r_time, mass, charge, "False", ms2=ms2))
            ms2 = []
        else:
            try:
                mz, i = sl[0].split('\t')
                ms2.append((float(mz), float(i)))
            except ValueError:
                continue
    return peaks


def read_msp(input_string, charge):
    peaks = []
    for spec in input_string.strip().split('\n\n'):
        ms2 = []
        inchikey = "False"
        r_time = 0
        for line in spec.split('\n'):
            sl = line.split(': ')
            sl[0] = sl[0].replace(' ', '').replace('/', '').upper()
            if sl[0] == "PRECURSORMZ":
                mass = sl[1]
            elif sl[0] == "NAME":
                name = sl[1]
            #elif sl[0] == "RETENTIONTIME":
                #r_time = sl[1]
            elif sl[0] == "IONMODE":
                charge = sl[1].capitalize()
            elif sl[0] == "INCHIKEY":
                inchikey = sl[1]
            elif line and line[0].isdigit():
                try:
                    row = re.split('[\t ]', line)
                    ms2.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue
        peaks.append(Peak(name, r_time, mass, charge, inchikey, ms2=ms2))
    return peaks


def read_mzXML(input_string, charge):
    peaks = []
    tree = ET.fromstring(input_string)
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
            peaks.append(Peak(name, r_time, mz, charge, "False"))
    return peaks


def dot_product(x, y, epsilon=0.01):
    """Calculate the dot_product of two spectra"""
    z = 0
    n_v1 = 0
    n_v2 = 0

    for int1, int2 in Utils.approximate_matches(x, y, epsilon):
        z += int1 * int2
        n_v1 += int1 * int1
        n_v2 += int2 * int2
    return z / (math.sqrt(n_v1) * math.sqrt(n_v2))


def jacquard(x, y, epsilon=0.01):
    """Calculate the Jacquard Index of two spectra"""
    intersect = 0
    for val1, val2 in Utils.approximate_matches(x, y, epsilon):
        if val1 and val2:
            intersect += 1
    return intersect/float((len(x)+len(y)-intersect))


class Peak:
    """A class holding information about an unknown peak"""
    def __init__(self, name, r_time, mz, charge, inchi_key, ms2=[]):
        self.name = name
        self.r_time = float(r_time)  # retention time
        self.mz = float(mz)  # mass to charge ratio
        self.charge = charge  # polarity of charge
        self.inchi_key = inchi_key  # the id of the peak if known, as an Inchikey
        self.ms2peaks = ms2
        self.isomers = []
        self.formulas = set()
        self.total_hits = 0
        self.native_hit = False
        self.min_steps = 99

    def __str__(self):
        return self.name

    def score_isomers(self, metric=dot_product, energy_level=1):
        """
        Calculates the cosign similarity score between the provided ms2 peak list and pre-calculated CFM-spectra and
        sorts the isomer list according to this metric.
        :param metric: The scoring metric to use for the spectra. Function must accept 2 lists of (mz,intensity) tuples
         and return a score. Defaults to dot_product.
        :type energy_level: function
        :param energy_level: The Fragmentation energy level to use. Ranges from 0-2. Defaults to 1
        :type energy_level: int
        :return:
        :rtype:
        """
        if not self.ms2peaks:
            raise ValueError('The ms2 peak list is empty')

        for i, hit in enumerate(self.isomers):
            if "CFM_spectra" in hit:
                hit_spec = hit["CFM_spectra"]['Energy_%s' % energy_level]
                self.isomers[i]['Spectral_score'] = round(metric(self.ms2peaks, hit_spec)*1000)
                del hit['CFM_spectra']
            else:
                self.isomers[i]['Spectral_score'] = None
        self.isomers.sort(key=lambda x: x['Spectral_score'], reverse=True)

######################################################################################
#                                   Main code                                        #
######################################################################################
if __name__ == '__main__':
    test = read_mgf("/Users/JGJeffryes/Documents/Research/Scripts/Python/database-scripts/Batch Adduct Query/Cultured Cells Replicate 1.mgf")
    for i in range(1, len(test)):
        pass

    sys.exit()
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
        data.unk_peaks = read_csv(open(unknowns_file).read())
    elif '.mgf' in unknowns_file:
        data.unk_peaks = read_mgf(open(unknowns_file).read())
    elif '.mzXML' in unknowns_file:
        data.unk_peaks = read_mzXML(open(unknowns_file).read())
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
    with open(data.name+'.pkl') as outfile:
        cPickle.dump(data, outfile)

    tend = time.time()
    print "BatchAdductQuery.py completed in %s seconds" %(tend-tstart)
