__author__ = 'JGJeffryes'

"""
This script reads a set of unknowns in a CSV, MGF, mzXML, or excel file and queries a mongo database for potential matches.
The matches are written back to a CSV, HTML doc, or printed to STDOUT
USAGE - BatchAdductQuery.py [options] run_name file_of_Unknowns
"""

import sys
from Utils import establish_db_client
from optparse import OptionParser
from HTMLPrinter import Printer
import xml.etree.ElementTree as ET
import numpy
import cPickle
import time
import string
import re


class Dataset():
    """A class containing all the information for a metabolomics data set"""
    def __init__(self, name, options):
        self.name = name
        self.options = options  # Object containing all options from the optparser
        dtype = numpy.dtype("S20, f8, f8")
        if hasattr(options, 'adduct_list'):
            all_pos_adducts = numpy.loadtxt('Positive Adducts full.txt', dtype=dtype)
            self.pos_adducts = numpy.array(filter(lambda x: x[0] in options.adduct_list, all_pos_adducts), dtype=dtype)
            all_neg_adducts = numpy.loadtxt('Negative Adducts full.txt', dtype=dtype)
            self.neg_adducts = numpy.array(filter(lambda x: x[0] in options.adduct_list, all_neg_adducts), dtype=dtype)
        if hasattr(options, 'positive_adduct_file'):
            self.pos_adducts = numpy.loadtxt(options.positive_adduct_file, dtype=dtype)
        if hasattr(options, 'negative_adduct_file'):
            self.neg_adducts = numpy.loadtxt(options.negative_adduct_file, dtype=dtype)
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
            precision = self.options.tolerance
        upper_bounds = potential_masses + precision
        lower_bounds = potential_masses - precision

        #search database for hits in the each adducts mass range that have no innate charge.
        hits = {}
        for i, adduct in enumerate(adducts):
            hits[adduct['f0']] = [x for x in db.compounds.find({"$and": [{"Mass": {"$gte": float(lower_bounds[i])}},
                                                            {"Mass": {"$lte": float(upper_bounds[i])}}, {'Charge': 0}]},
                                                            {"FP2": 0, "Reactant in": 0, "Product of": 0})]
            for compound in hits[adduct['f0']]:
                #Filters out halogens if the flag is enabled by moving to the next compound before the current compound
                # is counted or stored.
                if not self.options.halogens:
                    if re.search('F[^e]|Cl|Br', compound['Formula']):
                        continue

                #update the total hits for the peak
                peak.total_hits += 1

                #create a dictionary of formulas keyed by the adduct that produces them
                try:
                    if not compound['Formula'] in peak.formulas[adduct['f0']]:
                        peak.formulas[adduct['f0']].append(compound['Formula'])
                        peak.total_formulas += 1
                except KeyError:
                    peak.formulas[adduct['f0']] = [compound['Formula']]
                    peak.total_formulas += 1
                    #if the compound is not in the dictionary of compounds that are isomers of a formula, add it
                try:
                    if not compound in self.isomers[compound['Formula']]:
                        self.isomers[compound['Formula']].append(compound)
                except KeyError:
                    self.isomers[compound['Formula']] = [compound]

    def find_M_hits(self, peak, db):
        """This function finds compounds that have a natural positive charge (M+)"""
        if self.options.ppm:
            precision = (self.options.tolerance/100000.)*peak.mz
        else:
            precision = self.options.tolerance
        upper_bound = peak.mz + precision
        lower_bound = peak.mz - precision
        #search database for hits in the mass range that have an innate positive charge.
        hits = [x for x in db.compounds.find({"$and": [{"Mass": {"$gte": lower_bound}}, {"Mass": {"$lte": upper_bound}},
                                            {'Charge': 1}]}, {"FP2": 0, "Reactant in": 0, "Product of": 0})]
        for compound in hits:
            #Filters out halogens if the flag is enabled by moving to the next compound before the current compound
            # is counted or stored.
            if not self.options.halogens:
                if re.search('F[^e]|Cl|Br', compound['Formula']):
                    continue
            peak.total_hits += 1

            #create a dictionary of formulas keyed by the adduct that produces them
            try:
                if not compound['Formula'] in peak.formulas['M+']:
                    peak.formulas['M+'].append(compound['Formula'])
                    peak.total_formulas += 1
            except KeyError:
                peak.formulas['M+'] = [compound['Formula']]
                peak.total_formulas += 1

            #if the compound is not in the dictionary of compounds that are isomers of a formula, add it
            try:
                if not compound in self.isomers[compound['Formula']]:
                    self.isomers[compound['Formula']].append(compound)
            except KeyError:
                self.isomers[compound['Formula']] = [compound]

    def annotate_peaks(self, db):
        """ This function iterates the through the unknown peaks in the data set searches the database for compounds
            that match a peak m/z given the adducts permitted. Statistics on the annotated data set are printed"""
        for i, peak in enumerate(self.unk_peaks):
            if peak.charge == '+' or peak.charge == 'Positive' or peak.charge:
                self.find_db_hits(peak, db, self.pos_adducts)
                self.find_M_hits(peak, db)

            elif peak.charge == '-' or peak.charge == 'Negative' or not peak.charge:
                self.find_db_hits(peak, db, self.neg_adducts)

            else:
                raise ValueError("Invalid compound charge specification. Please use \"+\" or \"Positive\" for positive"
                                 " ions and \"-\" or \"Negative\"for negative ions")

            if peak.total_hits > 0:
                self.matched_peaks += 1
                self.total_hits += peak.total_hits
                self.total_formulas += len(peak.formulas)
            print "%s percent of peaks processed" % int(float(i)/float(len(self.unk_peaks)) * 100)

        print "Proposed matches for %s of %s peaks" % (self.matched_peaks, len(self.unk_peaks))
        try:
            print "Average hits per peak: %s" % (float(self.total_hits)/float(self.matched_peaks))
            print "Average formulas per peak: %s" % (float(self.total_formulas)/float(self.matched_peaks))
        except ZeroDivisionError:
            pass

    def cluster_peaks(self):
        """ This function looks for peaks at about the same retention time that appear to be related. This could either
            be because the peaks are different adducts of the same compound or potentally that they are isomers of the
            same formula. """
        tstart = time.time()
        r_times = set()
        tolerance = options.cluster/2  # the +/- value is half the total span of the cluster
        all_peaks = self.unk_peaks + self.known_peaks  # try to cluster known peaks with the unknown if possible
        #sort the peaks and iterate through the list. Time we look backward and forward at peaks and compare the
        # difference in retention times. if its in the tolerance add it to the group, if not we can stop looking in that
        # direction because the difference can only get bigger.
        all_peaks.sort(key=lambda x: x.r_time)
        for i, peak in enumerate(all_peaks):
            group = [peak]
            step = 1
            #look backward in the list
            while (i - step) >= 0:  # this makes sure we don't step off the index.
                if peak.r_time - all_peaks[i - step].r_time <= tolerance:
                    group.append(all_peaks[i - step])
                else:
                    break  # stop looking this direction
                step += 1
            #look forward in the list
            step = 1
            while (i + step) < len(all_peaks):
                spread = all_peaks[i + step].r_time - peak.r_time
                if spread <= tolerance:
                    group.append(all_peaks[i + step])
                else:
                    break
                step += 1
            #the above search is exhaustive but is prone to producing duplicate groups. Therefore we sort the group and
            # convert the list to an immutable type which can be hashed and rapidly checked against previously generated
            # groups
            if len(group) > 1:
                group.sort(key=lambda x: x.mz)
                group = tuple(group)
                if not group in r_times:
                    r_times.add(group)
        #regroup the peaks in each group in a dictionary with a formula proposed for that peak as the key.
        for group in r_times:
            cluster = {}
            for peak in group:
                for adduct in peak.formulas:
                    for formula in peak.formulas[adduct]:
                        try:
                            cluster[formula].append((peak.name, adduct))
                        except KeyError:
                            cluster[formula] = [(peak.name, adduct)]
            #can't really have a cluster of one peak can we? so we only store the formulas with more than one peak
            for formula in cluster:
                if len(cluster[formula]) > 1:
                    self.clusters.append((formula, cluster[formula]))
        tend = time.time()
        print "Clustering peaks took %s sec" % (tend-tstart)

    def print_ranked_isomers(self, formula):
        #print each isomer hit to stdout
        known_isomers = []
        native_isomers = []
        other_isomers = []
        for isomer in self.isomers[formula]:
            if isomer['InChI Key'].split('-')[0] in self.known_set:
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

    def write_ranked_isomers(self, formula):
        #write each isomer hit to csv file
        known_isomers = []
        native_isomers = []
        other_isomers = []
        for isomer in self.isomers[formula]:
            if isomer['InChI Key'].split('-')[0] in self.known_set:
                known_isomers.append(isomer)
            elif isomer['_id'] in self.native_set:
                native_isomers.append(isomer)
            else:
                other_isomers.append(isomer)
        if len(known_isomers) > 0:
            outfile.write("Matches with known/target compounds: %s\n" % len(known_isomers))
            for isomer in sort_NPLike(known_isomers):
                outfile.write(str(isomer)+'\n')
            outfile.writelines(known_isomers)
        if len(native_isomers) > 0:
            outfile.write("Matches with native compounds %s\n" % len(native_isomers))
            for isomer in sort_NPLike(native_isomers):
                outfile.write(str(isomer)+'\n')
        if len(other_isomers) > 0:
            outfile.write("Matches with promiscuity products %s\n" % len(other_isomers))
            for isomer in sort_NPLike(other_isomers):
                outfile.write(str(isomer)+'\n')

    def html_ranked_isomers(self, formula):
        #write each isomer to a html file
        labels = ['KEGG Code', 'Names', 'InChI Key', 'NP Likeness']
        known_isomers = []
        native_isomers = []
        other_isomers = []
        for isomer in self.isomers[formula]:
            if isomer['InChI Key'].split('-')[0] in self.known_set:
                known_isomers.append(isomer)
            elif isomer['_id'] in self.native_set:
                native_isomers.append(isomer)
            else:
                other_isomers.append(isomer)
        if len(known_isomers) > 0:
            outfile.write('<tr> <td colspan="2"> Matches with target/known compounds: %s</td></tr>\n' % len(known_isomers))
            printer.print_compound_html(sort_NPLike(known_isomers), outfile, labels=labels)
        if len(native_isomers) > 0:
            outfile.write('<tr> <td colspan="2"> Matches with native compounds: %s</td></tr>\n' % len(native_isomers))
            printer.print_compound_html(sort_NPLike(native_isomers), outfile, labels=labels)
        if len(other_isomers) > 0:
            outfile.write('<tr> <td colspan="2"> Matches with promiscuity products: %s</td></tr>\n' % len(other_isomers))
            printer.print_compound_html(sort_NPLike(other_isomers), outfile, labels=labels)


def get_modelSEED_comps(kb_db, db, models):
    seed_ids, _ids = set(), set()
    for x in models:
        for z in kb_db.models.find_one({'_id': x}, {'Compounds': 1})['Compounds']:
            seed_ids.add(z)
    for y in seed_ids:
        y = y.split('.')[1].zfill(5)
        for w in db.compounds.find({'DB_links.Model_SEED': y}, {'_id': 1}):
            _ids.add(w['_id'])
    return _ids


def sort_NPLike(dic_list):
    return sorted(dic_list, key=lambda x: float(x['NP Likeness']), reverse=True)


def validate():
    #this verifies that certain known peaks are returned. Should be rewritten as an external test
    found = 0
    for kpeak in data.known_peaks:
        right_key = kpeak.known
        for upeak in data.unk_peaks:
            if upeak.name == kpeak.name:
                for adduct in upeak.formulas:
                    for formula in upeak.formulas[adduct]:
                        for compound in data.isomers[formula]:
                            if compound["InChI Key"] == right_key:
                                print upeak.name
                                found += 1
    print found / float(len(data.known_peaks)) * 100


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
    def __init__(self, name, r_time, mz, charge, formula, id):
        self.name = name
        self.r_time = float(r_time)  # retention time
        self.mz = float(mz)  # mass to charge ratio
        self.charge = charge  # polarity of charge
        self.inchi_key = id  # the id of the peak if known, as an InChI Key
        self.formulas = formula  # a dictionary of all the formulas that match the peak for a given adduct as a key
        self.total_formulas = 0
        self.total_hits = 0

    def __str__(self):
        return self.name

    def print_formulas(self):
        #this sends the results to stdout in a "hopefully" human readable form"
        print "%s: %s formulas and %s structures" % (self.name, self.total_formulas, self.total_hits)
        print ''
        for adduct in self.formulas:
            num_form = len(self.formulas[adduct])
            if num_form > 0:
                print string.center(adduct, 25)
                print string.ljust("Formula", 18), string.center("Isomers", 7)
                for formula in self.formulas[adduct]:
                    print string.ljust(formula, 18), string.center(repr(len(data.isomers[formula])), 7)
                print ''

    def write_formulas_csv(self, outfile):
        #writes a comma separated value file that can be opened in excel
        outfile.write("%s,%s,%s\n" % (self.name, len(self.formulas), self.total_hits))
        outfile.write("\n")
        for adduct in self.formulas:
            num_form = len(self.formulas[adduct])
            if num_form > 0:
                outfile.write("%s\n" % adduct)
                outfile.write("Formula,Isomers\n")
                for formula in self.formulas[adduct]:
                    outfile.write("%s,%s\n" % (formula, len(data.isomers[formula])))
                    data.write_ranked_isomers(formula)
                outfile.write("\n")

    def write_formulas_html(self, outfile):
        #Writes the formulas as html documents
        outfile.write('<h2> %s</h2>' % self.name)
        outfile.write('<h3> Formulas: %s    Isomers: %s </h3>' % (len(self.formulas), self.total_hits))
        for adduct in self.formulas:
            num_form = len(self.formulas[adduct])
            if num_form > 0:
                outfile.write('<h3> %s</h3>\n' % adduct)
                outfile.write('<table border="1" cellspacing="0" cellpadding="6" style="font-size:14pt; text-align:center;">\n')
                outfile.write('\t<tr>\n\t\t<th>Formula</th>\n\t\t<th> Isomers</th>\n\t</tr>\n')
                for formula in self.formulas[adduct]:
                    outfile.write('\t<tr>\n\t\t<td><b> %s </b></td>\n\t\t<td> %s </td>\n\t</tr>\n'
                                  % (formula, len(data.isomers[formula])))
                    data.html_ranked_isomers(formula)
            outfile.write('</table>\n')
        outfile.write('<br>\n')


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
    parser.add_option("-o", "--outfile", dest="output_file", default="null",
                      help="The file name of the desired output file. Must be .csv or .pkl")
    parser.add_option("-d", "--database", dest="database", default="1GenKEGG",
                      help="The name of the database to search")
    parser.add_option("-m", "--modelSEED", dest="modelSEED", default="null",
                      help="The model SEED id of the organism for the data set")
    parser.add_option("-k", "--known", dest="known_file", default="null",
                      help="The path to a file containing known or targeted peaks")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float", default=0.002,
                      help="The m/z tolerance(precision) for peak searching")
    parser.add_option("-c", "--cluster", dest="cluster", type="float", default=-1, help="Sort peaks by retention time "
                      "and search for common formulas. The maximum spread for a cluster must be input as a float")
    parser.add_option("--ppm", dest="ppm", action="store_true", default=False,
                      help="Tolerance is in Parts Per Million")
    parser.add_option("-x", "--halogens", dest="halogens", action="store_true", default=False,
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
    elif ('.xls' or '.xlsx') in unknowns_file:
        data.unk_peaks = read_excel(unknowns_file)
    else:
        sys.exit("Unknown file type not recognised. Please use .mgf, .xlsx, .txt or .csv file")

    client = establish_db_client()
    db = client[options.database]
    if not db.compounds.count():
        sys.exit('No compounds in supplied database')

    if options.modelSEED != 'null':
        print "Loading native compounds"
        kbase_db = client['KBase']
        data.native_set = get_modelSEED_comps(kbase_db, db, [options.modelSEED])

    data.annotate_peaks(db)

    #if an outfile is specified write to it(faster) otherwise print to stdout

    num_files = 0
    if options.output_file != 'null':
        outfile = open(options.output_file, 'w')
    for i, peak in enumerate(data.unk_peaks):
        if options.output_file != 'null':
            if '.csv' in options.output_file:
                peak.write_formulas_csv(outfile)
            elif '.pkl' in options.output_file:
                cPickle.dump(peak, outfile, 2)
            elif '.html' in options.output_file:
                if i % 20 == 0:
                    outfile.close()
                    num_files += 1
                    outfile = open(str(num_files)+options.output_file, 'w')
                printer = Printer(db)
                peak.write_formulas_html(outfile)
            else:
                raise ValueError('Output file type not recognised. Please use .html, .csv, or .pkl extensions')
        else:
            peak.print_formulas()

    #if the cluster flag is true, do the clustering calculations
    if options.cluster >= 0:

        data.cluster_peaks()

        #print formulas and isomers for clusters
        for cluster in data.clusters:
            print cluster
            try:
                data.print_ranked_isomers(cluster[0])
            except:
                print "Could not retrieve matching isomers"
            print ''
            print ''

    tend = time.time()
    print "BatchAdductQuery.py completed in %s seconds" %(tend-tstart)
