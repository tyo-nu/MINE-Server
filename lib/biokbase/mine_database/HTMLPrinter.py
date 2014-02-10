__author__ = 'JGJeffryes'
"""
This is a library of HTML printing functions.
"""
import pymongo
import sys
import re


class Printer():
    def __init__(self, database, image_type = 'svg', image_source='http://lincolnpark.chem-eng.northwestern.edu/Smiles_dump/N@ME'):
        self.cimg_html = '\t\t<td><img align="center" height="200" src=%s.%s></td>\n' % (image_source, image_type)
        self.cstoic_html = '\t\t<td><p style="font-size:26pt;"> NUM</td>\n'
        self.clbl_html = '\t\t<td></td><td><p style="font-size:16pt; text-align:center;">L@BEL'
        self.plus_html = '\t\t<td><p style="font-size:26pt;"> + </td>\n'
        self.arrow_html = '\t\t<td><p style="font-size:26pt;"> &lt;==&gt;</td>\n'

        if isinstance(database, pymongo.database.Database):
            self.db = database
        else:
            try:
                client = pymongo.MongoClient()
            except:
                sys.exit("Failed to load database client. Please verify that mongod is running")
            self.db = client[database]

    def _print_comp_image(self, term, outfile):
        if term[0] != "1":
                outfile.write(re.sub('NUM', str(term[0]), self.cstoic_html))
        else:
            outfile.write('<td></td>')

        #special case to deal with H+
        if term[1] == 'X71306b6c4efe11bc7c485fbc71932f3deb14fa2c':
            outfile.write(re.sub('N@ME\.svg', "H+.png", self.cimg_html))
        #normal compound annotation
        else:
            outfile.write(re.sub('N@ME', str(term[1]), self.cimg_html))

    def _print_rxn_images(self, rxn, outfile):
        outfile.write('\t<tr>\n')
        for term in rxn['Reactants'][:-1]:
            self._print_comp_image(term, outfile)
            outfile.write(self.plus_html)
        self._print_comp_image(rxn['Reactants'][-1], outfile)
        outfile.write(self.arrow_html)
        for term in rxn['Products'][:-1]:
            self._print_comp_image(term, outfile)
            outfile.write(self.plus_html)
        self._print_comp_image(rxn['Products'][-1], outfile)
        outfile.write('\t</tr>\n')

    def _print_comp_label(self, _id, key, outfile):
        comp = self.db.compounds.find_one({'_id': _id}, {key: 1})
        try:
            outfile.write(re.sub('L@BEL', str(comp[key]), self.clbl_html))
        except (KeyError, TypeError):
            outfile.write('\t\t<td></td><td>')

    def _print_rxn_labels(self, rxn, key, outfile):
        outfile.write('\t<tr>\n')
        for term in rxn['Reactants']:
            self._print_comp_label(term[1], key, outfile)
            outfile.write('</td><td></td>\n')
        for term in rxn['Products'][:-1]:
            self._print_comp_label(term[1], key, outfile)
            outfile.write('</td><td></td>\n')
        self._print_comp_label(rxn['Products'][-1][1], key, outfile)
        outfile.write('\t</tr>\n')

    def print_rxn_html(self, reaction_id, outfile, rxn_labels=None, comp_labels=None):
        """This function prints a whole reaction with help from above"""
        rxn = self.db.reactions.find_one({"_id": reaction_id})
        if rxn_labels:
            outfile.write('<h3>%s</h3>' % '\t'.join([key + ': ' + rxn[key] for key in rxn_labels]))
        outfile.write('<table border="0" cellspacing="0" cellpadding="0" >\n')
        self._print_rxn_images(rxn, outfile)
        if comp_labels:
            self._print_rxn_labels(rxn, [key for key in comp_labels], outfile)
        outfile.write('</table>\n')

    def _write_compound_info(self, compound, labels, outfile):
        #print the compound image
        outfile.write('\t<tr>\n')
        outfile.write(re.sub('N@ME', str(compound['_id']), self.cimg_html))
        for key in labels:
            try:
                if isinstance(compound[key], list):
                    outfile.write('\t\t<td> %s </td>\n' % str(", ".join(compound[key])))
                else:
                    outfile.write('\t\t<td> %s </td>\n' % compound[key])
            except KeyError:
                outfile.write('\t\t<td> Not Found </td>\n')
        outfile.write('\t</tr>\n')

    def print_compound_html(self, compounds, outfile, labels=None):
        outfile.write('<table border="1" cellspacing="0" cellpadding="6"'
                      ' style="font-size:14pt; text-align:center;">\n')
        #print the heading block for the compounds
        outfile.write('\t<tr><th>Image</th>')
        for key in labels:
            outfile.write('<th>%s</th>' % key)
        outfile.write('</tr>\n')
        if isinstance(compounds, dict):
            self._write_compound_info(compounds, labels, outfile)
        elif isinstance(compounds, list):
            for compound in compounds:
                comp = self.db.compounds.find_one({"_id": compound})
                if comp:
                    self._write_compound_info(comp, labels, outfile)
                outfile.write('<tr>Compound not found</tr>\n')
        else:
            comp = self.db.compounds.find_one({'_id': compounds})
            if comp:
                self._write_compound_info(comp, labels, outfile)
            else:
                outfile.write('<tr>Compound not found</tr>\n')
        outfile.write('</table>\n')

if __name__ == '__main__':
    _id = sys.argv[1]
    if len(sys.argv) > 2:
        db_name = sys.argv[-1]
    else:
        db_name = '1GenEcoCyc'
    foo = Printer(db_name)
    outfile = open(_id + '.html', 'w')
    if _id[0] == 'R':
        foo.print_rxn_html(_id, outfile)
    if _id[0] == 'C':
        foo.print_compound_html()
