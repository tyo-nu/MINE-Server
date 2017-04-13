#BEGIN_HEADER
import pybel
import time
import Utils
import BatchAdductQuery
from ast import literal_eval
from minedatabase import databases, queries
import ast

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

search_projection = {'SMILES': 1, 'Formula': 1, 'MINE_id': 1, 'Names': 1, 'Inchikey': 1, 'Mass': 1, 'Sources': 1,
                     'Generation': 1, 'NP_likeness': 1}
#END_HEADER


class mineDatabaseServices:
    '''
    Module Name:
    mineDatabaseServices

    Module Description:
    =head1 mineDatabaseServices

=head2 SYNOPSIS

The MINE database is fundamentally composed of two different types of documents, which are represented by the Compound
and Reaction objects. Users can use text-matching queries to access these records directly or perform two types of more
advanced queries: Mass Adduct queries and pathway queries. Mass Adduct queries return a list of compounds that might
match the m/z of an unknown compound. Pathway queries return either the shortest path or all paths between two compounds
 in the database.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.models = []
        self.db_client = Utils.establish_db_client()
        self.kbase_db = self.db_client['KBase']
        self.keggdb = self.db_client['KEGGdb']
        self.legacy_dbs = {"CDMINE", "KEGGexp2", "YMDBexp2", "EcoCycexp2"}
        for model in self.kbase_db.models.find({}, {'Name': 1}):
            self.models.append((model['_id'], model['Name']))
        with open('./lib/biokbase/mine_database/Positive Adducts full.txt') as infile:
            self.pos_adducts = [line.split(' \t')[0] for line in infile if not line[0] == '#']
        with open('./lib/biokbase/mine_database/Negative Adducts full.txt') as infile:
            self.neg_adducts = [line.split(' \t')[0] for line in infile if not line[0] == '#']
        self.timeout = 10
        #END_CONSTRUCTOR
        pass

    def model_search(self, query):
        # self.ctx is set by the wsgi application class
        # return variables are: models
        #BEGIN model_search
        print('<Model Search: %s>' % query)
        db = self.db_client['KEGGdb']
        cursor = db.models.find({"$text": {"$search": query}}, {"score": {"$meta": "textScore"}, "_id": 1})
        models = [x["_id"] for x in cursor.sort([("score", {"$meta": "textScore"})])]
        #END model_search

        #At some point might do deeper type checking...
        if not isinstance(models, list):
            raise ValueError('Method model_search return value ' +
                             'models is not type list as required.')
        # return the results
        return [models]

    def quick_search(self, db, query):
        # self.ctx is set by the wsgi application class
        # return variables are: quick_search_results
        #BEGIN quick_search
        print("<Quick Search: DB=%s, Query=%s>" % (db, query))
        db = self.db_client[db]
        if db.name not in self.legacy_dbs:
            quick_search_results = queries.quick_search(db, query, search_projection)
        else:
            quick_search_results = Utils.quick_search(db, query, search_projection)
        #END quick_search

        #At some point might do deeper type checking...
        if not isinstance(quick_search_results, list):
            raise ValueError('Method quick_search return value ' +
                             'quick_search_results is not type list as required.')
        # return the results
        return [quick_search_results]

    def similarity_search(self, db, comp_structure, min_tc, fp_type, limit, parent_filter, reaction_filter):
        # self.ctx is set by the wsgi application class
        # return variables are: similarity_search_results
        #BEGIN similarity_search
        print("<Similarity Search: DB=%s, Structure=%s, MinTC=%s, FPType=%s, Limit=%s>" % (db, comp_structure,
                                                                                           min_tc, fp_type, limit))
        db = self.db_client[db]
        if db.name not in self.legacy_dbs:
            return [queries.similarity_search(db, comp_structure, min_tc, fp_type, limit, search_projection)]
        similarity_search_results = []
        fp_type = str(fp_type)
        if "\n" in comp_structure:
            mol = pybel.readstring('mol', str(comp_structure))
        else:
            mol = pybel.readstring('smi', str(comp_structure))
        query_fp = set(mol.calcfp(fp_type).bits)
        len_fp = len(query_fp)
        for x in db.compounds.find({"$and": [{"len_"+fp_type: {"$gte": min_tc*len_fp}},
                                   {"len_"+fp_type: {"$lte": len_fp/min_tc}}]},
                                   dict([(fp_type, 1)]+search_projection.items())):
            test_fp = set(x[fp_type])
            tc = len(query_fp & test_fp)/float(len(query_fp | test_fp))
            if tc >= min_tc:
                del x[fp_type]
                similarity_search_results.append(x)
                if len(similarity_search_results) == limit:
                    break
        similarity_search_results = Utils.score_compounds(db, similarity_search_results, parent_filter, parent_frac=.75, reaction_frac=.25)
        #END similarity_search

        #At some point might do deeper type checking...
        if not isinstance(similarity_search_results, list):
            raise ValueError('Method similarity_search return value ' +
                             'similarity_search_results is not type list as required.')
        # return the results
        return [similarity_search_results]

    def structure_search(self, db, input_format, comp_structure, parent_filter, reaction_filter):
        # self.ctx is set by the wsgi application class
        # return variables are: structure_search_results
        #BEGIN structure_search
        print("<Structure Search: DB=%s, Structure=%s, Format=%s>" % (db, comp_structure, input_format))
        db = self.db_client[db]
        if db.name not in self.legacy_dbs:
            return [queries.structure_search(db, comp_structure, search_projection)]
        mol = pybel.readstring(str(input_format), str(comp_structure))
        inchi_key = mol.write("inchikey").strip()
        # sure, we could look for a matching SMILES but this is faster
        structure_search_results = Utils.quick_search(db, inchi_key, search_projection)
        structure_search_results = Utils.score_compounds(db, structure_search_results, parent_filter, parent_frac=.75, reaction_frac=.25)
        #END structure_search

        #At some point might do deeper type checking...
        if not isinstance(structure_search_results, list):
            raise ValueError('Method structure_search return value ' +
                             'structure_search_results is not type list as required.')
        # return the results
        return [structure_search_results]

    def substructure_search(self, db, substructure, limit, parent_filter, reaction_filter):
        # self.ctx is set by the wsgi application class
        # return variables are: substructure_search_results
        #BEGIN substructure_search
        print("<Substructure Search: DB=%s, Structure=%s, Limit=%s>" % (db, substructure, limit))
        db = self.db_client[db]
        if db.name not in self.legacy_dbs:
            return [queries.substructure_search(db, substructure, limit, search_projection)]
        substructure_search_results = []
        if "\n" in substructure:
            query_mol = pybel.readstring('mol', str(substructure))
        else:
            query_mol = pybel.readstring('smi', str(substructure))
        query_fp = query_mol.calcfp("FP4").bits
        smarts = pybel.Smarts(query_mol.write('smi').strip())
        for x in db.compounds.find({"FP4": {"$all": query_fp}}, search_projection):
            try:
                if smarts.findall(pybel.readstring("smi", str(x["SMILES"]))):
                    substructure_search_results.append(x)
                    if len(substructure_search_results) == limit:
                        break
            except IOError:  # Smarts searches fail for generalized compounds so we just skip over them
                continue
        substructure_search_results = Utils.score_compounds(db, substructure_search_results, parent_filter, parent_frac=.75, reaction_frac=.25)
        #END substructure_search

        #At some point might do deeper type checking...
        if not isinstance(substructure_search_results, list):
            raise ValueError('Method substructure_search return value ' +
                             'substructure_search_results is not type list as required.')
        # return the results
        return [substructure_search_results]

    def database_query(self, db, mongo_query, parent_filter, reaction_filter):
        # self.ctx is set by the wsgi application class
        # return variables are: database_query_results
        #BEGIN database_query
        print("<Database Search: DB=%s, Query=%s>" % (db, mongo_query))
        db = self.db_client[db]
        database_query_results = queries.advanced_search(db, mongo_query, search_projection)
        database_query_results = Utils.score_compounds(db, database_query_results, parent_filter, parent_frac=.75, reaction_frac=.25)
        #END database_query

        #At some point might do deeper type checking...
        if not isinstance(database_query_results, list):
            raise ValueError('Method database_query return value ' +
                             'database_query_results is not type list as required.')
        # return the results
        return [database_query_results]

    def get_ids(self, db, collection, query):
        # self.ctx is set by the wsgi application class
        # return variables are: ids
        #BEGIN get_ids
        if db == 'admin':
            raise ValueError('Illegal query')

        db = self.db_client[db]
        if query:
            query = literal_eval(query)  # this transforms the string into a dictionary
        else:
            query = {}
        ids = [x['_id'] for x in db[collection].find(query)]
        #END get_ids

        #At some point might do deeper type checking...
        if not isinstance(ids, list):
            raise ValueError('Method get_ids return value ' +
                             'ids is not type list as required.')
        # return the results
        return [ids]

    def get_comps(self, db, ids):
        # self.ctx is set by the wsgi application class
        # return variables are: objects
        #BEGIN get_comps
        objects = []
        db = self.db_client[db]
        if len(ids) > 50:
            raise ValueError("List of compound ids is too long. get_comps is limited to 50 or fewer compounds.")
        for x in ids:
            if isinstance(x, int):
                meh = db.compounds.find_one({'MINE_id': x}, {"len_FP2": 0, "FP2": 0, "len_FP4": 0, "FP4": 0})
            else:
                meh = db.compounds.find_one({'_id': x}, {"len_FP2": 0, "FP2": 0, "len_FP4": 0, "FP4": 0})
            # New MINEs won't have this precomputed
            if meh and 'Reactant_in' not in meh and 'Product_of' not in meh:
                meh['Reactant_in'] = [x['_id'] for x in db.reactions.find({'Reactants.c_id': meh['_id']})]
                meh['Product_of'] = [x['_id'] for x in db.reactions.find({'Products.c_id': meh['_id']})]
            objects.append(meh)
        #END get_comps

        #At some point might do deeper type checking...
        if not isinstance(objects, list):
            raise ValueError('Method get_comps return value ' +
                             'objects is not type list as required.')
        # return the results
        return [objects]

    def get_rxns(self, db, ids):
        # self.ctx is set by the wsgi application class
        # return variables are: objects
        #BEGIN get_rxns
        objects = []
        db = self.db_client[db]
        for x in ids:
            objects.append(db.reactions.find_one({'_id': x}))
        #END get_rxns

        #At some point might do deeper type checking...
        if not isinstance(objects, list):
            raise ValueError('Method get_rxns return value ' +
                             'objects is not type list as required.')
        # return the results
        return [objects]

    def get_ops(self, db, operator_names):
        # self.ctx is set by the wsgi application class
        # return variables are: objects
        #BEGIN get_ops
        objects = []
        db = self.db_client[db]
        if not operator_names:
            objects = [x for x in db.operators.find()]
        else:
            objects = [db.operators.find_one({'$or': [{'_id': x}, {"Name": x}]}) for x in operator_names]
        #END get_ops

        #At some point might do deeper type checking...
        if not isinstance(objects, list):
            raise ValueError('Method get_ops return value ' +
                             'objects is not type list as required.')
        # return the results
        return [objects]

    def get_operator(self, db, operator_name):
        # self.ctx is set by the wsgi application class
        # return variables are: operator
        #BEGIN get_operator
        db = self.db_client[db]
        operator = db.operators.find_one({'$or': [{'_id': operator_name}, {"Name": operator_name}]})
        operator['Reaction_ids'] = db.reactions.find({"Operators": operator_name}).distinct("_id")
        #END get_operator

        #At some point might do deeper type checking...
        if not isinstance(operator, dict):
            raise ValueError('Method get_operator return value ' +
                             'operator is not type dict as required.')
        # return the results
        return [operator]

    def get_adducts(self):
        # self.ctx is set by the wsgi application class
        # return variables are: adducts
        #BEGIN get_adducts
        adducts = [self.pos_adducts, self.neg_adducts]
        #END get_adducts

        #At some point might do deeper type checking...
        if not isinstance(adducts, list):
            raise ValueError('Method get_adducts return value ' +
                             'adducts is not type list as required.')
        # return the results
        return [adducts]

    def ms_adduct_search(self, text, text_type, ms_params):
        # self.ctx is set by the wsgi application class
        # return variables are: ms_adduct_output
        #BEGIN ms_adduct_search
        print("<MS Adduct Search: TextType=%s, Text=%s, Parameters=%s>" % (text_type, text, ms_params))
        name = text_type+time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())
        if isinstance(ms_params, dict):
            ms_params = Struct(**ms_params)
        db = self.db_client[ms_params.db]
        ms_params.verbose = False
        dataset = BatchAdductQuery.Dataset(name, ms_params)
        ms_adduct_output = []
        if text_type == 'form':
            for mz in text.split('\n'):
                dataset.unk_peaks.append(BatchAdductQuery.Peak(mz, 0, float(mz), ms_params.charge, "False"))
        elif text_type == 'mgf':
            dataset.unk_peaks = BatchAdductQuery.read_mgf(text, ms_params.charge)
        elif text_type == 'mzXML':
            dataset.unk_peaks = BatchAdductQuery.read_mzXML(text, ms_params.charge)
        elif text_type == 'msp':
            dataset.unk_peaks = BatchAdductQuery.read_msp(text, ms_params.charge)
        else:
            raise IOError('%s files not supported' % text_type)
        dataset.native_set = BatchAdductQuery.get_KEGG_comps(db, self.keggdb, ms_params.models)
        dataset.annotate_peaks(db)
        for peak in dataset.unk_peaks:
            for hit in peak.isomers:
                if 'CFM_spectra' in hit:
                    del hit['CFM_spectra']
                ms_adduct_output.append(hit)
        if ms_params.models:
            ms_adduct_output = Utils.score_compounds(db, ms_adduct_output, ms_params.models[0], parent_frac=.75, reaction_frac=.25)
        #END ms_adduct_search

        #At some point might do deeper type checking...
        if not isinstance(ms_adduct_output, list):
            raise ValueError('Method ms_adduct_search return value ' +
                             'ms_adduct_output is not type list as required.')
        # return the results
        return [ms_adduct_output]

    def ms2_search(self, text, text_type, ms_params):
        # self.ctx is set by the wsgi application class
        # return variables are: ms_adduct_output
        #BEGIN ms2_search
        print("<MS Adduct Search: TextType=%s, Parameters=%s>" % (text_type, ms_params))
        name = text_type+time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())
        if isinstance(ms_params, dict):
            ms_params = Struct(**ms_params)
        db = self.db_client[ms_params.db]
        ms_params.verbose = False
        dataset = BatchAdductQuery.Dataset(name, ms_params)
        ms_adduct_output = []
        if text_type == 'form':
            split_form = [x.split() for x in text.strip().split('\n')]
            dataset.unk_peaks.append(BatchAdductQuery.Peak(split_form[0][0], 0, float(split_form[0][0]), ms_params.charge,
                                                           "False", ms2=[(float(mz), float(i)) for mz, i in split_form[1:]]))
        elif text_type == 'mgf':
            dataset.unk_peaks = BatchAdductQuery.read_mgf(text, ms_params.charge)
        elif text_type == 'mzXML':
            dataset.unk_peaks = BatchAdductQuery.read_mzXML(text, ms_params.charge)
        elif text_type == 'msp':
            dataset.unk_peaks = BatchAdductQuery.read_msp(text, ms_params.charge)
        else:
            raise IOError('%s files not supported' % text_type)
        dataset.native_set = BatchAdductQuery.get_KEGG_comps(db, self.keggdb, ms_params.models)
        dataset.annotate_peaks(db)
        for peak in dataset.unk_peaks:
            if ms_params.scoring_function == 'jaccard':
                if not ms_params.ppm:
                    peak.score_isomers(metric=BatchAdductQuery.jaccard, energy_level=ms_params.energy_level,
                                       tolerance=float(ms_params.tolerance)/1000)
                else:
                    peak.score_isomers(metric=BatchAdductQuery.jaccard, energy_level=ms_params.energy_level)
            else:
                if not ms_params.ppm:
                    peak.score_isomers(metric=BatchAdductQuery.dot_product, energy_level=ms_params.energy_level,
                                       tolerance=float(ms_params.tolerance)/1000)
                else:
                    peak.score_isomers(metric=BatchAdductQuery.dot_product, energy_level=ms_params.energy_level)
            for hit in peak.isomers:
                ms_adduct_output.append(hit)
            if ms_params.models:
                ms_adduct_output = Utils.score_compounds(db, ms_adduct_output, ms_params.models[0], parent_frac=.75, reaction_frac=.25)
        #END ms2_search

        #At some point might do deeper type checking...
        if not isinstance(ms_adduct_output, list):
            raise ValueError('Method ms2_search return value ' +
                             'ms_adduct_output is not type list as required.')
        # return the results
        return [ms_adduct_output]

    def spectra_download(self, db, mongo_query, parent_filter, putative, spec_type):
        # self.ctx is set by the wsgi application class
        # return variables are: spectral_library
        #BEGIN spectra_download
        def print_peaklist(peaklist):
            text = ["Num Peaks: %s" % len(peaklist)]
            for x in peaklist:
                text.append("%s %s" % (x[0], x[1]))
            text.append("")
            return text

        spectral_library = []
        msp_projection = {'MINE_id': 1, 'Names': 1, 'Mass': 1, 'Generation': 1,
                          'Inchikey': 1, 'Formula': 1, 'SMILES': 1,
                          'Sources': 1, 'Pos_CFM_spectra': 1,
                          'Neg_CFM_spectra': 1}
        db = self.db_client[db]
        if mongo_query:
            query_dict = ast.literal_eval(mongo_query)
        else:
            query_dict = {}
        if not putative:
            query_dict['Generation'] = 0
        if parent_filter:
            model = db.models.find_one({"_id": parent_filter})
            if not model:
                raise ValueError('Invalid Model specified')
            parents = model["Compound_ids"]
            query_dict['$or'] = [{'_id': {'$in': parents}},
                                 {'Sources.Compound': {'$in': parents}}]
        results = db.compounds.find(query_dict, msp_projection)

        for compound in results:
            # make header
            header = []
            if "Names" in compound and len(compound['Names']):
                header.append("Name: %s" % compound['Names'][0])
                for alt in compound['Names'][1:]:
                    header.append("Synonym: %s" % alt)
            for k, v in compound.items():
                if k not in {"Names", "Pos_CFM_spectra", "Neg_CFM_spectra"}:
                    header.append("%s: %s" % (k, v))
            header.append("Instrument: CFM-ID")

            # add peak lists
            if 'Pos_CFM_spectra' in compound:
                for energy, spec in compound['Pos_CFM_spectra'].items():
                    if not spec_type or (True, int(energy[:2])) in spec_type:
                        spectral_library += header
                        spectral_library += ["Ionization Mode: Positive",
                                             "Energy: %s" % energy]
                        spectral_library += print_peaklist(spec)

            if 'Neg_CFM_spectra' in compound:
                for energy, spec in compound['Neg_CFM_spectra'].items():
                    if not spec_type or (False, int(energy[:2])) in spec_type:
                        spectral_library += header
                        spectral_library += ["Ionization Mode: Negative",
                                             "Energy: %s" % energy]
                        spectral_library += print_peaklist(spec)
        spectral_library = "\n".join(spectral_library)
        #END spectra_download

        #At some point might do deeper type checking...
        if not isinstance(spectral_library, basestring):
            raise ValueError('Method spectra_download return value ' +
                             'spectral_library is not type basestring as required.')
        # return the results
        return [spectral_library]

    def pathway_search(self, db, start_comp, end_comp, len_limit, all_paths):
        # self.ctx is set by the wsgi application class
        # return variables are: pathway_query_results
        #BEGIN pathway_search
        pathway_query_results = ['Not yet implemented']
        #END pathway_search

        #At some point might do deeper type checking...
        if not isinstance(pathway_query_results, list):
            raise ValueError('Method pathway_search return value ' +
                             'pathway_query_results is not type list as required.')
        # return the results
        return [pathway_query_results]

meh = mineDatabaseServices(None)
meh.spectra_download('EcoCycexp2', False, False, False,
                                    [(True, 20), (False, 40)])