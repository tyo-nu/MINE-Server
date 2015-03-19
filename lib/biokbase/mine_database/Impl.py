#BEGIN_HEADER
import pybel
import time
import Utils
import BatchAdductQuery
from PathwaySearch import PathwaySearch
from ast import literal_eval


class Pathway_query_params():
    def __init__(self, db, start, end, length, all_path):
        self.db = db
        self.start_comp = start
        self.end_comp = end
        self.len_limit = length
        self.all_paths = all_path
        self.np_min = -3
        self.gibbs_cap = 100
        self.verbose = False


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

search_projection = {'SMILES': 1, 'Formula': 1, 'MINE_id': 1, 'Names': 1, 'Inchikey': 1, 'Mass': 1}
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
        db = self.db_client[db]
        quick_search_results = Utils.quick_search(db, query, search_projection)
        for x in quick_search_results:
            if not isinstance(x['_id'], unicode):
                x['_id'] = unicode(x['_id'])
        #END quick_search

        #At some point might do deeper type checking...
        if not isinstance(quick_search_results, list):
            raise ValueError('Method quick_search return value ' +
                             'quick_search_results is not type list as required.')
        # return the results
        return [quick_search_results]

    def similarity_search(self, db, comp_structure, min_tc, fp_type, limit):
        # self.ctx is set by the wsgi application class
        # return variables are: similarity_search_results
        #BEGIN similarity_search
        similarity_search_results = []
        fp_type = str(fp_type)
        db = self.db_client[db]
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

        #END similarity_search

        #At some point might do deeper type checking...
        if not isinstance(similarity_search_results, list):
            raise ValueError('Method similarity_search return value ' +
                             'similarity_search_results is not type list as required.')
        # return the results
        return [similarity_search_results]

    def structure_search(self, db, input_format, comp_structure):
        # self.ctx is set by the wsgi application class
        # return variables are: structure_search_results
        #BEGIN structure_search
        db = self.db_client[db]
        mol = pybel.readstring(str(input_format), str(comp_structure))
        inchi_key = mol.write("inchikey").strip()
        # sure, we could look for a matching SMILES but this is faster
        structure_search_results = Utils.quick_search(db, inchi_key, search_projection)
        #END structure_search

        #At some point might do deeper type checking...
        if not isinstance(structure_search_results, list):
            raise ValueError('Method structure_search return value ' +
                             'structure_search_results is not type list as required.')
        # return the results
        return [structure_search_results]

    def substructure_search(self, db, substructure, limit):
        # self.ctx is set by the wsgi application class
        # return variables are: substructure_search_results
        #BEGIN substructure_search
        substructure_search_results = []
        db = self.db_client[db]
        if "\n" in substructure:
            query_mol = pybel.readstring('mol', str(substructure))
        else:
            query_mol = pybel.readstring('smi', str(substructure))
        query_fp = query_mol.calcfp("FP4").bits
        smarts = pybel.Smarts(query_mol.write('smi').strip())
        for x in db.compounds.find({"FP4": {"$all": query_fp}}, search_projection):
            if smarts.findall(pybel.readstring("smi", str(x["SMILES"]))):
                substructure_search_results.append(x)
                if len(substructure_search_results) == limit:
                    break
        #END substructure_search

        #At some point might do deeper type checking...
        if not isinstance(substructure_search_results, list):
            raise ValueError('Method substructure_search return value ' +
                             'substructure_search_results is not type list as required.')
        # return the results
        return [substructure_search_results]

    def database_query(self, db, mongo_query):
        # self.ctx is set by the wsgi application class
        # return variables are: database_query_results
        #BEGIN database_query
        if db != 'admin':  # we don't want users poking around here
            db = self.db_client[db]
            query_dict = literal_eval(mongo_query)  # this transforms the string into a dictionary
            database_query_results = [x for x in db.compounds.find(query_dict, search_projection)]
        else:
            database_query_results = ['Illegal query']
        #END database_query

        #At some point might do deeper type checking...
        if not isinstance(database_query_results, list):
            raise ValueError('Method database_query return value ' +
                             'database_query_results is not type list as required.')
        # return the results
        return [database_query_results]

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
            if ('Reactant_in' in meh) and (len(meh['Reactant_in']) > 1000):
                meh['Reactant_in'] = meh['Reactant_in'][:1000]
            if ('Product_of' in meh) and len(meh['Product_of']) > 1000:
                meh['Product_of'] = meh['Product_of'][:1000]

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
        for x in operator_names:
            op = db.operators.find_one({'_id': x})
            #if op:
                #op['Reaction_ids'] = [x['_id'] for x in db.reactions.find({'Operators': x})]
            objects.append(op)
        #END get_ops

        #At some point might do deeper type checking...
        if not isinstance(objects, list):
            raise ValueError('Method get_ops return value ' +
                             'objects is not type list as required.')
        # return the results
        return [objects]

    def get_models(self):
        # self.ctx is set by the wsgi application class
        # return variables are: models
        #BEGIN get_models
        models = self.models
        #END get_models

        #At some point might do deeper type checking...
        if not isinstance(models, list):
            raise ValueError('Method get_models return value ' +
                             'models is not type list as required.')
        # return the results
        return [models]

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
        name = text_type+time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())
        if isinstance(ms_params, dict):
            ms_params = Struct(**ms_params)
        db = self.db_client[ms_params.db]
        ms_params.verbose = False
        dataset = BatchAdductQuery.Dataset(name, ms_params)
        ms_adduct_output = []
        if text_type == 'form':
            for mz in text.split('\n'):
                dataset.unk_peaks.append(BatchAdductQuery.Peak(mz, 0, float(mz), ms_params.charge, {}, "False"))
        else:
            raise IOError('%s files not supported' % text_type)
        dataset.native_set = BatchAdductQuery.get_KEGG_comps(db, self.keggdb, ms_params.models)
        dataset.annotate_peaks(db)
        for peak in dataset.unk_peaks:
            for adduct in peak.formulas:
                for formula in peak.formulas[adduct]:
                    for hit in dataset.isomers[formula]:
                        hit['peak_name'] = peak.name
                        hit['adduct'] = adduct
                        ms_adduct_output.append(hit)
        #ms_adduct_output.sort(key=lambda x: (x['steps_from_source'], -x['NP_likeness']))
        #END ms_adduct_search

        #At some point might do deeper type checking...
        if not isinstance(ms_adduct_output, list):
            raise ValueError('Method ms_adduct_search return value ' +
                             'ms_adduct_output is not type list as required.')
        # return the results
        return [ms_adduct_output]

    def mz_search(self, text, text_type, mz_params):
        # self.ctx is set by the wsgi application class
        # return variables are: batch_output
        #BEGIN mz_search
        name = text_type+time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())
        if isinstance(mz_params, dict):
            mz_params = Struct(**mz_params)
        db = self.db_client[mz_params.db]
        mz_params.verbose = False
        dataset = BatchAdductQuery.Dataset(name, mz_params)
        batch_output = []
        if text_type == 'form':
            for mz in text.split('\n'):
                dataset.unk_peaks.append(BatchAdductQuery.Peak(mz, 0, float(mz), mz_params.charge, {}, "False"))
        else:
            raise IOError('%s files not supported' % text_type)
        dataset.native_set = BatchAdductQuery.get_KEGG_comps(db, self.keggdb, mz_params.models)
        dataset.annotate_peaks(db)
        for peak in sorted(dataset.unk_peaks, key=lambda x: x.min_steps):
            peak.adducts = []
            for adduct in peak.formulas:
                for formula in peak.formulas[adduct]:
                    peak.adducts.append({'adduct': adduct, 'formula': formula,
                                         'isomers': sorted(dataset.isomers[formula], key=lambda x: x['steps_from_source'])})
            del peak.formulas, peak.inchi_key
            batch_output.append(peak.__dict__)
        #END mz_search

        #At some point might do deeper type checking...
        if not isinstance(batch_output, list):
            raise ValueError('Method mz_search return value ' +
                             'batch_output is not type list as required.')
        # return the results
        return [batch_output]

    def pathway_search(self, db, start_comp, end_comp, len_limit, all_paths):
        # self.ctx is set by the wsgi application class
        # return variables are: pathway_query_results
        #BEGIN pathway_search
        params = Pathway_query_params(db, start_comp, end_comp, len_limit, all_paths)
        pathsearch = PathwaySearch(params)
        if params.all_paths:
            pathway_query_results = pathsearch.dfs()
        else:
            pathway_query_results = pathsearch.bfs()
        #END pathway_search

        #At some point might do deeper type checking...
        if not isinstance(pathway_query_results, list):
            raise ValueError('Method pathway_search return value ' +
                             'pathway_query_results is not type list as required.')
        # return the results
        return [pathway_query_results]
