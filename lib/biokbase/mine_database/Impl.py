#BEGIN_HEADER
import pybel
import Utils
import BatchAdductQuery
from PathwaySearch import PathwaySearch


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


class Adduct_search_params():
    def __init__(self, db, tolerance, adducts, charge, models, ppm, halogens):
        self.db = db
        self.tolerance = tolerance
        self.adduct_list = adducts
        self.models = models
        self.ppm = ppm
        self.charge = charge
        self.halogens = halogens
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
        kbase_db = self.db_client['KBase']
        for model in kbase_db.models.find({}, {'Name': 1}):
            self.models.append((model['_id'], model['Name']))
        with open('/vol/model-prod/mine-server/lib/biokbase/mine_database/Positive Adducts full.txt') as infile:
            self.pos_adducts = [line.split('\t')[0] for line in infile if not line[0] == '#']
        with open('/vol/model-prod/mine-server/lib/biokbase/mine_database/Negative Adducts full.txt') as infile:
            self.neg_adducts = [line.split('\t')[0] for line in infile if not line[0] == '#']
        #END_CONSTRUCTOR

    def quick_search(self, db, query):
        # self.ctx is set by the wsgi application class
        # return variables are: quick_search_results
        #BEGIN quick_search
        mdb = self.db_client[db]
        quick_search_results = Utils.quick_search(mdb, query)
        #END quick_search

        #At some point might do deeper type checking...
        if not isinstance(quick_search_results, list):
            raise ValueError('Method quick_search return value ' +
                             'quick_search_results is not type list as required.')
        # return the results
        return [quick_search_results]

    def similarity_search(self, db, smiles, min_tc):
        # self.ctx is set by the wsgi application class
        # return variables are: similarity_search_results
        #BEGIN similarity_search
        similarity_search_results = []
        db = self.db_client[db]
        mol = pybel.readstring('smi', smiles)
        query_fp = set(mol.calcfp().bits)
        len_fp = len(query_fp)
        comps = [x for x in db.compounds.find({"$and": [{"len_FP2": {"$gte": min_tc*len_fp}},
                    {"len_FP2": {"$lte": len_fp/min_tc}}]},
                    {"FP2": 1, 'Mass': 1, 'Formula': 1, 'Inchi_key': 1, 'KEGG_code': 1, 'Names': 1})]
        for x in comps:
            test_fp = set(x['FP2'])
            tc = len(query_fp & test_fp)/float(len(query_fp | test_fp))
            if tc >= min_tc:
                del x['FP2']
                similarity_search_results.append(x)

        #END similarity_search

        #At some point might do deeper type checking...
        if not isinstance(similarity_search_results, list):
            raise ValueError('Method similarity_search return value ' +
                             'similarity_search_results is not type list as required.')
        # return the results
        return [similarity_search_results]

    def database_query(self, db, field, value, regex):
        # self.ctx is set by the wsgi application class
        # return variables are: database_query_results
        #BEGIN database_query
        if db != 'admin':
            db = self.db_client[db]
            if regex:
                database_query_results = [x for x in db.compounds.find({field: {'$regex': value}},
                                                 {'Mass': 1, 'Formula': 1, 'Inchi_key': 1, 'KEGG_code': 1, 'Names': 1})]
            else:
                database_query_results = [x for x in db.compounds.find({field: value},
                                                 {'Mass': 1, 'Formula': 1, 'Inchi_key': 1, 'KEGG_code': 1, 'Names': 1})]
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
        for x in ids:
            objects.append(db.compounds.find({'_id': x}, {'_id': 1}))
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
            objects.append(db.reactions.find({'_id': x}, {'_id': 1}))
        #END get_rxns

        #At some point might do deeper type checking...
        if not isinstance(objects, list):
            raise ValueError('Method get_rxns return value ' +
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

    def adduct_db_search(self, db, mz, tolerance, adduct_list, models, ppm, charge, halogens):
        # self.ctx is set by the wsgi application class
        # return variables are: output
        #BEGIN adduct_db_search
        output = []
        db = self.db_client[db]
        name = "Single Peak: " + str(mz)
        params = Adduct_search_params(db, tolerance, adduct_list, charge, models, ppm, halogens)
        dataset = BatchAdductQuery.Dataset(name, params)
        dataset.unk_peaks = [BatchAdductQuery.Peak(name, 0, mz, charge, {}, "False")]
        dataset.annotate_peaks(db)
        result = dataset.unk_peaks[0]
        for adduct in result.formulas:
            for formula in result.formulas[adduct]:
                output.append((adduct, formula, [x['_id'] for x in dataset.isomers[formula]]))
        #END adduct_db_search

        #At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method adduct_db_search return value ' +
                             'output is not type list as required.')
        # return the results
        return [output]

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
