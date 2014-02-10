#BEGIN_HEADER
from pymongo import MongoClient
import Utils
import BatchAdductQuery
from PathwaySearch import PathwaySearch
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
	self.db_client = MongoClient(host='branch')
        admin = self.db_client['admin']
        admin.authenticate('worker', 'bnice14bot')
        kbase_db = self.db_client['KBase']
        for model in kbase_db.models.find({}, {'Name': 1}):
            self.models.append((model['_id'], model['Name']))
        with open('/vol/model-prod/mine-server/lib/biokbase/mine_database/Positive Adducts full.txt') as infile:
            self.pos_adducts = [line.split('\t')[0] for line in infile if not line[0] == '#']
        with open('/vol/model-prod/mine-server/lib/biokbase/mine_database/Negative Adducts full.txt') as infile:
            self.neg_adducts = [line.split('\t')[0] for line in infile if not line[0] == '#']
        #END_CONSTRUCTOR
        pass

    def quick_search(self, db, query):
        # self.ctx is set by the wsgi application class
        # return variables are: quick_search_results
        #BEGIN quick_search
        db = self.db_client[db]
        quick_search_results = Utils.quick_search(db, query)
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
        #TODO implement similarity search
        similarity_search_results = ['Not Yet Implemented']
        #END similarity_search

        #At some point might do deeper type checking...
        if not isinstance(similarity_search_results, list):
            raise ValueError('Method similarity_search return value ' +
                             'similarity_search_results is not type list as required.')
        # return the results
        return [similarity_search_results]

    def database_query(self, params):
        # self.ctx is set by the wsgi application class
        # return variables are: database_query_results
        #BEGIN database_query
        if params.db != 'admin':
            db = self.db_client[params.db]
            if params.regex:
                database_query_results = [x['_id'] for x in db.compounds.find({params.field: {'$regex': params.value}},
                                            {'Mass': 1, 'Formula': 1, 'Inchi_key': 1, 'KEGG_code': 1, 'Names': 1})]
            else:
                database_query_results = [x['_id'] for x in db.compounds.find({params.field: params.value},
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
        adducts = (self.pos_adducts, self.neg_adducts)
        #END get_adducts

        #At some point might do deeper type checking...
        if not isinstance(adducts, list):
            raise ValueError('Method get_adducts return value ' +
                             'adducts is not type list as required.')
        # return the results
        return [adducts]

    def adduct_db_search(self, params):
        # self.ctx is set by the wsgi application class
        # return variables are: output
        #BEGIN adduct_db_search
        output = []
        db = self.db_client[params.db]
        name = "Single Peak: " + str(params.mz)
        dataset = BatchAdductQuery.Dataset(name, params)
        dataset.unk_peaks = [BatchAdductQuery.Peak(name, 0, params.mz, params.charge, {}, "False")]
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

    def pathway_search(self, pathway_query_params):
        # self.ctx is set by the wsgi application class
        # return variables are: pathway_query_results
        #BEGIN pathway_search
        setattr(pathway_query_params, 'verbose', False)
        pathsearch = PathwaySearch(pathway_query_params)
        if pathway_query_params.all_paths:
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
