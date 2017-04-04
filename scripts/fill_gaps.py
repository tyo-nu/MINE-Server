from lib.biokbase.mine_database.Client import mineDatabaseServices
services = mineDatabaseServices('http://bio-data-1.mcs.anl.gov/services/mine-database')
db = 'EcoCycexp2'
print("Unknown sources")
with open('no_source.tsv') as infile:
    for line in infile:
        sl = line.split('\t')
        if len(sl) == 5:
            inchi = sl[4].split('=')[1].split('-')[0]
            res = services.database_query(db, "{'Inchikey': {'$regex':'" + inchi +
                                          "'}, 'Product_of': {'$exists' : 1}}", "", "")
            if res:
                n_source = len(services.get_comps(db, [res[0]['_id']])[0]["Product_of"])
                #if "Sources" in res[0]:
                    #print(sl[0], res[0]['Sources'][0]['Operators'])
                #else:
                print('"%s";"%s"' % (sl[0], n_source))
            else:
                print('"%s";"0"' % sl[0])
print("\n\n")
print("Unknown sink")
with open('no_sink.tsv') as infile:
    for line in infile:
        sl = line.split('\t')
        if len(sl) == 5:
            n_sinks = 0
            inchi = sl[4].split('=')[1].split('-')[0]
            res = services.database_query(db, "{'Inchikey': {'$regex':'" + inchi +
                                          "'}, 'Reactant_in': {'$exists' : 1}}", "", "")
            if res:
                for rxn in services.get_comps(db, [res[0]['_id']])[0]["Reactant_in"]:
                    meh = services.database_query(db, "{'Generation': 0, 'Product_of':'" + rxn +"'}", "", "")
                    if meh:
                        n_sinks += 1
                        for comp in meh:
                            print(comp['Names'])
            print('"%s";"%s"' % (sl[0], n_sinks))