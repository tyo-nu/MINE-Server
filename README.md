Server code for the MINE Databases

## MINE-Server

This repository contains the RESTful API client for the MINE databases, written in Flask. The API also serves as the backend for the MINE website.

The following databases are featured on the MINE website:
TODO: update database list here
* E. coli - “EcoCycexp2”
* Yeast - “YMDBexp2”
* KEGG - “KEGGexp2”

An example in Python accessing the API. TODO: Update examples

	>>> from Client import mineDatabaseServices
	>>> service = mineDatabaseServices("http://modelseed.org/services/mine-database")
	>>> service.quick_search("EcoCycexp2","ATP")
	[{u'SMILES': u'OC1C(COP(=O)(OP(=O)(OP(=O)(O)O)O)O)OC(C1O)n1cnc2c1ncnc2N', u'NP_likeness': 1.0719, u'Generation': 0.0, u'MINE_id': 12815, u'Inchikey': u'ZKHQWZAMYRWXGA-UHFFFAOYSA-N', u'Mass': 506.99574548699997, u'Names': [u'ATP', u"Adenosine 5'-triphosphate"], u'Formula': u'C10H16N5O13P3', u'_id': u'C25a9fafebc1b08a0ae0fec015803771c73485a61'}]

* To learn more about the MINE databases, visit http://minedatabase.mcs.anl.gov/. 
* For API documentation, please refer to the [API specification](http://jamesjeffryes.github.io/MINE-Server/)
