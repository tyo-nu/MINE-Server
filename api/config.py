"""Configuration file. All settings are stored in Config class which gets
instantiated in __init__.py."""

import os

import minedatabase

try:
    from api.credentials import (MONGO_PASSWORD, MONGO_USERNAME, POSTGRES_PASSWORD,
                                 POSTGRES_USERNAME)
except ImportError:
    raise FileNotFoundError("MINE-Server/api/credentials.py not found.")

APP_DIR = os.path.abspath(os.path.dirname(__file__))
MINEDB_DIR = os.path.dirname(minedatabase.__file__)


class Config(object):
    """Class to configure important filepaths and settings."""
    # ------------------------------ Switches ------------------------------- #
    #: Whether to use thermo database - useful to switch off for development
    THERMO_ON = False

    # ------------------------------ Filepaths ------------------------------ #
    # Local filepaths are defined here

    #: Path to "Positive Adducts full.txt"
    POS_ADDUCT_PATH = os.path.join(MINEDB_DIR,
                                   'data/adducts/Positive Adducts full.txt')

    #: Path to "Negative Adducts full.txt"
    NEG_ADDUCT_PATH = os.path.join(MINEDB_DIR,
                                   'data/adducts/Negative Adducts full.txt')

    #: Path to directory with test data
    TEST_DATA_DIR = os.path.join(APP_DIR, '../tests/data')

    #: Path to operator images
    OP_IMG_DIR = os.path.join(APP_DIR, '../static/operator_images')

    # ------------------------------ Databases ------------------------------ #
    # Settings for interface with MongoDB and PostgreSQL

    #: URI to MINE MongoDB  # TODO: change this back to port 27017
    #MONGO_URI = f'mongodb://{MONGO_USERNAME}:{MONGO_PASSWORD}@localhost:27017/'
    MONGO_URI = f'mongodb://jrs9291:TurtleUkulele21!@minedatabase.ci.northwestern.edu:27017'

    #: URI to Thermo PostgreSQL DB
    POSTGRES_URI = f'postgres://{POSTGRES_USERNAME}:{POSTGRES_PASSWORD}@minedatabase.ci.northwestern.edu:5432/eq_compounds'

    #: Name of core compound database with spectra and fingerprint info
    CORE_DB_NAME = 'core'

    #: Name of compound external references database (e.g. pubchem IDs, KEGG IDs, etc.)
    REF_DB_NAME = 'compound_references'

    #: Name of KEGG database with models collection
    KEGG_DB_NAME = 'kegg'
