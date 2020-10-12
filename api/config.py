"""Configuration file. All settings are stored in Config class which gets
instantiated in __init__.py."""

import os
import minedatabase

try:
    from api.credentials import MONGO_USERNAME, MONGO_PASSWORD
except ImportError:
    raise FileNotFoundError("MINE-Server/api/credentials.py not found.")

APP_DIR = os.path.abspath(os.path.dirname(__file__))
MINEDB_DIR = os.path.dirname(minedatabase.__file__)


class Config(object):
    """Class to configure important filepaths and settings."""
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

    # ------------------------------- MongoDB ------------------------------- #
    # Settings for interface with MongoDB

    #: URI to MINE MongoDB
    MONGO_URI = f'mongodb://{MONGO_USERNAME}:{MONGO_PASSWORD}@localhost/'

    #: Name of KEGG database with models collection
    KEGG_DB_NAME = 'kegg'
