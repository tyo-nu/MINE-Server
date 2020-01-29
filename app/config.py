"""Configuration file. All settings are stored in Config class which gets
instantiated in __init__.py."""

import os

#: Path to 'app' directory
APP_DIR = os.path.abspath(os.path.dirname(__file__))


class Config(object):
    """Class to configure important filepaths and settings."""
    # ---------------------------- Main Settings ---------------------------- #
    # Main settings for MINE-Server API - things most likely to be changed

    #: Mode - options are 'production' or 'development'
    MODE = 'development'

    # ------------------------------ Filepaths ------------------------------ #
    # Local filepaths are defined here

    #: Path to "Positive Adducts full.txt"
    POS_ADDUCT_PATH = os.path.join(APP_DIR, 'data/Positive Adducts full.txt')

    #: Path to "Negative Adducts full.txt"
    NEG_ADDUCT_PATH = os.path.join(APP_DIR, 'data/Negative Adducts full.txt')

    # ------------------------------- MongoDB ------------------------------- #
    # Settings for interface with MongoDB

    #: URI to MINE MongoDB
    MONGO_URI = 'mongodb://localhost/'

    #: Name of KEGG database with models
    KEGG_DB_NAME = ''
